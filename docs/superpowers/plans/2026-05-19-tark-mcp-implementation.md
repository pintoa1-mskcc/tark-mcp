# TARK MCP Server Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a Python MCP server wrapping the Ensembl TARK REST API, exposing transcript, gene, sequence, diff, MANE, and release tools, including protein and CDS sequence comparison in `diff_transcripts`.

**Architecture:** `TarkClient` (httpx + pagination + retry + cache) → tool modules (business logic + model construction) → `server.py` (MCP stdio). All TARK API complexity (pagination, deduplication, 1-based coordinate conversion, nested response unpacking) is handled transparently. The `diff_transcripts` tool enriches structural diffs with CDS and protein sequences obtained by calling the diff endpoint plus the translation endpoint for each protein-coding transcript.

**Tech Stack:** Python 3.11+, `mcp` SDK, `httpx`, `pydantic` v2, `cachetools`, `tenacity`, `pytest`, `respx` (HTTP mocking), `pytest-asyncio`

---

## API Response Format Reference

These facts about the live TARK API inform the models and client code throughout:

- `GET /api/transcript/?stable_id=...&expand_all=true` returns paginated results
- `assembly` on transcript root is a nested object: `{"assembly_name": "GRCh38", ...}`
- `assembly` on exons, genes, translations is a plain string: `"GRCh38"`
- `sequence` is nested: `{"sequence": "ATCG...", "seq_checksum": "..."}`
- `exon_order` is the exon position field (not `order`)
- `transcript_release_set` is an array of objects with `release_date` strings
- Translations in transcript response have no `sequence` — use `GET /api/translation/?stable_id=<ENSP>&expand_all=true` for protein sequences
- Coordinates are 1-based → convert to 0-based: `model_start = api_start - 1`, `model_end = api_end`
- CDS offsets: `cds_start = len(five_prime_utr_seq)`, `cds_end = len(transcript_seq) - len(three_prime_utr_seq)`
- Diff endpoint: `GET /api/transcript/diff/?diff_me_stable_id=...&diff_with_stable_id=...` returns top-level `results` (diff flags object), `diff_me_transcript`, `diff_with_transcript`; `has_translation_seq_changed` is `null` for non-coding

---

## File Structure

```
tark-mcp/
├── pyproject.toml
├── src/
│   └── tark_mcp/
│       ├── __init__.py
│       ├── server.py          # MCP entry point, tool registration
│       ├── client.py          # TarkClient: HTTP, pagination, cache, retry
│       ├── models.py          # Pydantic models + coordinate normalization
│       └── tools/
│           ├── __init__.py
│           ├── releases.py    # get_releases
│           ├── transcripts.py # get_transcript, search_transcripts_by_region
│           ├── genes.py       # get_gene_transcripts
│           ├── sequences.py   # get_transcript_sequence, get_transcript_exons, get_protein_for_transcript
│           ├── mane.py        # get_mane_transcripts
│           └── diff.py        # diff_transcripts (with protein/CDS comparison)
└── tests/
    ├── conftest.py
    ├── unit/
    │   ├── test_models.py
    │   ├── test_client.py
    │   └── test_tools/
    │       ├── test_releases.py
    │       ├── test_transcripts.py
    │       ├── test_genes.py
    │       ├── test_sequences.py
    │       ├── test_mane.py
    │       └── test_diff.py
    └── integration/
        └── test_live_api.py
```

---

## Task 1: Project Scaffolding

**Files:**
- Create: `pyproject.toml`
- Create: `src/tark_mcp/__init__.py`
- Create: `src/tark_mcp/tools/__init__.py`
- Create: `tests/__init__.py`
- Create: `tests/unit/__init__.py`
- Create: `tests/unit/test_tools/__init__.py`
- Create: `tests/integration/__init__.py`
- Create: `tests/conftest.py`

- [ ] **Step 1: Create `pyproject.toml`**

```toml
[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[project]
name = "tark-mcp"
version = "0.1.0"
requires-python = ">=3.11"
dependencies = [
    "mcp>=1.0",
    "httpx>=0.27",
    "pydantic>=2.0",
    "cachetools>=5.0",
    "tenacity>=8.0",
]

[project.scripts]
tark-mcp = "tark_mcp.server:main"

[tool.hatch.build.targets.wheel]
packages = ["src/tark_mcp"]

[tool.pytest.ini_options]
asyncio_mode = "auto"
testpaths = ["tests"]

[project.optional-dependencies]
dev = [
    "pytest>=8.0",
    "pytest-asyncio>=0.23",
    "respx>=0.21",
    "httpx>=0.27",
]
```

- [ ] **Step 2: Create package skeleton**

```bash
mkdir -p src/tark_mcp/tools tests/unit/test_tools tests/integration
touch src/tark_mcp/__init__.py
touch src/tark_mcp/tools/__init__.py
touch tests/__init__.py
touch tests/unit/__init__.py
touch tests/unit/test_tools/__init__.py
touch tests/integration/__init__.py
```

- [ ] **Step 3: Create `tests/conftest.py`**

```python
import pytest
import httpx
import respx

BASE_URL = "https://tark.ensembl.org/api/"

# ---------------------------------------------------------------------------
# Minimal API fixture data — mirrors real TARK API response shapes
# ---------------------------------------------------------------------------

TRANSCRIPT_BRCA2_RAW = {
    "stable_id": "ENST00000380152",
    "stable_id_version": 7,
    "assembly": {"assembly_id": 1, "assembly_name": "GRCh38", "genome": 1, "session": 1},
    "loc_start": 32315475,   # 1-based → model stores 32315474
    "loc_end": 32400266,
    "loc_strand": 1,
    "loc_region": "13",
    "biotype": "protein_coding",
    "sequence": {"sequence": "ATCGATCGATCGATCGATCGATCGATCGATCG", "seq_checksum": "ABC"},
    "five_prime_utr_seq": "ATCG",
    "three_prime_utr_seq": "CG",
    "three_prime_utr_start": 32398771,
    "three_prime_utr_end": 32400266,
    "five_prime_utr_start": 32315475,
    "five_prime_utr_end": 32316460,
    "transcript_release_set": [
        {"assembly": "GRCh38", "shortname": "110", "description": "Ensembl release",
         "release_date": "2023-04-01", "source": "Ensembl"},
    ],
    "genes": [
        {
            "stable_id": "ENSG00000139618",
            "stable_id_version": 15,
            "assembly": "GRCh38",
            "loc_start": 32315475,
            "loc_end": 32400266,
            "loc_strand": 1,
            "loc_region": "13",
            "name": "BRCA2",
        }
    ],
    "translations": [
        {
            "stable_id": "ENSP00000369497",
            "stable_id_version": 3,
            "assembly": "GRCh38",
            "loc_start": 32316461,
            "loc_end": 32398770,
            "loc_strand": 1,
            "loc_region": "13",
        }
    ],
    "exons": [
        {
            "stable_id": "ENSE00001184784",
            "stable_id_version": 4,
            "assembly": "GRCh38",
            "loc_start": 32315475,
            "loc_end": 32315667,
            "loc_strand": 1,
            "loc_region": "13",
            "exon_order": 1,
            "transcript_stable_id": "ENST00000380152",
            "transcript_stable_id_version": 7,
        }
    ],
}

TRANSLATION_BRCA2_RAW = {
    "stable_id": "ENSP00000369497",
    "stable_id_version": 3,
    "assembly": {"assembly_id": 1, "assembly_name": "GRCh38", "genome": 1, "session": 1},
    "loc_start": 32316461,
    "loc_end": 32398770,
    "loc_strand": 1,
    "loc_region": "13",
    "sequence": {"sequence": "MPIGSKERP", "seq_checksum": "DEF"},
}

TRANSCRIPT_NONCODING_RAW = {
    "stable_id": "ENST00000614536",
    "stable_id_version": 1,
    "assembly": {"assembly_id": 1, "assembly_name": "GRCh38", "genome": 1, "session": 1},
    "loc_start": 54529654,
    "loc_end": 54529961,
    "loc_strand": 1,
    "loc_region": "19",
    "biotype": "unprocessed_pseudogene",
    "sequence": {"sequence": "GGCTTGTTCACA", "seq_checksum": "GHI"},
    "five_prime_utr_seq": None,
    "three_prime_utr_seq": None,
    "three_prime_utr_start": None,
    "three_prime_utr_end": None,
    "five_prime_utr_start": None,
    "five_prime_utr_end": None,
    "transcript_release_set": [
        {"assembly": "GRCh38", "shortname": "92", "description": "Ensembl release",
         "release_date": "2018-04-05", "source": "Ensembl"},
    ],
    "genes": [],
    "translations": [],
    "exons": [
        {
            "stable_id": "ENSE00003719241",
            "stable_id_version": 1,
            "assembly": "GRCh38",
            "loc_start": 54529654,
            "loc_end": 54529961,
            "loc_strand": 1,
            "loc_region": "19",
            "exon_order": 1,
            "transcript_stable_id": "ENST00000614536",
            "transcript_stable_id_version": 1,
        }
    ],
}

DIFF_RESPONSE_RAW = {
    "count": 1,
    "next": None,
    "previous": None,
    "results": {
        "diff_me_stable_id": "ENST00000380152",
        "diff_with_stable_id": "ENST00000614536",
        "diff_me_stable_id_version": 7,
        "diff_with_stable_id_version": 1,
        "diff_me_assembly": {"assembly_name": "GRCh38"},
        "diff_with_assembly": {"assembly_name": "GRCh38"},
        "has_stable_id_changed": True,
        "has_transcript_changed": True,
        "has_seq_changed": True,
        "has_exon_set_changed": True,
        "has_translation_stable_id_changed": True,
        "has_translation_seq_changed": None,
    },
    "diff_me_transcript": TRANSCRIPT_BRCA2_RAW,
    "diff_with_transcript": TRANSCRIPT_NONCODING_RAW,
}

RELEASE_LIST_RAW = [
    {
        "shortname": "110",
        "description": "Ensembl release 110",
        "release_date": "2023-04-01",
        "assembly": "GRCh38",
        "source": "Ensembl",
    }
]

MANE_LIST_RESPONSE_RAW = {
    "count": 1,
    "next": None,
    "previous": None,
    "results": [TRANSCRIPT_BRCA2_RAW],
}
```

- [ ] **Step 4: Install dependencies**

```bash
pip install -e ".[dev]"
```

Expected: `Successfully installed tark-mcp-0.1.0` (plus dependencies)

- [ ] **Step 5: Commit**

```bash
git add pyproject.toml src/ tests/
git commit -m "chore: scaffold project structure and conftest fixtures"
```

---

## Task 2: Pydantic Models

**Files:**
- Create: `src/tark_mcp/models.py`
- Create: `tests/unit/test_models.py`

- [ ] **Step 1: Write failing tests for models**

Create `tests/unit/test_models.py`:

```python
import pytest
from tark_mcp.models import (
    Exon, Gene, Translation, Transcript, TranscriptDiff, ExonDiff, Release
)


# ---------------------------------------------------------------------------
# Exon
# ---------------------------------------------------------------------------

def test_exon_coordinate_normalization():
    raw = {
        "stable_id": "ENSE00001184784",
        "stable_id_version": 4,
        "transcript_stable_id": "ENST00000380152",
        "transcript_stable_id_version": 7,
        "assembly": "GRCh38",
        "loc_start": 32315475,   # 1-based
        "loc_end": 32315667,
        "loc_strand": 1,
        "loc_region": "13",
        "exon_order": 1,
    }
    exon = Exon.model_validate(raw)
    assert exon.loc_start == 32315474  # 0-based
    assert exon.loc_end == 32315667    # unchanged (exclusive end)
    assert exon.order == 1
    assert exon.stable_id == "ENSE00001184784"
    assert exon.stable_id_version == 4
    assert exon.transcript_stable_id == "ENST00000380152"
    assert exon.transcript_stable_id_version == 7
    assert exon.assembly == "GRCh38"


# ---------------------------------------------------------------------------
# Gene
# ---------------------------------------------------------------------------

def test_gene_coordinate_normalization():
    raw = {
        "stable_id": "ENSG00000139618",
        "stable_id_version": 15,
        "assembly": "GRCh38",
        "loc_start": 32315475,
        "loc_end": 32400266,
        "loc_strand": 1,
        "loc_region": "13",
        "name": "BRCA2",
    }
    gene = Gene.model_validate(raw)
    assert gene.loc_start == 32315474
    assert gene.loc_end == 32400266
    assert gene.name == "BRCA2"


def test_gene_name_can_be_none():
    raw = {
        "stable_id": "ENSG00000274968",
        "stable_id_version": 1,
        "assembly": "GRCh38",
        "loc_start": 100,
        "loc_end": 200,
        "loc_strand": 1,
        "loc_region": "19",
        "name": None,
    }
    gene = Gene.model_validate(raw)
    assert gene.name is None


# ---------------------------------------------------------------------------
# Translation
# ---------------------------------------------------------------------------

def test_translation_with_sequence():
    raw = {
        "stable_id": "ENSP00000369497",
        "stable_id_version": 3,
        "transcript_stable_id": "ENST00000380152",
        "transcript_stable_id_version": 7,
        "assembly": "GRCh38",
        "sequence": {"sequence": "MPIGSKERP", "seq_checksum": "DEF"},
    }
    t = Translation.model_validate(raw)
    assert t.stable_id == "ENSP00000369497"
    assert t.sequence == "MPIGSKERP"


def test_translation_without_sequence():
    raw = {
        "stable_id": "ENSP00000369497",
        "stable_id_version": 3,
        "transcript_stable_id": "ENST00000380152",
        "transcript_stable_id_version": 7,
        "assembly": "GRCh38",
    }
    t = Translation.model_validate(raw)
    assert t.sequence is None


def test_translation_nested_assembly_is_ignored():
    """Translation stable_id_version and assembly come as plain strings in transcript context."""
    raw = {
        "stable_id": "ENSP00000369497",
        "stable_id_version": 3,
        "transcript_stable_id": "ENST00000380152",
        "transcript_stable_id_version": 7,
        "assembly": "GRCh38",
    }
    t = Translation.model_validate(raw)
    assert t.assembly == "GRCh38"


# ---------------------------------------------------------------------------
# Transcript
# ---------------------------------------------------------------------------

def test_transcript_coordinate_normalization():
    raw = {
        "stable_id": "ENST00000380152",
        "stable_id_version": 7,
        "assembly": {"assembly_name": "GRCh38"},
        "loc_start": 32315475,
        "loc_end": 32400266,
        "loc_strand": 1,
        "loc_region": "13",
        "biotype": "protein_coding",
        "sequence": {"sequence": "ATCGATCGATCGATCGATCGATCGATCGATCG"},
        "five_prime_utr_seq": "ATCG",
        "three_prime_utr_seq": "CG",
        "transcript_release_set": [
            {"release_date": "2023-04-01", "assembly": "GRCh38", "shortname": "110",
             "description": "Ensembl release", "source": "Ensembl"}
        ],
        "genes": [],
        "translations": [],
        "exons": [],
    }
    t = Transcript.model_validate(raw)
    assert t.assembly == "GRCh38"
    assert t.loc_start == 32315474  # 0-based
    assert t.loc_end == 32400266
    assert t.sequence == "ATCGATCGATCGATCGATCGATCGATCGATCG"


def test_transcript_cds_boundaries_computed_from_utrs():
    """cds_start = len(5'UTR), cds_end = len(seq) - len(3'UTR)"""
    raw = {
        "stable_id": "ENST00000380152",
        "stable_id_version": 7,
        "assembly": {"assembly_name": "GRCh38"},
        "loc_start": 1,
        "loc_end": 100,
        "loc_strand": 1,
        "loc_region": "13",
        "biotype": "protein_coding",
        "sequence": {"sequence": "AAAACCCCGGGG"},  # len=12
        "five_prime_utr_seq": "AAAA",               # len=4  → cds_start=4
        "three_prime_utr_seq": "GGGG",              # len=4  → cds_end=12-4=8
        "transcript_release_set": [],
        "genes": [],
        "translations": [],
        "exons": [],
    }
    t = Transcript.model_validate(raw)
    assert t.cds_start == 4
    assert t.cds_end == 8


def test_transcript_cds_none_when_no_utrs():
    raw = {
        "stable_id": "ENST00000614536",
        "stable_id_version": 1,
        "assembly": {"assembly_name": "GRCh38"},
        "loc_start": 1,
        "loc_end": 100,
        "loc_strand": 1,
        "loc_region": "19",
        "biotype": "unprocessed_pseudogene",
        "sequence": {"sequence": "GGCTTGTTCACA"},
        "five_prime_utr_seq": None,
        "three_prime_utr_seq": None,
        "transcript_release_set": [],
        "genes": [],
        "translations": [],
        "exons": [],
    }
    t = Transcript.model_validate(raw)
    assert t.cds_start is None
    assert t.cds_end is None


def test_transcript_latest_release_date():
    """The model keeps the most-recent release_date from transcript_release_set."""
    raw = {
        "stable_id": "ENST00000380152",
        "stable_id_version": 7,
        "assembly": {"assembly_name": "GRCh38"},
        "loc_start": 1, "loc_end": 100, "loc_strand": 1, "loc_region": "13",
        "biotype": "protein_coding",
        "sequence": {"sequence": "ATCG"},
        "five_prime_utr_seq": None, "three_prime_utr_seq": None,
        "transcript_release_set": [
            {"release_date": "2014-05-01", "assembly": "GRCh38", "shortname": "76",
             "description": "Ensembl release", "source": "Ensembl"},
            {"release_date": "2023-04-01", "assembly": "GRCh38", "shortname": "110",
             "description": "Ensembl release", "source": "Ensembl"},
        ],
        "genes": [], "translations": [], "exons": [],
    }
    t = Transcript.model_validate(raw)
    assert t.latest_release_date == "2023-04-01"


def test_transcript_exons_parsed_correctly():
    raw = {
        "stable_id": "ENST00000380152",
        "stable_id_version": 7,
        "assembly": {"assembly_name": "GRCh38"},
        "loc_start": 1, "loc_end": 1000, "loc_strand": 1, "loc_region": "13",
        "biotype": "protein_coding",
        "sequence": {"sequence": "ATCG"},
        "five_prime_utr_seq": None, "three_prime_utr_seq": None,
        "transcript_release_set": [],
        "genes": [],
        "translations": [],
        "exons": [
            {
                "stable_id": "ENSE00001184784",
                "stable_id_version": 4,
                "transcript_stable_id": "ENST00000380152",
                "transcript_stable_id_version": 7,
                "assembly": "GRCh38",
                "loc_start": 32315475,
                "loc_end": 32315667,
                "loc_strand": 1,
                "loc_region": "13",
                "exon_order": 1,
            }
        ],
    }
    t = Transcript.model_validate(raw)
    assert len(t.exons) == 1
    assert t.exons[0].order == 1
    assert t.exons[0].loc_start == 32315474  # 0-based


# ---------------------------------------------------------------------------
# Release
# ---------------------------------------------------------------------------

def test_release_parsed():
    raw = {"shortname": "110", "description": "Ensembl release 110",
           "release_date": "2023-04-01", "assembly": "GRCh38", "source": "Ensembl"}
    r = Release.model_validate(raw)
    assert r.shortname == "110"
    assert r.assembly == "GRCh38"
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest tests/unit/test_models.py -v 2>&1 | head -30
```

Expected: `ModuleNotFoundError: No module named 'tark_mcp'` (or similar import errors)

- [ ] **Step 3: Create `src/tark_mcp/models.py`**

```python
from __future__ import annotations
from pydantic import BaseModel, model_validator


def _normalize_start(v: int | None) -> int | None:
    """Convert 1-based inclusive start → 0-based."""
    return v - 1 if v is not None else None


class Exon(BaseModel):
    stable_id: str
    stable_id_version: int
    transcript_stable_id: str
    transcript_stable_id_version: int
    assembly: str
    order: int
    loc_region: str
    loc_start: int
    loc_end: int
    loc_strand: int

    @model_validator(mode="before")
    @classmethod
    def _normalize(cls, data: dict) -> dict:
        data = dict(data)
        if "loc_start" in data and data["loc_start"] is not None:
            data["loc_start"] = data["loc_start"] - 1
        if "exon_order" in data:
            data["order"] = data.pop("exon_order")
        return data


class Gene(BaseModel):
    stable_id: str
    stable_id_version: int
    name: str | None
    loc_region: str
    loc_start: int
    loc_end: int
    loc_strand: int
    assembly: str

    @model_validator(mode="before")
    @classmethod
    def _normalize(cls, data: dict) -> dict:
        data = dict(data)
        if "loc_start" in data and data["loc_start"] is not None:
            data["loc_start"] = data["loc_start"] - 1
        return data


class Translation(BaseModel):
    stable_id: str
    stable_id_version: int
    transcript_stable_id: str
    transcript_stable_id_version: int
    assembly: str
    sequence: str | None = None

    @model_validator(mode="before")
    @classmethod
    def _normalize(cls, data: dict) -> dict:
        data = dict(data)
        # Flatten nested assembly object if present
        if isinstance(data.get("assembly"), dict):
            data["assembly"] = data["assembly"]["assembly_name"]
        # Flatten nested sequence object if present
        seq = data.get("sequence")
        if isinstance(seq, dict):
            data["sequence"] = seq.get("sequence")
        return data


class Transcript(BaseModel):
    stable_id: str
    stable_id_version: int
    assembly: str
    biotype: str
    loc_region: str
    loc_start: int
    loc_end: int
    loc_strand: int
    cds_start: int | None
    cds_end: int | None
    exons: list[Exon]
    genes: list[Gene]
    translations: list[Translation]
    sequence: str | None
    latest_release_date: str | None = None

    @model_validator(mode="before")
    @classmethod
    def _normalize(cls, data: dict) -> dict:
        data = dict(data)

        # Flatten nested assembly object
        if isinstance(data.get("assembly"), dict):
            data["assembly"] = data["assembly"]["assembly_name"]

        # Convert 1-based start → 0-based
        if "loc_start" in data and data["loc_start"] is not None:
            data["loc_start"] = data["loc_start"] - 1

        # Flatten nested sequence object
        seq = data.get("sequence")
        if isinstance(seq, dict):
            data["sequence"] = seq.get("sequence")

        # Compute CDS boundaries from UTR sequences
        five_utr = data.get("five_prime_utr_seq")
        three_utr = data.get("three_prime_utr_seq")
        transcript_seq = data.get("sequence") if not isinstance(seq, dict) else (seq or {}).get("sequence")
        if isinstance(data.get("sequence"), str):
            transcript_seq = data["sequence"]

        if five_utr is not None and transcript_seq is not None:
            data["cds_start"] = len(five_utr)
            data["cds_end"] = len(transcript_seq) - (len(three_utr) if three_utr else 0)
        else:
            data.setdefault("cds_start", None)
            data.setdefault("cds_end", None)

        # Compute latest_release_date from transcript_release_set array
        release_set = data.get("transcript_release_set", [])
        if isinstance(release_set, list) and release_set:
            dates = [r["release_date"] for r in release_set if r.get("release_date")]
            data["latest_release_date"] = max(dates) if dates else None
        elif isinstance(release_set, dict):
            data["latest_release_date"] = release_set.get("release_date")
        else:
            data["latest_release_date"] = None

        return data


class Release(BaseModel):
    shortname: str
    description: str | None
    release_date: str
    assembly: str
    source: str


class ExonDiff(BaseModel):
    order: int
    change: str  # "added", "removed", "modified", "unchanged"
    ref_coords: tuple[int, int] | None
    candidate_coords: tuple[int, int] | None


class TranscriptDiff(BaseModel):
    reference_stable_id: str
    candidate_stable_id: str
    reference_assembly: str
    candidate_assembly: str
    biotype_changed: bool
    cds_changed: bool
    exon_count_changed: bool
    sequence_changed: bool
    exon_diffs: list[ExonDiff]
    reference_protein_coding: bool
    candidate_protein_coding: bool
    cds_sequence_changed: bool | None
    ref_cds_sequence: str | None
    candidate_cds_sequence: str | None
    protein_sequence_changed: bool | None
    ref_protein_sequence: str | None
    candidate_protein_sequence: str | None
```

- [ ] **Step 4: Fix the CDS computation — the validator runs before `sequence` is flattened**

The `_normalize` validator flattens `sequence` from a dict to a string in the same pass. The CDS computation must read the already-flattened string, not the original dict. The code above handles this by checking `isinstance(data.get("sequence"), str)` after the flatten. Verify logic is correct:

The flatten block sets `data["sequence"] = seq.get("sequence")` (a string). Then the CDS block reads `data["sequence"]` which is now a string. This is correct as written.

- [ ] **Step 5: Run tests to verify they pass**

```bash
pytest tests/unit/test_models.py -v
```

Expected: all tests PASS

- [ ] **Step 6: Commit**

```bash
git add src/tark_mcp/models.py tests/unit/test_models.py
git commit -m "feat: add Pydantic models with coordinate normalization and CDS computation"
```

---

## Task 3: TarkClient

**Files:**
- Create: `src/tark_mcp/client.py`
- Create: `tests/unit/test_client.py`

- [ ] **Step 1: Write failing tests for TarkClient**

Create `tests/unit/test_client.py`:

```python
import pytest
import httpx
import respx

from tark_mcp.client import TarkClient

BASE = "https://tark.ensembl.org/api/"


@pytest.fixture
def client():
    return TarkClient()


@respx.mock
@pytest.mark.asyncio
async def test_get_single_page(client):
    """Fetches a single page and returns results list."""
    respx.get(BASE + "release/nopagination/").mock(
        return_value=httpx.Response(200, json=[{"shortname": "110"}])
    )
    result = await client.get("release/nopagination/")
    assert result == [{"shortname": "110"}]


@respx.mock
@pytest.mark.asyncio
async def test_get_paginates_automatically(client):
    """Follows next links and aggregates all pages."""
    page1 = {"count": 2, "next": BASE + "transcript/?page=2", "previous": None,
             "results": [{"stable_id": "ENST000001"}]}
    page2 = {"count": 2, "next": None, "previous": BASE + "transcript/",
             "results": [{"stable_id": "ENST000002"}]}
    respx.get(BASE + "transcript/").mock(return_value=httpx.Response(200, json=page1))
    respx.get(BASE + "transcript/?page=2").mock(return_value=httpx.Response(200, json=page2))

    results = await client.get("transcript/")
    assert len(results) == 2
    assert results[0]["stable_id"] == "ENST000001"
    assert results[1]["stable_id"] == "ENST000002"


@respx.mock
@pytest.mark.asyncio
async def test_get_rewrites_http_to_https(client):
    """HTTPS is enforced even if next link is http://."""
    page1 = {"count": 2, "next": "http://tark.ensembl.org/api/transcript/?page=2",
             "previous": None, "results": [{"stable_id": "ENST000001"}]}
    page2 = {"count": 2, "next": None, "previous": None,
             "results": [{"stable_id": "ENST000002"}]}
    respx.get(BASE + "transcript/").mock(return_value=httpx.Response(200, json=page1))
    respx.get(BASE + "transcript/?page=2").mock(return_value=httpx.Response(200, json=page2))

    results = await client.get("transcript/")
    assert len(results) == 2


@respx.mock
@pytest.mark.asyncio
async def test_get_returns_dict_for_non_paginated_response(client):
    """When response is a plain dict (e.g. diff endpoint), return it directly."""
    payload = {"results": {"diff_me_stable_id": "X"}, "diff_me_transcript": {}}
    respx.get(BASE + "transcript/diff/").mock(return_value=httpx.Response(200, json=payload))

    result = await client.get_raw("transcript/diff/")
    assert result["results"]["diff_me_stable_id"] == "X"


@respx.mock
@pytest.mark.asyncio
async def test_404_returns_empty_list(client):
    respx.get(BASE + "transcript/").mock(return_value=httpx.Response(404, json={}))
    result = await client.get("transcript/")
    assert result == []


@respx.mock
@pytest.mark.asyncio
async def test_http_error_raises_mcp_error(client):
    from mcp.server.fastmcp import FastMCP
    respx.get(BASE + "transcript/").mock(return_value=httpx.Response(400, json={"detail": "bad"}))
    with pytest.raises(Exception, match="400"):
        await client.get("transcript/")


@respx.mock
@pytest.mark.asyncio
async def test_cache_hit_avoids_second_request(client):
    """Second identical request is served from cache."""
    route = respx.get(BASE + "release/nopagination/").mock(
        return_value=httpx.Response(200, json=[{"shortname": "110"}])
    )
    await client.get("release/nopagination/")
    await client.get("release/nopagination/")
    assert route.call_count == 1
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest tests/unit/test_client.py -v 2>&1 | head -20
```

Expected: `ImportError: cannot import name 'TarkClient'`

- [ ] **Step 3: Create `src/tark_mcp/client.py`**

```python
from __future__ import annotations
import os
from typing import Any

import httpx
from cachetools import TTLCache
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception

BASE_URL = os.environ.get("TARK_BASE_URL", "https://tark.ensembl.org/api/")
CACHE_TTL = int(os.environ.get("TARK_CACHE_TTL", "3600"))
REQUEST_TIMEOUT = int(os.environ.get("TARK_REQUEST_TIMEOUT", "30"))
MAX_RETRIES = int(os.environ.get("TARK_MAX_RETRIES", "3"))

_cache: TTLCache = TTLCache(maxsize=512, ttl=CACHE_TTL)


def _is_transient(exc: BaseException) -> bool:
    if isinstance(exc, httpx.HTTPStatusError):
        return exc.response.status_code >= 500
    return isinstance(exc, (httpx.NetworkError, httpx.TimeoutException))


class TarkClient:
    def __init__(self) -> None:
        self._http = httpx.AsyncClient(base_url=BASE_URL, timeout=REQUEST_TIMEOUT)

    @retry(
        stop=stop_after_attempt(MAX_RETRIES),
        wait=wait_exponential(multiplier=1, min=1, max=10),
        retry=retry_if_exception(_is_transient),
        reraise=True,
    )
    async def _fetch(self, url: str, params: dict | None = None) -> Any:
        cache_key = url + str(sorted((params or {}).items()))
        if cache_key in _cache:
            return _cache[cache_key]

        response = await self._http.get(url, params=params)

        if response.status_code == 404:
            return None
        if response.status_code >= 400:
            raise Exception(f"TARK API error {response.status_code}: {response.text[:200]}")

        data = response.json()
        _cache[cache_key] = data
        return data

    async def get(self, path: str, params: dict | None = None) -> list[dict]:
        """Fetch a paginated endpoint and return aggregated results list."""
        url = BASE_URL + path
        results: list[dict] = []
        while url:
            data = await self._fetch(url, params if url == BASE_URL + path else None)
            if data is None:
                break
            if isinstance(data, list):
                results.extend(data)
                break
            if "results" in data and isinstance(data["results"], list):
                results.extend(data["results"])
                next_url = data.get("next")
                if next_url:
                    # Enforce HTTPS
                    url = next_url.replace("http://", "https://")
                    params = None
                else:
                    break
            else:
                # Single-result or non-list results (shouldn't paginate further)
                results.append(data.get("results", data))
                break
        return results

    async def get_raw(self, path: str, params: dict | None = None) -> dict:
        """Fetch an endpoint that returns a plain dict (e.g. diff endpoint)."""
        url = BASE_URL + path
        data = await self._fetch(url, params)
        return data or {}

    async def aclose(self) -> None:
        await self._http.aclose()
```

- [ ] **Step 4: Run tests**

```bash
pytest tests/unit/test_client.py -v
```

Expected: all tests PASS (or debug any failures)

- [ ] **Step 5: Commit**

```bash
git add src/tark_mcp/client.py tests/unit/test_client.py
git commit -m "feat: add TarkClient with pagination, caching, retry"
```

---

## Task 4: Releases Tool

**Files:**
- Create: `src/tark_mcp/tools/releases.py`
- Create: `tests/unit/test_tools/test_releases.py`

- [ ] **Step 1: Write failing test**

Create `tests/unit/test_tools/test_releases.py`:

```python
import pytest
import httpx
import respx

from tark_mcp.client import TarkClient
from tark_mcp.tools.releases import get_releases
from tests.conftest import RELEASE_LIST_RAW

BASE = "https://tark.ensembl.org/api/"


@respx.mock
@pytest.mark.asyncio
async def test_get_releases_returns_release_list():
    client = TarkClient()
    respx.get(BASE + "release/nopagination/").mock(
        return_value=httpx.Response(200, json=RELEASE_LIST_RAW)
    )
    releases = await get_releases(client)
    assert len(releases) == 1
    assert releases[0].shortname == "110"
    assert releases[0].assembly == "GRCh38"
    assert releases[0].release_date == "2023-04-01"
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest tests/unit/test_tools/test_releases.py -v 2>&1 | head -10
```

Expected: `ImportError`

- [ ] **Step 3: Create `src/tark_mcp/tools/releases.py`**

```python
from __future__ import annotations
from tark_mcp.client import TarkClient
from tark_mcp.models import Release


async def get_releases(client: TarkClient) -> list[Release]:
    data = await client.get("release/nopagination/")
    return [Release.model_validate(r) for r in data]
```

- [ ] **Step 4: Run test to verify it passes**

```bash
pytest tests/unit/test_tools/test_releases.py -v
```

Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/tark_mcp/tools/releases.py tests/unit/test_tools/test_releases.py
git commit -m "feat: add get_releases tool"
```

---

## Task 5: Transcripts Tool

**Files:**
- Create: `src/tark_mcp/tools/transcripts.py`
- Create: `tests/unit/test_tools/test_transcripts.py`

- [ ] **Step 1: Write failing tests**

Create `tests/unit/test_tools/test_transcripts.py`:

```python
import pytest
import httpx
import respx

from tark_mcp.client import TarkClient
from tark_mcp.tools.transcripts import get_transcript, search_transcripts_by_region
from tests.conftest import TRANSCRIPT_BRCA2_RAW, TRANSCRIPT_NONCODING_RAW

BASE = "https://tark.ensembl.org/api/"

PAGINATED_TWO = {
    "count": 2,
    "next": None,
    "previous": None,
    "results": [
        TRANSCRIPT_BRCA2_RAW,
        {**TRANSCRIPT_BRCA2_RAW, "stable_id_version": 6,
         "transcript_release_set": [{"assembly": "GRCh38", "shortname": "109",
                                     "description": "Ensembl release",
                                     "release_date": "2022-09-01", "source": "Ensembl"}]},
    ],
}


@respx.mock
@pytest.mark.asyncio
async def test_get_transcript_returns_latest_version():
    """When multiple versions exist, deduplicate and return the latest."""
    client = TarkClient()
    respx.get(BASE + "transcript/").mock(
        return_value=httpx.Response(200, json=PAGINATED_TWO)
    )
    result = await get_transcript("ENST00000380152", client=client)
    assert result is not None
    assert result.stable_id == "ENST00000380152"
    # Should return version 7 (latest release_date 2023-04-01)
    assert result.stable_id_version == 7


@respx.mock
@pytest.mark.asyncio
async def test_get_transcript_strips_version_suffix():
    """ENST00000380152.7 → queries for stable_id=ENST00000380152."""
    client = TarkClient()
    route = respx.get(BASE + "transcript/").mock(
        return_value=httpx.Response(200, json={
            "count": 1, "next": None, "previous": None,
            "results": [TRANSCRIPT_BRCA2_RAW]
        })
    )
    await get_transcript("ENST00000380152.7", client=client)
    assert "stable_id=ENST00000380152" in str(route.calls[0].request.url)


@respx.mock
@pytest.mark.asyncio
async def test_get_transcript_not_found_returns_none():
    client = TarkClient()
    respx.get(BASE + "transcript/").mock(
        return_value=httpx.Response(200, json={"count": 0, "next": None,
                                               "previous": None, "results": []})
    )
    result = await get_transcript("ENST00000999999", client=client)
    assert result is None


@respx.mock
@pytest.mark.asyncio
async def test_get_transcript_assembly_both_returns_list():
    """assembly='both' fans out to GRCh37 and GRCh38, returns list."""
    client = TarkClient()
    t38 = {**TRANSCRIPT_BRCA2_RAW}
    t37 = {**TRANSCRIPT_BRCA2_RAW,
           "assembly": {"assembly_name": "GRCh37", "assembly_id": 2, "genome": 2, "session": 1}}
    respx.get(BASE + "transcript/").mock(side_effect=[
        httpx.Response(200, json={"count": 1, "next": None, "previous": None, "results": [t38]}),
        httpx.Response(200, json={"count": 1, "next": None, "previous": None, "results": [t37]}),
    ])
    result = await get_transcript("ENST00000380152", assembly="both", client=client)
    assert isinstance(result, list)
    assert len(result) == 2


@respx.mock
@pytest.mark.asyncio
async def test_search_transcripts_by_region():
    """0-based input is converted to 1-based for the API; chr prefix is stripped."""
    client = TarkClient()
    route = respx.get(BASE + "transcript/").mock(
        return_value=httpx.Response(200, json={
            "count": 1, "next": None, "previous": None,
            "results": [TRANSCRIPT_BRCA2_RAW]
        })
    )
    results = await search_transcripts_by_region(
        "chr13", start=32315474, end=32400266, client=client
    )
    assert len(results) == 1
    url = str(route.calls[0].request.url)
    # 0-based 32315474 → 1-based 32315475 in query
    assert "loc_start=32315475" in url
    assert "loc_region=13" in url   # chr prefix stripped
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest tests/unit/test_tools/test_transcripts.py -v 2>&1 | head -10
```

Expected: `ImportError`

- [ ] **Step 3: Create `src/tark_mcp/tools/transcripts.py`**

```python
from __future__ import annotations
import asyncio
from tark_mcp.client import TarkClient
from tark_mcp.models import Transcript


def _strip_version(stable_id: str) -> tuple[str, int | None]:
    """Split 'ENST00000380152.7' → ('ENST00000380152', 7). Returns None version if absent."""
    if "." in stable_id:
        parts = stable_id.rsplit(".", 1)
        try:
            return parts[0], int(parts[1])
        except ValueError:
            pass
    return stable_id, None


def _deduplicate(transcripts: list[Transcript]) -> list[Transcript]:
    """For each (assembly, stable_id, stable_id_version) keep the one with the most recent release."""
    best: dict[tuple, Transcript] = {}
    for t in transcripts:
        key = (t.assembly, t.stable_id, t.stable_id_version)
        existing = best.get(key)
        if existing is None or (t.latest_release_date or "") > (existing.latest_release_date or ""):
            best[key] = t
    return list(best.values())


async def _fetch_for_assembly(
    stable_id: str,
    version: int | None,
    assembly: str,
    client: TarkClient,
) -> list[Transcript]:
    params: dict = {"stable_id": stable_id, "expand_all": "true",
                    "assembly_name": assembly}
    if version is not None:
        params["stable_id_version"] = version
    data = await client.get("transcript/", params=params)
    transcripts = [Transcript.model_validate(r) for r in data]
    return _deduplicate(transcripts)


async def get_transcript(
    stable_id: str,
    assembly: str = "GRCh38",
    client: TarkClient | None = None,
) -> Transcript | list[Transcript] | None:
    if client is None:
        client = TarkClient()
    sid, version = _strip_version(stable_id)

    if assembly == "both":
        results = await asyncio.gather(
            _fetch_for_assembly(sid, version, "GRCh38", client),
            _fetch_for_assembly(sid, version, "GRCh37", client),
        )
        combined = results[0] + results[1]
        return combined if combined else None

    transcripts = await _fetch_for_assembly(sid, version, assembly, client)
    if not transcripts:
        return None
    # Return the single most-recently-released record for this assembly
    return max(transcripts, key=lambda t: t.latest_release_date or "")


async def search_transcripts_by_region(
    region: str,
    start: int,
    end: int,
    assembly: str = "GRCh38",
    client: TarkClient | None = None,
) -> list[Transcript]:
    if client is None:
        client = TarkClient()
    # Strip chr prefix
    loc_region = region.removeprefix("chr")
    # Convert 0-based → 1-based
    params = {
        "assembly_name": assembly,
        "loc_region": loc_region,
        "loc_start": start + 1,
        "loc_end": end,
        "expand": "transcript_release_set",
    }
    if assembly == "both":
        results = await asyncio.gather(
            client.get("transcript/", {**params, "assembly_name": "GRCh38"}),
            client.get("transcript/", {**params, "assembly_name": "GRCh37"}),
        )
        data = results[0] + results[1]
    else:
        data = await client.get("transcript/", params)

    return _deduplicate([Transcript.model_validate(r) for r in data])
```

- [ ] **Step 4: Run tests**

```bash
pytest tests/unit/test_tools/test_transcripts.py -v
```

Expected: all PASS

- [ ] **Step 5: Commit**

```bash
git add src/tark_mcp/tools/transcripts.py tests/unit/test_tools/test_transcripts.py
git commit -m "feat: add get_transcript and search_transcripts_by_region tools"
```

---

## Task 6: Genes Tool

**Files:**
- Create: `src/tark_mcp/tools/genes.py`
- Create: `tests/unit/test_tools/test_genes.py`

- [ ] **Step 1: Write failing tests**

Create `tests/unit/test_tools/test_genes.py`:

```python
import pytest
import httpx
import respx

from tark_mcp.client import TarkClient
from tark_mcp.tools.genes import get_gene_transcripts
from tests.conftest import TRANSCRIPT_BRCA2_RAW, TRANSCRIPT_NONCODING_RAW

BASE = "https://tark.ensembl.org/api/"

MIXED_ASSEMBLY_RESULTS = {
    "count": 2, "next": None, "previous": None,
    "results": [
        TRANSCRIPT_BRCA2_RAW,
        {**TRANSCRIPT_NONCODING_RAW,
         "assembly": {"assembly_name": "GRCh37", "assembly_id": 2, "genome": 2, "session": 1}},
    ],
}


@respx.mock
@pytest.mark.asyncio
async def test_get_gene_transcripts_filters_by_assembly():
    """Search endpoint returns mixed assemblies; client-side filters to GRCh38."""
    client = TarkClient()
    respx.get(BASE + "transcript/search/").mock(
        return_value=httpx.Response(200, json=MIXED_ASSEMBLY_RESULTS)
    )
    results = await get_gene_transcripts("BRCA2", assembly="GRCh38", client=client)
    assert all(t.assembly == "GRCh38" for t in results)
    assert len(results) == 1


@respx.mock
@pytest.mark.asyncio
async def test_get_gene_transcripts_assembly_both_returns_all():
    client = TarkClient()
    respx.get(BASE + "transcript/search/").mock(
        return_value=httpx.Response(200, json=MIXED_ASSEMBLY_RESULTS)
    )
    results = await get_gene_transcripts("BRCA2", assembly="both", client=client)
    assert len(results) == 2


@respx.mock
@pytest.mark.asyncio
async def test_get_gene_transcripts_passes_identifier_in_query():
    client = TarkClient()
    route = respx.get(BASE + "transcript/search/").mock(
        return_value=httpx.Response(200, json={"count": 0, "next": None, "previous": None,
                                               "results": []})
    )
    await get_gene_transcripts("ENSG00000139618", client=client)
    assert "identifier_field=ENSG00000139618" in str(route.calls[0].request.url)
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest tests/unit/test_tools/test_genes.py -v 2>&1 | head -10
```

Expected: `ImportError`

- [ ] **Step 3: Create `src/tark_mcp/tools/genes.py`**

```python
from __future__ import annotations
from tark_mcp.client import TarkClient
from tark_mcp.models import Transcript
from tark_mcp.tools.transcripts import _deduplicate


async def get_gene_transcripts(
    gene_identifier: str,
    assembly: str = "GRCh38",
    client: TarkClient | None = None,
) -> list[Transcript]:
    if client is None:
        client = TarkClient()
    params = {
        "identifier_field": gene_identifier,
        "expand": "exons,genes,sequence",
    }
    data = await client.get("transcript/search/", params)
    transcripts = [Transcript.model_validate(r) for r in data]
    if assembly != "both":
        transcripts = [t for t in transcripts if t.assembly == assembly]
    return _deduplicate(transcripts)
```

- [ ] **Step 4: Run tests**

```bash
pytest tests/unit/test_tools/test_genes.py -v
```

Expected: all PASS

- [ ] **Step 5: Commit**

```bash
git add src/tark_mcp/tools/genes.py tests/unit/test_tools/test_genes.py
git commit -m "feat: add get_gene_transcripts tool"
```

---

## Task 7: Sequences Tool

**Files:**
- Create: `src/tark_mcp/tools/sequences.py`
- Create: `tests/unit/test_tools/test_sequences.py`

- [ ] **Step 1: Write failing tests**

Create `tests/unit/test_tools/test_sequences.py`:

```python
import pytest
import httpx
import respx

from tark_mcp.client import TarkClient
from tark_mcp.tools.sequences import (
    get_transcript_sequence, get_transcript_exons, get_protein_for_transcript
)
from tests.conftest import TRANSCRIPT_BRCA2_RAW, TRANSCRIPT_NONCODING_RAW, TRANSLATION_BRCA2_RAW

BASE = "https://tark.ensembl.org/api/"

SINGLE_RESULT = {"count": 1, "next": None, "previous": None,
                 "results": [TRANSCRIPT_BRCA2_RAW]}
NONCODING_RESULT = {"count": 1, "next": None, "previous": None,
                    "results": [TRANSCRIPT_NONCODING_RAW]}


@respx.mock
@pytest.mark.asyncio
async def test_get_transcript_sequence_returns_sequence():
    client = TarkClient()
    respx.get(BASE + "transcript/").mock(return_value=httpx.Response(200, json=SINGLE_RESULT))
    result = await get_transcript_sequence("ENST00000380152", client=client)
    assert result is not None
    assert result["stable_id"] == "ENST00000380152"
    assert result["sequence"] == "ATCGATCGATCGATCGATCGATCGATCGATCG"
    assert result["assembly"] == "GRCh38"


@respx.mock
@pytest.mark.asyncio
async def test_get_transcript_sequence_not_found_returns_none():
    client = TarkClient()
    respx.get(BASE + "transcript/").mock(
        return_value=httpx.Response(200, json={"count": 0, "next": None,
                                               "previous": None, "results": []})
    )
    result = await get_transcript_sequence("ENST00000999999", client=client)
    assert result is None


@respx.mock
@pytest.mark.asyncio
async def test_get_transcript_exons_returns_ordered_exons():
    client = TarkClient()
    respx.get(BASE + "transcript/").mock(return_value=httpx.Response(200, json=SINGLE_RESULT))
    exons = await get_transcript_exons("ENST00000380152", client=client)
    assert len(exons) == 1
    assert exons[0].order == 1
    assert exons[0].stable_id == "ENSE00001184784"


@respx.mock
@pytest.mark.asyncio
async def test_get_transcript_exons_negative_strand_reversed():
    """Exons on negative strand are returned in reverse order."""
    client = TarkClient()
    neg_transcript = {
        **TRANSCRIPT_BRCA2_RAW,
        "loc_strand": -1,
        "exons": [
            {**TRANSCRIPT_BRCA2_RAW["exons"][0], "exon_order": 1},
            {**TRANSCRIPT_BRCA2_RAW["exons"][0], "stable_id": "ENSE00002",
             "exon_order": 2, "loc_start": 32315700, "loc_end": 32315800},
        ],
    }
    respx.get(BASE + "transcript/").mock(
        return_value=httpx.Response(200, json={
            "count": 1, "next": None, "previous": None, "results": [neg_transcript]
        })
    )
    exons = await get_transcript_exons("ENST00000380152", client=client)
    # Negative strand: returned highest order first
    assert exons[0].order == 2
    assert exons[1].order == 1


@respx.mock
@pytest.mark.asyncio
async def test_get_protein_for_transcript_returns_translation():
    client = TarkClient()
    respx.get(BASE + "transcript/").mock(return_value=httpx.Response(200, json=SINGLE_RESULT))
    result = await get_protein_for_transcript("ENST00000380152", client=client)
    assert result is not None
    assert result.stable_id == "ENSP00000369497"


@respx.mock
@pytest.mark.asyncio
async def test_get_protein_for_noncoding_returns_none():
    client = TarkClient()
    respx.get(BASE + "transcript/").mock(return_value=httpx.Response(200, json=NONCODING_RESULT))
    result = await get_protein_for_transcript("ENST00000614536", client=client)
    assert result is None
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest tests/unit/test_tools/test_sequences.py -v 2>&1 | head -10
```

Expected: `ImportError`

- [ ] **Step 3: Create `src/tark_mcp/tools/sequences.py`**

```python
from __future__ import annotations
from tark_mcp.client import TarkClient
from tark_mcp.models import Exon, Translation
from tark_mcp.tools.transcripts import get_transcript


async def get_transcript_sequence(
    stable_id: str,
    assembly: str = "GRCh38",
    client: TarkClient | None = None,
) -> dict | list[dict] | None:
    if client is None:
        client = TarkClient()
    t = await get_transcript(stable_id, assembly=assembly, client=client)
    if t is None:
        return None
    if isinstance(t, list):
        return [
            {"stable_id": x.stable_id, "stable_id_version": x.stable_id_version,
             "assembly": x.assembly, "sequence": x.sequence}
            for x in t
        ]
    return {"stable_id": t.stable_id, "stable_id_version": t.stable_id_version,
            "assembly": t.assembly, "sequence": t.sequence}


async def get_transcript_exons(
    stable_id: str,
    assembly: str = "GRCh38",
    client: TarkClient | None = None,
) -> list[Exon]:
    if client is None:
        client = TarkClient()
    t = await get_transcript(stable_id, assembly=assembly, client=client)
    if t is None:
        return []
    transcripts = t if isinstance(t, list) else [t]
    exons: list[Exon] = []
    for transcript in transcripts:
        ordered = sorted(transcript.exons, key=lambda e: e.order)
        if transcript.loc_strand == -1:
            ordered = list(reversed(ordered))
        exons.extend(ordered)
    return exons


async def get_protein_for_transcript(
    stable_id: str,
    assembly: str = "GRCh38",
    client: TarkClient | None = None,
) -> Translation | list[Translation] | None:
    if client is None:
        client = TarkClient()
    t = await get_transcript(stable_id, assembly=assembly, client=client)
    if t is None:
        return None
    if isinstance(t, list):
        results = [x.translations[0] if x.translations else None for x in t]
        return results
    return t.translations[0] if t.translations else None
```

- [ ] **Step 4: Run tests**

```bash
pytest tests/unit/test_tools/test_sequences.py -v
```

Expected: all PASS

- [ ] **Step 5: Commit**

```bash
git add src/tark_mcp/tools/sequences.py tests/unit/test_tools/test_sequences.py
git commit -m "feat: add get_transcript_sequence, get_transcript_exons, get_protein_for_transcript"
```

---

## Task 8: MANE Tool

**Files:**
- Create: `src/tark_mcp/tools/mane.py`
- Create: `tests/unit/test_tools/test_mane.py`

- [ ] **Step 1: Write failing tests**

Create `tests/unit/test_tools/test_mane.py`:

```python
import pytest
import httpx
import respx

from tark_mcp.client import TarkClient
from tark_mcp.tools.mane import get_mane_transcripts
from tests.conftest import TRANSCRIPT_BRCA2_RAW, MANE_LIST_RESPONSE_RAW

BASE = "https://tark.ensembl.org/api/"


@respx.mock
@pytest.mark.asyncio
async def test_get_mane_transcripts_returns_all():
    client = TarkClient()
    respx.get(BASE + "transcript/manelist/").mock(
        return_value=httpx.Response(200, json=MANE_LIST_RESPONSE_RAW)
    )
    results = await get_mane_transcripts(client=client)
    assert len(results) == 1
    assert results[0].stable_id == "ENST00000380152"


@respx.mock
@pytest.mark.asyncio
async def test_get_mane_transcripts_filters_by_gene_name():
    client = TarkClient()
    two_genes = {
        "count": 2, "next": None, "previous": None,
        "results": [
            TRANSCRIPT_BRCA2_RAW,
            {**TRANSCRIPT_BRCA2_RAW, "stable_id": "ENST00000999999",
             "genes": [{"stable_id": "ENSG00000012048", "stable_id_version": 1,
                        "assembly": "GRCh38", "loc_start": 100, "loc_end": 200,
                        "loc_strand": 1, "loc_region": "1", "name": "BRCA1"}]},
        ],
    }
    respx.get(BASE + "transcript/manelist/").mock(
        return_value=httpx.Response(200, json=two_genes)
    )
    results = await get_mane_transcripts(gene_identifier="BRCA2", client=client)
    assert len(results) == 1
    assert results[0].stable_id == "ENST00000380152"


@respx.mock
@pytest.mark.asyncio
async def test_get_mane_transcripts_filters_by_gene_stable_id():
    client = TarkClient()
    respx.get(BASE + "transcript/manelist/").mock(
        return_value=httpx.Response(200, json=MANE_LIST_RESPONSE_RAW)
    )
    results = await get_mane_transcripts(gene_identifier="ENSG00000139618", client=client)
    assert len(results) == 1
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest tests/unit/test_tools/test_mane.py -v 2>&1 | head -10
```

Expected: `ImportError`

- [ ] **Step 3: Create `src/tark_mcp/tools/mane.py`**

```python
from __future__ import annotations
from tark_mcp.client import TarkClient
from tark_mcp.models import Transcript


async def get_mane_transcripts(
    gene_identifier: str | None = None,
    client: TarkClient | None = None,
) -> list[Transcript]:
    if client is None:
        client = TarkClient()
    data = await client.get("transcript/manelist/")
    transcripts = [Transcript.model_validate(r) for r in data]
    if gene_identifier is None:
        return transcripts
    # Filter client-side by gene name or stable ID
    filtered = []
    for t in transcripts:
        for gene in t.genes:
            if (gene.name and gene.name.upper() == gene_identifier.upper()) or \
               gene.stable_id == gene_identifier:
                filtered.append(t)
                break
    return filtered
```

- [ ] **Step 4: Run tests**

```bash
pytest tests/unit/test_tools/test_mane.py -v
```

Expected: all PASS

- [ ] **Step 5: Commit**

```bash
git add src/tark_mcp/tools/mane.py tests/unit/test_tools/test_mane.py
git commit -m "feat: add get_mane_transcripts tool"
```

---

## Task 9: Diff Tool

**Files:**
- Create: `src/tark_mcp/tools/diff.py`
- Create: `tests/unit/test_tools/test_diff.py`

This is the most complex tool. It:
1. Calls `GET /api/transcript/diff/` to get structural diff and both transcript objects (with cDNA sequences + translation stable IDs)
2. For each transcript that has translations, calls `GET /api/translation/?stable_id=<ENSP>&expand_all=true` to get protein sequences
3. Computes `ExonDiff` list client-side by matching exons by order
4. Populates all `TranscriptDiff` fields including the new protein/CDS sequence fields

The diff endpoint response shape:
```
{
  "count": 1, "next": null, "previous": null,
  "results": { diff-flags object },
  "diff_me_transcript": { full transcript },
  "diff_with_transcript": { full transcript }
}
```

- [ ] **Step 1: Write failing tests**

Create `tests/unit/test_tools/test_diff.py`:

```python
import pytest
import httpx
import respx

from tark_mcp.client import TarkClient
from tark_mcp.tools.diff import diff_transcripts
from tests.conftest import (
    TRANSCRIPT_BRCA2_RAW, TRANSCRIPT_NONCODING_RAW,
    TRANSLATION_BRCA2_RAW, DIFF_RESPONSE_RAW
)

BASE = "https://tark.ensembl.org/api/"

DIFF_BOTH_CODING = {
    "count": 1, "next": None, "previous": None,
    "results": {
        "diff_me_stable_id": "ENST00000380152",
        "diff_with_stable_id": "ENST00000614536",
        "has_seq_changed": True,
        "has_exon_set_changed": True,
        "has_translation_seq_changed": True,
    },
    "diff_me_transcript": {
        **TRANSCRIPT_BRCA2_RAW,
        "exons": [
            {**TRANSCRIPT_BRCA2_RAW["exons"][0], "exon_order": 1},
            {**TRANSCRIPT_BRCA2_RAW["exons"][0], "stable_id": "ENSE00002",
             "exon_order": 2, "loc_start": 32316000, "loc_end": 32316500},
        ],
    },
    "diff_with_transcript": {
        **TRANSCRIPT_BRCA2_RAW,
        "stable_id": "ENST00000614536",
        "stable_id_version": 1,
        "biotype": "protein_coding",
        "sequence": {"sequence": "TTTTGGGGCCCCAAAA", "seq_checksum": "XYZ"},
        "five_prime_utr_seq": "TTTT",
        "three_prime_utr_seq": "AAAA",
        "translations": [
            {"stable_id": "ENSP00000999999", "stable_id_version": 1,
             "assembly": "GRCh38", "loc_start": 100, "loc_end": 200,
             "loc_strand": 1, "loc_region": "13",
             "transcript_stable_id": "ENST00000614536",
             "transcript_stable_id_version": 1}
        ],
        "exons": [
            {**TRANSCRIPT_BRCA2_RAW["exons"][0], "exon_order": 1},
        ],
    },
}

TRANSLATION_CANDIDATE_RAW = {
    "count": 1, "next": None, "previous": None,
    "results": [{
        "stable_id": "ENSP00000999999",
        "stable_id_version": 1,
        "assembly": {"assembly_name": "GRCh38", "assembly_id": 1, "genome": 1, "session": 1},
        "loc_start": 100, "loc_end": 200, "loc_strand": 1, "loc_region": "13",
        "sequence": {"sequence": "MVLSPAD", "seq_checksum": "ZZZ"},
    }]
}

TRANSLATION_REF_RESPONSE = {
    "count": 1, "next": None, "previous": None,
    "results": [{**TRANSLATION_BRCA2_RAW,
                 "assembly": {"assembly_name": "GRCh38", "assembly_id": 1,
                              "genome": 1, "session": 1}}]
}


@respx.mock
@pytest.mark.asyncio
async def test_diff_transcripts_coding_pair_populates_all_sequence_fields():
    """Both transcripts coding: all sequence fields populated, changed flags computed."""
    client = TarkClient()
    respx.get(BASE + "transcript/diff/").mock(
        return_value=httpx.Response(200, json=DIFF_BOTH_CODING)
    )
    respx.get(BASE + "translation/").mock(side_effect=[
        httpx.Response(200, json=TRANSLATION_REF_RESPONSE),
        httpx.Response(200, json=TRANSLATION_CANDIDATE_RAW),
    ])

    results = await diff_transcripts(
        ["ENST00000380152", "ENST00000614536"], client=client
    )
    assert len(results) == 1
    diff = results[0]

    assert diff.reference_stable_id == "ENST00000380152"
    assert diff.candidate_stable_id == "ENST00000614536"
    assert diff.reference_protein_coding is True
    assert diff.candidate_protein_coding is True

    # ref CDS: seq="ATCGATCGATCGATCGATCGATCGATCGATCG"(len=32), 5'UTR="ATCG"(4), 3'UTR="CG"(2)
    # cds_seq = seq[4:30]
    assert diff.ref_cds_sequence == "ATCGATCGATCGATCGATCGATCG"
    # candidate: seq="TTTTGGGGCCCCAAAA"(16), 5'UTR="TTTT"(4), 3'UTR="AAAA"(4)
    # cds_seq = seq[4:12]
    assert diff.candidate_cds_sequence == "GGGGCCCC"
    assert diff.cds_sequence_changed is True

    assert diff.ref_protein_sequence == "MPIGSKERP"
    assert diff.candidate_protein_sequence == "MVLSPAD"
    assert diff.protein_sequence_changed is True


@respx.mock
@pytest.mark.asyncio
async def test_diff_transcripts_noncoding_ref_sets_none_sentinels():
    """Non-coding reference: protein_coding=False, sequence comparison fields=None."""
    client = TarkClient()
    diff_response = {
        **DIFF_RESPONSE_RAW,
        "diff_me_transcript": TRANSCRIPT_NONCODING_RAW,
        "diff_with_transcript": TRANSCRIPT_NONCODING_RAW,
    }
    respx.get(BASE + "transcript/diff/").mock(
        return_value=httpx.Response(200, json=diff_response)
    )
    results = await diff_transcripts(
        ["ENST00000614536", "ENST00000614536"], client=client
    )
    diff = results[0]
    assert diff.reference_protein_coding is False
    assert diff.candidate_protein_coding is False
    assert diff.cds_sequence_changed is None
    assert diff.protein_sequence_changed is None
    assert diff.ref_cds_sequence is None
    assert diff.candidate_cds_sequence is None
    assert diff.ref_protein_sequence is None
    assert diff.candidate_protein_sequence is None


@respx.mock
@pytest.mark.asyncio
async def test_diff_transcripts_mixed_pair_sets_none_sentinels():
    """One coding, one non-coding: sequence comparison fields are None."""
    client = TarkClient()
    respx.get(BASE + "transcript/diff/").mock(
        return_value=httpx.Response(200, json=DIFF_RESPONSE_RAW)
    )
    # Only one translation fetch (for the coding ref)
    respx.get(BASE + "translation/").mock(
        return_value=httpx.Response(200, json=TRANSLATION_REF_RESPONSE)
    )
    results = await diff_transcripts(
        ["ENST00000380152", "ENST00000614536"], client=client
    )
    diff = results[0]
    assert diff.reference_protein_coding is True
    assert diff.candidate_protein_coding is False
    assert diff.cds_sequence_changed is None
    assert diff.protein_sequence_changed is None


@respx.mock
@pytest.mark.asyncio
async def test_diff_transcripts_requires_at_least_two_ids():
    client = TarkClient()
    with pytest.raises(Exception, match="At least 2 stable IDs"):
        await diff_transcripts(["ENST00000380152"], client=client)


@respx.mock
@pytest.mark.asyncio
async def test_diff_transcripts_exon_diffs_computed():
    """ExonDiff list is computed from exon lists of both transcripts."""
    client = TarkClient()
    respx.get(BASE + "transcript/diff/").mock(
        return_value=httpx.Response(200, json=DIFF_RESPONSE_RAW)
    )
    respx.get(BASE + "translation/").mock(
        return_value=httpx.Response(200, json=TRANSLATION_REF_RESPONSE)
    )
    results = await diff_transcripts(
        ["ENST00000380152", "ENST00000614536"], client=client
    )
    diff = results[0]
    assert len(diff.exon_diffs) >= 1
    assert all(d.change in ("added", "removed", "modified", "unchanged") for d in diff.exon_diffs)


@respx.mock
@pytest.mark.asyncio
async def test_diff_transcripts_multiple_pairs():
    """Three stable IDs → two (ref, candidate) pairs, both processed."""
    client = TarkClient()
    respx.get(BASE + "transcript/diff/").mock(
        return_value=httpx.Response(200, json=DIFF_RESPONSE_RAW)
    )
    respx.get(BASE + "translation/").mock(
        return_value=httpx.Response(200, json=TRANSLATION_REF_RESPONSE)
    )
    results = await diff_transcripts(
        ["ENST00000380152", "ENST00000614536", "ENST00000614536"], client=client
    )
    assert len(results) == 2
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest tests/unit/test_tools/test_diff.py -v 2>&1 | head -10
```

Expected: `ImportError`

- [ ] **Step 3: Create `src/tark_mcp/tools/diff.py`**

```python
from __future__ import annotations
import asyncio
from tark_mcp.client import TarkClient
from tark_mcp.models import Transcript, Translation, ExonDiff, TranscriptDiff


def _compute_exon_diffs(ref_exons: list, candidate_exons: list) -> list[ExonDiff]:
    """Compare two exon lists by order. Returns ExonDiff per exon."""
    ref_map = {e.order: e for e in ref_exons}
    cand_map = {e.order: e for e in candidate_exons}
    all_orders = sorted(set(ref_map) | set(cand_map))
    diffs = []
    for order in all_orders:
        ref = ref_map.get(order)
        cand = cand_map.get(order)
        if ref is None:
            change = "added"
            ref_coords = None
            cand_coords = (cand.loc_start, cand.loc_end)
        elif cand is None:
            change = "removed"
            ref_coords = (ref.loc_start, ref.loc_end)
            cand_coords = None
        elif (ref.loc_start, ref.loc_end) != (cand.loc_start, cand.loc_end):
            change = "modified"
            ref_coords = (ref.loc_start, ref.loc_end)
            cand_coords = (cand.loc_start, cand.loc_end)
        else:
            change = "unchanged"
            ref_coords = (ref.loc_start, ref.loc_end)
            cand_coords = (cand.loc_start, cand.loc_end)
        diffs.append(ExonDiff(order=order, change=change,
                               ref_coords=ref_coords, candidate_coords=cand_coords))
    return diffs


def _extract_cds_sequence(t: Transcript) -> str | None:
    """Slice CDS from transcript sequence using cds_start/cds_end offsets."""
    if t.sequence is None or t.cds_start is None or t.cds_end is None:
        return None
    return t.sequence[t.cds_start:t.cds_end]


async def _fetch_protein_sequence(
    translation_stable_id: str,
    assembly: str,
    client: TarkClient,
) -> str | None:
    """Fetch protein sequence from /api/translation/ endpoint."""
    data = await client.get("translation/", {
        "stable_id": translation_stable_id,
        "expand_all": "true",
        "assembly_name": assembly,
    })
    if not data:
        return None
    # Parse first result; sequence is nested {"sequence": "...", "seq_checksum": "..."}
    raw = data[0]
    seq = raw.get("sequence")
    if isinstance(seq, dict):
        return seq.get("sequence")
    return seq


async def _build_diff(
    ref: Transcript,
    candidate: Transcript,
    client: TarkClient,
) -> TranscriptDiff:
    ref_cds_seq = _extract_cds_sequence(ref)
    cand_cds_seq = _extract_cds_sequence(candidate)

    ref_protein_seq: str | None = None
    cand_protein_seq: str | None = None

    async def _get_prot(t: Transcript) -> str | None:
        if not t.translations:
            return None
        return await _fetch_protein_sequence(t.translations[0].stable_id, t.assembly, client)

    ref_protein_seq, cand_protein_seq = await asyncio.gather(
        _get_prot(ref), _get_prot(candidate)
    )

    ref_coding = ref_cds_seq is not None and ref.translations != []
    cand_coding = cand_cds_seq is not None and candidate.translations != []

    cds_changed: bool | None = None
    if ref_cds_seq is not None and cand_cds_seq is not None:
        cds_changed = ref_cds_seq != cand_cds_seq

    protein_changed: bool | None = None
    if ref_protein_seq is not None and cand_protein_seq is not None:
        protein_changed = ref_protein_seq != cand_protein_seq

    exon_diffs = _compute_exon_diffs(ref.exons, candidate.exons)

    return TranscriptDiff(
        reference_stable_id=ref.stable_id,
        candidate_stable_id=candidate.stable_id,
        reference_assembly=ref.assembly,
        candidate_assembly=candidate.assembly,
        biotype_changed=ref.biotype != candidate.biotype,
        cds_changed=(ref.cds_start, ref.cds_end) != (candidate.cds_start, candidate.cds_end),
        exon_count_changed=len(ref.exons) != len(candidate.exons),
        sequence_changed=ref.sequence != candidate.sequence,
        exon_diffs=exon_diffs,
        reference_protein_coding=ref_coding,
        candidate_protein_coding=cand_coding,
        cds_sequence_changed=cds_changed,
        ref_cds_sequence=ref_cds_seq,
        candidate_cds_sequence=cand_cds_seq,
        protein_sequence_changed=protein_changed,
        ref_protein_sequence=ref_protein_seq,
        candidate_protein_sequence=cand_protein_seq,
    )


async def _fetch_diff_pair(
    ref_stable_id: str,
    candidate_stable_id: str,
    ref_assembly: str,
    candidate_assembly: str,
    client: TarkClient,
) -> TranscriptDiff:
    raw = await client.get_raw(
        "transcript/diff/",
        params={
            "diff_me_stable_id": ref_stable_id,
            "diff_with_stable_id": candidate_stable_id,
        },
    )
    ref_data = raw.get("diff_me_transcript", {})
    cand_data = raw.get("diff_with_transcript", {})

    # Ensure assembly is set correctly from context if nested object is missing
    if isinstance(ref_data.get("assembly"), dict):
        pass  # model_validate will flatten it
    else:
        ref_data = {**ref_data, "assembly": {"assembly_name": ref_assembly}}

    if isinstance(cand_data.get("assembly"), dict):
        pass
    else:
        cand_data = {**cand_data, "assembly": {"assembly_name": candidate_assembly}}

    ref = Transcript.model_validate(ref_data)
    candidate = Transcript.model_validate(cand_data)
    return await _build_diff(ref, candidate, client)


async def diff_transcripts(
    stable_ids: list[str],
    assemblies: list[str] | None = None,
    client: TarkClient | None = None,
) -> list[TranscriptDiff]:
    if len(stable_ids) < 2:
        raise ValueError("At least 2 stable IDs required for diff")
    if client is None:
        client = TarkClient()

    # Resolve per-entry assembly: fill missing with "GRCh38"
    resolved_assemblies = list(assemblies or [])
    while len(resolved_assemblies) < len(stable_ids):
        resolved_assemblies.append("GRCh38")

    ref_id = stable_ids[0]
    ref_assembly = resolved_assemblies[0]

    pairs = [
        (ref_id, stable_ids[i], ref_assembly, resolved_assemblies[i])
        for i in range(1, len(stable_ids))
    ]

    results = await asyncio.gather(*[
        _fetch_diff_pair(ref, cand, ra, ca, client)
        for ref, cand, ra, ca in pairs
    ])
    return list(results)
```

- [ ] **Step 4: Run tests**

```bash
pytest tests/unit/test_tools/test_diff.py -v
```

Expected: all PASS. If any fail, debug and fix before continuing.

- [ ] **Step 5: Commit**

```bash
git add src/tark_mcp/tools/diff.py tests/unit/test_tools/test_diff.py
git commit -m "feat: add diff_transcripts tool with CDS and protein sequence comparison"
```

---

## Task 10: MCP Server

**Files:**
- Create: `src/tark_mcp/server.py`

- [ ] **Step 1: Create `src/tark_mcp/server.py`**

```python
from __future__ import annotations
from mcp.server.fastmcp import FastMCP

from tark_mcp.client import TarkClient
from tark_mcp.tools.releases import get_releases
from tark_mcp.tools.transcripts import get_transcript, search_transcripts_by_region
from tark_mcp.tools.genes import get_gene_transcripts
from tark_mcp.tools.sequences import (
    get_transcript_sequence, get_transcript_exons, get_protein_for_transcript
)
from tark_mcp.tools.mane import get_mane_transcripts
from tark_mcp.tools.diff import diff_transcripts

mcp = FastMCP("tark")
_client = TarkClient()


@mcp.tool()
async def tark_get_releases() -> list[dict]:
    """List all available TARK releases with metadata (short name, date, assembly, source)."""
    releases = await get_releases(_client)
    return [r.model_dump() for r in releases]


@mcp.tool()
async def tark_get_transcript(stable_id: str, assembly: str = "GRCh38") -> dict | list[dict] | None:
    """Retrieve a transcript by Ensembl stable ID with full exon structure, CDS boundaries, and genes.

    Args:
        stable_id: Ensembl transcript stable ID, e.g. 'ENST00000380152' or 'ENST00000380152.7'
        assembly: Genome build — 'GRCh37', 'GRCh38' (default), or 'both'
    """
    result = await get_transcript(stable_id, assembly=assembly, client=_client)
    if result is None:
        return None
    if isinstance(result, list):
        return [t.model_dump() for t in result]
    return result.model_dump()


@mcp.tool()
async def tark_search_transcripts_by_region(
    region: str, start: int, end: int, assembly: str = "GRCh38"
) -> list[dict]:
    """Find all transcripts overlapping a genomic region (0-based half-open coordinates).

    Args:
        region: Chromosome, e.g. '13' or 'chr13'
        start: 0-based start position (inclusive)
        end: 0-based end position (exclusive)
        assembly: 'GRCh37', 'GRCh38' (default), or 'both'
    """
    results = await search_transcripts_by_region(region, start, end, assembly=assembly, client=_client)
    return [t.model_dump() for t in results]


@mcp.tool()
async def tark_get_gene_transcripts(gene_identifier: str, assembly: str = "GRCh38") -> list[dict]:
    """Get all transcripts for a gene symbol or Ensembl gene ID.

    Args:
        gene_identifier: Gene symbol (e.g. 'BRCA2') or Ensembl gene ID (e.g. 'ENSG00000139618')
        assembly: 'GRCh37', 'GRCh38' (default), or 'both'
    """
    results = await get_gene_transcripts(gene_identifier, assembly=assembly, client=_client)
    return [t.model_dump() for t in results]


@mcp.tool()
async def tark_get_transcript_sequence(
    stable_id: str, assembly: str = "GRCh38"
) -> dict | list[dict] | None:
    """Fetch the cDNA sequence for a transcript.

    Args:
        stable_id: Ensembl transcript stable ID
        assembly: 'GRCh37', 'GRCh38' (default), or 'both'
    """
    return await get_transcript_sequence(stable_id, assembly=assembly, client=_client)


@mcp.tool()
async def tark_get_transcript_exons(
    stable_id: str, assembly: str = "GRCh38"
) -> list[dict]:
    """Return the ordered exon list for a transcript with 0-based genomic coordinates.

    Args:
        stable_id: Ensembl transcript stable ID
        assembly: 'GRCh37', 'GRCh38' (default), or 'both'
    """
    exons = await get_transcript_exons(stable_id, assembly=assembly, client=_client)
    return [e.model_dump() for e in exons]


@mcp.tool()
async def tark_get_protein_for_transcript(
    stable_id: str, assembly: str = "GRCh38"
) -> dict | list[dict | None] | None:
    """Return the protein (translation) stable ID and version for a transcript.

    Args:
        stable_id: Ensembl transcript stable ID
        assembly: 'GRCh37', 'GRCh38' (default), or 'both'
    """
    result = await get_protein_for_transcript(stable_id, assembly=assembly, client=_client)
    if result is None:
        return None
    if isinstance(result, list):
        return [t.model_dump() if t else None for t in result]
    return result.model_dump()


@mcp.tool()
async def tark_get_mane_transcripts(gene_identifier: str | None = None) -> list[dict]:
    """Return MANE Select and MANE Plus Clinical transcripts, optionally filtered by gene.

    Args:
        gene_identifier: Optional gene symbol or Ensembl gene ID to filter results
    """
    results = await get_mane_transcripts(gene_identifier=gene_identifier, client=_client)
    return [t.model_dump() for t in results]


@mcp.tool()
async def tark_diff_transcripts(
    stable_ids: list[str],
    assemblies: list[str] | None = None,
) -> list[dict]:
    """Compare transcripts against a reference. First ID is reference; all subsequent are candidates.

    Returns structural diff (exon-level, CDS boundaries, biotype) plus CDS nucleotide sequence
    comparison and protein (amino acid) sequence comparison. Non-coding transcripts produce
    None for sequence comparison fields.

    Args:
        stable_ids: List of ≥2 Ensembl transcript stable IDs; first is the reference
        assemblies: Optional per-entry assembly override list; defaults to 'GRCh38' for missing entries
    """
    results = await diff_transcripts(stable_ids, assemblies=assemblies, client=_client)
    return [d.model_dump() for d in results]


def main() -> None:
    mcp.run()


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Verify server imports cleanly**

```bash
python -c "from tark_mcp.server import mcp; print('OK')"
```

Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add src/tark_mcp/server.py
git commit -m "feat: add MCP server with all tool registrations"
```

---

## Task 11: Integration Tests

**Files:**
- Create: `tests/integration/test_live_api.py`

- [ ] **Step 1: Create `tests/integration/test_live_api.py`**

```python
"""Integration tests against the live TARK API.

Run with:  TARK_INTEGRATION=1 pytest tests/integration/ -v
Skipped by default.
"""
import os
import pytest

pytestmark = pytest.mark.skipif(
    os.environ.get("TARK_INTEGRATION") != "1",
    reason="Set TARK_INTEGRATION=1 to run live API tests"
)

from tark_mcp.client import TarkClient
from tark_mcp.tools.releases import get_releases
from tark_mcp.tools.transcripts import get_transcript, search_transcripts_by_region
from tark_mcp.tools.genes import get_gene_transcripts
from tark_mcp.tools.sequences import (
    get_transcript_sequence, get_transcript_exons, get_protein_for_transcript
)
from tark_mcp.tools.mane import get_mane_transcripts
from tark_mcp.tools.diff import diff_transcripts


@pytest.fixture(scope="module")
def client():
    return TarkClient()


@pytest.mark.asyncio
async def test_get_releases(client):
    releases = await get_releases(client)
    assert len(releases) > 0
    assert all(hasattr(r, "shortname") for r in releases)


@pytest.mark.asyncio
async def test_get_transcript_brca2(client):
    t = await get_transcript("ENST00000380152", client=client)
    assert t is not None
    assert t.stable_id == "ENST00000380152"
    assert t.assembly == "GRCh38"
    assert len(t.exons) > 0


@pytest.mark.asyncio
async def test_get_transcript_both_assemblies(client):
    result = await get_transcript("ENST00000380152", assembly="both", client=client)
    assert isinstance(result, list)
    assemblies = {t.assembly for t in result}
    # Should have at least GRCh38
    assert "GRCh38" in assemblies


@pytest.mark.asyncio
async def test_search_transcripts_by_region(client):
    # BRCA2 locus on chr13
    results = await search_transcripts_by_region("13", 32315474, 32400266, client=client)
    assert len(results) > 0
    stable_ids = {t.stable_id for t in results}
    assert "ENST00000380152" in stable_ids


@pytest.mark.asyncio
async def test_get_gene_transcripts_by_symbol(client):
    results = await get_gene_transcripts("BRCA2", client=client)
    assert len(results) > 0
    assert all(t.assembly == "GRCh38" for t in results)


@pytest.mark.asyncio
async def test_get_transcript_sequence(client):
    result = await get_transcript_sequence("ENST00000380152", client=client)
    assert result is not None
    assert result["sequence"] is not None
    assert len(result["sequence"]) > 100


@pytest.mark.asyncio
async def test_get_transcript_exons(client):
    exons = await get_transcript_exons("ENST00000380152", client=client)
    assert len(exons) > 0
    # Exons should be in ascending order (positive strand)
    orders = [e.order for e in exons]
    assert orders == sorted(orders)


@pytest.mark.asyncio
async def test_get_protein_for_transcript(client):
    result = await get_protein_for_transcript("ENST00000380152", client=client)
    assert result is not None
    assert result.stable_id.startswith("ENSP")


@pytest.mark.asyncio
async def test_get_mane_transcripts(client):
    results = await get_mane_transcripts(client=client)
    assert len(results) > 0


@pytest.mark.asyncio
async def test_get_mane_transcripts_filtered(client):
    results = await get_mane_transcripts(gene_identifier="BRCA2", client=client)
    assert len(results) > 0
    for t in results:
        gene_names = {g.name for g in t.genes}
        assert "BRCA2" in gene_names


@pytest.mark.asyncio
async def test_diff_transcripts_two_coding(client):
    # BRCA2 stable versions
    results = await diff_transcripts(
        ["ENST00000380152", "ENST00000544455"], client=client
    )
    assert len(results) == 1
    diff = results[0]
    assert diff.reference_stable_id == "ENST00000380152"
    assert diff.candidate_stable_id == "ENST00000544455"
    # Both are in protein_coding transcripts
    assert diff.reference_protein_coding is True


@pytest.mark.asyncio
async def test_diff_transcripts_protein_sequences_populated(client):
    results = await diff_transcripts(
        ["ENST00000380152", "ENST00000544455"], client=client
    )
    diff = results[0]
    if diff.reference_protein_coding and diff.candidate_protein_coding:
        assert diff.ref_protein_sequence is not None
        assert diff.candidate_protein_sequence is not None
        assert diff.protein_sequence_changed is not None
```

- [ ] **Step 2: Verify integration tests are skipped in normal runs**

```bash
pytest tests/integration/ -v
```

Expected: all tests SKIPPED with message "Set TARK_INTEGRATION=1 to run live API tests"

- [ ] **Step 3: Commit**

```bash
git add tests/integration/test_live_api.py
git commit -m "test: add integration tests against live TARK API"
```

---

## Task 12: Full Test Suite Pass

- [ ] **Step 1: Run the complete test suite**

```bash
pytest tests/unit/ -v
```

Expected: all unit tests PASS

- [ ] **Step 2: Fix any failures**

If tests fail, read the error output carefully and fix the implementation. Common issues:
- Model validator reading stale `data["sequence"]` before/after the flatten — check ordering in `Transcript._normalize`
- `_extract_cds_sequence` returning wrong slice — print `t.cds_start`, `t.cds_end`, `len(t.sequence)` to verify
- `_fetch_protein_sequence` not finding sequence in API response — check the `translation/` endpoint response parsing

- [ ] **Step 3: Final commit**

```bash
pytest tests/unit/ -v  # must be all green
git add -A
git commit -m "chore: all unit tests passing"
```

---

## Self-Review Checklist

- Spec coverage: all 9 MCP tools implemented ✓, protein/CDS sequence comparison in `diff_transcripts` ✓, `Translation.sequence` field ✓, `TranscriptDiff` 8 new fields ✓
- No TBDs or placeholders ✓
- Types consistent: `Transcript.exons` → `list[Exon]`, `Exon.order` (from `exon_order`), `Translation.sequence: str | None`, `TranscriptDiff.ref_cds_sequence: str | None` — all consistent across tasks ✓
- `_deduplicate` imported from `transcripts.py` in `genes.py` ✓
- `get_transcript` used by `sequences.py` tools (no duplication of fetch logic) ✓
