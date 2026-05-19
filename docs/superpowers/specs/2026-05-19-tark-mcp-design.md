# TARK MCP Server — Design Spec

**Date:** 2026-05-19  
**Status:** Approved

---

## Problem Statement

The Ensembl TARK (Transcript Archive) REST API at `https://tark.ensembl.org/api/` exposes a rich genomics dataset covering transcripts, genes, exons, sequences, and releases across GRCh37 and GRCh38 genome builds. No MCP server wrapping this API exists. This document specifies a Python MCP server that exposes TARK data through a semantic, LLM-friendly tool interface, abstracting away API pagination, duplicate records, coordinate systems, and checksum-based lookups.

---

## Approach

**Semantic Genomics Layer.** Tools are organized around genomics intent (e.g., "get a transcript", "find transcripts in a region") rather than raw API endpoints. All TARK API complexity — pagination, deduplication across releases, 1-based-to-0-based coordinate conversion, CDS boundary calculation from UTR sequences — is handled transparently by the server.

---

## Architecture

```
tark-mcp/
├── pyproject.toml
├── README.md
├── src/
│   └── tark_mcp/
│       ├── __init__.py
│       ├── server.py        # MCP server entry point (stdio transport)
│       ├── client.py        # TARK API HTTP client
│       ├── models.py        # Pydantic data models + coordinate normalization
│       └── tools/
│           ├── __init__.py
│           ├── transcripts.py
│           ├── genes.py
│           ├── sequences.py
│           ├── diff.py
│           ├── mane.py
│           └── releases.py
└── tests/
    ├── conftest.py
    ├── unit/
    │   ├── test_client.py
    │   ├── test_models.py
    │   └── test_tools/
    │       └── (one file per tool module)
    └── integration/
        └── test_live_api.py   # skipped unless TARK_INTEGRATION=1
```

### Layer Responsibilities

**`TarkClient` (`client.py`)**  
Low-level HTTP client built on `httpx`. Responsibilities:
- Base URL: `https://tark.ensembl.org/api/`
- Auto-pagination: follows `next` links, aggregates all pages
- HTTPS enforcement: rewrites any `http://` URLs to `https://`
- TTL caching: `cachetools.TTLCache`, default TTL 3600 s (configurable via `TARK_CACHE_TTL` env var)
- Retry: up to 3 attempts with exponential backoff (`tenacity`) on 5xx and network errors

**`models.py`**  
Pydantic v2 models for all TARK resource types. Coordinate normalization occurs at parse time:
- API 1-based inclusive coordinates → 0-based half-open `[start, end)` intervals
- CDS start/end calculated from `five_prime_utr_seq` / `three_prime_utr_seq` lengths
- Deduplication: for identical `(assembly, stable_id, stable_id_version)` tuples, the record with the most recent `transcript_release_set.release_date` is kept

**`server.py`**  
MCP server using the official `mcp` Python SDK (stdio transport). Registers all tools and wires them to the tool modules.

---

## Data Models

```python
class Exon(BaseModel):
    stable_id: str
    stable_id_version: int
    assembly: str             # which genome build this exon belongs to
    order: int                # 1-based exon number
    loc_region: str
    loc_start: int            # 0-based
    loc_end: int              # exclusive
    loc_strand: int           # 1 or -1

class Gene(BaseModel):
    stable_id: str
    stable_id_version: int
    name: str | None
    loc_region: str
    loc_start: int            # 0-based
    loc_end: int              # exclusive
    assembly: str

class Translation(BaseModel):
    stable_id: str
    stable_id_version: int

class Transcript(BaseModel):
    stable_id: str
    stable_id_version: int
    assembly: str
    biotype: str
    loc_region: str
    loc_start: int            # 0-based
    loc_end: int              # exclusive
    loc_strand: int
    cds_start: int | None     # 0-based offset into transcript sequence
    cds_end: int | None
    exons: list[Exon]
    genes: list[Gene]
    translations: list[Translation]
    sequence: str | None

class Release(BaseModel):
    shortname: str
    description: str | None
    release_date: str
    assembly: str
    source: str

class ExonDiff(BaseModel):
    order: int
    change: str           # "added", "removed", "modified", "unchanged"
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
```

---

## MCP Tools

All tools accept `assembly` defaulting to `"GRCh38"`. Accepted values: `"GRCh37"`, `"GRCh38"`, or `"both"`. When `"both"` is specified, the tool fans out requests to both assemblies in parallel and returns combined results, with each record's `assembly` field indicating its origin. For tools that return a single object (e.g., `get_transcript`), specifying `"both"` returns a list with up to two entries (one per assembly).

### `get_transcript`
Retrieve a single transcript by stable ID, with full exon structure, CDS boundaries, UTR sequences, associated genes, and protein accession.

- **Parameters:** `stable_id: str`, `assembly: str = "GRCh38"` (`"GRCh37"` | `"GRCh38"` | `"both"`)
- **Returns:** `Transcript | list[Transcript]` — single object when one assembly specified, list of up to 2 when `"both"`
- **API calls:** `GET /api/transcript/?stable_id=...&expand_all=true`
- **Notes:** Strips version suffix from stable_id if provided (e.g., `ENST00000614536.1` → id=`ENST00000614536`, version=`1`). When `assembly="both"`, requests are issued in parallel.

### `search_transcripts_by_region`
Find all transcripts whose genomic footprint overlaps a given region.

- **Parameters:** `region: str` (chromosome, e.g. `"1"`, `"chrX"`), `start: int` (0-based), `end: int` (exclusive), `assembly: str = "GRCh38"` (`"GRCh37"` | `"GRCh38"` | `"both"`)
- **Returns:** `list[Transcript]`
- **API calls:** `GET /api/transcript/?assembly_name=...&loc_region=...&loc_start=...&loc_end=...&expand=transcript_release_set`
- **Notes:** Converts 0-based input to 1-based for the API. Strips `chr` prefix from region if present.

### `get_gene_transcripts`
Retrieve all transcripts associated with a gene symbol or Ensembl gene ID.

- **Parameters:** `gene_identifier: str` (e.g. `"BRCA2"` or `"ENSG00000139618"`), `assembly: str = "GRCh38"` (`"GRCh37"` | `"GRCh38"` | `"both"`)
- **Returns:** `list[Transcript]`
- **API calls:** `GET /api/transcript/search/?identifier_field=...&expand=exons,genes,sequence`
- **Notes:** The search endpoint does not support server-side assembly filtering; results are filtered client-side by the `assembly` field on each returned transcript.

### `get_transcript_sequence`
Fetch the cDNA sequence for a transcript.

- **Parameters:** `stable_id: str`, `assembly: str = "GRCh38"` (`"GRCh37"` | `"GRCh38"` | `"both"`)
- **Returns:** `{ "stable_id": str, "assembly": str, "sequence": str }` or `list[...]` when `"both"`
- **API calls:** `GET /api/transcript/?stable_id=...&expand_all=true` (sequence is nested in the transcript response)

### `get_transcript_exons`
Return the ordered exon list for a transcript with 0-based genomic coordinates.

- **Parameters:** `stable_id: str`, `assembly: str = "GRCh38"` (`"GRCh37"` | `"GRCh38"` | `"both"`)
- **Returns:** `list[Exon]` (combined from both assemblies when `"both"`, each Exon carries its `assembly` field)
- **API calls:** `GET /api/transcript/?stable_id=...&expand_all=true`
- **Notes:** Exons are returned in transcript order (reversed for negative-strand transcripts).

### `get_protein_for_transcript`
Return the protein (translation) stable ID and version for a transcript.

- **Parameters:** `stable_id: str`, `assembly: str = "GRCh38"` (`"GRCh37"` | `"GRCh38"` | `"both"`)
- **Returns:** `Translation | list[Translation] | None`
- **API calls:** `GET /api/transcript/?stable_id=...&expand_all=true`

### `diff_transcripts`
Compare two or more transcripts against a reference transcript. The first stable ID is the reference; all subsequent IDs are compared to it. Returns a structured diff per pair including changes in exon count, exon coordinates, CDS start/end, biotype, and sequence.

- **Parameters:** `stable_ids: list[str]` (minimum 2), `assemblies: list[str] | None` (per-entry assembly override; if shorter than `stable_ids`, defaults to `"GRCh38"` for missing entries)
- **Returns:** `list[TranscriptDiff]` with one entry per (reference, candidate) pair
- **API calls:** `GET /api/transcript/diff/?diff_me_stable_id=...&diff_with_stable_id=...` (one call per pair, fanned out)

### `get_mane_transcripts`
Return the MANE Select and MANE Plus Clinical transcript list, optionally filtered by gene identifier.

- **Parameters:** `gene_identifier: str | None`
- **Returns:** `list[Transcript]`
- **API calls:** `GET /api/transcript/manelist/`
- **Notes:** If `gene_identifier` is supplied, the results are filtered client-side by gene name or stable ID.

### `get_releases`
List all available TARK releases with metadata (short name, date, assembly, source).

- **Parameters:** *(none)*
- **Returns:** `list[Release]`
- **API calls:** `GET /api/release/nopagination/`

---

## Error Handling

| Condition | Behaviour |
|-----------|-----------|
| HTTP 404 / empty results | Return empty list or `None`; no exception |
| HTTP 4xx (other) | Raise `McpError` with descriptive message including the parameters used |
| HTTP 5xx / network error | Retry up to 3× with exponential backoff; raise `McpError` after exhaustion |
| Invalid `assembly` value | Immediate `McpError` listing accepted values (`"GRCh37"`, `"GRCh38"`, `"both"`) before any HTTP call |
| `stable_ids` list < 2 in `diff_transcripts` | Immediate `McpError`: "At least 2 stable IDs required for diff" |

---

## Caching

- `TarkClient` uses `cachetools.TTLCache` keyed on the full request URL (including query string)
- Default TTL: 3600 s (1 hour), configurable via `TARK_CACHE_TTL` environment variable (seconds)
- Cache size: 512 entries maximum
- Cache is per-process (not shared across server restarts)

---

## Configuration

All configuration via environment variables:

| Variable | Default | Description |
|----------|---------|-------------|
| `TARK_BASE_URL` | `https://tark.ensembl.org/api/` | Override API base URL |
| `TARK_CACHE_TTL` | `3600` | Cache TTL in seconds |
| `TARK_REQUEST_TIMEOUT` | `30` | HTTP request timeout in seconds |
| `TARK_MAX_RETRIES` | `3` | Max retry attempts on transient errors |

---

## Testing Strategy

**Unit tests** (`tests/unit/`) — `pytest` + `respx` (mock `httpx`):
- `test_client.py`: pagination aggregation, HTTPS rewriting, deduplication logic, TTL cache hit/miss, retry behaviour
- `test_models.py`: coordinate conversion (1-based → 0-based), CDS boundary calculation, version-suffix stripping
- `test_tools/`: one file per tool module testing parameter validation, response shaping, error propagation

**Integration tests** (`tests/integration/test_live_api.py`) — one test per tool hitting the live TARK API. Marked with `pytest.mark.integration`; skipped unless `TARK_INTEGRATION=1` is set. Validates that real API responses deserialize correctly.

---

## Packaging & Installation

```toml
# pyproject.toml (excerpt)
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
```

### Claude Desktop config snippet
```json
{
  "mcpServers": {
    "tark": {
      "command": "tark-mcp",
      "env": {
        "TARK_CACHE_TTL": "3600"
      }
    }
  }
}
```

### VS Code MCP config snippet
```json
{
  "mcp": {
    "servers": {
      "tark": {
        "type": "stdio",
        "command": "tark-mcp"
      }
    }
  }
}
```
