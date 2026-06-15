# Diff Versioned and RefSeq Transcripts Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `tark_diff_transcripts` correctly handle versioned IDs (e.g. `NM_001128425.2`, `ENST00000710952.1`) and RefSeq IDs (`NM_...`) by resolving each transcript individually instead of relying on TARK's `/transcript/diff/` endpoint.

**Architecture:** Replace `_fetch_diff_pair` in `diff.py` to fetch each transcript individually via `/transcript/?stable_id=...&stable_id_version=...` (stripping the version suffix before calling the API), then compute the diff locally with the existing `_build_diff` function. This bypasses the TARK diff endpoint, which does not accept versioned IDs or RefSeq stable IDs.

**Tech Stack:** Python 3.11+, httpx, respx (test mocking), pytest-asyncio, Pydantic v2

---

### Task 1: Add `_resolve_transcript` helper and update `_fetch_diff_pair`

**Files:**
- Modify: `src/tark_mcp/tools/diff.py`

**Context on current code:**
- `_fetch_diff_pair` currently calls `client.get_raw("transcript/diff/", ...)` which only works for bare Ensembl IDs.
- `_build_diff(ref, candidate, client)` takes two `Transcript` objects and computes the full diff — this is unchanged.
- `_strip_version` already exists in `transcripts.py` and does the same thing we need here. We'll duplicate the logic inline rather than import to avoid coupling.

- [ ] **Step 1: Replace `_fetch_diff_pair` in `src/tark_mcp/tools/diff.py`**

Replace the current `_fetch_diff_pair` function entirely with the following. Keep everything else in the file unchanged:

```python
async def _resolve_transcript(
    stable_id: str,
    assembly: str,
    client: TarkClient,
) -> Transcript:
    """Fetch a single Transcript by stable ID, handling version suffixes and RefSeq IDs."""
    sid, version = stable_id, None
    if "." in stable_id:
        parts = stable_id.rsplit(".", 1)
        try:
            version = int(parts[1])
            sid = parts[0]
        except ValueError:
            pass  # not a version suffix (e.g. NM_001128425.something non-numeric)

    params: dict = {"stable_id": sid, "expand_all": "true", "assembly_name": assembly}
    if version is not None:
        params["stable_id_version"] = version

    data = await client.get("transcript/", params=params)
    if not data:
        raise ValueError(
            f"Transcript not found: {stable_id} (assembly={assembly})"
        )

    transcripts = [Transcript.model_validate(r) for r in data]

    if version is not None:
        matching = [t for t in transcripts if t.stable_id_version == version]
        if matching:
            return matching[0]

    return max(transcripts, key=lambda t: t.latest_release_date or "")


async def _fetch_diff_pair(
    ref_stable_id: str,
    candidate_stable_id: str,
    ref_assembly: str,
    candidate_assembly: str,
    client: TarkClient,
) -> TranscriptDiff:
    ref, candidate = await asyncio.gather(
        _resolve_transcript(ref_stable_id, ref_assembly, client),
        _resolve_transcript(candidate_stable_id, candidate_assembly, client),
    )
    return await _build_diff(ref, candidate, client)
```

Note: remove the old `import` of anything from the previous `_fetch_diff_pair` if it becomes unused (the `get_raw` import line is not a direct import — `client.get_raw` is a method — so no import to remove).

- [ ] **Step 2: Verify the file compiles**

```bash
cd /Users/adymun/tark-mcp && python -c "from tark_mcp.tools.diff import diff_transcripts; print('OK')"
```

Expected: `OK`

---

### Task 2: Update existing diff tests to match new fetch pattern

**Files:**
- Modify: `tests/unit/test_tools/test_diff.py`

**Context:** The existing tests mock `GET /transcript/diff/`. After Task 1, that endpoint is no longer called — instead, two calls to `GET /transcript/` are made (one per transcript). The mocks must be updated accordingly.

The existing fixture data in `conftest.py`:
- `TRANSCRIPT_BRCA2_RAW` — Ensembl coding transcript, `stable_id="ENST00000380152"`, `stable_id_version=7`
- `TRANSCRIPT_NONCODING_RAW` — non-coding transcript, `stable_id="ENST00000614536"`, `stable_id_version=1`
- `TRANSLATION_BRCA2_RAW` — protein translation for BRCA2
- `DIFF_RESPONSE_RAW` — old TARK diff response (no longer used)

The TARK `/transcript/` endpoint returns paginated results: `{"count": N, "next": null, "results": [...]}`.

- [ ] **Step 3: Rewrite all existing tests in `test_diff.py` to use `/transcript/` mocks**

Replace the entire contents of `tests/unit/test_tools/test_diff.py` with:

```python
import pytest
import httpx
import respx

from tark_mcp.client import TarkClient
from tark_mcp.tools.diff import diff_transcripts
from tests.conftest import (
    TRANSCRIPT_BRCA2_RAW, TRANSCRIPT_NONCODING_RAW,
    TRANSLATION_BRCA2_RAW,
)

BASE = "https://tark.ensembl.org/api/"

# Paginated wrapper used by client.get()
def _page(items):
    return {"count": len(items), "next": None, "previous": None, "results": items}

TRANSCRIPT_BRCA2_PAGE = _page([TRANSCRIPT_BRCA2_RAW])
TRANSCRIPT_NONCODING_PAGE = _page([TRANSCRIPT_NONCODING_RAW])

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

TRANSCRIPT_CODING_CANDIDATE_RAW = {
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
}
TRANSCRIPT_CODING_CANDIDATE_PAGE = _page([TRANSCRIPT_CODING_CANDIDATE_RAW])


@respx.mock
@pytest.mark.asyncio
async def test_diff_transcripts_coding_pair_populates_all_sequence_fields():
    """Both transcripts coding: all sequence fields populated, changed flags computed."""
    client = TarkClient()
    respx.get(BASE + "transcript/").mock(side_effect=[
        httpx.Response(200, json=TRANSCRIPT_BRCA2_PAGE),
        httpx.Response(200, json=TRANSCRIPT_CODING_CANDIDATE_PAGE),
    ])
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

    # ref CDS: seq="ATCGATCGATCGATCGATCGATCGATCGATCG"(32), 5'UTR="ATCG"(4), 3'UTR="CG"(2)
    assert diff.ref_cds_sequence == "ATCGATCGATCGATCGATCGATCGAT"
    # candidate: seq="TTTTGGGGCCCCAAAA"(16), 5'UTR="TTTT"(4), 3'UTR="AAAA"(4)
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
    respx.get(BASE + "transcript/").mock(side_effect=[
        httpx.Response(200, json=TRANSCRIPT_NONCODING_PAGE),
        httpx.Response(200, json=TRANSCRIPT_NONCODING_PAGE),
    ])
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
    respx.get(BASE + "transcript/").mock(side_effect=[
        httpx.Response(200, json=TRANSCRIPT_BRCA2_PAGE),
        httpx.Response(200, json=TRANSCRIPT_NONCODING_PAGE),
    ])
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
    respx.get(BASE + "transcript/").mock(side_effect=[
        httpx.Response(200, json=TRANSCRIPT_BRCA2_PAGE),
        httpx.Response(200, json=TRANSCRIPT_NONCODING_PAGE),
    ])
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
    # ref fetched once, then two candidates fetched once each
    # asyncio.gather fetches ref+cand1 in parallel, then ref+cand2 in parallel (4 total)
    respx.get(BASE + "transcript/").mock(side_effect=[
        httpx.Response(200, json=TRANSCRIPT_BRCA2_PAGE),
        httpx.Response(200, json=TRANSCRIPT_NONCODING_PAGE),
        httpx.Response(200, json=TRANSCRIPT_BRCA2_PAGE),
        httpx.Response(200, json=TRANSCRIPT_NONCODING_PAGE),
    ])
    respx.get(BASE + "translation/").mock(
        return_value=httpx.Response(200, json=TRANSLATION_REF_RESPONSE)
    )
    results = await diff_transcripts(
        ["ENST00000380152", "ENST00000614536", "ENST00000614536"], client=client
    )
    assert len(results) == 2


@respx.mock
@pytest.mark.asyncio
async def test_diff_transcripts_strips_version_suffix():
    """Versioned IDs like 'ENST00000380152.7' are stripped and stable_id_version passed separately."""
    client = TarkClient()
    # Both versioned IDs should result in /transcript/ calls (version stripped)
    respx.get(BASE + "transcript/").mock(side_effect=[
        httpx.Response(200, json=TRANSCRIPT_BRCA2_PAGE),
        httpx.Response(200, json=TRANSCRIPT_NONCODING_PAGE),
    ])
    respx.get(BASE + "translation/").mock(
        return_value=httpx.Response(200, json=TRANSLATION_REF_RESPONSE)
    )
    results = await diff_transcripts(
        ["ENST00000380152.7", "ENST00000614536.1"], client=client
    )
    assert len(results) == 1
    diff = results[0]
    # stable_id stored on Transcript is the bare ID (no version suffix)
    assert diff.reference_stable_id == "ENST00000380152"
    assert diff.candidate_stable_id == "ENST00000614536"


@respx.mock
@pytest.mark.asyncio
async def test_diff_transcripts_refseq_ids():
    """RefSeq stable IDs (NM_...) are accepted and fetched via /transcript/ endpoint."""
    # Build minimal RefSeq-like transcript fixtures
    refseq_v2_raw = {
        **TRANSCRIPT_BRCA2_RAW,
        "stable_id": "NM_001128425",
        "stable_id_version": 2,
    }
    refseq_v1_raw = {
        **TRANSCRIPT_BRCA2_RAW,
        "stable_id": "NM_001128425",
        "stable_id_version": 1,
        "sequence": {"sequence": "TTTTGGGGCCCCAAAA", "seq_checksum": "XYZ"},
        "five_prime_utr_seq": "TTTT",
        "three_prime_utr_seq": "AAAA",
        "translations": [],
    }

    client = TarkClient()
    respx.get(BASE + "transcript/").mock(side_effect=[
        httpx.Response(200, json=_page([refseq_v2_raw])),
        httpx.Response(200, json=_page([refseq_v1_raw])),
    ])
    respx.get(BASE + "translation/").mock(
        return_value=httpx.Response(200, json=TRANSLATION_REF_RESPONSE)
    )

    results = await diff_transcripts(
        ["NM_001128425.2", "NM_001128425.1"], client=client
    )
    assert len(results) == 1
    diff = results[0]
    assert diff.reference_stable_id == "NM_001128425"
    assert diff.candidate_stable_id == "NM_001128425"


@respx.mock
@pytest.mark.asyncio
async def test_diff_transcripts_mixed_refseq_ensembl():
    """RefSeq and Ensembl IDs can be diffed against each other."""
    refseq_raw = {
        **TRANSCRIPT_BRCA2_RAW,
        "stable_id": "NM_001128425",
        "stable_id_version": 2,
    }

    client = TarkClient()
    respx.get(BASE + "transcript/").mock(side_effect=[
        httpx.Response(200, json=_page([refseq_raw])),
        httpx.Response(200, json=TRANSCRIPT_NONCODING_PAGE),
    ])
    respx.get(BASE + "translation/").mock(
        return_value=httpx.Response(200, json=TRANSLATION_REF_RESPONSE)
    )

    results = await diff_transcripts(
        ["NM_001128425.2", "ENST00000614536"], client=client
    )
    assert len(results) == 1
    diff = results[0]
    assert diff.reference_stable_id == "NM_001128425"
    assert diff.candidate_stable_id == "ENST00000614536"


@respx.mock
@pytest.mark.asyncio
async def test_diff_transcripts_not_found_raises():
    """ValueError raised when a transcript cannot be found."""
    client = TarkClient()
    respx.get(BASE + "transcript/").mock(side_effect=[
        httpx.Response(200, json=_page([])),   # ref not found
        httpx.Response(200, json=TRANSCRIPT_NONCODING_PAGE),
    ])
    with pytest.raises(ValueError, match="Transcript not found"):
        await diff_transcripts(["NM_NOTREAL.1", "ENST00000614536"], client=client)
```

- [ ] **Step 4: Run all diff tests to verify they fail (expected at this point — implementation not done yet)**

```bash
cd /Users/adymun/tark-mcp && python -m pytest tests/unit/test_tools/test_diff.py -v 2>&1 | tail -30
```

Expected: Several FAILED tests because `_fetch_diff_pair` still calls `/transcript/diff/` and the mocks no longer set that up.

---

### Task 3: Implement the code change

**Files:**
- Modify: `src/tark_mcp/tools/diff.py`

- [ ] **Step 5: Apply the implementation from Task 1 Step 1**

Open `src/tark_mcp/tools/diff.py`. Replace the entire `_fetch_diff_pair` function (lines ~115–139) with:

```python
async def _resolve_transcript(
    stable_id: str,
    assembly: str,
    client: TarkClient,
) -> Transcript:
    """Fetch a single Transcript by stable ID, handling version suffixes and RefSeq IDs."""
    sid, version = stable_id, None
    if "." in stable_id:
        parts = stable_id.rsplit(".", 1)
        try:
            version = int(parts[1])
            sid = parts[0]
        except ValueError:
            pass

    params: dict = {"stable_id": sid, "expand_all": "true", "assembly_name": assembly}
    if version is not None:
        params["stable_id_version"] = version

    data = await client.get("transcript/", params=params)
    if not data:
        raise ValueError(
            f"Transcript not found: {stable_id} (assembly={assembly})"
        )

    transcripts = [Transcript.model_validate(r) for r in data]

    if version is not None:
        matching = [t for t in transcripts if t.stable_id_version == version]
        if matching:
            return matching[0]

    return max(transcripts, key=lambda t: t.latest_release_date or "")


async def _fetch_diff_pair(
    ref_stable_id: str,
    candidate_stable_id: str,
    ref_assembly: str,
    candidate_assembly: str,
    client: TarkClient,
) -> TranscriptDiff:
    ref, candidate = await asyncio.gather(
        _resolve_transcript(ref_stable_id, ref_assembly, client),
        _resolve_transcript(candidate_stable_id, candidate_assembly, client),
    )
    return await _build_diff(ref, candidate, client)
```

- [ ] **Step 6: Run all diff tests**

```bash
cd /Users/adymun/tark-mcp && python -m pytest tests/unit/test_tools/test_diff.py -v 2>&1 | tail -40
```

Expected: All tests PASS.

- [ ] **Step 7: Run full test suite to check for regressions**

```bash
cd /Users/adymun/tark-mcp && python -m pytest tests/unit/ -v 2>&1 | tail -40
```

Expected: All tests PASS (or the same set of pre-existing failures as before this change).

---

### Task 4: Update server.py docstring and commit

**Files:**
- Modify: `src/tark_mcp/server.py`

- [ ] **Step 8: Update the `tark_diff_transcripts` docstring**

In `src/tark_mcp/server.py`, replace the `stable_ids` arg description in `tark_diff_transcripts`:

Old:
```python
        stable_ids: List of ≥2 Ensembl transcript stable IDs; first is the reference
```

New:
```python
        stable_ids: List of ≥2 transcript stable IDs (Ensembl or RefSeq); first is the
            reference. Version suffixes are supported (e.g. 'ENST00000380152.7',
            'NM_001128425.2').
```

- [ ] **Step 9: Commit everything**

```bash
cd /Users/adymun/tark-mcp && git add src/tark_mcp/tools/diff.py src/tark_mcp/server.py tests/unit/test_tools/test_diff.py && git commit -m "feat: support versioned and RefSeq IDs in tark_diff_transcripts

Replace TARK /transcript/diff/ endpoint call with individual /transcript/
lookups so versioned IDs (e.g. NM_001128425.2, ENST00000380152.7) and
RefSeq stable IDs (NM_...) are resolved correctly before diffing.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

Expected: Commit succeeds with 3 files changed.
