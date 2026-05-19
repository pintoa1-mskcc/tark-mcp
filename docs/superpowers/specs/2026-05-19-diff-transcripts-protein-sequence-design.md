# diff_transcripts — Protein & CDS Sequence Comparison — Design Spec

**Date:** 2026-05-19  
**Status:** Approved  
**Extends:** `2026-05-19-tark-mcp-design.md`

---

## Problem Statement

The existing `diff_transcripts` tool compares transcript structure (exon count, exon coordinates, CDS boundaries, biotype, cDNA sequence flag) but does not compare the actual CDS nucleotide sequence or the translated amino acid (protein) sequence. Bioinformaticians need to know whether the coding sequence and the protein product changed, not just whether the coordinates shifted.

---

## Approach

Enrich the existing `TranscriptDiff` model with sequence comparison fields. For each diff pair, fetch the full transcript objects (with `expand_all=true`) in parallel alongside the existing structural diff API call. Derive CDS sequences from the transcript sequence using the already-computed `cds_start`/`cds_end` offsets. Obtain protein sequences from the `Translation.sequence` field returned by `expand_all=true`.

---

## Model Changes

### `Translation` model — add `sequence` field

```python
class Translation(BaseModel):
    stable_id: str
    stable_id_version: int
    transcript_stable_id: str
    transcript_stable_id_version: int
    assembly: str
    sequence: str | None    # amino acid sequence; None if not returned by API
```

### `TranscriptDiff` model — add eight fields

```python
class TranscriptDiff(BaseModel):
    # existing fields (unchanged)
    reference_stable_id: str
    candidate_stable_id: str
    reference_assembly: str
    candidate_assembly: str
    biotype_changed: bool
    cds_changed: bool
    exon_count_changed: bool
    sequence_changed: bool
    exon_diffs: list[ExonDiff]

    # new fields
    reference_protein_coding: bool        # True if reference has a Translation with a non-None sequence
    candidate_protein_coding: bool        # True if candidate has a Translation with a non-None sequence
    cds_sequence_changed: bool | None     # None if either transcript is non-coding
    ref_cds_sequence: str | None          # transcript.sequence[cds_start:cds_end]; None if non-coding or sequence absent
    candidate_cds_sequence: str | None
    protein_sequence_changed: bool | None # None if either transcript is non-coding
    ref_protein_sequence: str | None      # from Translation.sequence; None if non-coding or absent
    candidate_protein_sequence: str | None
```

**Derivation rules:**
- `ref_cds_sequence` = `ref_transcript.sequence[cds_start:cds_end]` where `cds_start` and `cds_end` are 0-based offsets into the transcript sequence. Set to `None` if `cds_start` is `None` or `transcript.sequence` is `None`.
- `reference_protein_coding` = `True` if the reference transcript has at least one `Translation` with a non-`None` `sequence`. `False` otherwise. When multiple translations exist (rare), use the first one with a non-`None` sequence.
- `cds_sequence_changed` = `None` if either `ref_cds_sequence` or `candidate_cds_sequence` is `None`; otherwise `ref_cds_sequence != candidate_cds_sequence`.
- `protein_sequence_changed` = `None` if either `ref_protein_sequence` or `candidate_protein_sequence` is `None`; otherwise `ref_protein_sequence != candidate_protein_sequence`.

---

## Tool Behavior Changes

### `diff_transcripts` — fetch strategy per pair

**Before:** one API call per pair.  
**After:** three concurrent requests per pair (all fanned out via `asyncio.gather`):

1. `GET /api/transcript/diff/?diff_me_stable_id=<ref>&diff_with_stable_id=<candidate>` — structural diff (existing)
2. `GET /api/transcript/?stable_id=<ref>&expand_all=true` — full ref transcript (sequence + translations)
3. `GET /api/transcript/?stable_id=<candidate>&expand_all=true` — full candidate transcript

Cross-pair parallelism is unchanged: all pairs are still dispatched concurrently.

After fetching, the tool:
1. Populates the new `TranscriptDiff` fields using the derivation rules above.
2. When `assembly="both"`, each (ref, candidate) combination uses the assembly-specific transcript fetched for that pair.

---

## Error Handling

No new error conditions are introduced:

| Condition | Behaviour |
|-----------|-----------|
| `expand_all=true` returns no sequence for a transcript | `cds_sequence` and `protein_sequence` fields are `None`; `protein_coding` is `False`; no error raised |
| `expand_all=true` returns a `Translation` but no sequence | `protein_coding` is `False`; `protein_sequence` is `None` |
| Diff API call fails | Existing retry + `McpError` behaviour unchanged |
| Transcript fetch (for sequences) fails | Propagate as `McpError`; the `None` sentinel distinguishes "not applicable" from "error" |

`None` values for sequence fields signal to the LLM that comparison was not possible (non-coding transcript or data unavailable), not that an error occurred.

---

## Testing

**Unit tests** (`tests/unit/test_tools/test_diff.py`) additions:
- Coding pair: verify `cds_sequence_changed`, `protein_sequence_changed`, and the actual sequence fields are populated correctly.
- Non-coding reference: verify `reference_protein_coding=False`, `cds_sequence_changed=None`, `protein_sequence_changed=None`.
- Non-coding candidate: same as above for candidate side.
- Mixed pair (one coding, one non-coding): verify `None` sentinels on both comparison fields.
- Parallel fetch: verify three concurrent requests are issued per pair (mock assertion).

**Unit tests** (`tests/unit/test_models.py`) additions:
- `Translation` with `sequence` populated deserializes correctly.
- `Translation` with no `sequence` field in API response → `sequence=None`.
