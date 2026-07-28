# Needleman-Wunsch Protein Alignment for `diff_transcripts`

## Problem

`diff_transcripts` currently reports whether two transcripts' protein sequences
differ (`protein_sequence_changed`) and returns the two raw sequences, but it
does not report *where* they differ. Comparing the raw strings naively (or by
eyeballing them) miscounts amino acid positions whenever an indel is present,
because downstream residues shift. Positions must instead be derived from a
proper global sequence alignment.

## Goals

- Perform a real Needleman-Wunsch global alignment between the reference and
  candidate protein sequences.
- Report each differing alignment column (substitution, insertion, or
  deletion) with correctly counted 1-based amino acid positions on each side.
- Keep this additive: do not change the meaning of existing `TranscriptDiff`
  fields (`protein_sequence_changed`, `ref_protein_sequence`,
  `candidate_protein_sequence`).

## Non-goals

- Biologically-weighted scoring (e.g. BLOSUM62, affine gap penalties). A
  simple match/mismatch/gap linear scoring scheme is sufficient and matches
  what was manually validated against the real TARK data
  (NM_005027.3 vs ENST00000222254.8 → R234S, P313S).
- CDS/nucleotide alignment. This applies only to protein sequences.
- Multiple sequence alignment (only pairwise ref vs candidate).

## Design

### New module: `src/tark_mcp/tools/alignment.py`

A generic, pure, unit-testable global alignment function:

```python
def needleman_wunsch(
    seq_a: str,
    seq_b: str,
    match: int = 2,
    mismatch: int = -1,
    gap: int = -2,
) -> tuple[str, str]:
    """Global (Needleman-Wunsch) alignment of two sequences.

    Returns the two aligned sequences, each the same length, with '-'
    inserted at gap positions.
    """
```

Implementation: standard O(n*m) DP matrix with traceback, linear gap penalty.
No external dependencies (no biopython) — this keeps the tool's dependency
footprint unchanged.

### New model: `AminoAcidDiff` (in `models.py`)

Mirrors the existing `ExonDiff` pattern:

```python
class AminoAcidDiff(BaseModel):
    ref_position: int | None        # 1-based; None for insertion columns
    candidate_position: int | None  # 1-based; None for deletion columns
    change: str                     # "substitution" | "insertion" | "deletion"
    ref_residue: str | None
    candidate_residue: str | None
```

Unlike `ExonDiff` (which lists every exon, changed or not), `AminoAcidDiff`
entries are only emitted for **differing** alignment columns — identical
columns are omitted, since protein sequences can be hundreds of residues long
and most columns are typically unchanged.

### `TranscriptDiff` changes (in `models.py`)

Add one new field:

```python
protein_diffs: list[AminoAcidDiff] | None
```

- `None` when either side's protein sequence is unavailable (e.g. one or both
  transcripts are non-coding) — same condition currently used to leave
  `protein_sequence_changed` as `None`.
- Otherwise a list, empty if the aligned sequences are identical.

### `diff.py` changes

New helper, analogous to `_compute_exon_diffs`:

```python
def _compute_protein_diffs(ref_seq: str, cand_seq: str) -> list[AminoAcidDiff]:
    aligned_ref, aligned_cand = needleman_wunsch(ref_seq, cand_seq)
    ref_pos = 0
    cand_pos = 0
    diffs = []
    for r, c in zip(aligned_ref, aligned_cand):
        if r != "-":
            ref_pos += 1
        if c != "-":
            cand_pos += 1
        if r != c:
            if r == "-":
                change = "insertion"
            elif c == "-":
                change = "deletion"
            else:
                change = "substitution"
            diffs.append(AminoAcidDiff(
                ref_position=ref_pos if r != "-" else None,
                candidate_position=cand_pos if c != "-" else None,
                change=change,
                ref_residue=None if r == "-" else r,
                candidate_residue=None if c == "-" else c,
            ))
    return diffs
```

In `_build_diff`, after computing `ref_protein_seq` / `cand_protein_seq`:

```python
protein_diffs = None
if ref_protein_seq is not None and cand_protein_seq is not None:
    protein_diffs = _compute_protein_diffs(ref_protein_seq, cand_protein_seq)
```

Pass `protein_diffs=protein_diffs` into the returned `TranscriptDiff`.

## Testing

1. **`tests/unit/test_tools/test_alignment.py`** (new): unit tests for
   `needleman_wunsch` covering identical sequences, a single substitution, a
   single insertion, a single deletion, and a mixed case.
2. **`tests/unit/test_tools/test_diff.py`**: 
   - Unit test for `_compute_protein_diffs` directly.
   - Regression test using the real-world case validated in this session:
     NM_005027.3 vs ENST00000222254.8 protein sequences → expects exactly two
     diffs, `R234S` (substitution at ref_position=234, candidate_position=234)
     and `P313S` (substitution at ref_position=313, candidate_position=313).
   - Update any existing fixtures/assertions in this file that construct or
     assert on `TranscriptDiff` to account for the new `protein_diffs` field.

## Backward compatibility

Purely additive: existing fields and their semantics are unchanged. Existing
consumers of `TranscriptDiff` that don't read `protein_diffs` are unaffected.
