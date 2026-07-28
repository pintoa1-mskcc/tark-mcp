# Needleman-Wunsch Protein Alignment Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a proper Needleman-Wunsch global alignment between reference and candidate protein sequences in `diff_transcripts`, reporting per-position amino acid differences with correctly counted 1-based positions.

**Architecture:** A new pure, dependency-free `alignment.py` module implements generic NW alignment. `diff.py` gains a helper that walks the aligned strings to produce `AminoAcidDiff` entries only for differing columns, using the same position-counting pattern already used by `_compute_exon_diffs`. A new `protein_diffs` field is added to `TranscriptDiff`, populated only when both sequences are available.

**Tech Stack:** Python, pydantic (models), pytest + pytest-asyncio + respx (tests). No new dependencies.

**Spec:** `docs/superpowers/specs/2026-07-28-nw-protein-alignment-design.md`

---

### Task 1: Needleman-Wunsch alignment module

**Files:**
- Create: `src/tark_mcp/tools/alignment.py`
- Test: `tests/unit/test_tools/test_alignment.py`

- [ ] **Step 1: Write the failing tests**

Create `tests/unit/test_tools/test_alignment.py`:

```python
from tark_mcp.tools.alignment import needleman_wunsch


def test_identical_sequences_align_with_no_gaps():
    aligned_a, aligned_b = needleman_wunsch("MATS", "MATS")
    assert aligned_a == "MATS"
    assert aligned_b == "MATS"


def test_single_substitution_aligns_without_gaps():
    aligned_a, aligned_b = needleman_wunsch("MAT", "MST")
    assert aligned_a == "MAT"
    assert aligned_b == "MST"


def test_insertion_in_second_sequence_produces_gap_in_first():
    aligned_a, aligned_b = needleman_wunsch("AC", "ABC")
    assert aligned_a == "A-C"
    assert aligned_b == "ABC"


def test_deletion_in_second_sequence_produces_gap_in_second():
    aligned_a, aligned_b = needleman_wunsch("ABC", "AC")
    assert aligned_a == "ABC"
    assert aligned_b == "A-C"


def test_aligned_sequences_are_always_equal_length():
    aligned_a, aligned_b = needleman_wunsch("MPIGSKERP", "MVLSPAD")
    assert len(aligned_a) == len(aligned_b)


def test_empty_sequence_aligns_as_all_gaps():
    aligned_a, aligned_b = needleman_wunsch("", "ABC")
    assert aligned_a == "---"
    assert aligned_b == "ABC"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/unit/test_tools/test_alignment.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'tark_mcp.tools.alignment'`

- [ ] **Step 3: Implement the alignment module**

Create `src/tark_mcp/tools/alignment.py`:

```python
from __future__ import annotations


def needleman_wunsch(
    seq_a: str,
    seq_b: str,
    match: int = 2,
    mismatch: int = -1,
    gap: int = -2,
) -> tuple[str, str]:
    """Global (Needleman-Wunsch) alignment of two sequences.

    Returns a pair of equal-length strings representing the alignment of
    seq_a and seq_b, with '-' inserted at gap positions. Uses a linear gap
    penalty and a simple match/mismatch scoring scheme (no substitution
    matrix).
    """
    n, m = len(seq_a), len(seq_b)

    # score[i][j] = best alignment score for seq_a[:i] vs seq_b[:j]
    score = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        score[i][0] = i * gap
    for j in range(1, m + 1):
        score[0][j] = j * gap

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = match if seq_a[i - 1] == seq_b[j - 1] else mismatch
            score[i][j] = max(
                score[i - 1][j - 1] + s,
                score[i - 1][j] + gap,
                score[i][j - 1] + gap,
            )

    aligned_a: list[str] = []
    aligned_b: list[str] = []
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            s = match if seq_a[i - 1] == seq_b[j - 1] else mismatch
            if score[i][j] == score[i - 1][j - 1] + s:
                aligned_a.append(seq_a[i - 1])
                aligned_b.append(seq_b[j - 1])
                i -= 1
                j -= 1
                continue
        if i > 0 and score[i][j] == score[i - 1][j] + gap:
            aligned_a.append(seq_a[i - 1])
            aligned_b.append("-")
            i -= 1
            continue
        aligned_a.append("-")
        aligned_b.append(seq_b[j - 1])
        j -= 1

    aligned_a.reverse()
    aligned_b.reverse()
    return "".join(aligned_a), "".join(aligned_b)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/unit/test_tools/test_alignment.py -v`
Expected: PASS (6 passed)

- [ ] **Step 5: Commit**

```bash
git add src/tark_mcp/tools/alignment.py tests/unit/test_tools/test_alignment.py
git commit -m "feat: add Needleman-Wunsch global alignment helper"
```

---

### Task 2: `AminoAcidDiff` model and `TranscriptDiff.protein_diffs` field

**Files:**
- Modify: `src/tark_mcp/models.py`

- [ ] **Step 1: Add the `AminoAcidDiff` model**

In `src/tark_mcp/models.py`, right after the existing `ExonDiff` class definition:

```python
class ExonDiff(BaseModel):
    order: int
    change: str  # "added", "removed", "modified", "unchanged"
    ref_coords: tuple[int, int] | None
    candidate_coords: tuple[int, int] | None


class AminoAcidDiff(BaseModel):
    ref_position: int | None        # 1-based; None for insertion columns
    candidate_position: int | None  # 1-based; None for deletion columns
    change: str  # "substitution", "insertion", "deletion"
    ref_residue: str | None
    candidate_residue: str | None
```

- [ ] **Step 2: Add `protein_diffs` field to `TranscriptDiff`**

In `src/tark_mcp/models.py`, add the new field to the end of `TranscriptDiff`:

```python
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
    protein_diffs: list[AminoAcidDiff] | None
```

- [ ] **Step 3: Verify the module still imports cleanly**

Run: `python -c "from tark_mcp.models import AminoAcidDiff, TranscriptDiff; print('ok')"`
Expected: `ok`

- [ ] **Step 4: Commit**

```bash
git add src/tark_mcp/models.py
git commit -m "feat: add AminoAcidDiff model and TranscriptDiff.protein_diffs field"
```

---

### Task 3: `_compute_protein_diffs` helper in `diff.py`

**Files:**
- Modify: `src/tark_mcp/tools/diff.py`
- Test: `tests/unit/test_tools/test_diff.py`

- [ ] **Step 1: Write the failing unit test**

Add to `tests/unit/test_tools/test_diff.py` (near the top, after imports — add `_compute_protein_diffs` to the import from `tark_mcp.tools.diff`):

```python
from tark_mcp.tools.diff import _compute_protein_diffs, _resolve_transcript, diff_transcripts
```

Then add this test function (anywhere after the imports, e.g. right before the first `@respx.mock` test):

```python
def test_compute_protein_diffs_reports_substitution_with_aligned_positions():
    """A single substitution between equal-length sequences is reported at the
    same 1-based position on both sides."""
    diffs = _compute_protein_diffs("MAT", "MST")
    assert len(diffs) == 1
    d = diffs[0]
    assert d.change == "substitution"
    assert d.ref_position == 2
    assert d.candidate_position == 2
    assert d.ref_residue == "A"
    assert d.candidate_residue == "S"


def test_compute_protein_diffs_reports_deletion_with_none_candidate_position():
    """A residue present only in the reference is a deletion; candidate_position
    is None since there's no corresponding aligned residue on that side."""
    diffs = _compute_protein_diffs("ABC", "AC")
    assert len(diffs) == 1
    d = diffs[0]
    assert d.change == "deletion"
    assert d.ref_position == 2
    assert d.candidate_position is None
    assert d.ref_residue == "B"
    assert d.candidate_residue is None


def test_compute_protein_diffs_reports_insertion_with_none_ref_position():
    """A residue present only in the candidate is an insertion; ref_position is
    None since there's no corresponding aligned residue on that side."""
    diffs = _compute_protein_diffs("AC", "ABC")
    assert len(diffs) == 1
    d = diffs[0]
    assert d.change == "insertion"
    assert d.ref_position is None
    assert d.candidate_position == 2
    assert d.ref_residue is None
    assert d.candidate_residue == "B"


def test_compute_protein_diffs_identical_sequences_returns_empty_list():
    assert _compute_protein_diffs("MATS", "MATS") == []


def test_compute_protein_diffs_real_world_regression_r234s_p313s():
    """Regression test: NM_005027.3 vs ENST00000222254.8 (GRCh37) protein
    sequences, validated manually against the live TARK API. Expect exactly
    two substitutions: R234S and P313S, with no indels."""
    ref_seq = (
        "MAGPEGFQYRALYPFRRERPEDLELLPGDVLVVSRAALQALGVAEGGERCPQSVGWMPGLNERTRQRGDF"
        "PGTYVEFLGPVALARPGPRPRGPRPLPARPRDGAPEPGLTLPDLPEQFSPPDVAPPLLVKLVEAIERTGL"
        "DSESHYRPELPAPRTDWSLSDVDQWDTAALADGIKSFLLALPAPLVTPEASAEARRALREAAGPVGPALE"
        "PPTLPLHRALTLRFLLQHLGRVARRAPALGPAVRALGATFGPLLLRAPPPPSSPPPGGAPDGSEPSPDFP"
        "ALLVEKLLQEHLEEQEVAPPALPPKPPKAKPAPTVLANGGSPPSLQDAEWYWGDISREEVNEKLRDTPDG"
        "TFLVRDASSKIQGEYTLTLRKGGNNKLIKVFHRDGHYGFSEPLTFCSVVDLINHYRHESLAQYNAKLDTR"
        "LLYPVSKYQQDQIVKEDSVEAVGAQLKVYHQQYQDKSREYDQLYEEYTRTSQELQMKRTAIEAFNETIKI"
        "FEEQGQTQEKCSKEYLERFRREGNEKEMQRILLNSERLKSRIAEIHESRTKLEQQLRAQASDNREIDKRM"
        "NSLKPDLMQLRKIRDQYLVWLTQKGARQKKINEWLGIKNETEDQYALMEDEDDLPHHEERTWYVGKINR"
        "TQAEEMLSGKRDGTFLIRESSQRGCYACSVVVDGDTKHCVIYRTATGFGFAEPYNLYGSLKELVLHYQHA"
        "SLVQHNDALTVTLAHPVRAPGPGPPPAAR"
    )
    cand_seq = (
        "MAGPEGFQYRALYPFRRERPEDLELLPGDVLVVSRAALQALGVAEGGERCPQSVGWMPGLNERTRQRGDF"
        "PGTYVEFLGPVALARPGPRPRGPRPLPARPRDGAPEPGLTLPDLPEQFSPPDVAPPLLVKLVEAIERTGL"
        "DSESHYRPELPAPRTDWSLSDVDQWDTAALADGIKSFLLALPAPLVTPEASAEARRALREAAGPVGPALE"
        "PPTLPLHRALTLRFLLQHLGRVASRAPALGPAVRALGATFGPLLLRAPPPPSSPPPGGAPDGSEPSPDFP"
        "ALLVEKLLQEHLEEQEVAPPALPPKPPKAKPASTVLANGGSPPSLQDAEWYWGDISREEVNEKLRDTPDG"
        "TFLVRDASSKIQGEYTLTLRKGGNNKLIKVFHRDGHYGFSEPLTFCSVVDLINHYRHESLAQYNAKLDTR"
        "LLYPVSKYQQDQIVKEDSVEAVGAQLKVYHQQYQDKSREYDQLYEEYTRTSQELQMKRTAIEAFNETIKI"
        "FEEQGQTQEKCSKEYLERFRREGNEKEMQRILLNSERLKSRIAEIHESRTKLEQQLRAQASDNREIDKRM"
        "NSLKPDLMQLRKIRDQYLVWLTQKGARQKKINEWLGIKNETEDQYALMEDEDDLPHHEERTWYVGKINR"
        "TQAEEMLSGKRDGTFLIRESSQRGCYACSVVVDGDTKHCVIYRTATGFGFAEPYNLYGSLKELVLHYQHA"
        "SLVQHNDALTVTLAHPVRAPGPGPPPAAR"
    )
    diffs = _compute_protein_diffs(ref_seq, cand_seq)
    assert len(diffs) == 2

    r234s = diffs[0]
    assert r234s.change == "substitution"
    assert r234s.ref_position == 234
    assert r234s.candidate_position == 234
    assert r234s.ref_residue == "R"
    assert r234s.candidate_residue == "S"

    p313s = diffs[1]
    assert p313s.change == "substitution"
    assert p313s.ref_position == 313
    assert p313s.candidate_position == 313
    assert p313s.ref_residue == "P"
    assert p313s.candidate_residue == "S"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/unit/test_tools/test_diff.py -k compute_protein_diffs -v`
Expected: FAIL with `ImportError: cannot import name '_compute_protein_diffs'`

- [ ] **Step 3: Implement `_compute_protein_diffs`**

In `src/tark_mcp/tools/diff.py`, add the import and the new helper function right after `_compute_exon_diffs`:

```python
from tark_mcp.client import TarkClient
from tark_mcp.models import AminoAcidDiff, ExonDiff, Transcript, TranscriptDiff
from tark_mcp.tools.alignment import needleman_wunsch
```

(This replaces the existing `from tark_mcp.client import TarkClient` and `from tark_mcp.models import Transcript, ExonDiff, TranscriptDiff` import lines at the top of the file — merge them into the two lines above, alphabetized.)

```python
def _compute_protein_diffs(ref_seq: str, cand_seq: str) -> list[AminoAcidDiff]:
    """Align two protein sequences with Needleman-Wunsch and return one
    AminoAcidDiff per differing column, with correctly counted 1-based
    positions on each side. Identical columns are omitted."""
    aligned_ref, aligned_cand = needleman_wunsch(ref_seq, cand_seq)

    ref_pos = 0
    cand_pos = 0
    diffs: list[AminoAcidDiff] = []
    for r, c in zip(aligned_ref, aligned_cand):
        if r != "-":
            ref_pos += 1
        if c != "-":
            cand_pos += 1
        if r == c:
            continue
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

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/unit/test_tools/test_diff.py -k compute_protein_diffs -v`
Expected: PASS (5 passed)

- [ ] **Step 5: Commit**

```bash
git add src/tark_mcp/tools/diff.py tests/unit/test_tools/test_diff.py
git commit -m "feat: compute amino-acid-level protein diffs via NW alignment"
```

---

### Task 4: Wire `protein_diffs` into `_build_diff` and update existing tests

**Files:**
- Modify: `src/tark_mcp/tools/diff.py`
- Modify: `tests/unit/test_tools/test_diff.py`

- [ ] **Step 1: Update the existing coding-pair test to assert on `protein_diffs`**

In `tests/unit/test_tools/test_diff.py`, in `test_diff_transcripts_coding_pair_populates_all_sequence_fields`, add these assertions right after the existing `assert diff.protein_sequence_changed is True` line:

```python
    assert diff.protein_sequence_changed is True

    # ref="MPIGSKERP" (9aa), candidate="MVLSPAD" (7aa) aligned via NW
    # (match=2, mismatch=-1, gap=-2) gives: "MPIGSKERP" / "M-VLS-PAD"
    assert [d.change for d in diff.protein_diffs] == [
        "deletion", "substitution", "substitution",
        "deletion", "substitution", "substitution", "substitution",
    ]
    assert diff.protein_diffs[0].ref_position == 2
    assert diff.protein_diffs[0].candidate_position is None
    assert diff.protein_diffs[0].ref_residue == "P"
    assert diff.protein_diffs[1].ref_position == 3
    assert diff.protein_diffs[1].candidate_position == 2
    assert diff.protein_diffs[1].ref_residue == "I"
    assert diff.protein_diffs[1].candidate_residue == "V"
```

- [ ] **Step 2: Update the non-coding test to assert `protein_diffs is None`**

In `test_diff_transcripts_noncoding_ref_sets_none_sentinels`, add after `assert diff.ref_protein_sequence is None`:

```python
    assert diff.ref_protein_sequence is None
    assert diff.candidate_protein_sequence is None
    assert diff.protein_diffs is None
```

(Remove the duplicate `assert diff.candidate_protein_sequence is None` that already exists right below — keep only one copy.)

- [ ] **Step 3: Update the mixed-pair test to assert `protein_diffs is None`**

In `test_diff_transcripts_mixed_pair_sets_none_sentinels`, add after `assert diff.protein_sequence_changed is None`:

```python
    assert diff.protein_sequence_changed is None
    assert diff.protein_diffs is None
```

- [ ] **Step 4: Run the full diff test file to verify these new assertions fail**

Run: `pytest tests/unit/test_tools/test_diff.py -v`
Expected: FAIL — `AttributeError: 'TranscriptDiff' object has no attribute 'protein_diffs'` (or similar, since `_build_diff` doesn't set it yet)

- [ ] **Step 5: Wire `protein_diffs` into `_build_diff`**

In `src/tark_mcp/tools/diff.py`, in `_build_diff`, after the block that computes `protein_changed`:

```python
    protein_changed: bool | None = None
    if ref_protein_seq is not None and cand_protein_seq is not None:
        protein_changed = ref_protein_seq != cand_protein_seq

    protein_diffs: list[AminoAcidDiff] | None = None
    if ref_protein_seq is not None and cand_protein_seq is not None:
        protein_diffs = _compute_protein_diffs(ref_protein_seq, cand_protein_seq)
```

Then add `protein_diffs=protein_diffs` to the `TranscriptDiff(...)` construction at the end of `_build_diff`, right after `candidate_protein_sequence=cand_protein_seq,`:

```python
        protein_sequence_changed=protein_changed,
        ref_protein_sequence=ref_protein_seq,
        candidate_protein_sequence=cand_protein_seq,
        protein_diffs=protein_diffs,
    )
```

- [ ] **Step 6: Run the full diff test file to verify it passes**

Run: `pytest tests/unit/test_tools/test_diff.py -v`
Expected: PASS (all tests in file pass)

- [ ] **Step 7: Commit**

```bash
git add src/tark_mcp/tools/diff.py tests/unit/test_tools/test_diff.py
git commit -m "feat: populate TranscriptDiff.protein_diffs from NW alignment"
```

---

### Task 5: Full test suite and manual verification

**Files:** None (verification only)

- [ ] **Step 1: Run the full unit test suite**

Run: `pytest tests/unit/ -v`
Expected: All tests PASS, no failures or errors.

- [ ] **Step 2: Verify MCP tool docstring/behavior is unaffected**

Run: `grep -n "diff_transcripts" src/tark_mcp/server.py`
Confirm the `@mcp.tool()` wrapper for `diff_transcripts` in `server.py` needs no changes (it should already just pass through to `tools/diff.py`'s `diff_transcripts` and return the model — no server.py changes are expected since `protein_diffs` is just a new field on the existing return model).

- [ ] **Step 3: Manual smoke test against the live API (optional, requires network)**

Run:
```bash
python -c "
import asyncio
from tark_mcp.tools.diff import diff_transcripts

async def main():
    results = await diff_transcripts(
        stable_ids=['NM_005027.3', 'ENST00000222254.8'],
        assemblies=['GRCh37', 'GRCh37'],
    )
    for d in results[0].protein_diffs:
        print(d)

asyncio.run(main())
"
```
Expected: exactly two `AminoAcidDiff` entries printed — one at position 234 (R→S) and one at position 313 (P→S).

- [ ] **Step 4: Final commit if any cleanup was needed**

If steps 1-3 required no code changes, no commit is needed here. If any issues were found and fixed, commit with:

```bash
git add -A
git commit -m "fix: address issues found during full-suite verification"
```
