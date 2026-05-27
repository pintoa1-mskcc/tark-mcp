# Transcript Release Version Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add release version information (e.g., "Ensembl v110") to transcript query responses alongside the existing release date.

**Architecture:** Extend the Transcript model's `_normalize()` method to extract and format version information from the `transcript_release_set` array that's already present in API responses. No changes to API calls or data fetching logic.

**Tech Stack:** Python 3.11+, Pydantic, pytest

---

## File Structure

**Modifications:**
- `src/tark_mcp/models.py` - Add field and extraction logic to Transcript model
- `src/tark_mcp/tools/versions.py` - Include version in history output
- `tests/unit/test_models.py` - Add tests for version extraction
- `CUSTOM_QUERIES.md` - Update examples to show new field

**No new files needed** - this is a pure extension of existing functionality.

---

## Task 1: Add Unit Tests for Release Version Extraction

**Files:**
- Modify: `tests/unit/test_models.py`

- [ ] **Step 1: Write test for single release version**

Add to `tests/unit/test_models.py` at the end of the file:

```python
def test_transcript_latest_release_version_single():
    """Test extraction of single release version."""
    data = {
        "stable_id": "ENST00000380152",
        "stable_id_version": 7,
        "assembly": "GRCh38",
        "loc_start": 100,
        "loc_end": 200,
        "loc_region": "13",
        "loc_strand": 1,
        "biotype": "protein_coding",
        "genes": [],
        "exons": [],
        "translations": [],
        "transcript_release_set": [
            {"source": "Ensembl", "shortname": "110", "release_date": "2023-04-01"}
        ],
    }
    t = Transcript.model_validate(data)
    assert t.latest_release_version == "Ensembl v110"
```

- [ ] **Step 2: Write test for multiple release versions**

Add below the previous test:

```python
def test_transcript_latest_release_version_multiple():
    """Test concatenation of multiple release versions."""
    data = {
        "stable_id": "ENST00000380152",
        "stable_id_version": 7,
        "assembly": "GRCh38",
        "loc_start": 100,
        "loc_end": 200,
        "loc_region": "13",
        "loc_strand": 1,
        "biotype": "protein_coding",
        "genes": [],
        "exons": [],
        "translations": [],
        "transcript_release_set": [
            {"source": "Ensembl", "shortname": "109", "release_date": "2022-12-01"},
            {"source": "Ensembl", "shortname": "110", "release_date": "2023-04-01"}
        ],
    }
    t = Transcript.model_validate(data)
    assert t.latest_release_version == "Ensembl v109, Ensembl v110"
```

- [ ] **Step 3: Write test for empty release set**

Add below the previous test:

```python
def test_transcript_latest_release_version_empty():
    """Test handling of empty release set."""
    data = {
        "stable_id": "ENST00000380152",
        "stable_id_version": 7,
        "assembly": "GRCh38",
        "loc_start": 100,
        "loc_end": 200,
        "loc_region": "13",
        "loc_strand": 1,
        "biotype": "protein_coding",
        "genes": [],
        "exons": [],
        "translations": [],
        "transcript_release_set": [],
    }
    t = Transcript.model_validate(data)
    assert t.latest_release_version is None
```

- [ ] **Step 4: Write test for dict release set**

Add below the previous test:

```python
def test_transcript_latest_release_version_dict():
    """Test handling of dict instead of list."""
    data = {
        "stable_id": "ENST00000380152",
        "stable_id_version": 7,
        "assembly": "GRCh38",
        "loc_start": 100,
        "loc_end": 200,
        "loc_region": "13",
        "loc_strand": 1,
        "biotype": "protein_coding",
        "genes": [],
        "exons": [],
        "translations": [],
        "transcript_release_set": {
            "source": "Ensembl",
            "shortname": "110",
            "release_date": "2023-04-01"
        },
    }
    t = Transcript.model_validate(data)
    assert t.latest_release_version == "Ensembl v110"
```

- [ ] **Step 5: Write test for missing fields in release set**

Add below the previous test:

```python
def test_transcript_latest_release_version_missing_fields():
    """Test handling of releases with missing source or shortname."""
    data = {
        "stable_id": "ENST00000380152",
        "stable_id_version": 7,
        "assembly": "GRCh38",
        "loc_start": 100,
        "loc_end": 200,
        "loc_region": "13",
        "loc_strand": 1,
        "biotype": "protein_coding",
        "genes": [],
        "exons": [],
        "translations": [],
        "transcript_release_set": [
            {"shortname": "110"},  # missing source
            {"source": "Ensembl"},  # missing shortname
        ],
    }
    t = Transcript.model_validate(data)
    assert t.latest_release_version is None
```

- [ ] **Step 6: Run tests to verify they fail**

Run:
```bash
cd /Users/adymun/tark-mcp
pytest tests/unit/test_models.py::test_transcript_latest_release_version_single -v
pytest tests/unit/test_models.py::test_transcript_latest_release_version_multiple -v
pytest tests/unit/test_models.py::test_transcript_latest_release_version_empty -v
pytest tests/unit/test_models.py::test_transcript_latest_release_version_dict -v
pytest tests/unit/test_models.py::test_transcript_latest_release_version_missing_fields -v
```

Expected: All tests FAIL with AttributeError: 'Transcript' object has no attribute 'latest_release_version'

- [ ] **Step 7: Commit the tests**

```bash
git add tests/unit/test_models.py
git commit -m "test: add tests for transcript release version extraction"
```

---

## Task 2: Add Release Version Field to Transcript Model

**Files:**
- Modify: `src/tark_mcp/models.py:67-83` (Transcript class field declarations)

- [ ] **Step 1: Add latest_release_version field to Transcript model**

In `src/tark_mcp/models.py`, add the new field after `latest_release_date` (around line 83):

```python
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
    latest_release_version: str | None = None  # NEW LINE
```

- [ ] **Step 2: Run tests to verify field exists but is None**

Run:
```bash
pytest tests/unit/test_models.py::test_transcript_latest_release_version_empty -v
```

Expected: PASS (empty release_set should result in None, which is the default)

Other tests should still FAIL because extraction logic isn't implemented yet.

- [ ] **Step 3: Commit the field addition**

```bash
git add src/tark_mcp/models.py
git commit -m "feat: add latest_release_version field to Transcript model"
```

---

## Task 3: Implement Release Version Extraction Logic

**Files:**
- Modify: `src/tark_mcp/models.py:114-122` (Transcript._normalize method)

- [ ] **Step 1: Add version extraction logic in _normalize method**

In `src/tark_mcp/models.py`, in the `Transcript._normalize()` method, add the following code immediately after the `latest_release_date` computation (after line 122):

```python
        # Compute latest_release_date from transcript_release_set array
        release_set = data.get("transcript_release_set", [])
        if isinstance(release_set, list) and release_set:
            dates = [r["release_date"] for r in release_set if r.get("release_date")]
            data["latest_release_date"] = max(dates) if dates else None
        elif isinstance(release_set, dict):
            data["latest_release_date"] = release_set.get("release_date")
        else:
            data["latest_release_date"] = None

        # Compute latest_release_version from transcript_release_set array
        release_set = data.get("transcript_release_set", [])
        if isinstance(release_set, list) and release_set:
            versions = [
                f"{r['source']} v{r['shortname']}"
                for r in release_set
                if r.get('source') and r.get('shortname')
            ]
            data["latest_release_version"] = ", ".join(versions) if versions else None
        elif isinstance(release_set, dict):
            source = release_set.get('source')
            shortname = release_set.get('shortname')
            if source and shortname:
                data["latest_release_version"] = f"{source} v{shortname}"
            else:
                data["latest_release_version"] = None
        else:
            data["latest_release_version"] = None
```

- [ ] **Step 2: Run all release version tests**

Run:
```bash
pytest tests/unit/test_models.py -k "release_version" -v
```

Expected: All 5 tests PASS

- [ ] **Step 3: Run all model tests to ensure nothing broke**

Run:
```bash
pytest tests/unit/test_models.py -v
```

Expected: All tests PASS

- [ ] **Step 4: Commit the implementation**

```bash
git add src/tark_mcp/models.py
git commit -m "feat: implement release version extraction in Transcript model

Extract source and shortname from transcript_release_set and format
as 'Source vVersion' (e.g., 'Ensembl v110'). Multiple releases are
comma-separated."
```

---

## Task 4: Update Version History Tool

**Files:**
- Modify: `src/tark_mcp/tools/versions.py:143-153` (version summary construction)

- [ ] **Step 1: Add release_version to version summaries**

In `src/tark_mcp/tools/versions.py`, in the `get_transcript_version_history()` function, add `release_version` to the version summary dictionary (around line 143):

```python
        version_summaries.append({
            "version": t.stable_id_version,
            "release_date": t.latest_release_date,
            "release_version": t.latest_release_version,  # NEW LINE
            "location": f"{t.loc_region}:{t.loc_start}-{t.loc_end}",
            "genomic_span": genomic_span,
            "transcript_length": len(t.sequence) if t.sequence else None,
            "exon_count": len(t.exons),
            "cds_length": cds_length,
            "protein_id": f"{t.translations[0].stable_id}.{t.translations[0].stable_id_version}" if t.translations else None,
            "protein_length": protein_length,
        })
```

- [ ] **Step 2: Run existing version tool tests**

Run:
```bash
pytest tests/unit/test_tools/ -k "version" -v
```

Expected: All tests PASS (or skip if no version tool tests exist)

- [ ] **Step 3: Commit the version tool update**

```bash
git add src/tark_mcp/tools/versions.py
git commit -m "feat: include release version in version history output"
```

---

## Task 5: Update Documentation

**Files:**
- Modify: `CUSTOM_QUERIES.md:22-26` (example output)

- [ ] **Step 1: Update example in CUSTOM_QUERIES.md**

In `CUSTOM_QUERIES.md`, update the example code (around line 22) to show the new field:

```python
from tark_mcp.tools.versions import get_transcript_all_versions

# Get all versions
versions = await get_transcript_all_versions("ENST00000263967", "GRCh38")

for v in versions:
    print(f"Version {v.stable_id_version}: {v.loc_start}-{v.loc_end}")
    print(f"  Release: {v.latest_release_version}")
    print(f"  Released: {v.latest_release_date}")
```

- [ ] **Step 2: Update version history example**

In `CUSTOM_QUERIES.md`, update the version history example (around line 36) to show the new field:

```python
from tark_mcp.tools.versions import get_transcript_version_history

history = await get_transcript_version_history("ENST00000263967", "GRCh38")

print(f"Found {history['version_count']} versions")

for ver in history['versions']:
    print(f"Version {ver['version']}")
    print(f"  Release: {ver['release_version']}")
    print(f"  Date: {ver['release_date']}")

for change in history['changes']:
    print(f"v{change['from_version']} → v{change['to_version']}:")
    print(f"  Transcript length: {change['transcript_length_delta']:+} bp")
    print(f"  CDS length: {change['cds_length_delta']:+} bp")
```

- [ ] **Step 3: Commit documentation updates**

```bash
git add CUSTOM_QUERIES.md
git commit -m "docs: update examples to show release version field"
```

---

## Task 6: Run Full Test Suite

**Files:**
- None (verification only)

- [ ] **Step 1: Run all unit tests**

Run:
```bash
pytest tests/unit/ -v
```

Expected: All tests PASS

- [ ] **Step 2: Run all integration tests**

Run:
```bash
pytest tests/integration/ -v
```

Expected: All tests PASS (the new field should appear in responses automatically)

- [ ] **Step 3: Run the full test suite**

Run:
```bash
pytest -v
```

Expected: All tests PASS

---

## Task 7: Manual Verification (Optional)

**Files:**
- None (manual testing)

- [ ] **Step 1: Test with example script**

Run the existing example script to see the new field in action:

```bash
cd /Users/adymun/tark-mcp
python examples/query_all_versions.py
```

Expected: Output should include lines like:
```
Release: Ensembl v110
```

- [ ] **Step 2: Test MCP server interactively (if desired)**

Start the MCP server and query a transcript to verify the field appears:

```bash
cd /Users/adymun/tark-mcp
# Start server and test with MCP client
```

Expected: JSON responses include `"latest_release_version": "Ensembl v110"`

---

## Success Criteria

All items should be true after completing this plan:

- [ ] `Transcript` model has `latest_release_version` field
- [ ] Single releases format as "Ensembl v110"
- [ ] Multiple releases concatenate as "Ensembl v109, Ensembl v110"
- [ ] Empty/missing release data results in `None`
- [ ] Version history tool includes `release_version` in output
- [ ] All unit tests pass
- [ ] All integration tests pass
- [ ] Documentation shows the new field in examples
- [ ] All changes committed with clear commit messages

---

## Rollback Plan

If issues are discovered, rollback is simple since this is purely additive:

```bash
git log --oneline -10  # Find commits before this feature
git reset --hard <commit-hash>
```

The feature is backward compatible, so partial rollback (reverting individual commits) is also safe.
