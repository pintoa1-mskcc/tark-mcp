# Design: Add Release Version Information to Transcript Queries

**Date:** 2026-05-27  
**Status:** Approved

## Problem Statement

Currently, when transcript data is queried from the TARK API, only the release date (e.g., "2023-04-01") is returned. Users need to also see the release version information (e.g., "Ensembl v110") to better understand which annotation release the data comes from.

## Goals

1. Expose release version information (e.g., "Ensembl v110") in transcript query responses
2. Support multiple releases when a transcript appears in more than one (e.g., "Ensembl v109, Ensembl v110")
3. Maintain full backward compatibility with existing code
4. Follow existing patterns in the codebase

## Non-Goals

- Changing the TARK API calls or query parameters
- Modifying the structure of existing fields
- Adding new MCP tool functions

## Architecture

### Data Flow

The TARK API already returns release information in the `transcript_release_set` array within transcript responses. The Transcript model currently extracts only the `latest_release_date` from this array. We'll extend the extraction logic to also capture version information.

```
TARK API Response
  └─ transcript_release_set: [
       {source: "Ensembl", shortname: "110", release_date: "2023-04-01"},
       {source: "Ensembl", shortname: "109", release_date: "2022-12-01"}
     ]
        ↓
Transcript._normalize() extracts:
  - latest_release_date: "2023-04-01" (max date - existing)
  - latest_release_version: "Ensembl v109, Ensembl v110" (all versions - new)
        ↓
Transcript model fields populated
```

### Key Principle

This change only extracts and formats data already present in the API response. No new API calls, no changes to query parameters, no modifications to data fetching logic.

## Data Model Changes

### Transcript Model

Add a new optional field to the `Transcript` class in `src/tark_mcp/models.py`:

```python
latest_release_version: str | None = None
```

**Field Semantics:**
- Contains formatted version string(s) from `transcript_release_set`
- Format: `"{source} v{shortname}"` where source is the database (e.g., "Ensembl") and shortname is the version number (e.g., "110")
- Examples:
  - Single release: `"Ensembl v110"`
  - Multiple releases: `"Ensembl v109, Ensembl v110"`
  - Multiple sources: `"Ensembl v110, GENCODE v42"` (if that occurs in practice)
- Value is `None` when `transcript_release_set` is empty or missing

### Version Extraction Logic

**Location:** `src/tark_mcp/models.py`, in the `Transcript._normalize()` method, immediately after the `latest_release_date` extraction (around line 122).

**Algorithm:**
1. Retrieve `transcript_release_set` from the API response data
2. Handle three cases:
   - **List of releases:** Extract source and shortname from each, format as `"{source} v{shortname}"`, join with ", "
   - **Single release (dict):** Extract and format as single version string
   - **Empty/missing:** Set to `None`
3. Assign the result to `data["latest_release_version"]`

**Implementation:**
```python
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

**Error Handling:**
- Use `.get()` for safe dictionary access
- Skip releases with missing `source` or `shortname` fields
- Never raise exceptions - worst case is `None` value
- Empty version list results in `None` (not empty string)

## Component Changes

### Primary: `src/tark_mcp/models.py`

**Changes:**
1. Add field declaration to `Transcript` class (line ~83)
2. Add extraction logic in `_normalize()` method (lines ~123-137)

**Impact:**
- All MCP tools that return Transcript objects will automatically include the new field
- No changes needed to tool implementations

### Secondary: `src/tark_mcp/tools/versions.py`

**Changes:**
Update `get_transcript_version_history()` to include release version in the version summaries.

**Location:** Line ~145 in the version summary dictionary construction

**Add:**
```python
"release_version": t.latest_release_version,
```

**Output format:**
```json
{
  "stable_id": "ENST00000263967",
  "assembly": "GRCh38",
  "version_count": 3,
  "versions": [
    {
      "version": 7,
      "release_date": "2023-04-01",
      "release_version": "Ensembl v110",
      "location": "13:32315508-32400266",
      ...
    }
  ]
}
```

## API Surface Changes

### MCP Tool Responses

All MCP tools that return transcript data will now include the `latest_release_version` field in their JSON responses. No changes to function signatures or parameters.

**Affected Tools:**
- `tark_get_transcript()` - returns single Transcript or list
- `tark_search_transcripts_by_region()` - returns list of Transcripts
- `tark_get_gene_transcripts()` - returns list of Transcripts
- Any other tool using the Transcript model

**Example Response:**
```json
{
  "stable_id": "ENST00000380152",
  "stable_id_version": 7,
  "assembly": "GRCh38",
  "latest_release_date": "2023-04-01",
  "latest_release_version": "Ensembl v110",
  "biotype": "protein_coding",
  ...
}
```

## Edge Cases

| Case | Input | Output |
|------|-------|--------|
| Single release | `[{source: "Ensembl", shortname: "110"}]` | `"Ensembl v110"` |
| Multiple releases | `[{source: "Ensembl", shortname: "109"}, {source: "Ensembl", shortname: "110"}]` | `"Ensembl v109, Ensembl v110"` |
| Empty array | `[]` | `None` |
| Dict instead of list | `{source: "Ensembl", shortname: "110"}` | `"Ensembl v110"` |
| Missing source | `[{shortname: "110"}]` | `None` |
| Missing shortname | `[{source: "Ensembl"}]` | `None` |
| Missing field entirely | (field not present) | `None` |
| Mixed sources | `[{source: "Ensembl", shortname: "110"}, {source: "GENCODE", shortname: "42"}]` | `"Ensembl v110, GENCODE v42"` |

## Testing Strategy

### Unit Tests

Add to `tests/unit/test_models.py`:

```python
def test_transcript_latest_release_version_single():
    """Test extraction of single release version."""
    data = {
        # ... minimal transcript fields ...
        "transcript_release_set": [
            {"source": "Ensembl", "shortname": "110", "release_date": "2023-04-01"}
        ],
    }
    t = Transcript.model_validate(data)
    assert t.latest_release_version == "Ensembl v110"

def test_transcript_latest_release_version_multiple():
    """Test concatenation of multiple release versions."""
    data = {
        # ... minimal transcript fields ...
        "transcript_release_set": [
            {"source": "Ensembl", "shortname": "109", "release_date": "2022-12-01"},
            {"source": "Ensembl", "shortname": "110", "release_date": "2023-04-01"}
        ],
    }
    t = Transcript.model_validate(data)
    assert t.latest_release_version == "Ensembl v109, Ensembl v110"

def test_transcript_latest_release_version_empty():
    """Test handling of empty release set."""
    data = {
        # ... minimal transcript fields ...
        "transcript_release_set": [],
    }
    t = Transcript.model_validate(data)
    assert t.latest_release_version is None

def test_transcript_latest_release_version_dict():
    """Test handling of dict instead of list."""
    data = {
        # ... minimal transcript fields ...
        "transcript_release_set": {
            "source": "Ensembl", 
            "shortname": "110", 
            "release_date": "2023-04-01"
        },
    }
    t = Transcript.model_validate(data)
    assert t.latest_release_version == "Ensembl v110"
```

### Integration Tests

Existing integration tests in `tests/integration/test_live_api.py` should pass without modification, demonstrating backward compatibility. The new field will simply appear in responses.

## Documentation Updates

Update `CUSTOM_QUERIES.md` to show the new field in examples:

```python
from tark_mcp.tools.versions import get_transcript_all_versions

versions = await get_transcript_all_versions("ENST00000263967", "GRCh38")

for v in versions:
    print(f"Version {v.stable_id_version}")
    print(f"  Release: {v.latest_release_version}")  # NEW
    print(f"  Released: {v.latest_release_date}")
```

## Backward Compatibility

✅ **Fully backward compatible:**

- New field is optional with default value of `None`
- No changes to existing fields or their types
- No changes to function signatures or parameters
- No changes to API query logic
- All existing code continues to work without modification
- Existing tests pass unchanged

## Implementation Plan

1. **Modify Transcript model** (`src/tark_mcp/models.py`)
   - Add `latest_release_version` field declaration
   - Add extraction logic in `_normalize()` method

2. **Update version history tool** (`src/tark_mcp/tools/versions.py`)
   - Add `release_version` to version summary dictionaries

3. **Add unit tests** (`tests/unit/test_models.py`)
   - Test single release version formatting
   - Test multiple release version concatenation
   - Test empty release set handling
   - Test dict vs list handling

4. **Update documentation** (`CUSTOM_QUERIES.md`)
   - Add examples showing the new field

5. **Run test suite**
   - Verify all existing tests pass
   - Verify new tests pass
   - Run integration tests if available

## Success Criteria

- [ ] Transcript model includes `latest_release_version` field
- [ ] Field correctly formats single releases as "Ensembl v110"
- [ ] Field correctly concatenates multiple releases as "Ensembl v109, Ensembl v110"
- [ ] Field is `None` when release data is missing
- [ ] Version history tool includes release version in output
- [ ] All existing tests pass
- [ ] New unit tests cover edge cases
- [ ] Documentation updated with examples
