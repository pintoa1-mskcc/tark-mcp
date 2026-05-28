# Custom TARK API Queries

This document describes custom query functions built on top of tark-mcp that provide additional functionality beyond the standard MCP tools.

## Version History Queries

**Location:** `src/tark_mcp/tools/versions.py`

The standard `get_transcript()` tool returns only the **latest version** of a transcript. The custom version query functions provide access to **all historical versions** stored in TARK.

### Functions

#### `get_transcript_all_versions(stable_id, assembly="GRCh38", client=None)`

Returns all versions of a transcript without deduplication.

```python
from tark_mcp.tools.versions import get_transcript_all_versions

# Get all versions
versions = await get_transcript_all_versions("ENST00000263967", "GRCh38")

for v in versions:
    print(f"Version {v.stable_id_version}: {v.loc_start}-{v.loc_end}")
    print(f"  Release: {v.latest_release_version}")
    print(f"  Released: {v.latest_release_date}")
```

#### `get_transcript_version_history(stable_id, assembly="GRCh38", client=None)`

Returns detailed version history with statistics and change tracking.

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

### Example Usage

See `examples/query_all_versions.py` for a complete working example:

```bash
cd /Users/adymun/tark-mcp
python examples/query_all_versions.py
```

### Why Custom Queries?

The standard TARK MCP tools follow the design spec requirement to deduplicate results and return only the latest version. This is correct behavior for most use cases, but researchers studying annotation history need access to all versions.

The custom queries:
- Call the TARK API directly via `TarkClient`
- Parse responses without deduplication
- Provide version-to-version change tracking
- Support both GRCh37 and GRCh38 assemblies

### Implementation Details

The key difference is skipping the deduplication step:

```python
# Standard tool (deduplicates)
def _deduplicate(transcripts: list[Transcript]) -> list[Transcript]:
    best: dict[tuple, Transcript] = {}
    for t in transcripts:
        key = (t.assembly, t.stable_id, t.stable_id_version)
        if existing is None or t.latest_release_date > existing.latest_release_date:
            best[key] = t
    return list(best.values())

# Custom query (no deduplication)
transcripts = [Transcript.model_validate(r) for r in api_data]
return sorted(transcripts, key=lambda t: t.stable_id_version)
```

## Adding Your Own Custom Queries

To add new custom queries:

1. Create a new module in `src/tark_mcp/tools/`
2. Import `TarkClient` and relevant models
3. Use `client.get(endpoint, params)` to call the API directly
4. Parse responses with Pydantic models
5. Add examples to `examples/`

Example template:

```python
from tark_mcp.client import TarkClient
from tark_mcp.models import Transcript

async def my_custom_query(stable_id: str, client: TarkClient | None = None):
    if client is None:
        client = TarkClient()
    
    # Call API directly
    data = await client.get("transcript/", params={"stable_id": stable_id})
    
    # Process results
    transcripts = [Transcript.model_validate(r) for r in data]
    
    # Your custom logic here
    return transcripts
```

## Available API Endpoints

The TARK API provides these endpoints (see design spec for details):

- `transcript/` - Transcript queries
- `transcript/search/` - Search by gene name/ID
- `transcript/diff/` - Compare transcripts
- `transcript/manelist/` - MANE transcripts
- `release/nopagination/` - Release metadata

All endpoints support `expand` and `expand_all` parameters to include nested data.
