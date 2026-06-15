# tark-mcp

An MCP server that exposes [TARK (Transcript ARKive)](https://tark.ensembl.org/) Ensembl transcript data as tools for AI coding assistants.

## Tools

| Tool | Description |
|------|-------------|
| `tark_get_releases` | List available TARK releases |
| `tark_get_transcript` | Fetch a transcript by Ensembl stable ID |
| `tark_search_transcripts_by_region` | Find transcripts overlapping a genomic region |
| `tark_get_gene_transcripts` | All transcripts for a gene symbol or Ensembl gene ID |
| `tark_get_transcript_sequence` | cDNA sequence for a transcript |
| `tark_get_transcript_exons` | Ordered exon list with genomic coordinates |
| `tark_get_protein_for_transcript` | Protein/translation stable ID for a transcript |
| `tark_get_mane_transcripts` | MANE Select and Plus Clinical transcripts |
| `tark_diff_transcripts` | Structural + CDS + protein diff between transcripts |

## Installation

### 1. Install the MCP server

```bash
pip install -e .
```

### 2. Register with Copilot CLI

Add to `~/.copilot/mcp.json`:

```json
{
  "mcpServers": {
    "tark": {
      "command": "tark-mcp",
      "args": []
    }
  }
}
```

### 3. Install the Copilot CLI skill (optional but recommended)

The `plugin/` directory contains a Copilot CLI skill that teaches agents how to use the TARK tools effectively.

```bash
# From the repo root — symlink so skill stays in sync with the repo
ln -s "$(pwd)/plugin" ~/.copilot/installed-plugins/tark-mcp
```

Or copy for a stable install:

```bash
cp -r plugin ~/.copilot/installed-plugins/tark-mcp
```

Once installed, agents will automatically use the skill when querying transcript data.

## Configuration

Environment variables (all optional):

| Variable | Default | Description |
|----------|---------|-------------|
| `TARK_BASE_URL` | `https://tark.ensembl.org/api/` | TARK API base URL |
| `TARK_CACHE_TTL` | `3600` | Cache TTL in seconds |
| `TARK_REQUEST_TIMEOUT` | `30` | HTTP timeout in seconds |
| `TARK_MAX_RETRIES` | `3` | Max retries on transient errors |

## Development

```bash
pip install -e ".[dev]"
pytest
```

## Custom Queries

For queries beyond the MCP tools (e.g. full version history for a transcript), see [CUSTOM_QUERIES.md](CUSTOM_QUERIES.md).
