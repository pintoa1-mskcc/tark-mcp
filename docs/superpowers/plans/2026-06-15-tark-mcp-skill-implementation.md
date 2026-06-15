# tark-mcp Skill Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Create a Copilot CLI plugin containing a workflow-oriented skill for the tark-mcp MCP server, and install it locally.

**Architecture:** A `plugin/` directory in the repo root contains a `plugin.json` manifest and a `skills/tark-mcp/SKILL.md` file. The plugin is symlinked into `~/.copilot/installed-plugins/tark-mcp/` so Copilot CLI discovers the skill automatically. A `README.md` documents setup.

**Tech Stack:** Markdown, JSON, Copilot CLI plugin format (same as superpowers-marketplace)

---

### Task 1: Create plugin manifest

**Files:**
- Create: `plugin/plugin.json`

- [ ] **Step 1: Create the plugin directory and manifest**

Create `plugin/plugin.json` with this exact content:

```json
{
  "name": "tark-mcp",
  "displayName": "TARK MCP",
  "description": "Workflow skill for querying Ensembl transcript data via the TARK MCP server",
  "version": "0.1.0",
  "author": {
    "name": "tark-mcp contributors"
  },
  "license": "MIT",
  "keywords": [
    "genomics",
    "ensembl",
    "transcripts",
    "bioinformatics",
    "mcp"
  ],
  "skills": "./skills/"
}
```

- [ ] **Step 2: Verify the file was created**

Run: `cat plugin/plugin.json`

Expected: JSON printed to stdout with `"name": "tark-mcp"`.

- [ ] **Step 3: Commit**

```bash
git add plugin/plugin.json
git commit -m "feat: add Copilot CLI plugin manifest"
```

---

### Task 2: Create the skill file

**Files:**
- Create: `plugin/skills/tark-mcp/SKILL.md`

- [ ] **Step 1: Create the skills directory**

```bash
mkdir -p plugin/skills/tark-mcp
```

- [ ] **Step 2: Create the SKILL.md**

Create `plugin/skills/tark-mcp/SKILL.md` with this exact content:

````markdown
---
name: tark-mcp
description: Use when querying Ensembl transcript data via TARK — transcripts by stable ID or gene, MANE Select/Plus Clinical, CDS and protein sequences, exon structure, genomic locus queries, GRCh37/GRCh38 assembly comparison, transcript structural diffs
---

# TARK MCP Skill

Use this skill when a user asks about Ensembl transcripts, genes, MANE transcripts, CDS/protein sequences, exon structure, genomic loci, or cross-assembly/cross-version comparison.

## Quick Tool Index

| Tool | What it does |
|------|-------------|
| `tark_get_releases` | List TARK releases with assembly, date, and source metadata |
| `tark_get_transcript` | Fetch a single transcript by Ensembl stable ID (with exons, CDS, genes, translations) |
| `tark_search_transcripts_by_region` | Find all transcripts overlapping a genomic locus (0-based half-open coords) |
| `tark_get_gene_transcripts` | All transcripts for a gene symbol (e.g. `BRCA2`) or Ensembl gene ID |
| `tark_get_transcript_sequence` | cDNA sequence for a transcript |
| `tark_get_transcript_exons` | Ordered exon list with 0-based genomic coordinates |
| `tark_get_protein_for_transcript` | Protein/translation stable ID and version for a transcript |
| `tark_get_mane_transcripts` | MANE Select and MANE Plus Clinical transcripts, optionally filtered by gene |
| `tark_diff_transcripts` | Structural + CDS + protein diff; first ID is the reference |

---

## Key Concepts

- **Coordinates are 0-based half-open.** The TARK API uses 1-based internally; tark-mcp converts on the way out. All `loc_start`/`loc_end` values you receive are 0-based half-open.
- **`assembly` parameter** accepts `"GRCh38"` (default), `"GRCh37"`, or `"both"`. Passing `"both"` returns a combined list; each record carries its own `assembly` field.
- **Versioned stable IDs** are accepted by most tools (e.g. `ENST00000380152.7`). The version is used as a filter if provided.
- **`tark_get_transcript` deduplicates** — it returns the single most-recently-released record per assembly. Use the `versions.py` custom queries (not MCP tools) if you need all historical versions.
- **`tark_diff_transcripts`: first ID is always the reference.** All subsequent IDs are compared against it.

---

## Workflow Recipes

### 1. Find the canonical transcript for a gene

**Use case:** "What is the canonical/MANE transcript for TP53?"

```
tark_get_mane_transcripts(gene_identifier="TP53")
```

Returns MANE Select and MANE Plus Clinical entries for the gene. If no MANE transcript exists (e.g. non-human genes, older annotations), fall back to:

```
tark_get_gene_transcripts(gene_identifier="TP53", assembly="GRCh38")
```

Then filter the results by `biotype == "protein_coding"` and pick the longest or highest-versioned transcript.

---

### 2. Look up a specific transcript

**Use case:** "Get me everything about ENST00000380152"

```
tark_get_transcript(stable_id="ENST00000380152", assembly="GRCh38")
```

The result includes exons, genes, CDS boundaries, and translation IDs. Chain additional tools as needed:

```
tark_get_transcript_sequence(stable_id="ENST00000380152")   # cDNA
tark_get_transcript_exons(stable_id="ENST00000380152")      # exon coords
tark_get_protein_for_transcript(stable_id="ENST00000380152") # translation ID
```

Pass a versioned ID (`ENST00000380152.7`) to pin to a specific version.

---

### 3. Investigate a genomic locus

**Use case:** "What transcripts overlap chr13:32,315,480-32,400,268?"

Strip `chr` prefix from the chromosome, convert to 0-based half-open if needed:

```
tark_search_transcripts_by_region(
    region="13",
    start=32315479,   # 0-based
    end=32400268,     # exclusive
    assembly="GRCh38"
)
```

Results include all overlapping transcripts with full exon structure.

---

### 4. Compare transcript versions across Ensembl releases

**Use case:** "What changed between ENST00000380152 version 6 and version 7?"

```
tark_diff_transcripts(
    stable_ids=["ENST00000380152.6", "ENST00000380152.7"]
)
```

First ID is the reference. Check the result fields:
- `exon_diffs` — per-exon change classification (`added`, `removed`, `modified`, `unchanged`)
- `cds_sequence_changed` — whether the coding sequence changed
- `protein_sequence_changed` — whether the amino acid sequence changed
- `biotype_changed`, `exon_count_changed` — structural flags

---

### 5. Cross-assembly comparison (GRCh37 vs GRCh38)

**Use case:** "How does ENST00000380152 differ between builds?"

First confirm the transcript exists in both builds:

```
tark_get_transcript(stable_id="ENST00000380152", assembly="both")
```

Then diff across assemblies — pass the same base ID but override assemblies:

```
tark_diff_transcripts(
    stable_ids=["ENST00000380152", "ENST00000380152"],
    assemblies=["GRCh37", "GRCh38"]
)
```

Interpret `exon_diffs` and coordinate fields relative to each build's reference sequence.

---

### 6. CDS / protein-focused query

**Use case:** "What is the protein sequence for ENST00000380152?"

`tark_get_transcript` returns `cds_start`, `cds_end`, and the full `sequence`. Slice the CDS:

```python
# cds_start and cds_end are 0-based offsets into the transcript sequence
cds_seq = transcript["sequence"][transcript["cds_start"]:transcript["cds_end"]]
```

To get the translation stable ID:

```
tark_get_protein_for_transcript(stable_id="ENST00000380152")
```

For the actual amino acid sequence, use `tark_diff_transcripts` (it fetches protein sequences internally and returns them in `ref_protein_sequence` / `candidate_protein_sequence`). Non-coding transcripts return `None` for these fields — always check `reference_protein_coding` before interpreting.

---

## Common Pitfalls

- **`tark_get_mane_transcripts` is slow for bulk use.** It fetches all MANE entries then filters client-side. Fine for a single gene; avoid calling it in a loop across many genes — use `tark_get_gene_transcripts` per gene instead.
- **Strip `chr` prefix for region searches.** `tark_search_transcripts_by_region` requires a numeric/plain chromosome string. `"chr13"` will return no results; use `"13"`.
- **`tark_diff_transcripts` passes IDs directly to the TARK `/diff/` endpoint without stripping versions.** If results are unexpected, try passing base stable IDs (without `.version` suffix).
- **Non-coding transcripts return `None` for CDS/protein fields in diffs.** Always check `reference_protein_coding` and `candidate_protein_coding` before interpreting `cds_sequence_changed` or `protein_sequence_changed`.
- **`assembly="both"` returns a flat list.** Check the `assembly` field on each record to identify which build it came from.
````

- [ ] **Step 3: Verify the file was created**

Run: `head -5 plugin/skills/tark-mcp/SKILL.md`

Expected output:
```
---
name: tark-mcp
description: Use when querying Ensembl transcript data via TARK ...
---
```

- [ ] **Step 4: Commit**

```bash
git add plugin/skills/tark-mcp/SKILL.md
git commit -m "feat: add tark-mcp workflow skill"
```

---

### Task 3: Install the plugin locally

**Files:**
- No repo files modified — creates a symlink in `~/.copilot/installed-plugins/`

- [ ] **Step 1: Symlink the plugin directory**

Run from the repo root:

```bash
ln -s "$(pwd)/plugin" ~/.copilot/installed-plugins/tark-mcp
```

- [ ] **Step 2: Verify the symlink was created**

```bash
ls -la ~/.copilot/installed-plugins/tark-mcp
```

Expected: a symlink pointing to the repo's `plugin/` directory.

- [ ] **Step 3: Verify the skill is discoverable**

```bash
cat ~/.copilot/installed-plugins/tark-mcp/plugin.json
cat ~/.copilot/installed-plugins/tark-mcp/skills/tark-mcp/SKILL.md | head -10
```

Expected: both files readable through the symlink.

---

### Task 4: Create README with setup instructions

**Files:**
- Create: `README.md`

- [ ] **Step 1: Create README.md**

Create `README.md` with this content:

```markdown
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
```

- [ ] **Step 2: Verify the file renders correctly**

```bash
head -5 README.md
```

Expected: `# tark-mcp`

- [ ] **Step 3: Commit**

```bash
git add README.md
git commit -m "docs: add README with installation and usage instructions"
```

---

## Done

At this point:
- `plugin/plugin.json` and `plugin/skills/tark-mcp/SKILL.md` exist in the repo
- The plugin is symlinked to `~/.copilot/installed-plugins/tark-mcp/`
- `README.md` documents installation
- All changes are committed

Copilot CLI will discover and load the `tark-mcp` skill automatically on next invocation.
