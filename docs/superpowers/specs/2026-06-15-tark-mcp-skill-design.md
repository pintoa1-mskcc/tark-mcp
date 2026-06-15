# tark-mcp Skill Design

**Date:** 2026-06-15  
**Status:** Approved  
**Audience:** Bioinformaticians / genomics researchers  

---

## Problem

The tark-mcp server exposes 9 Ensembl transcript tools via MCP. AI agents working in genomics contexts have no structured guidance on which tool to use for which task, how to chain tools together, or what pitfalls to avoid. Without a skill, agents default to guessing or querying inefficiently.

---

## Goal

Create a Copilot CLI skill that teaches agents to work with TARK data task-first — organized around what bioinformaticians actually do, not around API endpoints.

---

## Architecture

Two new files added to the tark-mcp repo:

```
plugin/
  plugin.json            ← Copilot CLI plugin manifest
  skills/
    tark-mcp/
      SKILL.md           ← skill content
```

The `plugin.json` follows the Copilot CLI plugin manifest format:
- `name`: `tark-mcp`
- `skills`: `"./skills/"`

Installation: symlink or copy `plugin/` into `~/.copilot/installed-plugins/tark-mcp/`. A note in `pyproject.toml` or `README` documents the install command.

The skill is also committed to the repo so it travels with the MCP server and can be installed by anyone who clones the repo.

---

## Skill Content Structure

### Frontmatter
- `name: tark-mcp`
- `description`: trigger phrase covering transcripts, genes, Ensembl IDs, MANE, assemblies, CDS, protein, exons, genomic locus queries

### Quick Tool Index
One-line summary of all 9 tools for fast scanning:
- `tark_get_releases` — list TARK releases with assembly/date metadata
- `tark_get_transcript` — fetch single transcript by Ensembl stable ID
- `tark_search_transcripts_by_region` — find transcripts overlapping a locus (0-based half-open)
- `tark_get_gene_transcripts` — all transcripts for a gene symbol or Ensembl gene ID
- `tark_get_transcript_sequence` — cDNA sequence for a transcript
- `tark_get_transcript_exons` — ordered exon list with genomic coordinates
- `tark_get_protein_for_transcript` — protein/translation stable ID and version
- `tark_get_mane_transcripts` — MANE Select and Plus Clinical transcripts, optionally filtered by gene
- `tark_diff_transcripts` — structural + CDS + protein diff; first ID is reference

### Key Concepts
Domain-assumed, kept brief:
- Coordinates returned are **0-based half-open** (TARK API is 1-based; conversion is handled internally)
- `assembly: "both"` returns results for GRCh37 and GRCh38 combined in one list; each record carries its assembly field
- `stable_id` accepts versioned IDs (`ENST00000380152.7`); version is used as a filter if provided
- `diff_transcripts`: the **first ID is always the reference**; all others are candidates compared against it
- `tark_get_transcript` returns the single most-recently-released record by default (deduplication applied)

### Workflow Recipes

Six task-oriented workflows, each showing which tools to call and in what order:

1. **Find the canonical transcript for a gene**  
   Use `tark_get_mane_transcripts(gene_identifier=...)` for MANE Select/Plus Clinical. If MANE is not required, use `tark_get_gene_transcripts` and filter by biotype.

2. **Look up a specific transcript**  
   Call `tark_get_transcript(stable_id, assembly)`. Chain `tark_get_transcript_sequence` for cDNA, `tark_get_transcript_exons` for exon structure, `tark_get_protein_for_transcript` for the translation ID.

3. **Investigate a genomic locus**  
   Call `tark_search_transcripts_by_region(region, start, end)`. Strip `chr` prefix from chromosome. Coordinates are 0-based half-open.

4. **Compare transcript versions across Ensembl releases**  
   Call `tark_diff_transcripts(["ENST00000XXXXXX.6", "ENST00000XXXXXX.7"])`. First ID is reference. Check `exon_diffs`, `cds_sequence_changed`, `protein_sequence_changed` in results.

5. **Cross-assembly comparison (GRCh37 vs GRCh38)**  
   Call `tark_get_transcript(stable_id, assembly="both")` to confirm the transcript exists in both builds. Then call `tark_diff_transcripts(["ENST...", "ENST..."], assemblies=["GRCh37", "GRCh38"])` for a structural diff across builds.

6. **CDS / protein-focused query**  
   Call `tark_get_transcript` (includes `cds_start`, `cds_end`, `sequence`). Then `tark_get_protein_for_transcript` for the translation stable ID. For the actual amino acid sequence use `tark_diff_transcripts` which fetches protein sequences internally, or query the TARK API translation endpoint directly via `TarkClient`.

### Common Pitfalls
- `tark_get_mane_transcripts` fetches **all** MANE entries then filters client-side — acceptable for single-gene queries, avoid in loops
- `tark_search_transcripts_by_region` requires a numeric chromosome string — strip the `chr` prefix (`"chr13"` → `"13"`)
- `tark_diff_transcripts` passes IDs directly to the TARK `/diff/` endpoint without stripping versions; if results are unexpected, try passing base stable IDs (without `.version` suffix)
- Non-coding transcripts return `None` for CDS/protein fields in diff results — always check `reference_protein_coding` / `candidate_protein_coding` before interpreting those fields

---

## Installation

```bash
# From the tark-mcp repo root
ln -s "$(pwd)/plugin" ~/.copilot/installed-plugins/tark-mcp
```

Or copy instead of symlink for a stable install:
```bash
cp -r plugin ~/.copilot/installed-plugins/tark-mcp
```

---

## Out of Scope

- Custom query functions in `versions.py` (not exposed as MCP tools) — documented separately in `CUSTOM_QUERIES.md`
- Adding new MCP tools
- CI/CD for the plugin
