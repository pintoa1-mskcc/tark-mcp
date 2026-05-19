from __future__ import annotations
from mcp.server.fastmcp import FastMCP

from tark_mcp.client import TarkClient
from tark_mcp.tools.releases import get_releases
from tark_mcp.tools.transcripts import get_transcript, search_transcripts_by_region
from tark_mcp.tools.genes import get_gene_transcripts
from tark_mcp.tools.sequences import (
    get_transcript_sequence, get_transcript_exons, get_protein_for_transcript
)
from tark_mcp.tools.mane import get_mane_transcripts
from tark_mcp.tools.diff import diff_transcripts

mcp = FastMCP("tark")
_client = TarkClient()


@mcp.tool()
async def tark_get_releases() -> list[dict]:
    """List all available TARK releases with metadata (short name, date, assembly, source)."""
    releases = await get_releases(_client)
    return [r.model_dump() for r in releases]


@mcp.tool()
async def tark_get_transcript(stable_id: str, assembly: str = "GRCh38") -> dict | list[dict] | None:
    """Retrieve a transcript by Ensembl stable ID with full exon structure, CDS boundaries, and genes.

    Args:
        stable_id: Ensembl transcript stable ID, e.g. 'ENST00000380152' or 'ENST00000380152.7'
        assembly: Genome build — 'GRCh37', 'GRCh38' (default), or 'both'
    """
    result = await get_transcript(stable_id, assembly=assembly, client=_client)
    if result is None:
        return None
    if isinstance(result, list):
        return [t.model_dump() for t in result]
    return result.model_dump()


@mcp.tool()
async def tark_search_transcripts_by_region(
    region: str, start: int, end: int, assembly: str = "GRCh38"
) -> list[dict]:
    """Find all transcripts overlapping a genomic region (0-based half-open coordinates).

    Args:
        region: Chromosome, e.g. '13' or 'chr13'
        start: 0-based start position (inclusive)
        end: 0-based end position (exclusive)
        assembly: 'GRCh37', 'GRCh38' (default), or 'both'
    """
    results = await search_transcripts_by_region(region, start, end, assembly=assembly, client=_client)
    return [t.model_dump() for t in results]


@mcp.tool()
async def tark_get_gene_transcripts(gene_identifier: str, assembly: str = "GRCh38") -> list[dict]:
    """Get all transcripts for a gene symbol or Ensembl gene ID.

    Args:
        gene_identifier: Gene symbol (e.g. 'BRCA2') or Ensembl gene ID (e.g. 'ENSG00000139618')
        assembly: 'GRCh37', 'GRCh38' (default), or 'both'
    """
    results = await get_gene_transcripts(gene_identifier, assembly=assembly, client=_client)
    return [t.model_dump() for t in results]


@mcp.tool()
async def tark_get_transcript_sequence(
    stable_id: str, assembly: str = "GRCh38"
) -> dict | list[dict] | None:
    """Fetch the cDNA sequence for a transcript.

    Args:
        stable_id: Ensembl transcript stable ID
        assembly: 'GRCh37', 'GRCh38' (default), or 'both'
    """
    return await get_transcript_sequence(stable_id, assembly=assembly, client=_client)


@mcp.tool()
async def tark_get_transcript_exons(
    stable_id: str, assembly: str = "GRCh38"
) -> list[dict]:
    """Return the ordered exon list for a transcript with 0-based genomic coordinates.

    Args:
        stable_id: Ensembl transcript stable ID
        assembly: 'GRCh37', 'GRCh38' (default), or 'both'
    """
    exons = await get_transcript_exons(stable_id, assembly=assembly, client=_client)
    return [e.model_dump() for e in exons]


@mcp.tool()
async def tark_get_protein_for_transcript(
    stable_id: str, assembly: str = "GRCh38"
) -> dict | list[dict | None] | None:
    """Return the protein (translation) stable ID and version for a transcript.

    Args:
        stable_id: Ensembl transcript stable ID
        assembly: 'GRCh37', 'GRCh38' (default), or 'both'
    """
    result = await get_protein_for_transcript(stable_id, assembly=assembly, client=_client)
    if result is None:
        return None
    if isinstance(result, list):
        return [t.model_dump() if t else None for t in result]
    return result.model_dump()


@mcp.tool()
async def tark_get_mane_transcripts(gene_identifier: str | None = None) -> list[dict]:
    """Return MANE Select and MANE Plus Clinical transcripts, optionally filtered by gene.

    Args:
        gene_identifier: Optional gene symbol or Ensembl gene ID to filter results
    """
    results = await get_mane_transcripts(gene_identifier=gene_identifier, client=_client)
    return [t.model_dump() for t in results]


@mcp.tool()
async def tark_diff_transcripts(
    stable_ids: list[str],
    assemblies: list[str] | None = None,
) -> list[dict]:
    """Compare transcripts against a reference. First ID is reference; all subsequent are candidates.

    Returns structural diff (exon-level, CDS boundaries, biotype) plus CDS nucleotide sequence
    comparison and protein (amino acid) sequence comparison. Non-coding transcripts produce
    None for sequence comparison fields.

    Args:
        stable_ids: List of ≥2 Ensembl transcript stable IDs; first is the reference
        assemblies: Optional per-entry assembly override list; defaults to 'GRCh38' for missing entries
    """
    results = await diff_transcripts(stable_ids, assemblies=assemblies, client=_client)
    return [d.model_dump() for d in results]


def main() -> None:
    mcp.run()


if __name__ == "__main__":
    main()
