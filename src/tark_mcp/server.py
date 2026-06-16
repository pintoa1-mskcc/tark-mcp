from __future__ import annotations
import asyncio
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
from tark_mcp.tools.formatters import format_transcripts_table

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
async def tark_get_transcripts(
    stable_ids: list[str],
    assemblies: list[str] | None = None,
) -> str:
    """Retrieve multiple transcripts in a single call and return a formatted summary table.

    Columns: Query, Assembly, Stable ID, Ver, Exons, 5'UTR, 3'UTR, CDS (bp),
             AA Len, First Release, Latest Release, Release Date, MANE.

    Args:
        stable_ids: List of transcript stable IDs (Ensembl or RefSeq). Version suffixes supported,
            e.g. ['ENST00000380152.7', 'NM_001128425.2']
        assemblies: Optional per-entry assembly override list ('GRCh37' or 'GRCh38').
            Defaults to 'GRCh38' for any missing entries.
    """
    resolved = list(assemblies or [])
    while len(resolved) < len(stable_ids):
        resolved.append("GRCh38")

    transcript_results, mane_raw = await asyncio.gather(
        asyncio.gather(*[
            get_transcript(sid, assembly=asm, client=_client)
            for sid, asm in zip(stable_ids, resolved)
        ]),
        _client.get("transcript/manelist/"),
    )

    # Build lookup: stable_id (no version) → normalised MANE type label
    mane_lookup: dict[str, str] = {}
    for entry in mane_raw:
        raw_type = (entry.get("mane_type") or "").upper()
        if "PLUS CLINICAL" in raw_type:
            label = "MANE Plus Clinical"
        elif "SELECT" in raw_type:
            label = "MANE Select"
        else:
            label = raw_type
        for key in ("ens_stable_id", "refseq_stable_id"):
            sid = entry.get(key)
            if sid:
                mane_lookup[sid] = label

    dicts: list[dict | list[dict] | None] = []
    for result in transcript_results:
        if result is None:
            dicts.append(None)
        elif isinstance(result, list):
            dicts.append([t.model_dump() for t in result])
        else:
            dicts.append(result.model_dump())

    return format_transcripts_table(stable_ids, resolved, dicts, mane_lookup=mane_lookup)


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
        stable_ids: List of ≥2 transcript stable IDs (Ensembl or RefSeq); first is the
            reference. Version suffixes are supported (e.g. 'ENST00000380152.7',
            'NM_001128425.2').
        assemblies: Optional per-entry assembly override list; defaults to 'GRCh38' for missing entries
    """
    results = await diff_transcripts(stable_ids, assemblies=assemblies, client=_client)
    return [d.model_dump() for d in results]


def main() -> None:
    mcp.run()


if __name__ == "__main__":
    main()
