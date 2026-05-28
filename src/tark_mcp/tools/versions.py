"""Tool for querying all versions of a transcript."""
from __future__ import annotations
import asyncio
from tark_mcp.client import TarkClient
from tark_mcp.models import Transcript


async def get_transcript_all_versions(
    stable_id: str,
    assembly: str = "GRCh38",
    client: TarkClient | None = None,
) -> list[Transcript]:
    """Get all versions of a transcript without deduplication.
    
    Unlike get_transcript which returns only the latest version per assembly,
    this function returns ALL versions that have ever existed in TARK.
    
    Args:
        stable_id: Ensembl transcript stable ID (without version suffix)
        assembly: 'GRCh37', 'GRCh38' (default), or 'both'
        client: Optional TarkClient instance
    
    Returns:
        List of all Transcript versions, sorted by stable_id_version
    
    Example:
        # Get all versions of ENST00000263967 in GRCh38
        versions = await get_transcript_all_versions("ENST00000263967", "GRCh38")
        
        # Get versions from both assemblies
        versions = await get_transcript_all_versions("ENST00000263967", "both")
    """
    if client is None:
        client = TarkClient()
    
    # Strip version suffix if provided
    if "." in stable_id:
        stable_id = stable_id.rsplit(".", 1)[0]
    
    if assembly == "both":
        # Fetch from both assemblies in parallel
        results = await asyncio.gather(
            _fetch_all_versions(stable_id, "GRCh38", client),
            _fetch_all_versions(stable_id, "GRCh37", client),
        )
        combined = results[0] + results[1]
        # Sort by assembly, then version
        combined.sort(key=lambda t: (t.assembly, t.stable_id_version))
        return combined
    else:
        transcripts = await _fetch_all_versions(stable_id, assembly, client)
        transcripts.sort(key=lambda t: t.stable_id_version)
        return transcripts


async def _fetch_all_versions(
    stable_id: str,
    assembly: str,
    client: TarkClient,
) -> list[Transcript]:
    """Fetch all versions for a single assembly without deduplication."""
    data = await client.get(
        "transcript/",
        params={
            "stable_id": stable_id,
            "assembly_name": assembly,
            "expand_all": "true",
        },
    )
    
    # Parse all transcripts WITHOUT deduplication
    # (the standard tools deduplicate to keep only latest per release)
    transcripts = [Transcript.model_validate(r) for r in data]
    return transcripts


async def get_transcript_version_history(
    stable_id: str,
    assembly: str = "GRCh38",
    client: TarkClient | None = None,
) -> dict:
    """Get version history with summary statistics.
    
    Args:
        stable_id: Ensembl transcript stable ID
        assembly: 'GRCh37', 'GRCh38' (default), or 'both'
        client: Optional TarkClient instance
    
    Returns:
        Dictionary with version history and comparison data:
        {
            "stable_id": str,
            "assembly": str,
            "version_count": int,
            "versions": [
                {
                    "version": int,
                    "release_date": str,
                    "release_version": str | None,
                    "location": str,
                    "genomic_span": int,
                    "transcript_length": int,
                    "exon_count": int,
                    "cds_length": int | None,
                    "protein_id": str | None,
                    "protein_length": int | None,
                },
                ...
            ],
            "changes": [
                {
                    "from_version": int,
                    "to_version": int,
                    "genomic_span_delta": int,
                    "transcript_length_delta": int,
                    "exon_count_delta": int,
                    "cds_length_delta": int | None,
                },
                ...
            ]
        }
    """
    if client is None:
        client = TarkClient()
    
    versions = await get_transcript_all_versions(stable_id, assembly, client)
    
    if not versions:
        return {
            "stable_id": stable_id,
            "assembly": assembly,
            "version_count": 0,
            "versions": [],
            "changes": [],
        }
    
    # Build version summaries
    version_summaries = []
    for t in versions:
        genomic_span = t.loc_end - t.loc_start
        cds_length = (t.cds_end - t.cds_start) if t.cds_start is not None else None
        protein_length = (cds_length // 3) if cds_length is not None else None
        
        version_summaries.append({
            "version": t.stable_id_version,
            "release_date": t.latest_release_date,
            "release_version": t.latest_release_version,
            "location": f"{t.loc_region}:{t.loc_start}-{t.loc_end}",
            "genomic_span": genomic_span,
            "transcript_length": len(t.sequence) if t.sequence else None,
            "exon_count": len(t.exons),
            "cds_length": cds_length,
            "protein_id": f"{t.translations[0].stable_id}.{t.translations[0].stable_id_version}" if t.translations else None,
            "protein_length": protein_length,
        })
    
    # Build change deltas between consecutive versions
    changes = []
    for i in range(len(versions) - 1):
        prev = version_summaries[i]
        curr = version_summaries[i + 1]
        
        changes.append({
            "from_version": prev["version"],
            "to_version": curr["version"],
            "genomic_span_delta": curr["genomic_span"] - prev["genomic_span"],
            "transcript_length_delta": (curr["transcript_length"] - prev["transcript_length"]) 
                                        if curr["transcript_length"] and prev["transcript_length"] else None,
            "exon_count_delta": curr["exon_count"] - prev["exon_count"],
            "cds_length_delta": (curr["cds_length"] - prev["cds_length"]) 
                                if curr["cds_length"] is not None and prev["cds_length"] is not None else None,
        })
    
    return {
        "stable_id": versions[0].stable_id,
        "assembly": assembly,
        "version_count": len(versions),
        "versions": version_summaries,
        "changes": changes,
    }
