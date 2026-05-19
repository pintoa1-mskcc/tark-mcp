from __future__ import annotations
from tark_mcp.client import TarkClient
from tark_mcp.models import Exon, Translation
from tark_mcp.tools.transcripts import get_transcript


async def get_transcript_sequence(
    stable_id: str,
    assembly: str = "GRCh38",
    client: TarkClient | None = None,
) -> dict | list[dict] | None:
    if client is None:
        client = TarkClient()
    t = await get_transcript(stable_id, assembly=assembly, client=client)
    if t is None:
        return None
    if isinstance(t, list):
        return [
            {"stable_id": x.stable_id, "stable_id_version": x.stable_id_version,
             "assembly": x.assembly, "sequence": x.sequence}
            for x in t
        ]
    return {"stable_id": t.stable_id, "stable_id_version": t.stable_id_version,
            "assembly": t.assembly, "sequence": t.sequence}


async def get_transcript_exons(
    stable_id: str,
    assembly: str = "GRCh38",
    client: TarkClient | None = None,
) -> list[Exon]:
    if client is None:
        client = TarkClient()
    t = await get_transcript(stable_id, assembly=assembly, client=client)
    if t is None:
        return []
    transcripts = t if isinstance(t, list) else [t]
    exons: list[Exon] = []
    for transcript in transcripts:
        ordered = sorted(transcript.exons, key=lambda e: e.order)
        if transcript.loc_strand == -1:
            ordered = list(reversed(ordered))
        exons.extend(ordered)
    return exons


async def get_protein_for_transcript(
    stable_id: str,
    assembly: str = "GRCh38",
    client: TarkClient | None = None,
) -> Translation | list[Translation] | None:
    if client is None:
        client = TarkClient()
    t = await get_transcript(stable_id, assembly=assembly, client=client)
    if t is None:
        return None
    if isinstance(t, list):
        results = [x.translations[0] if x.translations else None for x in t]
        return results
    return t.translations[0] if t.translations else None
