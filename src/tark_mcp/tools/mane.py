from __future__ import annotations
from tark_mcp.client import TarkClient
from tark_mcp.models import ManeTranscript


async def get_mane_transcripts(
    gene_identifier: str | None = None,
    client: TarkClient | None = None,
) -> list[ManeTranscript]:
    if client is None:
        client = TarkClient()
    data = await client.get("transcript/manelist/")
    entries = [ManeTranscript.model_validate(r) for r in data]
    if gene_identifier is None:
        return entries
    # The endpoint carries gene symbols only (no Ensembl gene IDs), so filter by symbol.
    return [
        e for e in entries
        if e.gene_name and e.gene_name.upper() == gene_identifier.upper()
    ]
