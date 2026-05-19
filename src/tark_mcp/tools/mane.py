from __future__ import annotations
from tark_mcp.client import TarkClient
from tark_mcp.models import Transcript


async def get_mane_transcripts(
    gene_identifier: str | None = None,
    client: TarkClient | None = None,
) -> list[Transcript]:
    if client is None:
        client = TarkClient()
    data = await client.get("transcript/manelist/")
    transcripts = [Transcript.model_validate(r) for r in data]
    if gene_identifier is None:
        return transcripts
    # Filter client-side by gene name or stable ID
    filtered = []
    for t in transcripts:
        for gene in t.genes:
            if (gene.name and gene.name.upper() == gene_identifier.upper()) or \
               gene.stable_id == gene_identifier:
                filtered.append(t)
                break
    return filtered
