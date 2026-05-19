from __future__ import annotations
from tark_mcp.client import TarkClient
from tark_mcp.models import Transcript
from tark_mcp.tools.transcripts import _deduplicate


async def get_gene_transcripts(
    gene_identifier: str,
    assembly: str = "GRCh38",
    client: TarkClient | None = None,
) -> list[Transcript]:
    if client is None:
        client = TarkClient()
    params = {
        "identifier_field": gene_identifier,
        "expand": "exons,genes,sequence",
    }
    data = await client.get("transcript/search/", params)
    transcripts = [Transcript.model_validate(r) for r in data]
    if assembly != "both":
        transcripts = [t for t in transcripts if t.assembly == assembly]
    return _deduplicate(transcripts)
