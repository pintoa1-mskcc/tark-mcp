import pytest
import httpx
import respx

from tark_mcp.client import TarkClient
from tark_mcp.tools.genes import get_gene_transcripts
from tests.conftest import TRANSCRIPT_BRCA2_RAW, TRANSCRIPT_NONCODING_RAW

BASE = "https://tark.ensembl.org/api/"

MIXED_ASSEMBLY_RESULTS = {
    "count": 2, "next": None, "previous": None,
    "results": [
        TRANSCRIPT_BRCA2_RAW,
        {**TRANSCRIPT_NONCODING_RAW,
         "assembly": {"assembly_name": "GRCh37", "assembly_id": 2, "genome": 2, "session": 1}},
    ],
}


@respx.mock
@pytest.mark.asyncio
async def test_get_gene_transcripts_filters_by_assembly():
    """Search endpoint returns mixed assemblies; client-side filters to GRCh38."""
    client = TarkClient()
    respx.get(BASE + "transcript/search/").mock(
        return_value=httpx.Response(200, json=MIXED_ASSEMBLY_RESULTS)
    )
    results = await get_gene_transcripts("BRCA2", assembly="GRCh38", client=client)
    assert all(t.assembly == "GRCh38" for t in results)
    assert len(results) == 1


@respx.mock
@pytest.mark.asyncio
async def test_get_gene_transcripts_assembly_both_returns_all():
    client = TarkClient()
    respx.get(BASE + "transcript/search/").mock(
        return_value=httpx.Response(200, json=MIXED_ASSEMBLY_RESULTS)
    )
    results = await get_gene_transcripts("BRCA2", assembly="both", client=client)
    assert len(results) == 2


@respx.mock
@pytest.mark.asyncio
async def test_get_gene_transcripts_passes_identifier_in_query():
    client = TarkClient()
    route = respx.get(BASE + "transcript/search/").mock(
        return_value=httpx.Response(200, json={"count": 0, "next": None, "previous": None,
                                               "results": []})
    )
    await get_gene_transcripts("ENSG00000139618", client=client)
    assert "identifier_field=ENSG00000139618" in str(route.calls[0].request.url)
