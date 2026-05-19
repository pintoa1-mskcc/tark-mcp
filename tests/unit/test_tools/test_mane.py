import pytest
import httpx
import respx

from tark_mcp.client import TarkClient
from tark_mcp.tools.mane import get_mane_transcripts
from tests.conftest import TRANSCRIPT_BRCA2_RAW, MANE_LIST_RESPONSE_RAW

BASE = "https://tark.ensembl.org/api/"


@respx.mock
@pytest.mark.asyncio
async def test_get_mane_transcripts_returns_all():
    client = TarkClient()
    respx.get(BASE + "transcript/manelist/").mock(
        return_value=httpx.Response(200, json=MANE_LIST_RESPONSE_RAW)
    )
    results = await get_mane_transcripts(client=client)
    assert len(results) == 1
    assert results[0].stable_id == "ENST00000380152"


@respx.mock
@pytest.mark.asyncio
async def test_get_mane_transcripts_filters_by_gene_name():
    client = TarkClient()
    two_genes = {
        "count": 2, "next": None, "previous": None,
        "results": [
            TRANSCRIPT_BRCA2_RAW,
            {**TRANSCRIPT_BRCA2_RAW, "stable_id": "ENST00000999999",
             "genes": [{"stable_id": "ENSG00000012048", "stable_id_version": 1,
                        "assembly": "GRCh38", "loc_start": 100, "loc_end": 200,
                        "loc_strand": 1, "loc_region": "1", "name": "BRCA1"}]},
        ],
    }
    respx.get(BASE + "transcript/manelist/").mock(
        return_value=httpx.Response(200, json=two_genes)
    )
    results = await get_mane_transcripts(gene_identifier="BRCA2", client=client)
    assert len(results) == 1
    assert results[0].stable_id == "ENST00000380152"


@respx.mock
@pytest.mark.asyncio
async def test_get_mane_transcripts_filters_by_gene_stable_id():
    client = TarkClient()
    respx.get(BASE + "transcript/manelist/").mock(
        return_value=httpx.Response(200, json=MANE_LIST_RESPONSE_RAW)
    )
    results = await get_mane_transcripts(gene_identifier="ENSG00000139618", client=client)
    assert len(results) == 1
