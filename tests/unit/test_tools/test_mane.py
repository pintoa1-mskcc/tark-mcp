import pytest
import httpx
import respx

from tark_mcp.client import TarkClient
from tark_mcp.tools.mane import get_mane_transcripts
from tests.conftest import MANE_LIST_RESPONSE_RAW

BASE = "https://tark.ensembl.org/api/"


@respx.mock
@pytest.mark.asyncio
async def test_get_mane_transcripts_returns_all():
    client = TarkClient()
    respx.get(BASE + "transcript/manelist/").mock(
        return_value=httpx.Response(200, json=MANE_LIST_RESPONSE_RAW)
    )
    results = await get_mane_transcripts(client=client)
    assert len(results) == 4


@respx.mock
@pytest.mark.asyncio
async def test_get_mane_transcripts_parses_raw_entry():
    client = TarkClient()
    respx.get(BASE + "transcript/manelist/").mock(
        return_value=httpx.Response(200, json=MANE_LIST_RESPONSE_RAW)
    )
    results = await get_mane_transcripts(gene_identifier="DAXX", client=client)
    assert len(results) == 1
    m = results[0]
    assert m.ensembl_id == "ENST00000374542.10"
    assert m.refseq_id == "NM_001141969.2"
    assert m.mane_type == "MANE Select"
    assert m.gene_name == "DAXX"


@respx.mock
@pytest.mark.asyncio
async def test_get_mane_transcripts_filters_by_gene_name_case_insensitive():
    client = TarkClient()
    respx.get(BASE + "transcript/manelist/").mock(
        return_value=httpx.Response(200, json=MANE_LIST_RESPONSE_RAW)
    )
    results = await get_mane_transcripts(gene_identifier="abcc8", client=client)
    assert {m.mane_type for m in results} == {"MANE Select", "MANE Plus Clinical"}


@respx.mock
@pytest.mark.asyncio
async def test_get_mane_transcripts_unknown_gene_returns_empty():
    client = TarkClient()
    respx.get(BASE + "transcript/manelist/").mock(
        return_value=httpx.Response(200, json=MANE_LIST_RESPONSE_RAW)
    )
    results = await get_mane_transcripts(gene_identifier="NOTAGENE", client=client)
    assert results == []
