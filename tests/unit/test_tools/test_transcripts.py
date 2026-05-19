import pytest
import httpx
import respx

from tark_mcp.client import TarkClient
from tark_mcp import client as client_module
from tark_mcp.tools.transcripts import get_transcript, search_transcripts_by_region
from tests.conftest import TRANSCRIPT_BRCA2_RAW, TRANSCRIPT_NONCODING_RAW

BASE = "https://tark.ensembl.org/api/"

PAGINATED_TWO = {
    "count": 2,
    "next": None,
    "previous": None,
    "results": [
        TRANSCRIPT_BRCA2_RAW,
        {**TRANSCRIPT_BRCA2_RAW, "stable_id_version": 6,
         "transcript_release_set": [{"assembly": "GRCh38", "shortname": "109",
                                     "description": "Ensembl release",
                                     "release_date": "2022-09-01", "source": "Ensembl"}]},
    ],
}


@pytest.fixture
def client():
    client_module._cache.clear()
    return TarkClient()


@respx.mock
@pytest.mark.asyncio
async def test_get_transcript_returns_latest_version(client):
    """When multiple versions exist, deduplicate and return the latest."""
    respx.get(BASE + "transcript/").mock(
        return_value=httpx.Response(200, json=PAGINATED_TWO)
    )
    result = await get_transcript("ENST00000380152", client=client)
    assert result is not None
    assert result.stable_id == "ENST00000380152"
    # Should return version 7 (latest release_date 2023-04-01)
    assert result.stable_id_version == 7


@respx.mock
@pytest.mark.asyncio
async def test_get_transcript_strips_version_suffix(client):
    """ENST00000380152.7 → queries for stable_id=ENST00000380152."""
    route = respx.get(BASE + "transcript/").mock(
        return_value=httpx.Response(200, json={
            "count": 1, "next": None, "previous": None,
            "results": [TRANSCRIPT_BRCA2_RAW]
        })
    )
    await get_transcript("ENST00000380152.7", client=client)
    assert "stable_id=ENST00000380152" in str(route.calls[0].request.url)


@respx.mock
@pytest.mark.asyncio
async def test_get_transcript_not_found_returns_none(client):
    respx.get(BASE + "transcript/").mock(
        return_value=httpx.Response(200, json={"count": 0, "next": None,
                                               "previous": None, "results": []})
    )
    result = await get_transcript("ENST00000999999", client=client)
    assert result is None


@respx.mock
@pytest.mark.asyncio
async def test_get_transcript_assembly_both_returns_list(client):
    """assembly='both' fans out to GRCh37 and GRCh38, returns list."""
    t38 = {**TRANSCRIPT_BRCA2_RAW}
    t37 = {**TRANSCRIPT_BRCA2_RAW,
           "assembly": {"assembly_name": "GRCh37", "assembly_id": 2, "genome": 2, "session": 1}}
    respx.get(BASE + "transcript/").mock(side_effect=[
        httpx.Response(200, json={"count": 1, "next": None, "previous": None, "results": [t38]}),
        httpx.Response(200, json={"count": 1, "next": None, "previous": None, "results": [t37]}),
    ])
    result = await get_transcript("ENST00000380152", assembly="both", client=client)
    assert isinstance(result, list)
    assert len(result) == 2


@respx.mock
@pytest.mark.asyncio
async def test_search_transcripts_by_region(client):
    """0-based input is converted to 1-based for the API; chr prefix is stripped."""
    route = respx.get(BASE + "transcript/").mock(
        return_value=httpx.Response(200, json={
            "count": 1, "next": None, "previous": None,
            "results": [TRANSCRIPT_BRCA2_RAW]
        })
    )
    results = await search_transcripts_by_region(
        "chr13", start=32315474, end=32400266, client=client
    )
    assert len(results) == 1
    url = str(route.calls[0].request.url)
    # 0-based 32315474 → 1-based 32315475 in query
    assert "loc_start=32315475" in url
    assert "loc_region=13" in url   # chr prefix stripped
