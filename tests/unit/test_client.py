import pytest
import httpx
import respx

from tark_mcp.client import TarkClient
from tark_mcp import client as client_module

BASE = "https://tark.ensembl.org/api/"


@pytest.fixture
def client():
    client_module._cache.clear()
    return TarkClient()


@respx.mock
@pytest.mark.asyncio
async def test_get_single_page(client):
    """Fetches a single page and returns results list."""
    respx.get(BASE + "release/nopagination/").mock(
        return_value=httpx.Response(200, json=[{"shortname": "110"}])
    )
    result = await client.get("release/nopagination/")
    assert result == [{"shortname": "110"}]


@respx.mock
@pytest.mark.asyncio
async def test_get_paginates_automatically(client):
    """Follows next links and aggregates all pages."""
    page1 = {"count": 2, "next": BASE + "transcript/?page=2", "previous": None,
             "results": [{"stable_id": "ENST000001"}]}
    page2 = {"count": 2, "next": None, "previous": BASE + "transcript/",
             "results": [{"stable_id": "ENST000002"}]}
    # Register specific route first so respx matches it before the general route
    respx.get(BASE + "transcript/?page=2").mock(return_value=httpx.Response(200, json=page2))
    respx.get(BASE + "transcript/").mock(return_value=httpx.Response(200, json=page1))

    results = await client.get("transcript/")
    assert len(results) == 2
    assert results[0]["stable_id"] == "ENST000001"
    assert results[1]["stable_id"] == "ENST000002"


@respx.mock
@pytest.mark.asyncio
async def test_get_rewrites_http_to_https(client):
    """HTTPS is enforced even if next link is http://."""
    page1 = {"count": 2, "next": "http://tark.ensembl.org/api/transcript/?page=2",
             "previous": None, "results": [{"stable_id": "ENST000001"}]}
    page2 = {"count": 2, "next": None, "previous": None,
             "results": [{"stable_id": "ENST000002"}]}
    # Register specific route first so respx matches it before the general route
    respx.get(BASE + "transcript/?page=2").mock(return_value=httpx.Response(200, json=page2))
    respx.get(BASE + "transcript/").mock(return_value=httpx.Response(200, json=page1))

    results = await client.get("transcript/")
    assert len(results) == 2


@respx.mock
@pytest.mark.asyncio
async def test_get_returns_dict_for_non_paginated_response(client):
    """When response is a plain dict (e.g. diff endpoint), return it directly."""
    payload = {"results": {"diff_me_stable_id": "X"}, "diff_me_transcript": {}}
    respx.get(BASE + "transcript/diff/").mock(return_value=httpx.Response(200, json=payload))

    result = await client.get_raw("transcript/diff/")
    assert result["results"]["diff_me_stable_id"] == "X"


@respx.mock
@pytest.mark.asyncio
async def test_404_returns_empty_list(client):
    respx.get(BASE + "transcript/").mock(return_value=httpx.Response(404, json={}))
    result = await client.get("transcript/")
    assert result == []


@respx.mock
@pytest.mark.asyncio
async def test_http_error_raises_exception(client):
    respx.get(BASE + "transcript/").mock(return_value=httpx.Response(400, json={"detail": "bad"}))
    with pytest.raises(Exception, match="400"):
        await client.get("transcript/")


@respx.mock
@pytest.mark.asyncio
async def test_cache_hit_avoids_second_request(client):
    """Second identical request is served from cache."""
    route = respx.get(BASE + "release/nopagination/").mock(
        return_value=httpx.Response(200, json=[{"shortname": "110"}])
    )
    await client.get("release/nopagination/")
    await client.get("release/nopagination/")
    assert route.call_count == 1
