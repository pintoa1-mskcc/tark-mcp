import pytest
import httpx
import respx

from tark_mcp.client import TarkClient
from tark_mcp.tools.releases import get_releases
from tests.conftest import RELEASE_LIST_RAW

BASE = "https://tark.ensembl.org/api/"


@respx.mock
@pytest.mark.asyncio
async def test_get_releases_returns_release_list():
    client = TarkClient()
    respx.get(BASE + "release/nopagination/").mock(
        return_value=httpx.Response(200, json=RELEASE_LIST_RAW)
    )
    releases = await get_releases(client)
    assert len(releases) == 1
    assert releases[0].shortname == "110"
    assert releases[0].assembly == "GRCh38"
    assert releases[0].release_date == "2023-04-01"
