from __future__ import annotations
from tark_mcp.client import TarkClient
from tark_mcp.models import Release


async def get_releases(client: TarkClient) -> list[Release]:
    data = await client.get("release/nopagination/")
    return [Release.model_validate(r) for r in data]
