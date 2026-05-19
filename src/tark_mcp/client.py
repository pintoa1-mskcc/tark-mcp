from __future__ import annotations
import os
from typing import Any

import httpx
from cachetools import TTLCache
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception

BASE_URL = os.environ.get("TARK_BASE_URL", "https://tark.ensembl.org/api/")
CACHE_TTL = int(os.environ.get("TARK_CACHE_TTL", "3600"))
REQUEST_TIMEOUT = int(os.environ.get("TARK_REQUEST_TIMEOUT", "30"))
MAX_RETRIES = int(os.environ.get("TARK_MAX_RETRIES", "3"))

_cache: TTLCache = TTLCache(maxsize=512, ttl=CACHE_TTL)


def _is_transient(exc: BaseException) -> bool:
    if isinstance(exc, httpx.HTTPStatusError):
        return exc.response.status_code >= 500
    return isinstance(exc, (httpx.NetworkError, httpx.TimeoutException))


class TarkClient:
    def __init__(self) -> None:
        self._http = httpx.AsyncClient(base_url=BASE_URL, timeout=REQUEST_TIMEOUT)

    @retry(
        stop=stop_after_attempt(MAX_RETRIES),
        wait=wait_exponential(multiplier=1, min=1, max=10),
        retry=retry_if_exception(_is_transient),
        reraise=True,
    )
    async def _fetch(self, url: str, params: dict | None = None) -> Any:
        cache_key = url + str(sorted((params or {}).items()))
        if cache_key in _cache:
            return _cache[cache_key]

        response = await self._http.get(url, params=params)

        if response.status_code == 404:
            return None
        if response.status_code >= 400:
            raise Exception(f"TARK API error {response.status_code}: {response.text[:200]}")

        data = response.json()
        _cache[cache_key] = data
        return data

    async def get(self, path: str, params: dict | None = None) -> list[dict]:
        """Fetch a paginated endpoint and return aggregated results list."""
        url = BASE_URL + path
        results: list[dict] = []
        while url:
            data = await self._fetch(url, params if url == BASE_URL + path else None)
            if data is None:
                break
            if isinstance(data, list):
                results.extend(data)
                break
            if "results" in data and isinstance(data["results"], list):
                results.extend(data["results"])
                next_url = data.get("next")
                if next_url:
                    url = next_url.replace("http://", "https://")
                    params = None
                else:
                    break
            else:
                results.append(data.get("results", data))
                break
        return results

    async def get_raw(self, path: str, params: dict | None = None) -> dict:
        """Fetch an endpoint that returns a plain dict (e.g. diff endpoint)."""
        url = BASE_URL + path
        data = await self._fetch(url, params)
        return data or {}

    async def aclose(self) -> None:
        await self._http.aclose()
