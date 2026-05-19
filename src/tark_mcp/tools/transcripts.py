from __future__ import annotations
import asyncio
from tark_mcp.client import TarkClient
from tark_mcp.models import Transcript


def _strip_version(stable_id: str) -> tuple[str, int | None]:
    """Split 'ENST00000380152.7' → ('ENST00000380152', 7). Returns None version if absent."""
    if "." in stable_id:
        parts = stable_id.rsplit(".", 1)
        try:
            return parts[0], int(parts[1])
        except ValueError:
            pass
    return stable_id, None


def _deduplicate(transcripts: list[Transcript]) -> list[Transcript]:
    """For each (assembly, stable_id, stable_id_version) keep the one with the most recent release."""
    best: dict[tuple, Transcript] = {}
    for t in transcripts:
        key = (t.assembly, t.stable_id, t.stable_id_version)
        existing = best.get(key)
        if existing is None or (t.latest_release_date or "") > (existing.latest_release_date or ""):
            best[key] = t
    return list(best.values())


async def _fetch_for_assembly(
    stable_id: str,
    version: int | None,
    assembly: str,
    client: TarkClient,
) -> list[Transcript]:
    params: dict = {"stable_id": stable_id, "expand_all": "true",
                    "assembly_name": assembly}
    if version is not None:
        params["stable_id_version"] = version
    data = await client.get("transcript/", params=params)
    transcripts = [Transcript.model_validate(r) for r in data]
    return _deduplicate(transcripts)


async def get_transcript(
    stable_id: str,
    assembly: str = "GRCh38",
    client: TarkClient | None = None,
) -> Transcript | list[Transcript] | None:
    if client is None:
        client = TarkClient()
    sid, version = _strip_version(stable_id)

    if assembly == "both":
        results = await asyncio.gather(
            _fetch_for_assembly(sid, version, "GRCh38", client),
            _fetch_for_assembly(sid, version, "GRCh37", client),
        )
        combined = results[0] + results[1]
        return combined if combined else None

    transcripts = await _fetch_for_assembly(sid, version, assembly, client)
    if not transcripts:
        return None
    # Return the single most-recently-released record for this assembly
    return max(transcripts, key=lambda t: t.latest_release_date or "")


async def search_transcripts_by_region(
    region: str,
    start: int,
    end: int,
    assembly: str = "GRCh38",
    client: TarkClient | None = None,
) -> list[Transcript]:
    if client is None:
        client = TarkClient()
    # Strip chr prefix
    loc_region = region.removeprefix("chr")
    # Convert 0-based → 1-based
    params = {
        "assembly_name": assembly,
        "loc_region": loc_region,
        "loc_start": start + 1,
        "loc_end": end,
        "expand": "transcript_release_set",
    }
    if assembly == "both":
        results = await asyncio.gather(
            client.get("transcript/", {**params, "assembly_name": "GRCh38"}),
            client.get("transcript/", {**params, "assembly_name": "GRCh37"}),
        )
        data = results[0] + results[1]
    else:
        data = await client.get("transcript/", params)

    return _deduplicate([Transcript.model_validate(r) for r in data])
