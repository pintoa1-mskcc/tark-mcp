import pytest
import httpx
import respx

from tark_mcp.client import TarkClient
from tark_mcp.models import Transcript
from tark_mcp.tools.versions import (
    get_transcript_all_versions,
    get_transcript_version_history,
    describe_divergent_earlier_version,
)
from tests.conftest import TRANSCRIPT_BRCA2_RAW

BASE = "https://tark.ensembl.org/api/"


@pytest.fixture
def client():
    return TarkClient()


def _make_version(version: int, shortname: str, release_date: str, extra_exon: bool = False) -> dict:
    raw = {
        **TRANSCRIPT_BRCA2_RAW,
        "stable_id_version": version,
        "transcript_release_set": [{
            "assembly": "GRCh38", "shortname": shortname, "description": "Ensembl release",
            "release_date": release_date, "source": "Ensembl",
        }],
    }
    if extra_exon:
        raw = {**raw, "exons": raw["exons"] + [{
            "stable_id": "ENSE00009999999", "stable_id_version": 1, "assembly": "GRCh38",
            "loc_start": 32400300, "loc_end": 32400400, "loc_strand": 1, "loc_region": "13",
            "exon_order": 2, "transcript_stable_id": "ENST00000380152",
            "transcript_stable_id_version": version,
        }]}
    return raw


@respx.mock
@pytest.mark.asyncio
async def test_get_transcript_all_versions_returns_every_version(client):
    v1 = _make_version(1, "108", "2020-01-01")
    v2 = _make_version(2, "116", "2026-01-01")
    respx.get(BASE + "transcript/").mock(
        return_value=httpx.Response(200, json={"count": 2, "next": None, "previous": None, "results": [v2, v1]})
    )
    versions = await get_transcript_all_versions("ENST00000380152", "GRCh38", client)
    assert [v.stable_id_version for v in versions] == [1, 2]  # sorted ascending


@respx.mock
@pytest.mark.asyncio
async def test_describe_divergent_earlier_version_flags_differing_content(client):
    """Reproduces the ENST00000706094 case: an earlier version (first release v108) has
    a different exon/CDS signature than the version get_transcript() auto-resolved to
    (v2, first release v116) — this must be surfaced, not silently hidden."""
    v1 = _make_version(1, "108", "2020-01-01")
    v2 = _make_version(2, "116", "2026-01-01", extra_exon=True)
    respx.get(BASE + "transcript/").mock(
        return_value=httpx.Response(200, json={"count": 2, "next": None, "previous": None, "results": [v1, v2]})
    )
    resolved = Transcript.model_validate(v2)
    note = await describe_divergent_earlier_version("ENST00000380152", resolved, "GRCh38", client)
    assert "v1" in note
    assert "108" in note


@respx.mock
@pytest.mark.asyncio
async def test_describe_divergent_earlier_version_silent_when_content_matches(client):
    """Versions that only differ in UTR trims (identical exon/CDS signature) — the common
    case — should not be flagged."""
    v1 = _make_version(1, "108", "2020-01-01")
    v2 = _make_version(2, "116", "2026-01-01")
    respx.get(BASE + "transcript/").mock(
        return_value=httpx.Response(200, json={"count": 2, "next": None, "previous": None, "results": [v1, v2]})
    )
    resolved = Transcript.model_validate(v2)
    note = await describe_divergent_earlier_version("ENST00000380152", resolved, "GRCh38", client)
    assert note == ""


@respx.mock
@pytest.mark.asyncio
async def test_describe_divergent_earlier_version_skips_refseq(client):
    """Only Ensembl ENST IDs are checked — RefSeq version semantics differ."""
    resolved = Transcript.model_validate(_make_version(1, "108", "2020-01-01"))
    note = await describe_divergent_earlier_version("NM_001141970", resolved, "GRCh38", client)
    assert note == ""


@respx.mock
@pytest.mark.asyncio
async def test_get_transcript_version_history_empty(client):
    respx.get(BASE + "transcript/").mock(
        return_value=httpx.Response(200, json={"count": 0, "next": None, "previous": None, "results": []})
    )
    history = await get_transcript_version_history("ENST00000999999", "GRCh38", client)
    assert history["version_count"] == 0
    assert history["versions"] == []
