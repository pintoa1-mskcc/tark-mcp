import pytest
import httpx
import respx

from tark_mcp.client import TarkClient
from tark_mcp.tools.sequences import (
    get_transcript_sequence, get_transcript_exons, get_protein_for_transcript
)
from tests.conftest import TRANSCRIPT_BRCA2_RAW, TRANSCRIPT_NONCODING_RAW, TRANSLATION_BRCA2_RAW

BASE = "https://tark.ensembl.org/api/"

SINGLE_RESULT = {"count": 1, "next": None, "previous": None,
                 "results": [TRANSCRIPT_BRCA2_RAW]}
NONCODING_RESULT = {"count": 1, "next": None, "previous": None,
                    "results": [TRANSCRIPT_NONCODING_RAW]}


@respx.mock
@pytest.mark.asyncio
async def test_get_transcript_sequence_returns_sequence():
    client = TarkClient()
    respx.get(BASE + "transcript/").mock(return_value=httpx.Response(200, json=SINGLE_RESULT))
    result = await get_transcript_sequence("ENST00000380152", client=client)
    assert result is not None
    assert result["stable_id"] == "ENST00000380152"
    assert result["sequence"] == "ATCGATCGATCGATCGATCGATCGATCGATCG"
    assert result["assembly"] == "GRCh38"


@respx.mock
@pytest.mark.asyncio
async def test_get_transcript_sequence_not_found_returns_none():
    client = TarkClient()
    respx.get(BASE + "transcript/").mock(
        return_value=httpx.Response(200, json={"count": 0, "next": None,
                                               "previous": None, "results": []})
    )
    result = await get_transcript_sequence("ENST00000999999", client=client)
    assert result is None


@respx.mock
@pytest.mark.asyncio
async def test_get_transcript_exons_returns_ordered_exons():
    client = TarkClient()
    respx.get(BASE + "transcript/").mock(return_value=httpx.Response(200, json=SINGLE_RESULT))
    exons = await get_transcript_exons("ENST00000380152", client=client)
    assert len(exons) == 1
    assert exons[0].order == 1
    assert exons[0].stable_id == "ENSE00001184784"


@respx.mock
@pytest.mark.asyncio
async def test_get_transcript_exons_negative_strand_reversed():
    """Exons on negative strand are returned in reverse order."""
    client = TarkClient()
    neg_transcript = {
        **TRANSCRIPT_BRCA2_RAW,
        "loc_strand": -1,
        "exons": [
            {**TRANSCRIPT_BRCA2_RAW["exons"][0], "exon_order": 1},
            {**TRANSCRIPT_BRCA2_RAW["exons"][0], "stable_id": "ENSE00002",
             "exon_order": 2, "loc_start": 32315700, "loc_end": 32315800},
        ],
    }
    respx.get(BASE + "transcript/").mock(
        return_value=httpx.Response(200, json={
            "count": 1, "next": None, "previous": None, "results": [neg_transcript]
        })
    )
    exons = await get_transcript_exons("ENST00000380152", client=client)
    # Negative strand: returned highest order first
    assert exons[0].order == 2
    assert exons[1].order == 1


@respx.mock
@pytest.mark.asyncio
async def test_get_protein_for_transcript_returns_translation():
    client = TarkClient()
    respx.get(BASE + "transcript/").mock(return_value=httpx.Response(200, json=SINGLE_RESULT))
    result = await get_protein_for_transcript("ENST00000380152", client=client)
    assert result is not None
    assert result.stable_id == "ENSP00000369497"


@respx.mock
@pytest.mark.asyncio
async def test_get_protein_for_noncoding_returns_none():
    client = TarkClient()
    respx.get(BASE + "transcript/").mock(return_value=httpx.Response(200, json=NONCODING_RESULT))
    result = await get_protein_for_transcript("ENST00000614536", client=client)
    assert result is None
