"""Integration tests against the live TARK API.

Run with:  TARK_INTEGRATION=1 pytest tests/integration/ -v
Skipped by default.
"""
import os
import pytest

pytestmark = pytest.mark.skipif(
    os.environ.get("TARK_INTEGRATION") != "1",
    reason="Set TARK_INTEGRATION=1 to run live API tests"
)

from tark_mcp.client import TarkClient
from tark_mcp.tools.releases import get_releases
from tark_mcp.tools.transcripts import get_transcript, search_transcripts_by_region
from tark_mcp.tools.genes import get_gene_transcripts
from tark_mcp.tools.sequences import (
    get_transcript_sequence, get_transcript_exons, get_protein_for_transcript
)
from tark_mcp.tools.mane import get_mane_transcripts
from tark_mcp.tools.diff import diff_transcripts


@pytest.fixture(scope="module")
def client():
    return TarkClient()


@pytest.mark.asyncio
async def test_get_releases(client):
    releases = await get_releases(client)
    assert len(releases) > 0
    assert all(hasattr(r, "shortname") for r in releases)


@pytest.mark.asyncio
async def test_get_transcript_brca2(client):
    t = await get_transcript("ENST00000380152", client=client)
    assert t is not None
    assert t.stable_id == "ENST00000380152"
    assert t.assembly == "GRCh38"
    assert len(t.exons) > 0


@pytest.mark.asyncio
async def test_get_transcript_both_assemblies(client):
    result = await get_transcript("ENST00000380152", assembly="both", client=client)
    assert isinstance(result, list)
    assemblies = {t.assembly for t in result}
    # Should have at least GRCh38
    assert "GRCh38" in assemblies


@pytest.mark.asyncio
async def test_search_transcripts_by_region(client):
    # BRCA2 locus on chr13
    results = await search_transcripts_by_region("13", 32315474, 32400266, client=client)
    assert len(results) > 0
    stable_ids = {t.stable_id for t in results}
    assert "ENST00000380152" in stable_ids


@pytest.mark.asyncio
async def test_get_gene_transcripts_by_symbol(client):
    results = await get_gene_transcripts("BRCA2", client=client)
    assert len(results) > 0
    assert all(t.assembly == "GRCh38" for t in results)


@pytest.mark.asyncio
async def test_get_transcript_sequence(client):
    result = await get_transcript_sequence("ENST00000380152", client=client)
    assert result is not None
    assert result["sequence"] is not None
    assert len(result["sequence"]) > 100


@pytest.mark.asyncio
async def test_get_transcript_exons(client):
    exons = await get_transcript_exons("ENST00000380152", client=client)
    assert len(exons) > 0
    # Exons should be in ascending order (positive strand)
    orders = [e.order for e in exons]
    assert orders == sorted(orders)


@pytest.mark.asyncio
async def test_get_protein_for_transcript(client):
    result = await get_protein_for_transcript("ENST00000380152", client=client)
    assert result is not None
    assert result.stable_id.startswith("ENSP")


@pytest.mark.asyncio
async def test_get_mane_transcripts(client):
    results = await get_mane_transcripts(client=client)
    assert len(results) > 0


@pytest.mark.asyncio
async def test_get_mane_transcripts_filtered(client):
    results = await get_mane_transcripts(gene_identifier="BRCA2", client=client)
    assert len(results) > 0
    for t in results:
        gene_names = {g.name for g in t.genes}
        assert "BRCA2" in gene_names


@pytest.mark.asyncio
async def test_diff_transcripts_two_coding(client):
    # BRCA2 stable versions
    results = await diff_transcripts(
        ["ENST00000380152", "ENST00000544455"], client=client
    )
    assert len(results) == 1
    diff = results[0]
    assert diff.reference_stable_id == "ENST00000380152"
    assert diff.candidate_stable_id == "ENST00000544455"
    # Both are in protein_coding transcripts
    assert diff.reference_protein_coding is True


@pytest.mark.asyncio
async def test_diff_transcripts_protein_sequences_populated(client):
    results = await diff_transcripts(
        ["ENST00000380152", "ENST00000544455"], client=client
    )
    diff = results[0]
    if diff.reference_protein_coding and diff.candidate_protein_coding:
        assert diff.ref_protein_sequence is not None
        assert diff.candidate_protein_sequence is not None
        assert diff.protein_sequence_changed is not None
