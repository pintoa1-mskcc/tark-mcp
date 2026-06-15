import pytest
import httpx
import respx

from tark_mcp.client import TarkClient
from tark_mcp.tools.diff import _resolve_transcript, diff_transcripts
from tests.conftest import (
    TRANSCRIPT_BRCA2_RAW, TRANSCRIPT_NONCODING_RAW,
    TRANSLATION_BRCA2_RAW,
)

BASE = "https://tark.ensembl.org/api/"


# Paginated wrapper used by client.get()
def _page(items):
    return {"count": len(items), "next": None, "previous": None, "results": items}


TRANSCRIPT_BRCA2_PAGE = _page([TRANSCRIPT_BRCA2_RAW])
TRANSCRIPT_NONCODING_PAGE = _page([TRANSCRIPT_NONCODING_RAW])

TRANSLATION_CANDIDATE_RAW = {
    "count": 1, "next": None, "previous": None,
    "results": [{
        "stable_id": "ENSP00000999999",
        "stable_id_version": 1,
        "assembly": {"assembly_name": "GRCh38", "assembly_id": 1, "genome": 1, "session": 1},
        "loc_start": 100, "loc_end": 200, "loc_strand": 1, "loc_region": "13",
        "sequence": {"sequence": "MVLSPAD", "seq_checksum": "ZZZ"},
    }]
}

TRANSLATION_REF_RESPONSE = {
    "count": 1, "next": None, "previous": None,
    "results": [{**TRANSLATION_BRCA2_RAW,
                 "assembly": {"assembly_name": "GRCh38", "assembly_id": 1,
                              "genome": 1, "session": 1}}]
}

TRANSCRIPT_CODING_CANDIDATE_RAW = {
    **TRANSCRIPT_BRCA2_RAW,
    "stable_id": "ENST00000614536",
    "stable_id_version": 1,
    "biotype": "protein_coding",
    "sequence": {"sequence": "TTTTGGGGCCCCAAAA", "seq_checksum": "XYZ"},
    "five_prime_utr_seq": "TTTT",
    "three_prime_utr_seq": "AAAA",
    "translations": [
        {"stable_id": "ENSP00000999999", "stable_id_version": 1,
         "assembly": "GRCh38", "loc_start": 100, "loc_end": 200,
         "loc_strand": 1, "loc_region": "13",
         "transcript_stable_id": "ENST00000614536",
         "transcript_stable_id_version": 1}
    ],
    "exons": [
        {**TRANSCRIPT_BRCA2_RAW["exons"][0], "exon_order": 1},
    ],
}
TRANSCRIPT_CODING_CANDIDATE_PAGE = _page([TRANSCRIPT_CODING_CANDIDATE_RAW])


@respx.mock
@pytest.mark.asyncio
async def test_resolve_transcript_uses_version_tiebreaker_when_release_dates_missing():
    """Unversioned lookup prefers highest stable_id_version when release dates are missing."""
    refseq_v2_raw = {
        **TRANSCRIPT_BRCA2_RAW,
        "stable_id": "NM_001128425",
        "stable_id_version": 2,
        "transcript_release_set": [],
    }
    refseq_v1_raw = {
        **TRANSCRIPT_BRCA2_RAW,
        "stable_id": "NM_001128425",
        "stable_id_version": 1,
        "transcript_release_set": [],
    }

    client = TarkClient()
    respx.get(BASE + "transcript/").mock(
        return_value=httpx.Response(200, json=_page([refseq_v1_raw, refseq_v2_raw]))
    )

    transcript = await _resolve_transcript("NM_001128425", "GRCh38", client)

    assert transcript.stable_id_version == 2


@respx.mock
@pytest.mark.asyncio
async def test_diff_transcripts_coding_pair_populates_all_sequence_fields():
    """Both transcripts coding: all sequence fields populated, changed flags computed."""
    client = TarkClient()
    respx.get(BASE + "transcript/").mock(side_effect=[
        httpx.Response(200, json=TRANSCRIPT_BRCA2_PAGE),
        httpx.Response(200, json=TRANSCRIPT_CODING_CANDIDATE_PAGE),
    ])
    respx.get(BASE + "translation/").mock(side_effect=[
        httpx.Response(200, json=TRANSLATION_REF_RESPONSE),
        httpx.Response(200, json=TRANSLATION_CANDIDATE_RAW),
    ])

    results = await diff_transcripts(
        ["ENST00000380152", "ENST00000614536"], client=client
    )
    assert len(results) == 1
    diff = results[0]

    assert diff.reference_stable_id == "ENST00000380152"
    assert diff.candidate_stable_id == "ENST00000614536"
    assert diff.reference_protein_coding is True
    assert diff.candidate_protein_coding is True

    # ref CDS: seq="ATCGATCGATCGATCGATCGATCGATCGATCG"(32), 5'UTR="ATCG"(4), 3'UTR="CG"(2)
    assert diff.ref_cds_sequence == "ATCGATCGATCGATCGATCGATCGAT"
    # candidate: seq="TTTTGGGGCCCCAAAA"(16), 5'UTR="TTTT"(4), 3'UTR="AAAA"(4)
    assert diff.candidate_cds_sequence == "GGGGCCCC"
    assert diff.cds_sequence_changed is True

    assert diff.ref_protein_sequence == "MPIGSKERP"
    assert diff.candidate_protein_sequence == "MVLSPAD"
    assert diff.protein_sequence_changed is True


@respx.mock
@pytest.mark.asyncio
async def test_diff_transcripts_noncoding_ref_sets_none_sentinels():
    """Non-coding reference: protein_coding=False, sequence comparison fields=None."""
    client = TarkClient()
    respx.get(BASE + "transcript/").mock(side_effect=[
        httpx.Response(200, json=TRANSCRIPT_NONCODING_PAGE),
        httpx.Response(200, json=TRANSCRIPT_NONCODING_PAGE),
    ])
    results = await diff_transcripts(
        ["ENST00000614536", "ENST00000614536"], client=client
    )
    diff = results[0]
    assert diff.reference_protein_coding is False
    assert diff.candidate_protein_coding is False
    assert diff.cds_sequence_changed is None
    assert diff.protein_sequence_changed is None
    assert diff.ref_cds_sequence is None
    assert diff.candidate_cds_sequence is None
    assert diff.ref_protein_sequence is None
    assert diff.candidate_protein_sequence is None


@respx.mock
@pytest.mark.asyncio
async def test_diff_transcripts_mixed_pair_sets_none_sentinels():
    """One coding, one non-coding: sequence comparison fields are None."""
    client = TarkClient()
    respx.get(BASE + "transcript/").mock(side_effect=[
        httpx.Response(200, json=TRANSCRIPT_BRCA2_PAGE),
        httpx.Response(200, json=TRANSCRIPT_NONCODING_PAGE),
    ])
    respx.get(BASE + "translation/").mock(
        return_value=httpx.Response(200, json=TRANSLATION_REF_RESPONSE)
    )
    results = await diff_transcripts(
        ["ENST00000380152", "ENST00000614536"], client=client
    )
    diff = results[0]
    assert diff.reference_protein_coding is True
    assert diff.candidate_protein_coding is False
    assert diff.cds_sequence_changed is None
    assert diff.protein_sequence_changed is None


@respx.mock
@pytest.mark.asyncio
async def test_diff_transcripts_requires_at_least_two_ids():
    client = TarkClient()
    with pytest.raises(Exception, match="At least 2 stable IDs"):
        await diff_transcripts(["ENST00000380152"], client=client)


@respx.mock
@pytest.mark.asyncio
async def test_diff_transcripts_exon_diffs_computed():
    """ExonDiff list is computed from exon lists of both transcripts."""
    client = TarkClient()
    respx.get(BASE + "transcript/").mock(side_effect=[
        httpx.Response(200, json=TRANSCRIPT_BRCA2_PAGE),
        httpx.Response(200, json=TRANSCRIPT_NONCODING_PAGE),
    ])
    respx.get(BASE + "translation/").mock(
        return_value=httpx.Response(200, json=TRANSLATION_REF_RESPONSE)
    )
    results = await diff_transcripts(
        ["ENST00000380152", "ENST00000614536"], client=client
    )
    diff = results[0]
    assert len(diff.exon_diffs) >= 1
    assert all(d.change in ("added", "removed", "modified", "unchanged") for d in diff.exon_diffs)


@respx.mock
@pytest.mark.asyncio
async def test_diff_transcripts_multiple_pairs():
    """Three stable IDs → two (ref, candidate) pairs, both processed."""
    second_candidate_raw = {
        **TRANSCRIPT_NONCODING_RAW,
        "stable_id": "ENST00000999999",
        "stable_id_version": 1,
    }
    second_candidate_page = _page([second_candidate_raw])

    client = TarkClient()
    # ref is resolved once before the gather; candidates fetched concurrently: 3 HTTP calls total
    respx.get(BASE + "transcript/").mock(side_effect=[
        httpx.Response(200, json=TRANSCRIPT_BRCA2_PAGE),     # ref (resolved once)
        httpx.Response(200, json=TRANSCRIPT_NONCODING_PAGE), # cand1
        httpx.Response(200, json=second_candidate_page),     # cand2
    ])
    respx.get(BASE + "translation/").mock(
        return_value=httpx.Response(200, json=TRANSLATION_REF_RESPONSE)
    )
    results = await diff_transcripts(
        ["ENST00000380152", "ENST00000614536", "ENST00000999999"], client=client
    )
    assert len(results) == 2
    assert results[0].candidate_stable_id == "ENST00000614536"
    assert results[1].candidate_stable_id == "ENST00000999999"


@respx.mock
@pytest.mark.asyncio
async def test_diff_transcripts_strips_version_suffix():
    """Versioned IDs like 'ENST00000380152.7' are stripped and stable_id_version passed separately."""
    client = TarkClient()
    respx.get(BASE + "transcript/").mock(side_effect=[
        httpx.Response(200, json=TRANSCRIPT_BRCA2_PAGE),
        httpx.Response(200, json=TRANSCRIPT_NONCODING_PAGE),
    ])
    respx.get(BASE + "translation/").mock(
        return_value=httpx.Response(200, json=TRANSLATION_REF_RESPONSE)
    )
    results = await diff_transcripts(
        ["ENST00000380152.7", "ENST00000614536.1"], client=client
    )
    assert len(results) == 1
    diff = results[0]
    # stable_id stored on Transcript is the bare ID (no version suffix)
    assert diff.reference_stable_id == "ENST00000380152"
    assert diff.candidate_stable_id == "ENST00000614536"


@respx.mock
@pytest.mark.asyncio
async def test_diff_transcripts_refseq_ids():
    """RefSeq stable IDs (NM_...) are accepted and fetched via /transcript/ endpoint."""
    refseq_v2_raw = {
        **TRANSCRIPT_BRCA2_RAW,
        "stable_id": "NM_001128425",
        "stable_id_version": 2,
    }
    refseq_v1_raw = {
        **TRANSCRIPT_BRCA2_RAW,
        "stable_id": "NM_001128425",
        "stable_id_version": 1,
        "sequence": {"sequence": "TTTTGGGGCCCCAAAA", "seq_checksum": "XYZ"},
        "five_prime_utr_seq": "TTTT",
        "three_prime_utr_seq": "AAAA",
        "translations": [],
    }

    client = TarkClient()
    respx.get(BASE + "transcript/").mock(side_effect=[
        httpx.Response(200, json=_page([refseq_v2_raw])),
        httpx.Response(200, json=_page([refseq_v1_raw])),
    ])
    respx.get(BASE + "translation/").mock(
        return_value=httpx.Response(200, json=TRANSLATION_REF_RESPONSE)
    )

    results = await diff_transcripts(
        ["NM_001128425.2", "NM_001128425.1"], client=client
    )
    assert len(results) == 1
    diff = results[0]
    assert diff.reference_stable_id == "NM_001128425"
    assert diff.candidate_stable_id == "NM_001128425"


@respx.mock
@pytest.mark.asyncio
async def test_diff_transcripts_mixed_refseq_ensembl():
    """RefSeq and Ensembl IDs can be diffed against each other."""
    refseq_raw = {
        **TRANSCRIPT_BRCA2_RAW,
        "stable_id": "NM_001128425",
        "stable_id_version": 2,
    }

    client = TarkClient()
    respx.get(BASE + "transcript/").mock(side_effect=[
        httpx.Response(200, json=_page([refseq_raw])),
        httpx.Response(200, json=TRANSCRIPT_NONCODING_PAGE),
    ])
    respx.get(BASE + "translation/").mock(
        return_value=httpx.Response(200, json=TRANSLATION_REF_RESPONSE)
    )

    results = await diff_transcripts(
        ["NM_001128425.2", "ENST00000614536"], client=client
    )
    assert len(results) == 1
    diff = results[0]
    assert diff.reference_stable_id == "NM_001128425"
    assert diff.candidate_stable_id == "ENST00000614536"


@respx.mock
@pytest.mark.asyncio
async def test_diff_transcripts_not_found_raises():
    """ValueError raised when a transcript cannot be found."""
    client = TarkClient()
    respx.get(BASE + "transcript/").mock(side_effect=[
        httpx.Response(200, json=_page([])),
        httpx.Response(200, json=TRANSCRIPT_NONCODING_PAGE),
    ])
    with pytest.raises(ValueError, match="Transcript not found"):
        await diff_transcripts(["NM_NOTREAL.1", "ENST00000614536"], client=client)


@respx.mock
@pytest.mark.asyncio
async def test_diff_transcripts_version_mismatch_raises():
    """ValueError raised when API returns transcripts but not the requested version."""
    mismatched_version_raw = {
        **TRANSCRIPT_BRCA2_RAW,
        "stable_id_version": 8,
    }

    client = TarkClient()
    respx.get(BASE + "transcript/").mock(side_effect=[
        httpx.Response(200, json=_page([mismatched_version_raw])),
        httpx.Response(200, json=TRANSCRIPT_NONCODING_PAGE),
    ])
    with pytest.raises(ValueError, match="Transcript not found"):
        await diff_transcripts(["ENST00000380152.7", "ENST00000614536"], client=client)
