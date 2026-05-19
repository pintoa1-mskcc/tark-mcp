import pytest
import httpx
import respx

from tark_mcp.client import TarkClient
from tark_mcp.tools.diff import diff_transcripts
from tests.conftest import (
    TRANSCRIPT_BRCA2_RAW, TRANSCRIPT_NONCODING_RAW,
    TRANSLATION_BRCA2_RAW, DIFF_RESPONSE_RAW
)

BASE = "https://tark.ensembl.org/api/"

DIFF_BOTH_CODING = {
    "count": 1, "next": None, "previous": None,
    "results": {
        "diff_me_stable_id": "ENST00000380152",
        "diff_with_stable_id": "ENST00000614536",
        "has_seq_changed": True,
        "has_exon_set_changed": True,
        "has_translation_seq_changed": True,
    },
    "diff_me_transcript": {
        **TRANSCRIPT_BRCA2_RAW,
        "exons": [
            {**TRANSCRIPT_BRCA2_RAW["exons"][0], "exon_order": 1},
            {**TRANSCRIPT_BRCA2_RAW["exons"][0], "stable_id": "ENSE00002",
             "exon_order": 2, "loc_start": 32316000, "loc_end": 32316500},
        ],
    },
    "diff_with_transcript": {
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
    },
}

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


@respx.mock
@pytest.mark.asyncio
async def test_diff_transcripts_coding_pair_populates_all_sequence_fields():
    """Both transcripts coding: all sequence fields populated, changed flags computed."""
    client = TarkClient()
    respx.get(BASE + "transcript/diff/").mock(
        return_value=httpx.Response(200, json=DIFF_BOTH_CODING)
    )
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

    # ref CDS: seq="ATCGATCGATCGATCGATCGATCGATCGATCG"(len=32), 5'UTR="ATCG"(4), 3'UTR="CG"(2)
    # cds_seq = seq[4:30] = 26 chars
    assert diff.ref_cds_sequence == "ATCGATCGATCGATCGATCGATCGAT"
    # candidate: seq="TTTTGGGGCCCCAAAA"(16), 5'UTR="TTTT"(4), 3'UTR="AAAA"(4)
    # cds_seq = seq[4:12]
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
    diff_response = {
        **DIFF_RESPONSE_RAW,
        "diff_me_transcript": TRANSCRIPT_NONCODING_RAW,
        "diff_with_transcript": TRANSCRIPT_NONCODING_RAW,
    }
    respx.get(BASE + "transcript/diff/").mock(
        return_value=httpx.Response(200, json=diff_response)
    )
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
    respx.get(BASE + "transcript/diff/").mock(
        return_value=httpx.Response(200, json=DIFF_RESPONSE_RAW)
    )
    # Only one translation fetch (for the coding ref)
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
    respx.get(BASE + "transcript/diff/").mock(
        return_value=httpx.Response(200, json=DIFF_RESPONSE_RAW)
    )
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
    client = TarkClient()
    respx.get(BASE + "transcript/diff/").mock(
        return_value=httpx.Response(200, json=DIFF_RESPONSE_RAW)
    )
    respx.get(BASE + "translation/").mock(
        return_value=httpx.Response(200, json=TRANSLATION_REF_RESPONSE)
    )
    results = await diff_transcripts(
        ["ENST00000380152", "ENST00000614536", "ENST00000614536"], client=client
    )
    assert len(results) == 2
