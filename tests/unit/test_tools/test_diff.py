import pytest
import httpx
import respx

from tark_mcp.client import TarkClient
from tark_mcp.tools.diff import _compute_protein_diffs, _resolve_transcript, diff_transcripts
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


def test_compute_protein_diffs_reports_substitution_with_aligned_positions():
    """A single substitution between equal-length sequences is reported at the
    same 1-based position on both sides."""
    diffs = _compute_protein_diffs("MAT", "MST")
    assert len(diffs) == 1
    d = diffs[0]
    assert d.change == "substitution"
    assert d.ref_position == 2
    assert d.candidate_position == 2
    assert d.ref_residue == "A"
    assert d.candidate_residue == "S"


def test_compute_protein_diffs_reports_deletion_with_none_candidate_position():
    """A residue present only in the reference is a deletion; candidate_position
    is None since there's no corresponding aligned residue on that side."""
    diffs = _compute_protein_diffs("ABC", "AC")
    assert len(diffs) == 1
    d = diffs[0]
    assert d.change == "deletion"
    assert d.ref_position == 2
    assert d.candidate_position is None
    assert d.ref_residue == "B"
    assert d.candidate_residue is None


def test_compute_protein_diffs_reports_insertion_with_none_ref_position():
    """A residue present only in the candidate is an insertion; ref_position is
    None since there's no corresponding aligned residue on that side."""
    diffs = _compute_protein_diffs("AC", "ABC")
    assert len(diffs) == 1
    d = diffs[0]
    assert d.change == "insertion"
    assert d.ref_position is None
    assert d.candidate_position == 2
    assert d.ref_residue is None
    assert d.candidate_residue == "B"


def test_compute_protein_diffs_identical_sequences_returns_empty_list():
    assert _compute_protein_diffs("MATS", "MATS") == []


def test_compute_protein_diffs_real_world_regression_r234s_p313s():
    """Regression test: NM_005027.3 vs ENST00000222254.8 (GRCh37) protein
    sequences, validated manually against the live TARK API. Expect exactly
    two substitutions: R234S and P313S, with no indels."""
    ref_seq = (
        "MAGPEGFQYRALYPFRRERPEDLELLPGDVLVVSRAALQALGVAEGGERCPQSVGWMPGLNERTRQRGDF"
        "PGTYVEFLGPVALARPGPRPRGPRPLPARPRDGAPEPGLTLPDLPEQFSPPDVAPPLLVKLVEAIERTGL"
        "DSESHYRPELPAPRTDWSLSDVDQWDTAALADGIKSFLLALPAPLVTPEASAEARRALREAAGPVGPALE"
        "PPTLPLHRALTLRFLLQHLGRVARRAPALGPAVRALGATFGPLLLRAPPPPSSPPPGGAPDGSEPSPDFP"
        "ALLVEKLLQEHLEEQEVAPPALPPKPPKAKPAPTVLANGGSPPSLQDAEWYWGDISREEVNEKLRDTPDG"
        "TFLVRDASSKIQGEYTLTLRKGGNNKLIKVFHRDGHYGFSEPLTFCSVVDLINHYRHESLAQYNAKLDTR"
        "LLYPVSKYQQDQIVKEDSVEAVGAQLKVYHQQYQDKSREYDQLYEEYTRTSQELQMKRTAIEAFNETIKI"
        "FEEQGQTQEKCSKEYLERFRREGNEKEMQRILLNSERLKSRIAEIHESRTKLEQQLRAQASDNREIDKRM"
        "NSLKPDLMQLRKIRDQYLVWLTQKGARQKKINEWLGIKNETEDQYALMEDEDDLPHHEERTWYVGKINR"
        "TQAEEMLSGKRDGTFLIRESSQRGCYACSVVVDGDTKHCVIYRTATGFGFAEPYNLYGSLKELVLHYQHA"
        "SLVQHNDALTVTLAHPVRAPGPGPPPAAR"
    )
    cand_seq = (
        "MAGPEGFQYRALYPFRRERPEDLELLPGDVLVVSRAALQALGVAEGGERCPQSVGWMPGLNERTRQRGDF"
        "PGTYVEFLGPVALARPGPRPRGPRPLPARPRDGAPEPGLTLPDLPEQFSPPDVAPPLLVKLVEAIERTGL"
        "DSESHYRPELPAPRTDWSLSDVDQWDTAALADGIKSFLLALPAPLVTPEASAEARRALREAAGPVGPALE"
        "PPTLPLHRALTLRFLLQHLGRVASRAPALGPAVRALGATFGPLLLRAPPPPSSPPPGGAPDGSEPSPDFP"
        "ALLVEKLLQEHLEEQEVAPPALPPKPPKAKPASTVLANGGSPPSLQDAEWYWGDISREEVNEKLRDTPDG"
        "TFLVRDASSKIQGEYTLTLRKGGNNKLIKVFHRDGHYGFSEPLTFCSVVDLINHYRHESLAQYNAKLDTR"
        "LLYPVSKYQQDQIVKEDSVEAVGAQLKVYHQQYQDKSREYDQLYEEYTRTSQELQMKRTAIEAFNETIKI"
        "FEEQGQTQEKCSKEYLERFRREGNEKEMQRILLNSERLKSRIAEIHESRTKLEQQLRAQASDNREIDKRM"
        "NSLKPDLMQLRKIRDQYLVWLTQKGARQKKINEWLGIKNETEDQYALMEDEDDLPHHEERTWYVGKINR"
        "TQAEEMLSGKRDGTFLIRESSQRGCYACSVVVDGDTKHCVIYRTATGFGFAEPYNLYGSLKELVLHYQHA"
        "SLVQHNDALTVTLAHPVRAPGPGPPPAAR"
    )
    diffs = _compute_protein_diffs(ref_seq, cand_seq)
    assert len(diffs) == 2

    r234s = diffs[0]
    assert r234s.change == "substitution"
    assert r234s.ref_position == 234
    assert r234s.candidate_position == 234
    assert r234s.ref_residue == "R"
    assert r234s.candidate_residue == "S"

    p313s = diffs[1]
    assert p313s.change == "substitution"
    assert p313s.ref_position == 313
    assert p313s.candidate_position == 313
    assert p313s.ref_residue == "P"
    assert p313s.candidate_residue == "S"


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

    # ref="MPIGSKERP" (9aa), candidate="MVLSPAD" (7aa) aligned via NW
    # (match=2, mismatch=-1, gap=-2) gives: "MPIGSKERP" / "M-VLS-PAD"
    assert [d.change for d in diff.protein_diffs] == [
        "deletion", "substitution", "substitution",
        "deletion", "substitution", "substitution", "substitution",
    ]
    assert diff.protein_diffs[0].ref_position == 2
    assert diff.protein_diffs[0].candidate_position is None
    assert diff.protein_diffs[0].ref_residue == "P"
    assert diff.protein_diffs[1].ref_position == 3
    assert diff.protein_diffs[1].candidate_position == 2
    assert diff.protein_diffs[1].ref_residue == "I"
    assert diff.protein_diffs[1].candidate_residue == "V"


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
    assert diff.protein_diffs is None


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
    assert diff.protein_diffs is None


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
