import pytest
import httpx
import respx

from tark_mcp import client as client_module

BASE_URL = "https://tark.ensembl.org/api/"


@pytest.fixture(autouse=True)
def clear_tark_cache():
    """Clear the module-level TarkClient cache before each test to prevent cross-test pollution."""
    client_module._cache.clear()



# ---------------------------------------------------------------------------
# Minimal API fixture data — mirrors real TARK API response shapes
# ---------------------------------------------------------------------------

TRANSCRIPT_BRCA2_RAW = {
    "stable_id": "ENST00000380152",
    "stable_id_version": 7,
    "assembly": {"assembly_id": 1, "assembly_name": "GRCh38", "genome": 1, "session": 1},
    "loc_start": 32315475,   # 1-based → model stores 32315474
    "loc_end": 32400266,
    "loc_strand": 1,
    "loc_region": "13",
    "biotype": "protein_coding",
    "sequence": {"sequence": "ATCGATCGATCGATCGATCGATCGATCGATCG", "seq_checksum": "ABC"},
    "five_prime_utr_seq": "ATCG",
    "three_prime_utr_seq": "CG",
    "three_prime_utr_start": 32398771,
    "three_prime_utr_end": 32400266,
    "five_prime_utr_start": 32315475,
    "five_prime_utr_end": 32316460,
    "transcript_release_set": [
        {"assembly": "GRCh38", "shortname": "110", "description": "Ensembl release",
         "release_date": "2023-04-01", "source": "Ensembl"},
    ],
    "genes": [
        {
            "stable_id": "ENSG00000139618",
            "stable_id_version": 15,
            "assembly": "GRCh38",
            "loc_start": 32315475,
            "loc_end": 32400266,
            "loc_strand": 1,
            "loc_region": "13",
            "name": "BRCA2",
        }
    ],
    "translations": [
        {
            "stable_id": "ENSP00000369497",
            "stable_id_version": 3,
            "assembly": "GRCh38",
            "loc_start": 32316461,
            "loc_end": 32398770,
            "loc_strand": 1,
            "loc_region": "13",
            "transcript_stable_id": "ENST00000380152",
            "transcript_stable_id_version": 7,
        }
    ],
    "exons": [
        {
            "stable_id": "ENSE00001184784",
            "stable_id_version": 4,
            "assembly": "GRCh38",
            "loc_start": 32315475,
            "loc_end": 32315667,
            "loc_strand": 1,
            "loc_region": "13",
            "exon_order": 1,
            "transcript_stable_id": "ENST00000380152",
            "transcript_stable_id_version": 7,
        }
    ],
}

TRANSLATION_BRCA2_RAW = {
    "stable_id": "ENSP00000369497",
    "stable_id_version": 3,
    "assembly": {"assembly_id": 1, "assembly_name": "GRCh38", "genome": 1, "session": 1},
    "loc_start": 32316461,
    "loc_end": 32398770,
    "loc_strand": 1,
    "loc_region": "13",
    "sequence": {"sequence": "MPIGSKERP", "seq_checksum": "DEF"},
}

TRANSCRIPT_NONCODING_RAW = {
    "stable_id": "ENST00000614536",
    "stable_id_version": 1,
    "assembly": {"assembly_id": 1, "assembly_name": "GRCh38", "genome": 1, "session": 1},
    "loc_start": 54529654,
    "loc_end": 54529961,
    "loc_strand": 1,
    "loc_region": "19",
    "biotype": "unprocessed_pseudogene",
    "sequence": {"sequence": "GGCTTGTTCACA", "seq_checksum": "GHI"},
    "five_prime_utr_seq": None,
    "three_prime_utr_seq": None,
    "three_prime_utr_start": None,
    "three_prime_utr_end": None,
    "five_prime_utr_start": None,
    "five_prime_utr_end": None,
    "transcript_release_set": [
        {"assembly": "GRCh38", "shortname": "92", "description": "Ensembl release",
         "release_date": "2018-04-05", "source": "Ensembl"},
    ],
    "genes": [],
    "translations": [],
    "exons": [
        {
            "stable_id": "ENSE00003719241",
            "stable_id_version": 1,
            "assembly": "GRCh38",
            "loc_start": 54529654,
            "loc_end": 54529961,
            "loc_strand": 1,
            "loc_region": "19",
            "exon_order": 1,
            "transcript_stable_id": "ENST00000614536",
            "transcript_stable_id_version": 1,
        }
    ],
}

DIFF_RESPONSE_RAW = {
    "count": 1,
    "next": None,
    "previous": None,
    "results": {
        "diff_me_stable_id": "ENST00000380152",
        "diff_with_stable_id": "ENST00000614536",
        "diff_me_stable_id_version": 7,
        "diff_with_stable_id_version": 1,
        "diff_me_assembly": {"assembly_name": "GRCh38"},
        "diff_with_assembly": {"assembly_name": "GRCh38"},
        "has_stable_id_changed": True,
        "has_transcript_changed": True,
        "has_seq_changed": True,
        "has_exon_set_changed": True,
        "has_translation_stable_id_changed": True,
        "has_translation_seq_changed": None,
    },
    "diff_me_transcript": TRANSCRIPT_BRCA2_RAW,
    "diff_with_transcript": TRANSCRIPT_NONCODING_RAW,
}

RELEASE_LIST_RAW = [
    {
        "shortname": "110",
        "description": "Ensembl release 110",
        "release_date": "2023-04-01",
        "assembly": "GRCh38",
        "source": "Ensembl",
    }
]

MANE_LIST_RESPONSE_RAW = {
    "count": 1,
    "next": None,
    "previous": None,
    "results": [TRANSCRIPT_BRCA2_RAW],
}
