import pytest
from tark_mcp.models import (
    Exon, Gene, Translation, Transcript, TranscriptDiff, ExonDiff, Release
)


def test_exon_coordinate_normalization():
    raw = {
        "stable_id": "ENSE00001184784",
        "stable_id_version": 4,
        "transcript_stable_id": "ENST00000380152",
        "transcript_stable_id_version": 7,
        "assembly": "GRCh38",
        "loc_start": 32315475,   # 1-based
        "loc_end": 32315667,
        "loc_strand": 1,
        "loc_region": "13",
        "exon_order": 1,
    }
    exon = Exon.model_validate(raw)
    assert exon.loc_start == 32315474  # 0-based
    assert exon.loc_end == 32315667    # unchanged (exclusive end)
    assert exon.order == 1
    assert exon.stable_id == "ENSE00001184784"
    assert exon.stable_id_version == 4
    assert exon.transcript_stable_id == "ENST00000380152"
    assert exon.transcript_stable_id_version == 7
    assert exon.assembly == "GRCh38"


def test_gene_coordinate_normalization():
    raw = {
        "stable_id": "ENSG00000139618",
        "stable_id_version": 15,
        "assembly": "GRCh38",
        "loc_start": 32315475,
        "loc_end": 32400266,
        "loc_strand": 1,
        "loc_region": "13",
        "name": "BRCA2",
    }
    gene = Gene.model_validate(raw)
    assert gene.loc_start == 32315474
    assert gene.loc_end == 32400266
    assert gene.name == "BRCA2"


def test_gene_name_can_be_none():
    raw = {
        "stable_id": "ENSG00000274968",
        "stable_id_version": 1,
        "assembly": "GRCh38",
        "loc_start": 100,
        "loc_end": 200,
        "loc_strand": 1,
        "loc_region": "19",
        "name": None,
    }
    gene = Gene.model_validate(raw)
    assert gene.name is None


def test_translation_with_sequence():
    raw = {
        "stable_id": "ENSP00000369497",
        "stable_id_version": 3,
        "transcript_stable_id": "ENST00000380152",
        "transcript_stable_id_version": 7,
        "assembly": "GRCh38",
        "sequence": {"sequence": "MPIGSKERP", "seq_checksum": "DEF"},
    }
    t = Translation.model_validate(raw)
    assert t.stable_id == "ENSP00000369497"
    assert t.sequence == "MPIGSKERP"


def test_translation_without_sequence():
    raw = {
        "stable_id": "ENSP00000369497",
        "stable_id_version": 3,
        "transcript_stable_id": "ENST00000380152",
        "transcript_stable_id_version": 7,
        "assembly": "GRCh38",
    }
    t = Translation.model_validate(raw)
    assert t.sequence is None


def test_translation_nested_assembly_is_ignored():
    raw = {
        "stable_id": "ENSP00000369497",
        "stable_id_version": 3,
        "transcript_stable_id": "ENST00000380152",
        "transcript_stable_id_version": 7,
        "assembly": "GRCh38",
    }
    t = Translation.model_validate(raw)
    assert t.assembly == "GRCh38"


def test_transcript_coordinate_normalization():
    raw = {
        "stable_id": "ENST00000380152",
        "stable_id_version": 7,
        "assembly": {"assembly_name": "GRCh38"},
        "loc_start": 32315475,
        "loc_end": 32400266,
        "loc_strand": 1,
        "loc_region": "13",
        "biotype": "protein_coding",
        "sequence": {"sequence": "ATCGATCGATCGATCGATCGATCGATCGATCG"},
        "five_prime_utr_seq": "ATCG",
        "three_prime_utr_seq": "CG",
        "transcript_release_set": [
            {"release_date": "2023-04-01", "assembly": "GRCh38", "shortname": "110",
             "description": "Ensembl release", "source": "Ensembl"}
        ],
        "genes": [],
        "translations": [],
        "exons": [],
    }
    t = Transcript.model_validate(raw)
    assert t.assembly == "GRCh38"
    assert t.loc_start == 32315474  # 0-based
    assert t.loc_end == 32400266
    assert t.sequence == "ATCGATCGATCGATCGATCGATCGATCGATCG"


def test_transcript_cds_boundaries_computed_from_utrs():
    """cds_start = len(5'UTR), cds_end = len(seq) - len(3'UTR)"""
    raw = {
        "stable_id": "ENST00000380152",
        "stable_id_version": 7,
        "assembly": {"assembly_name": "GRCh38"},
        "loc_start": 1,
        "loc_end": 100,
        "loc_strand": 1,
        "loc_region": "13",
        "biotype": "protein_coding",
        "sequence": {"sequence": "AAAACCCCGGGG"},  # len=12
        "five_prime_utr_seq": "AAAA",               # len=4  → cds_start=4
        "three_prime_utr_seq": "GGGG",              # len=4  → cds_end=12-4=8
        "transcript_release_set": [],
        "genes": [],
        "translations": [],
        "exons": [],
    }
    t = Transcript.model_validate(raw)
    assert t.cds_start == 4
    assert t.cds_end == 8


def test_transcript_cds_none_when_no_utrs():
    raw = {
        "stable_id": "ENST00000614536",
        "stable_id_version": 1,
        "assembly": {"assembly_name": "GRCh38"},
        "loc_start": 1,
        "loc_end": 100,
        "loc_strand": 1,
        "loc_region": "19",
        "biotype": "unprocessed_pseudogene",
        "sequence": {"sequence": "GGCTTGTTCACA"},
        "five_prime_utr_seq": None,
        "three_prime_utr_seq": None,
        "transcript_release_set": [],
        "genes": [],
        "translations": [],
        "exons": [],
    }
    t = Transcript.model_validate(raw)
    assert t.cds_start is None
    assert t.cds_end is None


def test_transcript_latest_release_date():
    raw = {
        "stable_id": "ENST00000380152",
        "stable_id_version": 7,
        "assembly": {"assembly_name": "GRCh38"},
        "loc_start": 1, "loc_end": 100, "loc_strand": 1, "loc_region": "13",
        "biotype": "protein_coding",
        "sequence": {"sequence": "ATCG"},
        "five_prime_utr_seq": None, "three_prime_utr_seq": None,
        "transcript_release_set": [
            {"release_date": "2014-05-01", "assembly": "GRCh38", "shortname": "76",
             "description": "Ensembl release", "source": "Ensembl"},
            {"release_date": "2023-04-01", "assembly": "GRCh38", "shortname": "110",
             "description": "Ensembl release", "source": "Ensembl"},
        ],
        "genes": [], "translations": [], "exons": [],
    }
    t = Transcript.model_validate(raw)
    assert t.latest_release_date == "2023-04-01"


def test_transcript_exons_parsed_correctly():
    raw = {
        "stable_id": "ENST00000380152",
        "stable_id_version": 7,
        "assembly": {"assembly_name": "GRCh38"},
        "loc_start": 1, "loc_end": 1000, "loc_strand": 1, "loc_region": "13",
        "biotype": "protein_coding",
        "sequence": {"sequence": "ATCG"},
        "five_prime_utr_seq": None, "three_prime_utr_seq": None,
        "transcript_release_set": [],
        "genes": [],
        "translations": [],
        "exons": [
            {
                "stable_id": "ENSE00001184784",
                "stable_id_version": 4,
                "transcript_stable_id": "ENST00000380152",
                "transcript_stable_id_version": 7,
                "assembly": "GRCh38",
                "loc_start": 32315475,
                "loc_end": 32315667,
                "loc_strand": 1,
                "loc_region": "13",
                "exon_order": 1,
            }
        ],
    }
    t = Transcript.model_validate(raw)
    assert len(t.exons) == 1
    assert t.exons[0].order == 1
    assert t.exons[0].loc_start == 32315474  # 0-based


def test_release_parsed():
    raw = {"shortname": "110", "description": "Ensembl release 110",
           "release_date": "2023-04-01", "assembly": "GRCh38", "source": "Ensembl"}
    r = Release.model_validate(raw)
    assert r.shortname == "110"
    assert r.assembly == "GRCh38"
