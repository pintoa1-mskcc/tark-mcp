from __future__ import annotations
from pydantic import BaseModel, model_validator


class Exon(BaseModel):
    stable_id: str
    stable_id_version: int
    transcript_stable_id: str | None = None
    transcript_stable_id_version: int | None = None
    assembly: str
    order: int
    loc_region: str
    loc_start: int
    loc_end: int
    loc_strand: int

    @model_validator(mode="before")
    @classmethod
    def _normalize(cls, data: dict) -> dict:
        data = dict(data)
        if "loc_start" in data and data["loc_start"] is not None:
            data["loc_start"] = data["loc_start"] - 1
        if "exon_order" in data:
            data["order"] = data.pop("exon_order")
        return data


class Gene(BaseModel):
    stable_id: str
    stable_id_version: int
    name: str | None
    loc_region: str
    loc_start: int
    loc_end: int
    loc_strand: int
    assembly: str

    @model_validator(mode="before")
    @classmethod
    def _normalize(cls, data: dict) -> dict:
        data = dict(data)
        if "loc_start" in data and data["loc_start"] is not None:
            data["loc_start"] = data["loc_start"] - 1
        return data


class Translation(BaseModel):
    stable_id: str
    stable_id_version: int
    transcript_stable_id: str | None = None
    transcript_stable_id_version: int | None = None
    assembly: str
    sequence: str | None = None

    @model_validator(mode="before")
    @classmethod
    def _normalize(cls, data: dict) -> dict:
        data = dict(data)
        if isinstance(data.get("assembly"), dict):
            data["assembly"] = data["assembly"]["assembly_name"]
        seq = data.get("sequence")
        if isinstance(seq, dict):
            data["sequence"] = seq.get("sequence")
        return data


class Transcript(BaseModel):
    stable_id: str
    stable_id_version: int
    assembly: str
    biotype: str
    loc_region: str
    loc_start: int
    loc_end: int
    loc_strand: int
    cds_start: int | None
    cds_end: int | None
    five_prime_utr_length: int | None = None
    three_prime_utr_length: int | None = None
    cds_seq: str | None = None
    exons: list[Exon]
    genes: list[Gene]
    translations: list[Translation]
    sequence: str | None
    latest_release_date: str | None = None
    latest_release_version: str | None = None

    @model_validator(mode="before")
    @classmethod
    def _normalize(cls, data: dict) -> dict:
        data = dict(data)

        # Flatten nested assembly object
        if isinstance(data.get("assembly"), dict):
            data["assembly"] = data["assembly"]["assembly_name"]

        # Convert 1-based start → 0-based
        if "loc_start" in data and data["loc_start"] is not None:
            data["loc_start"] = data["loc_start"] - 1

        # Flatten nested sequence object FIRST, then use for CDS computation
        seq = data.get("sequence")
        if isinstance(seq, dict):
            data["sequence"] = seq.get("sequence")

        # Compute CDS boundaries from UTR sequences
        five_utr = data.get("five_prime_utr_seq")
        three_utr = data.get("three_prime_utr_seq")
        transcript_seq = data.get("sequence")  # now a string (or None)

        if five_utr is not None and transcript_seq is not None:
            data["cds_start"] = len(five_utr)
            data["cds_end"] = len(transcript_seq) - (len(three_utr) if three_utr else 0)
            data["five_prime_utr_length"] = len(five_utr)
            data["three_prime_utr_length"] = len(three_utr) if three_utr else 0
        else:
            data.setdefault("cds_start", None)
            data.setdefault("cds_end", None)

        # Pull cds_seq and UTR lengths directly from TARK's cds_info when available
        cds_info = data.get("cds_info")
        if isinstance(cds_info, dict):
            data["cds_seq"] = cds_info.get("cds_seq")
            if cds_info.get("five_prime_utr_length") is not None:
                data["five_prime_utr_length"] = cds_info["five_prime_utr_length"]
            if cds_info.get("three_prime_utr_length") is not None:
                data["three_prime_utr_length"] = cds_info["three_prime_utr_length"]

        # Compute latest_release_date from transcript_release_set array
        release_set = data.get("transcript_release_set", [])
        if isinstance(release_set, list) and release_set:
            dates = [r["release_date"] for r in release_set if r.get("release_date")]
            data["latest_release_date"] = max(dates) if dates else None
        elif isinstance(release_set, dict):
            data["latest_release_date"] = release_set.get("release_date")
        else:
            data["latest_release_date"] = None

        # Compute latest_release_version from transcript_release_set array
        release_set = data.get("transcript_release_set", [])
        if isinstance(release_set, list) and release_set:
            versions = [
                f"{r['source']} v{r['shortname']}"
                for r in release_set
                if r.get('source') and r.get('shortname')
            ]
            data["latest_release_version"] = ", ".join(versions) if versions else None
        elif isinstance(release_set, dict):
            source = release_set.get('source')
            shortname = release_set.get('shortname')
            if source and shortname:
                data["latest_release_version"] = f"{source} v{shortname}"
            else:
                data["latest_release_version"] = None
        else:
            data["latest_release_version"] = None

        return data


class Release(BaseModel):
    shortname: str
    description: str | None
    release_date: str
    assembly: str
    source: str


class ExonDiff(BaseModel):
    order: int
    change: str  # "added", "removed", "modified", "unchanged"
    ref_coords: tuple[int, int] | None
    candidate_coords: tuple[int, int] | None


class TranscriptDiff(BaseModel):
    reference_stable_id: str
    candidate_stable_id: str
    reference_assembly: str
    candidate_assembly: str
    biotype_changed: bool
    cds_changed: bool
    exon_count_changed: bool
    sequence_changed: bool
    exon_diffs: list[ExonDiff]
    reference_protein_coding: bool
    candidate_protein_coding: bool
    cds_sequence_changed: bool | None
    ref_cds_sequence: str | None
    candidate_cds_sequence: str | None
    protein_sequence_changed: bool | None
    ref_protein_sequence: str | None
    candidate_protein_sequence: str | None
