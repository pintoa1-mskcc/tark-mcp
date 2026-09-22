from tark_mcp.tools.formatters import format_transcripts_table
from tests.conftest import TRANSCRIPT_BRCA2_RAW


def _dump(raw: dict) -> dict:
    from tark_mcp.models import Transcript
    return Transcript.model_validate(raw).model_dump()


def test_note_column_present_and_populated():
    table = format_transcripts_table(
        ["ENST00000380152"], ["GRCh38"], [_dump(TRANSCRIPT_BRCA2_RAW)],
        notes=["v1 differs (first release 108)"],
    )
    assert "Note" in table
    assert "v1 differs (first release 108)" in table


def test_note_column_defaults_to_empty_when_not_provided():
    table = format_transcripts_table(["ENST00000380152"], ["GRCh38"], [_dump(TRANSCRIPT_BRCA2_RAW)])
    header, _, row = table.splitlines()
    assert "Note" in header
    assert row.rstrip().endswith("2023-04-01")  # MANE + Note columns both blank, trimmed off


def test_not_found_row_still_includes_note():
    table = format_transcripts_table(
        ["ENST00000999999"], ["GRCh38"], [None], notes=["some warning"],
    )
    assert "NOT FOUND" in table
    assert "some warning" in table


def test_mane_column_requires_matching_version():
    """MANE pairs are version-specific: ENST00000380152.8 <-> NM_000059.4. A record for
    another version of the same stable ID (e.g. NM_000059.3, or the fixture's .7) is not MANE."""
    lookup = {"ENST00000380152.8": "MANE Select", "NM_000059.4": "MANE Select"}
    v7 = _dump(TRANSCRIPT_BRCA2_RAW)
    v8 = {**v7, "stable_id_version": 8}
    table = format_transcripts_table(
        ["ENST00000380152.7", "ENST00000380152.8"], ["GRCh38", "GRCh38"], [v7, v8],
        mane_lookup=lookup,
    )
    _, _, row_v7, row_v8 = table.splitlines()
    assert "MANE Select" not in row_v7
    assert "MANE Select" in row_v8
