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
