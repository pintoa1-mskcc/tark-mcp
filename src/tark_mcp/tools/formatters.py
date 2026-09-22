from __future__ import annotations


def _fmt_row(row: list, widths: list[int]) -> str:
    return "  ".join(str(v).ljust(w) for v, w in zip(row, widths))


def format_transcripts_table(
    stable_ids: list[str],
    assemblies: list[str],
    results: list[dict | list[dict] | None],
    mane_lookup: dict[str, str] | None = None,
    notes: list[str] | None = None,
) -> str:
    """Format tark_get_transcripts results as a human-readable summary table.

    Columns: Query, Assembly, Stable ID, Ver, Exons, 5'UTR, 3'UTR, CDS (bp),
             AA Len, First Release, Latest Release, Release Date, MANE, Note.

    mane_lookup: optional dict mapping versioned stable ID → MANE type string,
                 e.g. {'ENST00000380152.8': 'MANE Select', 'NM_000059.4': 'MANE Select'}.
                 Versioned because MANE pairs one specific version of each transcript.
    notes: optional per-row warning strings (e.g. flagging that an unversioned query
           resolved to the latest version while an earlier, materially different
           version also exists — see versions.py's describe_divergent_earlier_version).
    """
    COL_HEADERS = [
        "Query", "Assembly", "Stable ID", "Ver", "Exons",
        "5'UTR", "3'UTR", "CDS (bp)", "AA Len",
        "First Release", "Latest Release", "Release Date", "MANE", "Note",
    ]
    COL_WIDTHS = [24, 10, 20, 5, 7, 8, 8, 10, 8, 22, 22, 14, 20, 40]

    rows: list[list] = []
    for i, (query, assembly, result) in enumerate(zip(stable_ids, assemblies, results)):
        note = notes[i] if notes and i < len(notes) else ""

        info: dict | None = None
        if isinstance(result, list):
            info = result[0] if result else None
        else:
            info = result

        if info is None:
            rows.append([query, assembly, "NOT FOUND", "", "", "", "", "", "", "", "", "", "", note])
            continue

        exon_count = len(info.get("exons") or [])
        cds_seq = info.get("cds_seq") or ""
        cds_len: int | str = len(cds_seq) if cds_seq else "N/A"
        aa_len: int | str = (len(cds_seq) // 3) - 1 if cds_seq else "N/A"

        rel_str = info.get("latest_release_version") or ""
        parts = [p.strip() for p in rel_str.split(",") if p.strip()]
        first_rel = parts[0] if parts else "N/A"
        last_rel = parts[-1] if parts else "N/A"
        release_date = info.get("latest_release_date") or "N/A"

        mane_status = ""
        if mane_lookup:
            sid = f"{info.get('stable_id', '')}.{info.get('stable_id_version', '')}"
            mane_status = mane_lookup.get(sid, "")

        rows.append([
            query,
            assembly,
            info.get("stable_id", ""),
            info.get("stable_id_version", ""),
            exon_count,
            info.get("five_prime_utr_length", "N/A"),
            info.get("three_prime_utr_length", "N/A"),
            cds_len,
            aa_len,
            first_rel,
            last_rel,
            release_date,
            mane_status,
            note,
        ])

    header_line = _fmt_row(COL_HEADERS, COL_WIDTHS)
    separator = "-" * (sum(COL_WIDTHS) + 2 * len(COL_WIDTHS))
    data_lines = [_fmt_row(row, COL_WIDTHS) for row in rows]
    return "\n".join([header_line, separator, *data_lines])
