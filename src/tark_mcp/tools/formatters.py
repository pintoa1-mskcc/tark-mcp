from __future__ import annotations


def _fmt_row(row: list, widths: list[int]) -> str:
    return "  ".join(str(v).ljust(w) for v, w in zip(row, widths))


def format_transcripts_table(
    stable_ids: list[str],
    assemblies: list[str],
    results: list[dict | list[dict] | None],
    mane_lookup: dict[str, str] | None = None,
) -> str:
    """Format tark_get_transcripts results as a human-readable summary table.

    Columns: Query, Assembly, Stable ID, Ver, Exons, 5'UTR, 3'UTR, CDS (bp),
             AA Len, First Release, Latest Release, Release Date, MANE.

    mane_lookup: optional dict mapping stable ID (no version) → MANE type string,
                 e.g. {'ENST00000380152': 'MANE SELECT', 'NM_024852': 'MANE SELECT'}
    """
    COL_HEADERS = [
        "Query", "Assembly", "Stable ID", "Ver", "Exons",
        "5'UTR", "3'UTR", "CDS (bp)", "AA Len",
        "First Release", "Latest Release", "Release Date", "MANE",
    ]
    COL_WIDTHS = [24, 10, 20, 5, 7, 8, 8, 10, 8, 22, 22, 14, 20]

    rows: list[list] = []
    for query, assembly, result in zip(stable_ids, assemblies, results):
        info: dict | None = None
        if isinstance(result, list):
            info = result[0] if result else None
        else:
            info = result

        if info is None:
            rows.append([query, assembly, "NOT FOUND", "", "", "", "", "", "", "", "", "", ""])
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
            sid = info.get("stable_id", "")
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
        ])

    header_line = _fmt_row(COL_HEADERS, COL_WIDTHS)
    separator = "-" * (sum(COL_WIDTHS) + 2 * len(COL_WIDTHS))
    data_lines = [_fmt_row(row, COL_WIDTHS) for row in rows]
    return "\n".join([header_line, separator, *data_lines])
