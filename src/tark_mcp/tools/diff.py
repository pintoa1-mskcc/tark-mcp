from __future__ import annotations
import asyncio
from tark_mcp.client import TarkClient
from tark_mcp.models import Transcript, ExonDiff, TranscriptDiff


def _compute_exon_diffs(ref_exons: list, candidate_exons: list) -> list[ExonDiff]:
    """Compare two exon lists by order. Returns ExonDiff per exon."""
    ref_map = {e.order: e for e in ref_exons}
    cand_map = {e.order: e for e in candidate_exons}
    all_orders = sorted(set(ref_map) | set(cand_map))
    diffs = []
    for order in all_orders:
        ref = ref_map.get(order)
        cand = cand_map.get(order)
        if ref is None:
            change = "added"
            ref_coords = None
            cand_coords = (cand.loc_start, cand.loc_end)
        elif cand is None:
            change = "removed"
            ref_coords = (ref.loc_start, ref.loc_end)
            cand_coords = None
        elif (ref.loc_start, ref.loc_end) != (cand.loc_start, cand.loc_end):
            change = "modified"
            ref_coords = (ref.loc_start, ref.loc_end)
            cand_coords = (cand.loc_start, cand.loc_end)
        else:
            change = "unchanged"
            ref_coords = (ref.loc_start, ref.loc_end)
            cand_coords = (cand.loc_start, cand.loc_end)
        diffs.append(ExonDiff(order=order, change=change,
                               ref_coords=ref_coords, candidate_coords=cand_coords))
    return diffs


def _extract_cds_sequence(t: Transcript) -> str | None:
    """Slice CDS from transcript sequence using cds_start/cds_end offsets."""
    if t.sequence is None or t.cds_start is None or t.cds_end is None:
        return None
    return t.sequence[t.cds_start:t.cds_end]


async def _fetch_protein_sequence(
    translation_stable_id: str,
    assembly: str,
    client: TarkClient,
) -> str | None:
    """Fetch protein sequence from /api/translation/ endpoint."""
    data = await client.get("translation/", {
        "stable_id": translation_stable_id,
        "expand_all": "true",
        "assembly_name": assembly,
    })
    if not data:
        return None
    raw = data[0]
    seq = raw.get("sequence")
    if isinstance(seq, dict):
        return seq.get("sequence")
    return seq


async def _build_diff(
    ref: Transcript,
    candidate: Transcript,
    client: TarkClient,
) -> TranscriptDiff:
    ref_cds_seq = _extract_cds_sequence(ref)
    cand_cds_seq = _extract_cds_sequence(candidate)

    ref_coding = ref_cds_seq is not None and bool(ref.translations)
    cand_coding = cand_cds_seq is not None and bool(candidate.translations)

    async def _get_prot(t: Transcript) -> str | None:
        if not t.translations:
            return None
        return await _fetch_protein_sequence(t.translations[0].stable_id, t.assembly, client)

    ref_protein_seq, cand_protein_seq = await asyncio.gather(
        _get_prot(ref), _get_prot(candidate)
    )

    cds_changed: bool | None = None
    if ref_cds_seq is not None and cand_cds_seq is not None:
        cds_changed = ref_cds_seq != cand_cds_seq

    protein_changed: bool | None = None
    if ref_protein_seq is not None and cand_protein_seq is not None:
        protein_changed = ref_protein_seq != cand_protein_seq

    exon_diffs = _compute_exon_diffs(ref.exons, candidate.exons)

    return TranscriptDiff(
        reference_stable_id=ref.stable_id,
        candidate_stable_id=candidate.stable_id,
        reference_assembly=ref.assembly,
        candidate_assembly=candidate.assembly,
        biotype_changed=ref.biotype != candidate.biotype,
        cds_changed=(ref.cds_start, ref.cds_end) != (candidate.cds_start, candidate.cds_end),
        exon_count_changed=len(ref.exons) != len(candidate.exons),
        sequence_changed=ref.sequence != candidate.sequence,
        exon_diffs=exon_diffs,
        reference_protein_coding=ref_coding,
        candidate_protein_coding=cand_coding,
        cds_sequence_changed=cds_changed,
        ref_cds_sequence=ref_cds_seq,
        candidate_cds_sequence=cand_cds_seq,
        protein_sequence_changed=protein_changed,
        ref_protein_sequence=ref_protein_seq,
        candidate_protein_sequence=cand_protein_seq,
    )


async def _fetch_diff_pair(
    ref_stable_id: str,
    candidate_stable_id: str,
    ref_assembly: str,
    candidate_assembly: str,
    client: TarkClient,
) -> TranscriptDiff:
    raw = await client.get_raw(
        "transcript/diff/",
        params={
            "diff_me_stable_id": ref_stable_id,
            "diff_with_stable_id": candidate_stable_id,
        },
    )
    ref_data = raw.get("diff_me_transcript", {})
    cand_data = raw.get("diff_with_transcript", {})

    if not isinstance(ref_data.get("assembly"), dict):
        ref_data = {**ref_data, "assembly": {"assembly_name": ref_assembly}}
    if not isinstance(cand_data.get("assembly"), dict):
        cand_data = {**cand_data, "assembly": {"assembly_name": candidate_assembly}}

    ref = Transcript.model_validate(ref_data)
    candidate = Transcript.model_validate(cand_data)
    return await _build_diff(ref, candidate, client)


async def diff_transcripts(
    stable_ids: list[str],
    assemblies: list[str] | None = None,
    client: TarkClient | None = None,
) -> list[TranscriptDiff]:
    if len(stable_ids) < 2:
        raise ValueError("At least 2 stable IDs required for diff")
    if client is None:
        client = TarkClient()

    resolved_assemblies = list(assemblies or [])
    while len(resolved_assemblies) < len(stable_ids):
        resolved_assemblies.append("GRCh38")

    ref_id = stable_ids[0]
    ref_assembly = resolved_assemblies[0]

    pairs = [
        (ref_id, stable_ids[i], ref_assembly, resolved_assemblies[i])
        for i in range(1, len(stable_ids))
    ]

    results = await asyncio.gather(*[
        _fetch_diff_pair(ref, cand, ra, ca, client)
        for ref, cand, ra, ca in pairs
    ])
    return list(results)
