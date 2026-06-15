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


async def _resolve_transcript(
    stable_id: str,
    assembly: str,
    client: TarkClient,
) -> Transcript:
    """Fetch a single Transcript by stable ID, handling version suffixes and RefSeq IDs."""
    sid, version = stable_id, None
    if "." in stable_id:
        parts = stable_id.rsplit(".", 1)
        try:
            version = int(parts[1])
            sid = parts[0]
        except ValueError:
            pass

    params: dict = {"stable_id": sid, "expand_all": "true", "assembly_name": assembly}
    if version is not None:
        params["stable_id_version"] = version

    data = await client.get("transcript/", params=params)
    if not data:
        raise ValueError(
            f"Transcript not found: {stable_id} (assembly={assembly})"
        )

    transcripts = [Transcript.model_validate(r) for r in data]

    if version is not None:
        matching = [t for t in transcripts if t.stable_id_version == version]
        if matching:
            return matching[0]
        raise ValueError(
            f"Transcript not found: {stable_id} (assembly={assembly})"
        )

    return max(
        transcripts,
        key=lambda t: (t.latest_release_date or "", t.stable_id_version),
    )


async def _fetch_diff_pair(
    ref: Transcript,
    candidate_stable_id: str,
    candidate_assembly: str,
    client: TarkClient,
) -> TranscriptDiff:
    candidate = await _resolve_transcript(candidate_stable_id, candidate_assembly, client)
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

    # Resolve the reference transcript once so all candidate pairs share it
    # without duplicate fetches (concurrent pairs would otherwise race on the cache).
    ref_transcript = await _resolve_transcript(ref_id, ref_assembly, client)

    results = await asyncio.gather(*[
        _fetch_diff_pair(ref_transcript, stable_ids[i], resolved_assemblies[i], client)
        for i in range(1, len(stable_ids))
    ])
    return list(results)
