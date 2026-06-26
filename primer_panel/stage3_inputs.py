"""Helpers for building Stage 3 in-silico PCR inputs."""

from __future__ import annotations

from .writers import PrimerRecord, TargetRecord


def build_stage3_inputs(
    target_records: list[TargetRecord],
    primer_records: list[PrimerRecord],
) -> tuple[list[dict], dict[str, tuple[str, int, int]]]:
    """Build isPcr primer batch and expected target coordinates.

    ``expected_start`` / ``expected_end`` are genomic product coordinates
    computed from template-relative Primer3 coordinates.
    """
    targets_by_id = {record.target_id: record for record in target_records}
    expected_coords = {
        record.target_id: (
            record.extended_chrom,
            record.extended_start,
            record.extended_end,
        )
        for record in target_records
    }

    primer_batch: list[dict] = []
    for primer in primer_records:
        if primer.primer3_status != "ok":
            continue

        target = targets_by_id.get(primer.target_id)
        if target is None:
            continue

        primer_batch.append({
            "name": f"{primer.target_id}_rank{primer.primer_rank}",
            "fwd": primer.forward_primer,
            "rev": primer.reverse_primer,
            "expected_chrom": target.template_chrom,
            "expected_start": target.template_start + primer.primer_left_start,
            "expected_end": target.template_start + primer.primer_right_start + 1,
        })

    return primer_batch, expected_coords


def build_stage3_inputs_from_target_coords(
    target_coords: dict[str, dict],
    primer_records: list[PrimerRecord],
) -> tuple[list[dict], dict[str, tuple[str, int, int]]]:
    """Build Stage 3 inputs from target_summary.tsv coordinate dicts."""
    expected_coords = {
        target_id: (
            coords["extended_chrom"],
            coords["extended_start"],
            coords["extended_end"],
        )
        for target_id, coords in target_coords.items()
    }

    primer_batch: list[dict] = []
    for primer in primer_records:
        if primer.primer3_status != "ok":
            continue

        target = target_coords.get(primer.target_id, {})

        template_start = target.get("template_start", 0)
        primer_batch.append({
            "name": f"{primer.target_id}_rank{primer.primer_rank}",
            "fwd": primer.forward_primer,
            "rev": primer.reverse_primer,
            "expected_chrom": target.get("template_chrom", ""),
            "expected_start": template_start + primer.primer_left_start,
            "expected_end": template_start + primer.primer_right_start + 1,
        })

    return primer_batch, expected_coords
