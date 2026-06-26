"""Adapter: bridge cds_handler output → primer-target-planner → v1 Target.

Converts CdsRequiredInterval → RequiredInterval, calls plan_targets(),
converts TargetWindow → target_grouper.Target for build_records() compatibility.
"""

from __future__ import annotations

import logging
from typing import Any

import primer_target_planner as ptp

from .cds_handler import CdsRequiredInterval
from .config import PipelineConfig

logger = logging.getLogger(__name__)


def _interval_to_required(ri: CdsRequiredInterval) -> ptp.RequiredInterval:
    """Convert CdsRequiredInterval to primer-target-planner RequiredInterval."""
    exon_nums = [c.cds_exon_number for c in ri.cds_exons]
    exon_ids = [c.cds_exon_id for c in ri.cds_exons]
    exon_coords = [(c.start, c.end) for c in ri.cds_exons]

    # Use first exon number as ID (stable, human-readable)
    if len(exon_nums) == 1:
        interval_id = f"cds{exon_nums[0]}"
    elif len(exon_nums) > 1:
        interval_id = f"cds{exon_nums[0]}_{exon_nums[-1]}"
    else:
        interval_id = f"iv_{ri.start}"

    return ptp.RequiredInterval(
        id=interval_id,
        start=ri.start,
        end=ri.end,
        metadata={
            "exon_numbers": exon_nums,
            "exon_ids": exon_ids,
            "exon_coords": exon_coords,
            "chrom": ri.chrom,
            "strand": ri.strand,
        },
    )


def _target_to_v1(tw: ptp.TargetWindow, interval_map: dict[str, CdsRequiredInterval],
                   chrom: str, strand: int) -> Any:
    """Convert TargetWindow to Target (for build_records)."""
    from .writers import Target

    # Collect CDS info from covered intervals
    all_exon_nums: list[int] = []
    all_exon_ids: list[str] = []
    all_exon_coords: list[tuple[int, int]] = []

    for iv_id in tw.covered_ids:
        ri = interval_map.get(iv_id)
        if ri is not None:
            for c in ri.cds_exons:
                all_exon_nums.append(c.cds_exon_number)
                all_exon_ids.append(c.cds_exon_id)
                all_exon_coords.append((c.start, c.end))

    # Deduplicate and sort
    seen_nums: set[int] = set()
    unique_nums: list[int] = []
    for n in all_exon_nums:
        if n not in seen_nums:
            seen_nums.add(n)
            unique_nums.append(n)
    unique_nums.sort()

    return Target(
        chrom=chrom,
        start=tw.start,
        end=tw.end,
        strand=strand,
        needs_review=False,
        status="ok",
        tiled=(tw.planning_mode == "tiled"),
        cds_exon_numbers=unique_nums,
        cds_exon_ids=all_exon_ids,
        cds_exon_coords=all_exon_coords,
    )


def plan_targets_with_external_planner(
    required_intervals: list[CdsRequiredInterval],
    cfg: PipelineConfig,
    gene_start: int | None = None,
    gene_end: int | None = None,
) -> list:
    """Plan targets using primer-target-planner package.

    Args:
        required_intervals: Output from cds_handler.build_required_intervals()
        cfg: Pipeline config (product_min, product_max, tile_overlap)
        gene_start: Gene-level start (from Ensembl gene span). If None, uses min of intervals.
        gene_end: Gene-level end (from Ensembl gene span). If None, uses max of intervals.

    Returns:
        List of v1 Target objects (from target_grouper module).
    """
    if not required_intervals:
        return []

    # Convert to planner format
    planner_intervals = [_interval_to_required(ri) for ri in required_intervals]
    interval_map = {iv.id: ri for iv, ri in zip(planner_intervals, required_intervals)}

    # Determine strand from first interval
    strand = required_intervals[0].strand
    strand_char = "+" if strand == 1 else "-"

    # Create planner config
    planner_cfg = ptp.PlannerConfig(
        product_min=cfg.product_min,
        product_max=cfg.product_max,
        strand=strand_char,
        tile_overlap=cfg.tile_overlap,
    )

    # Use gene-level bounds (from Ensembl gene span) for correct terminal_reverse logic
    bounds = ptp.PlanningBounds(
        start=gene_start if gene_start is not None else min(ri.start for ri in required_intervals),
        end=gene_end if gene_end is not None else max(ri.end for ri in required_intervals),
    )

    # Run planner
    logger.info("Running primer-target-planner (strand=%s, product_min=%d, product_max=%d)",
                strand_char, cfg.product_min, cfg.product_max)
    windows = ptp.plan_targets(planner_intervals, planner_cfg, bounds)

    logger.info("Planner produced %d targets", len(windows))

    # Convert to v1 Targets
    chrom = required_intervals[0].chrom
    targets = []
    for tw in windows:
        t = _target_to_v1(tw, interval_map, chrom, strand)
        targets.append(t)
        logger.debug("  %s: %d-%d (%dbp) mode=%s reason=%s",
                     tw.anchor_id, tw.start, tw.end, tw.length,
                     tw.planning_mode, tw.reason)

    return targets
