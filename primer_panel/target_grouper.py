"""Target grouping: CDS-based packing with product-size constraints.

The grouper produces targets whose start/end represent the *required_region*
(CDS exons + cds_buffer).  Short required_regions are NOT expanded here —
downstream in writers.py, they are extended toward gene interior to
product_min (→ extended_target) and then padded with primer_flank
(→ design_template) before being given to Primer3.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .config import PipelineConfig
from .cds_handler import CdsRequiredInterval
from .ensembl_client import CdsExon


@dataclass
class Target:
    """A PCR target region covering one or more CDS required intervals.

    ``start`` / ``end`` define the **required_region** — the genomic interval
    that must be covered by any valid PCR product.  Downstream in writers.py,
    short required_regions are extended toward gene interior (→ extended_target)
    and then padded with primer_flank (→ design_template) for Primer3.
    """

    chrom: str
    start: int
    end: int
    strand: int
    needs_review: bool = False
    status: str = "ok"
    tiled: bool = False                     # True if created by tiling an oversized interval
    cds_exon_numbers: list[int] = field(default_factory=list)
    cds_exon_ids: list[str] = field(default_factory=list)
    cds_exon_coords: list[tuple[int, int]] = field(default_factory=list)  # genomic (start, end) per CDS exon

    @property
    def length(self) -> int:
        return self.end - self.start


# ──────────────────────────────────────────────────────────────────────
# Gap finding (for splitting oversized groups)
# ──────────────────────────────────────────────────────────────────────

def _find_largest_gap(
    intervals: list[CdsRequiredInterval],
    indices: list[int],
) -> tuple[int, int]:
    """Find the largest gap between consecutive intervals in *indices*.

    Returns (position_in_indices_list, gap_size).
    ``position`` is the index where the gap sits (between indices[pos] and indices[pos+1]).
    Returns (-1, 0) if fewer than 2 intervals.
    """
    if len(indices) < 2:
        return -1, 0

    max_gap = 0
    max_pos = -1
    for i in range(len(indices) - 1):
        a = intervals[indices[i]]
        b = intervals[indices[i + 1]]
        gap = b.start - a.end
        if gap > max_gap:
            max_gap = gap
            max_pos = i
    return max_pos, max_gap


# ──────────────────────────────────────────────────────────────────────
# CDS exon assignment helper
# ──────────────────────────────────────────────────────────────────────

def _assign_cds_exons(
    target: Target,
    required_intervals: list[CdsRequiredInterval],
    interval_indices: list[int],
) -> None:
    """Populate target.cds_exon_numbers, cds_exon_ids, cds_exon_coords.

    Uses (number, id) pairs to preserve correspondence.
    cds_exon_coords stores genomic (start, end) for each CDS exon so that
    downstream code can validate primer coverage.
    """
    seen: set[int] = set()
    triples: list[tuple[int, str, tuple[int, int]]] = []
    for idx in interval_indices:
        ri = required_intervals[idx]
        for cex in ri.cds_exons:
            if cex.cds_exon_number not in seen:
                seen.add(cex.cds_exon_number)
                triples.append((cex.cds_exon_number, cex.cds_exon_id, (cex.start, cex.end)))
    triples.sort(key=lambda t: t[0])
    target.cds_exon_numbers = [t[0] for t in triples]
    target.cds_exon_ids = [t[1] for t in triples]
    target.cds_exon_coords = [t[2] for t in triples]


# ──────────────────────────────────────────────────────────────────────
# Classification
# ──────────────────────────────────────────────────────────────────────

def _classify_target(t: Target, cfg: PipelineConfig) -> None:
    """Assign status based on required_region length vs product-size constraints.

    - required_length > product_max  → needs_review_exceeds_product_max
    - tiled targets                  → tiled (already flagged via t.tiled)
    - required_length <= product_max → ok
    - required_length < product_min  → ok (short is fine; Primer3 uses flanks)
    """
    if t.needs_review:
        return
    length = t.length
    if length > cfg.product_max:
        t.needs_review = True
        t.status = "needs_review_exceeds_product_max"
    elif t.tiled:
        t.status = "tiled"
    else:
        t.status = "ok"


# ──────────────────────────────────────────────────────────────────────
# Main grouping
# ──────────────────────────────────────────────────────────────────────

def group_targets(
    required_intervals: list[CdsRequiredInterval],
    cfg: PipelineConfig,
) -> list[Target]:
    """Group CDS required intervals into PCR targets.

    Algorithm:
      1. Greedy pack adjacent intervals while combined span ≤ product_max.
      2. For each group:
         a. If span > product_max  → split at the largest CDS gap.
         b. Otherwise → keep as-is (no expansion).
      3. Single intervals > product_max → tile with overlap.
    """
    if not required_intervals:
        return []

    targets: list[Target] = []

    # --- Step 1: Greedy pack adjacent intervals ---
    groups: list[list[int]] = []  # each group is a list of interval indices
    cur_indices: list[int] = []

    for idx, ri in enumerate(required_intervals):
        if not cur_indices:
            cur_indices = [idx]
            continue

        # Check if this interval is on the same chromosome
        prev_ri = required_intervals[cur_indices[-1]]
        if ri.chrom != prev_ri.chrom:
            groups.append(cur_indices)
            cur_indices = [idx]
            continue

        # Check combined span
        group_start = required_intervals[cur_indices[0]].start
        combined_end = max(required_intervals[cur_indices[-1]].end, ri.end)
        combined_span = combined_end - group_start

        if combined_span <= cfg.product_max:
            cur_indices.append(idx)
        else:
            groups.append(cur_indices)
            cur_indices = [idx]

    if cur_indices:
        groups.append(cur_indices)

    # --- Step 2: Process each group ---
    for group_indices in groups:
        _process_group(group_indices, required_intervals, cfg, targets)

    # --- Step 3: Assign CDS exons to targets ---
    # Use overlap-based assignment so tiled targets also get CDS exon numbers.
    for t in targets:
        covered_indices = []
        for idx, ri in enumerate(required_intervals):
            if ri.chrom != t.chrom:
                continue
            if max(ri.start, t.start) < min(ri.end, t.end):
                covered_indices.append(idx)
        _assign_cds_exons(t, required_intervals, covered_indices)

    return targets


def _process_group(
    group_indices: list[int],
    required_intervals: list[CdsRequiredInterval],
    cfg: PipelineConfig,
    targets: list[Target],
) -> None:
    """Process a group of adjacent required intervals into targets.

    No expansion is performed — the target's start/end equals the required
    region bounds.  Short required_regions are valid; Primer3 will use the
    flanking design_template to find primers.
    """

    group_start = required_intervals[group_indices[0]].start
    group_end = required_intervals[group_indices[-1]].end
    group_span = group_end - group_start
    chrom = required_intervals[group_indices[0]].chrom
    strand = required_intervals[group_indices[0]].strand

    # Case 1: Single interval exceeding product_max → tile
    if len(group_indices) == 1 and group_span > cfg.product_max:
        _tile_interval(group_indices[0], required_intervals, cfg, targets)
        return

    # Case 2: Group exceeds product_max → split at largest gap
    if group_span > cfg.product_max:
        pos, gap = _find_largest_gap(required_intervals, group_indices)
        if pos >= 0 and gap > 0:
            left_indices = group_indices[:pos + 1]
            right_indices = group_indices[pos + 1:]
            _process_group(left_indices, required_intervals, cfg, targets)
            _process_group(right_indices, required_intervals, cfg, targets)
            return
        else:
            # Cannot split (single interval or no gap) → tile
            if len(group_indices) == 1:
                _tile_interval(group_indices[0], required_intervals, cfg, targets)
            else:
                _tile_group(group_indices, required_intervals, cfg, targets)
            return

    # Case 3: Group span within product_max — keep as-is, no expansion
    t = Target(
        chrom=chrom,
        start=group_start,
        end=group_end,
        strand=strand,
    )
    _classify_target(t, cfg)
    targets.append(t)


def _tile_interval(
    interval_idx: int,
    required_intervals: list[CdsRequiredInterval],
    cfg: PipelineConfig,
    targets: list[Target],
) -> None:
    """Tile a single required interval that exceeds product_max."""
    ri = required_intervals[interval_idx]
    overlap = cfg.tile_overlap
    step = cfg.product_max - overlap

    pos = ri.start
    while pos < ri.end:
        end = min(pos + cfg.product_max, ri.end)

        t = Target(
            chrom=ri.chrom,
            start=pos,
            end=end,
            strand=ri.strand,
            tiled=True,
        )
        _classify_target(t, cfg)
        targets.append(t)

        # Advance by step (product_max - overlap)
        if end >= ri.end:
            break
        pos += step
        # Ensure we don't get stuck
        if pos <= targets[-1].start and end < ri.end:
            pos = end - overlap + 1


def _tile_group(
    group_indices: list[int],
    required_intervals: list[CdsRequiredInterval],
    cfg: PipelineConfig,
    targets: list[Target],
) -> None:
    """Tile a group of intervals that cannot be split at gaps.

    Treats the entire group span as one region to be tiled.
    """
    group_start = required_intervals[group_indices[0]].start
    group_end = required_intervals[group_indices[-1]].end
    chrom = required_intervals[group_indices[0]].chrom
    strand = required_intervals[group_indices[0]].strand
    overlap = cfg.tile_overlap
    step = cfg.product_max - overlap

    pos = group_start
    while pos < group_end:
        end = min(pos + cfg.product_max, group_end)

        t = Target(
            chrom=chrom,
            start=pos,
            end=end,
            strand=strand,
            tiled=True,
        )
        _classify_target(t, cfg)
        targets.append(t)

        if end >= group_end:
            break
        pos += step
        if pos <= targets[-1].start and end < group_end:
            pos = end - overlap + 1
