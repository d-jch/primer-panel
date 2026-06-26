"""Target grouping v2: strand-aware sliding-window strategy.

Groups CDS required_regions into PCR targets:

  +strand: sweep rightward from 5' end (leftmost genomic region).
  -strand: sweep leftward from 3' end (rightmost genomic region).

For each expansion step, try product_min first; if it doesn't cover the
next required_region, fall back to product_max.

This module is self-contained and does not depend on the v1 grouper.
"""

from __future__ import annotations

from dataclasses import dataclass, field


# ──────────────────────────────────────────────────────────────────────
# Data classes
# ──────────────────────────────────────────────────────────────────────

@dataclass(frozen=True)
class RequiredRegion:
    """A single CDS exon or merged CDS interval that must be covered."""

    chrom: str
    start: int
    end: int
    strand: int = 1
    cds_exon_numbers: list[int] = field(default_factory=list)
    cds_exon_ids: list[str] = field(default_factory=list)
    cds_exon_coords: list[tuple[int, int]] = field(default_factory=list)

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass
class TargetGroup:
    """A group of required_regions forming one PCR target."""

    chrom: str
    strand: int
    required_start: int
    required_end: int
    extended_start: int
    extended_end: int
    required_regions: list[RequiredRegion] = field(default_factory=list)
    tiled: bool = False
    needs_review: bool = False
    status: str = "ok"

    @property
    def required_length(self) -> int:
        return self.required_end - self.required_start

    @property
    def extended_length(self) -> int:
        return self.extended_end - self.extended_start

    @property
    def cds_exon_numbers(self) -> list[int]:
        nums: list[int] = []
        for r in self.required_regions:
            nums.extend(r.cds_exon_numbers)
        return sorted(set(nums))

    @property
    def cds_exon_ids(self) -> list[str]:
        ids: list[str] = []
        for r in self.required_regions:
            ids.extend(r.cds_exon_ids)
        return ids

    @property
    def cds_exon_coords(self) -> list[tuple[int, int]]:
        coords: list[tuple[int, int]] = []
        for r in self.required_regions:
            coords.extend(r.cds_exon_coords)
        return coords


# ──────────────────────────────────────────────────────────────────────
# Tiling
# ──────────────────────────────────────────────────────────────────────

def _tile_region(
    region: RequiredRegion,
    product_max: int,
    tile_overlap: int = 200,
) -> list[TargetGroup]:
    """Tile a single required_region that exceeds product_max."""
    targets: list[TargetGroup] = []
    step = product_max - tile_overlap
    pos = region.start

    while pos < region.end:
        end = min(pos + product_max, region.end)
        targets.append(TargetGroup(
            chrom=region.chrom,
            strand=region.strand,
            required_start=pos,
            required_end=end,
            extended_start=pos,
            extended_end=end,
            required_regions=[region],
            tiled=True,
            status="tiled",
        ))
        if end >= region.end:
            break
        pos += step

    return targets


# ──────────────────────────────────────────────────────────────────────
# +strand sweep (left → right)
# ──────────────────────────────────────────────────────────────────────

def _sweep_plus(
    regions: list[RequiredRegion],
    product_min: int,
    product_max: int,
    tile_overlap: int,
) -> list[TargetGroup]:
    """Sweep left-to-right for +strand genes.

    Regions must be sorted by start ascending.
    Anchor is fixed at the first region's start.  Window does not shift
    as regions are absorbed (no snowball effect).
    """
    targets: list[TargetGroup] = []
    n = len(regions)
    i = 0

    while i < n:
        r = regions[i]

        if r.length > product_max:
            targets.extend(_tile_region(r, product_max, tile_overlap))
            i += 1
            continue

        anchor = r.start
        window_min = anchor + product_min
        window_max = anchor + product_max

        # Absorb regions whose end <= window_max (fixed anchor, no snowball)
        j = i
        while j + 1 < n and regions[j + 1].end <= window_max:
            j += 1

        if j == i:
            # No region absorbed
            if i == n - 1 and n > 1:
                # Last region → backward sweep
                _add_backward_plus(targets, regions, i, product_min, product_max)
                break
            else:
                # Single region, extend to product_min
                _add_target(targets, [r], anchor, anchor + product_min)
                i += 1
                continue

        # Determine extension: product_min if all fit within window_min,
        # product_max otherwise
        all_within_min = all(regions[k].end <= window_min for k in range(i, j + 1))
        ext_end = window_min if all_within_min else window_max

        _add_target(targets, regions[i:j + 1], anchor, ext_end)

        if j == n - 1:
            break
        i = j + 1

    return targets


def _add_backward_plus(
    targets: list[TargetGroup],
    regions: list[RequiredRegion],
    idx: int,
    product_min: int,
    product_max: int,
) -> None:
    """Last region, +strand: extend leftward."""
    region = regions[idx]
    end = region.end

    start_min = end - product_min
    k_min = _count_backward(regions, idx - 1, start_min)

    start_max = end - product_max
    k_max = _count_backward(regions, idx - 1, start_max)

    if k_min >= 1:
        k = idx - k_min
        ext_start = start_min
    elif k_max >= 1:
        k = idx - k_max
        ext_start = start_max
    else:
        k = idx
        ext_start = max(0, end - product_min)

    ext_start = max(0, ext_start)
    ext_end = end
    if ext_end - ext_start < product_min:
        ext_end = ext_start + product_min

    _add_target(targets, regions[k:idx + 1], ext_start, ext_end)


def _count_backward(
    regions: list[RequiredRegion],
    start_idx: int,
    window_start: int,
) -> int:
    """Count consecutive regions backward whose start >= window_start."""
    count = 0
    for j in range(start_idx, -1, -1):
        if regions[j].start >= window_start:
            count += 1
        else:
            break
    return count


# ──────────────────────────────────────────────────────────────────────
# -strand sweep (right → left)
# ──────────────────────────────────────────────────────────────────────

def _sweep_minus(
    regions: list[RequiredRegion],
    product_min: int,
    product_max: int,
    tile_overlap: int,
) -> list[TargetGroup]:
    """Sweep right-to-left for -strand genes.

    Regions must be sorted by start ascending (genomic order).
    Processing starts from the rightmost region (3' end) and extends
    leftward (toward gene interior).

    Anchor is fixed at the first region's end (genomic rightmost).
    Window does not shift as regions are absorbed (no snowball effect).

    Absorb: region.end >= anchor - product_min (reachable by min extension)
    Terminate: region.end < anchor - product_max (unreachable even by max extension)
    Gap zone: absorb with product_max extension

    Returns targets in genomic order (left to right).
    """
    n = len(regions)
    raw_targets: list[TargetGroup] = []

    i = n - 1  # start from rightmost
    while i >= 0:
        r = regions[i]

        if r.length > product_max:
            raw_targets.extend(_tile_region(r, product_max, tile_overlap))
            i -= 1
            continue

        # Anchor fixed at this region's end (rightmost edge)
        anchor = r.end
        window_min = anchor - product_min
        window_max = anchor - product_max

        # Absorb regions to the left whose end >= window_max
        # (reachable by at least product_max extension)
        j = i
        while j > 0 and regions[j - 1].end >= window_max:
            j -= 1

        if j == i:
            # No region absorbed
            if i == 0 and n > 1:
                # First region (leftmost) → forward sweep
                _add_forward_minus(raw_targets, regions, i, product_min, product_max)
                break
            else:
                # Single region, extend leftward to product_min
                ext_start = max(0, anchor - product_min)
                _add_target(raw_targets, [r], ext_start, anchor)
                i -= 1
                continue

        # Determine extension
        absorbed = regions[j:i + 1]
        group_left = min(r.start for r in absorbed)

        # If all absorbed regions fit within window_min → use product_min
        # Otherwise → use product_max
        all_within_min = all(regions[k].end >= window_min for k in range(j, i + 1))
        ext_start = max(0, group_left - product_min) if all_within_min else max(0, anchor - product_max)

        _add_target(raw_targets, absorbed, ext_start, anchor)
        i = j - 1

    # Reverse to genomic order
    raw_targets.reverse()
    return raw_targets


def _add_forward_minus(
    targets: list[TargetGroup],
    regions: list[RequiredRegion],
    idx: int,
    product_min: int,
    product_max: int,
) -> None:
    """First region (leftmost), -strand: extend rightward."""
    region = regions[idx]
    start = region.start

    end_min = start + product_min
    k_min = 0
    for j in range(idx + 1, len(regions)):
        if regions[j].end <= end_min:
            k_min += 1
        else:
            break

    end_max = start + product_max
    k_max = 0
    for j in range(idx + 1, len(regions)):
        if regions[j].end <= end_max:
            k_max += 1
        else:
            break

    if k_min >= 1:
        k = idx + k_min
        ext_end = end_min
    elif k_max >= 1:
        k = idx + k_max
        ext_end = end_max
    else:
        k = idx
        ext_end = start + product_min

    _add_target(targets, regions[idx:k + 1], start, ext_end)


# ──────────────────────────────────────────────────────────────────────
# Shared helpers
# ──────────────────────────────────────────────────────────────────────

def _add_target(
    targets: list[TargetGroup],
    absorbed: list[RequiredRegion],
    ext_start: int,
    ext_end: int,
) -> None:
    """Add a target from absorbed regions with extended boundaries."""
    req_start = min(r.start for r in absorbed)
    req_end = max(r.end for r in absorbed)

    needs_review = False
    status = "ok"
    if ext_start > req_start:
        needs_review = True
        status = "extended_does_not_cover_required_start"
    elif ext_end < req_end:
        needs_review = True
        status = "extended_does_not_cover_required_end"

    targets.append(TargetGroup(
        chrom=absorbed[0].chrom,
        strand=absorbed[0].strand,
        required_start=req_start,
        required_end=req_end,
        extended_start=ext_start,
        extended_end=ext_end,
        required_regions=absorbed,
        needs_review=needs_review,
        status=status,
    ))


# ──────────────────────────────────────────────────────────────────────
# Main entry point
# ──────────────────────────────────────────────────────────────────────

def group_targets(
    required_regions: list[RequiredRegion],
    product_min: int,
    product_max: int,
    tile_overlap: int = 200,
) -> list[TargetGroup]:
    """Group required_regions into PCR targets using strand-aware sliding-window.

    +strand: sweep rightward from 5' end.
    -strand: sweep leftward from 3' end.

    Args:
        required_regions: Sorted list of RequiredRegion (by chrom, start).
        product_min: Minimum extended target length (bp).
        product_max: Maximum extended target length (bp).
        tile_overlap: Overlap for tiled targets (bp).

    Returns:
        List of TargetGroup in genomic order.
    """
    if not required_regions:
        return []

    # Validate sorting
    for i in range(1, len(required_regions)):
        prev = required_regions[i - 1]
        curr = required_regions[i]
        if curr.chrom < prev.chrom or (curr.chrom == prev.chrom and curr.start < prev.start):
            raise ValueError(
                f"required_regions must be sorted by (chrom, start); "
                f"got {prev.chrom}:{prev.start} before {curr.chrom}:{curr.start}"
            )

    # Split by chromosome
    all_targets: list[TargetGroup] = []
    chrom_groups: list[list[RequiredRegion]] = []
    current: list[RequiredRegion] = []

    for r in required_regions:
        if current and r.chrom != current[0].chrom:
            chrom_groups.append(current)
            current = []
        current.append(r)
    if current:
        chrom_groups.append(current)

    for chrom_regions in chrom_groups:
        strand = chrom_regions[0].strand

        if strand == -1:
            targets = _sweep_minus(chrom_regions, product_min, product_max, tile_overlap)
        else:
            targets = _sweep_plus(chrom_regions, product_min, product_max, tile_overlap)

        all_targets.extend(targets)

    return all_targets


# ──────────────────────────────────────────────────────────────────────
# Conversion to v1 Target (for backward compatibility)
# ──────────────────────────────────────────────────────────────────────

def to_v1_targets(groups: list[TargetGroup]) -> list:
    """Convert TargetGroup list to v1 Target list (from target_grouper module)."""
    from .target_grouper import Target

    v1_targets: list[Target] = []
    for g in groups:
        t = Target(
            chrom=g.chrom,
            start=g.required_start,
            end=g.required_end,
            strand=g.strand,
            needs_review=g.needs_review,
            status=g.status,
            tiled=g.tiled,
            cds_exon_numbers=g.cds_exon_numbers,
            cds_exon_ids=g.cds_exon_ids,
            cds_exon_coords=g.cds_exon_coords,
        )
        v1_targets.append(t)
    return v1_targets
