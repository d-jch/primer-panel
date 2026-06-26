"""Tests for target_grouper_v2 sliding-window grouping strategy."""

import pytest
from primer_panel.target_grouper_v2 import (
    RequiredRegion,
    TargetGroup,
    group_targets,
    to_v1_targets,
)


def _r(start: int, end: int, chrom: str = "chr1", strand: int = 1,
       exon_nums: list[int] | None = None) -> RequiredRegion:
    """Helper to create a RequiredRegion."""
    return RequiredRegion(
        chrom=chrom,
        start=start,
        end=end,
        strand=strand,
        cds_exon_numbers=exon_nums or [],
    )


# ──────────────────────────────────────────────────────────────────────
# Basic grouping
# ──────────────────────────────────────────────────────────────────────

class TestBasicGrouping:
    """Core grouping scenarios."""

    def test_single_region(self):
        """Single region → one target, extended to product_min."""
        regions = [_r(1000, 1200)]
        targets = group_targets(regions, product_min=2700, product_max=3300)
        assert len(targets) == 1
        t = targets[0]
        assert t.required_start == 1000
        assert t.required_end == 1200
        assert t.extended_length == 2700
        assert t.extended_start == 1000
        assert t.extended_end == 3700

    def test_two_regions_close(self):
        """Two close regions merged by product_min."""
        # R1: 1000-1200, R2: 2000-2200
        # product_min=2700 → window_end = 1000+2700 = 3700
        # R2.end=2200 ≤ 3700 → covered by product_min
        regions = [_r(1000, 1200), _r(2000, 2200)]
        targets = group_targets(regions, product_min=2700, product_max=3300)
        assert len(targets) == 1
        t = targets[0]
        assert t.required_start == 1000
        assert t.required_end == 2200
        assert len(t.required_regions) == 2
        assert t.extended_length == 2700

    def test_two_regions_far(self):
        """Two far regions → separate targets."""
        # R1: 1000-1200, R2: 10000-10200
        # product_min window = 3700, doesn't reach R2
        # product_max window = 4300, doesn't reach R2
        regions = [_r(1000, 1200), _r(10000, 10200)]
        targets = group_targets(regions, product_min=2700, product_max=3300)
        assert len(targets) == 2

    def test_three_regions_min_covers_two(self):
        """product_min covers first two, not the third."""
        # R1: 1000-1200, R2: 2000-2200, R3: 8000-8200
        # product_min window=3700 covers R2 (end=2200≤3700), not R3 (end=8200)
        # product_max window=4300 also doesn't reach R3
        regions = [_r(1000, 1200), _r(2000, 2200), _r(8000, 8200)]
        targets = group_targets(regions, product_min=2700, product_max=3300)
        assert len(targets) == 2
        assert targets[0].required_end == 2200
        assert targets[1].required_start == 8000


# ──────────────────────────────────────────────────────────────────────
# product_min vs product_max fallback
# ──────────────────────────────────────────────────────────────────────

class TestMinMaxFallback:
    """product_min/product_max fallback logic."""

    def test_min_insufficient_max_sufficient(self):
        """product_min doesn't cover next region, product_max does."""
        # R1: 1000-1200, R2: 3500-3700
        # product_min=2700 → window_end=3700, R2.end=3700 → covered!
        # Actually R2.end=3700 == 3700, so product_min covers it.
        # Let me adjust: R2: 3500-3800 → R2.end=3800 > 3700 (min fails)
        # product_max=3300 → window_end=4300, R2.end=3800 ≤ 4300 (max works)
        regions = [_r(1000, 1200), _r(3500, 3800)]
        targets = group_targets(regions, product_min=2700, product_max=3300)
        assert len(targets) == 1
        assert targets[0].extended_length == 3300
        assert len(targets[0].required_regions) == 2

    def test_min_sufficient_uses_min(self):
        """product_min covers next region → use product_min, not product_max."""
        regions = [_r(1000, 1200), _r(2000, 2200)]
        targets = group_targets(regions, product_min=2700, product_max=3300)
        assert len(targets) == 1
        assert targets[0].extended_length == 2700  # not 3300


# ──────────────────────────────────────────────────────────────────────
# Last region backward expansion
# ──────────────────────────────────────────────────────────────────────

class TestBackwardExpansion:
    """Last uncovered region extends leftward."""

    def test_last_region_backward(self):
        """Last region not covered → extend leftward."""
        # R1: 1000-1200, R2: 50000-50200
        # Forward from R1: product_min=2700 → end=3700, doesn't reach R2
        # R2 is last uncovered → extend leftward from R2.end
        regions = [_r(1000, 1200), _r(50000, 50200)]
        targets = group_targets(regions, product_min=2700, product_max=3300)
        assert len(targets) == 2
        t2 = targets[1]
        assert t2.required_start == 50000
        assert t2.required_end == 50200
        assert t2.extended_end == 50200
        assert t2.extended_length == 2700
        assert t2.extended_start == 50200 - 2700

    def test_single_region_near_chr0_clamp(self):
        """Single region near chr0 → clamp start to 0, compensate rightward."""
        regions = [_r(100, 300)]
        targets = group_targets(regions, product_min=2700, product_max=3300)
        assert len(targets) == 1
        t = targets[0]
        assert t.extended_start >= 0
        assert t.extended_length == 2700

    def test_last_region_backward_absorbs_previous(self):
        """Last region backward expansion absorbs a previous region."""
        # R1: 48000-48200, R2: 50000-50200
        # Forward from R1: product_min=2700 → end=50700
        # R2.end=50200 ≤ 50700 → R2 is covered by forward expansion!
        # So this should be 1 target.
        regions = [_r(48000, 48200), _r(50000, 50200)]
        targets = group_targets(regions, product_min=2700, product_max=3300)
        assert len(targets) == 1


# ──────────────────────────────────────────────────────────────────────
# Tiling
# ──────────────────────────────────────────────────────────────────────

class TestTiling:
    """Oversized single regions are tiled."""

    def test_single_region_exceeds_product_max(self):
        """Region > product_max → tiled."""
        regions = [_r(1000, 5000)]  # 4000bp > 3300
        targets = group_targets(regions, product_min=2700, product_max=3300)
        assert len(targets) > 1
        assert all(t.tiled for t in targets)
        assert all(t.status == "tiled" for t in targets)

    def test_tile_overlap(self):
        """Tiled targets should have overlap."""
        regions = [_r(0, 10000)]
        targets = group_targets(regions, product_min=2700, product_max=3300, tile_overlap=200)
        for i in range(1, len(targets)):
            prev_end = targets[i - 1].required_end
            curr_start = targets[i].required_start
            assert curr_start < prev_end  # overlap exists


# ──────────────────────────────────────────────────────────────────────
# QC
# ──────────────────────────────────────────────────────────────────────

class TestQC:
    """Quality checks on targets."""

    def test_extended_covers_required(self):
        """extended_target must cover required_region."""
        regions = [_r(1000, 1200), _r(2000, 2200)]
        targets = group_targets(regions, product_min=2700, product_max=3300)
        for t in targets:
            assert t.extended_start <= t.required_start
            assert t.extended_end >= t.required_end

    def test_cds_exon_numbers_preserved(self):
        """CDS exon numbers must be preserved through grouping."""
        regions = [
            _r(1000, 1200, exon_nums=[1]),
            _r(2000, 2200, exon_nums=[2]),
        ]
        targets = group_targets(regions, product_min=2700, product_max=3300)
        assert targets[0].cds_exon_numbers == [1, 2]


# ──────────────────────────────────────────────────────────────────────
# Sorting validation
# ──────────────────────────────────────────────────────────────────────

class TestSorting:
    """Input must be sorted."""

    def test_unsorted_raises(self):
        regions = [_r(2000, 2200), _r(1000, 1200)]
        with pytest.raises(ValueError, match="sorted"):
            group_targets(regions, product_min=2700, product_max=3300)


# ──────────────────────────────────────────────────────────────────────
# v1 backward compatibility
# ──────────────────────────────────────────────────────────────────────

class TestV1Compat:
    """Conversion to v1 Target format."""

    def test_to_v1_targets(self):
        regions = [_r(1000, 1200, exon_nums=[1])]
        groups = group_targets(regions, product_min=2700, product_max=3300)
        v1 = to_v1_targets(groups)
        assert len(v1) == 1
        assert v1[0].start == 1000
        assert v1[0].end == 1200
        assert v1[0].cds_exon_numbers == [1]


# ──────────────────────────────────────────────────────────────────────
# Edge cases
# ──────────────────────────────────────────────────────────────────────

class TestEdgeCases:
    """Edge cases."""

    def test_empty_input(self):
        assert group_targets([], product_min=2700, product_max=3300) == []

    def test_region_equals_product_min(self):
        """Region exactly equals product_min → one target, no extension needed."""
        regions = [_r(1000, 3700)]  # 2700bp = product_min
        targets = group_targets(regions, product_min=2700, product_max=3300)
        assert len(targets) == 1
        assert targets[0].extended_length == 2700

    def test_region_at_chrom_start(self):
        """Region starts at 0 → clamp."""
        regions = [_r(0, 200)]
        targets = group_targets(regions, product_min=2700, product_max=3300)
        assert len(targets) == 1
        assert targets[0].extended_start >= 0


# ──────────────────────────────────────────────────────────────────────
# Strand-aware grouping
# ──────────────────────────────────────────────────────────────────────

class TestStrandAware:
    """-strand genes should sweep from right (3') to left (5')."""

    def test_minus_strand_sweeps_from_right(self):
        """-strand: start from rightmost region, extend leftward."""
        # Simulates FTH1: 3 intervals, -strand
        # cds1 (rightmost): 6000-6200
        # cds2 (middle):    4000-4200
        # cds3 (leftmost):  2000-2200
        regions = [
            _r(2000, 2200, strand=-1, exon_nums=[3]),
            _r(4000, 4200, strand=-1, exon_nums=[2]),
            _r(6000, 6200, strand=-1, exon_nums=[1]),
        ]
        targets = group_targets(regions, product_min=2700, product_max=3300)
        # -strand: reverse → [cds1(6000), cds2(4000), cds3(2000)]
        # sweep from cds1: window [6000, 6000+2700]=[6000,8700]
        # cds2.end=4200 ≤ 8700 ✓, cds3.end=2200 ≤ 8700 ✓
        # → all 3 in one target
        assert len(targets) == 1
        t = targets[0]
        assert t.required_start == 2000
        assert t.required_end == 6200
        assert t.strand == -1
        assert t.cds_exon_numbers == [1, 2, 3]

    def test_minus_strand_genomic_order_preserved(self):
        """-strand targets should be in genomic order (left to right)."""
        regions = [
            _r(2000, 2200, strand=-1, exon_nums=[3]),
            _r(4000, 4200, strand=-1, exon_nums=[2]),
            _r(6000, 6200, strand=-1, exon_nums=[1]),
        ]
        targets = group_targets(regions, product_min=2700, product_max=3300)
        # All in one target, genomic order
        assert targets[0].required_start == 2000
        assert targets[0].required_end == 6200
        # Required regions within target should be in genomic order
        starts = [r.start for r in targets[0].required_regions]
        assert starts == sorted(starts)

    def test_minus_strand_wide_gene_splits_correctly(self):
        """-strand gene wider than product_max splits correctly."""
        # cds1 (right): 10000-10200
        # cds2 (left):   1000-1200
        regions = [
            _r(1000, 1200, strand=-1, exon_nums=[2]),
            _r(10000, 10200, strand=-1, exon_nums=[1]),
        ]
        targets = group_targets(regions, product_min=2700, product_max=3300)
        # -strand: reverse → [cds1(10000), cds2(1000)]
        # sweep from cds1: window [10000, 12700]
        # cds2.end=1200 ≤ 12700 ✓ → both in one target
        # BUT required span = 10200-1000 = 9200bp > product_max
        # So they should be separate targets
        assert len(targets) == 2
        # Genomic order: left target first
        assert targets[0].required_start == 1000
        assert targets[1].required_start == 10000

    def test_minus_strand_fth1_scenario(self):
        """FTH1-like scenario: -strand, all CDS within product_min span."""
        # Simulates FTH1 without buffer
        # cds4+cds3 (left): 61964726-61965112
        # cds2 (middle):    61965368-61965515
        # cds1 (right):     61967311-61967425
        regions = [
            _r(61964726, 61965112, strand=-1, exon_nums=[3, 4]),
            _r(61965368, 61965515, strand=-1, exon_nums=[2]),
            _r(61967311, 61967425, strand=-1, exon_nums=[1]),
        ]
        targets = group_targets(regions, product_min=2700, product_max=3300)
        # -strand: reverse → [cds1(61967311), cds2(61965368), cds3+4(61964726)]
        # sweep from cds1: window [61967311, 61967311+2700]=[61967311, 6200011]
        # cds2.end=61965515 ≤ 6200011 ✓
        # cds3+4.end=61965112 ≤ 6200011 ✓
        # → all in one target
        assert len(targets) == 1
        t = targets[0]
        assert t.cds_exon_numbers == [1, 2, 3, 4]
        assert t.required_start == 61964726
        assert t.required_end == 61967425
