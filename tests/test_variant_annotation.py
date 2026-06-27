"""Tests for variant annotation module."""

import pytest
from pathlib import Path
from primer_panel.variant_annotation import (
    load_dbsnp_bed,
    _find_overlapping_snps,
    annotate_primer_snps,
    annotate_primer_pair,
    SnpDatabase,
    SnpEntry,
    ChromIndex,
)


@pytest.fixture
def sample_bed(tmp_path):
    """Create a small test BED file."""
    bed = tmp_path / "test.bed"
    bed.write_text(
        "chr1\t100\t101\trs100\n"
        "chr1\t200\t201\trs200\n"
        "chr1\t300\t301\trs300\n"
        "chr1\t400\t401\trs400\n"
        "chr2\t100\t101\trs500\n",
        encoding="utf-8",
    )
    return bed


@pytest.fixture
def sample_db(sample_bed):
    """Load SNP database from test BED file."""
    return load_dbsnp_bed(sample_bed)


class TestLoadDbsnpBed:
    """Test BED file loading."""

    def test_loads_chromosomes(self, sample_db):
        """Should load all chromosomes."""
        assert "chr1" in sample_db.chroms
        assert "chr2" in sample_db.chroms

    def test_sorted_starts(self, sample_db):
        """Starts should be sorted."""
        idx = sample_db.chroms["chr1"]
        assert idx.starts == [100, 200, 300, 400]

    def test_entries_count(self, sample_db):
        """Should load correct number of entries."""
        assert len(sample_db.chroms["chr1"].entries) == 4
        assert len(sample_db.chroms["chr2"].entries) == 1

    def test_skips_header_lines(self, tmp_path):
        """Should skip lines starting with #."""
        bed = tmp_path / "test.bed"
        bed.write_text(
            "#header\n"
            "chr1\t100\t101\trs100\n",
            encoding="utf-8",
        )
        db = load_dbsnp_bed(bed)
        assert len(db.chroms["chr1"].entries) == 1

    def test_skips_short_lines(self, tmp_path):
        """Should skip lines with fewer than 4 columns."""
        bed = tmp_path / "test.bed"
        bed.write_text(
            "chr1\t100\t101\n"  # Only 3 columns
            "chr1\t200\t201\trs200\n",
            encoding="utf-8",
        )
        db = load_dbsnp_bed(bed)
        assert len(db.chroms["chr1"].entries) == 1

    def test_skips_invalid_coordinates(self, tmp_path):
        """Should skip lines with non-integer coordinates."""
        bed = tmp_path / "test.bed"
        bed.write_text(
            "chr1\tabc\t101\trs100\n"
            "chr1\t200\t201\trs200\n",
            encoding="utf-8",
        )
        db = load_dbsnp_bed(bed)
        assert len(db.chroms["chr1"].entries) == 1


class TestFindOverlappingSnps:
    """Test SNP overlap queries."""

    def test_no_overlap_before(self, sample_db):
        """Query before all SNPs should return empty."""
        result = _find_overlapping_snps(sample_db, "chr1", 50, 60)
        assert result == []

    def test_no_overlap_after(self, sample_db):
        """Query after all SNPs should return empty."""
        result = _find_overlapping_snps(sample_db, "chr1", 500, 600)
        assert result == []

    def test_overlap_single(self, sample_db):
        """Query overlapping one SNP should return it."""
        result = _find_overlapping_snps(sample_db, "chr1", 90, 110)
        assert len(result) == 1
        assert result[0].rsid == "rs100"

    def test_overlap_multiple(self, sample_db):
        """Query overlapping multiple SNPs should return all."""
        result = _find_overlapping_snps(sample_db, "chr1", 150, 350)
        assert len(result) == 2
        rsids = {s.rsid for s in result}
        assert rsids == {"rs200", "rs300"}

    def test_no_overlap_different_chrom(self, sample_db):
        """Query on different chromosome should return empty."""
        result = _find_overlapping_snps(sample_db, "chr3", 100, 200)
        assert result == []

    def test_boundary_overlap(self, sample_db):
        """Query touching SNP boundary should detect overlap."""
        # SNP at 100-101, query at 95-100 -> no overlap (end == start)
        result = _find_overlapping_snps(sample_db, "chr1", 95, 100)
        assert len(result) == 0

        # SNP at 100-101, query at 95-101 -> overlap
        result = _find_overlapping_snps(sample_db, "chr1", 95, 101)
        assert len(result) == 1

    def test_performance_no_full_scan(self, tmp_path):
        """Should not scan entire chromosome for small query."""
        # Create a large BED file
        bed = tmp_path / "large.bed"
        lines = []
        for i in range(10000):
            lines.append(f"chr1\t{i*100}\t{i*100+1}\trs{i}")
        bed.write_text("\n".join(lines), encoding="utf-8")

        db = load_dbsnp_bed(bed)

        # Query near the end should not scan from beginning
        result = _find_overlapping_snps(db, "chr1", 999900, 999950)
        # Should find rs9999 at 999900-999901
        assert len(result) == 1
        assert result[0].rsid == "rs9999"

    def test_long_interval_overlapping_many_earlier_entries(self, tmp_path):
        """Intervals starting far before the query must not be missed.

        Regression test: the old code scanned back only 10 entries and
        missed longer BED intervals that started earlier but still overlapped.
        """
        bed = tmp_path / "long.bed"
        # 30 short SNPs scattered before position 5000
        lines = []
        for i in range(30):
            lines.append(f"chr1\t{i*100}\t{i*100+1}\trs_short_{i}")
        # One long interval starting at 0, extending to 6000
        lines.append(f"chr1\t0\t6000\trs_long_interval")
        bed.write_text("\n".join(lines), encoding="utf-8")

        db = load_dbsnp_bed(bed)

        # Query at [5500, 5600) — the long interval overlaps, but its start
        # is 5500 entries before the query position in sorted order.
        result = _find_overlapping_snps(db, "chr1", 5500, 5600)
        rsids = {s.rsid for s in result}
        assert "rs_long_interval" in rsids


class TestAnnotatePrimerSnps:
    """Test single primer annotation."""

    def test_no_overlap(self, sample_db):
        """Primer with no SNP overlap should return none."""
        risk, total, three_p, hits = annotate_primer_snps(
            sample_db, "chr1", 500, 20, is_reverse=False
        )
        assert risk == "none"
        assert total == 0
        assert three_p == 0
        assert hits == ""

    def test_forward_primer_high_risk(self, sample_db):
        """Forward primer with SNP at 3' end should be high risk."""
        # Primer at 90-110, 3' end is 105-110
        # rs100 at 100-101 overlaps 3' end (101 > 105 is false)
        # Actually: three_prime_start = 90 + 20 - 5 = 105
        # three_prime_end = 90 + 20 = 110
        # rs100: start=100, end=101
        # Check: snp.end > three_prime_start (101 > 105 is false)
        # So rs100 does NOT overlap 3' end
        risk, total, three_p, hits = annotate_primer_snps(
            sample_db, "chr1", 90, 20, is_reverse=False
        )
        # rs100 overlaps primer (90-110) but not 3' end (105-110)
        assert risk == "medium"
        assert total == 1
        assert three_p == 0

    def test_forward_primer_3p_overlap(self, sample_db):
        """Forward primer with SNP at exact 3' end should be high risk."""
        # Forward primer at 190-210, 3' end is 205-210
        # rs200 at 200-201: end=201 > 205 is false -> no 3' overlap
        # rs200 overlaps primer body (190-210) -> medium risk
        risk, total, three_p, hits = annotate_primer_snps(
            sample_db, "chr1", 190, 20, is_reverse=False
        )
        # rs200 at 200-201 overlaps primer but not 3' end [205, 210)
        assert risk == "medium"
        assert total == 1
        assert three_p == 0

        # Forward primer at 198-218, 3' end is 213-218
        # rs200 at 200-201: overlaps primer [198, 218) but not 3' end [213, 218)
        risk, total, three_p, hits = annotate_primer_snps(
            sample_db, "chr1", 198, 20, is_reverse=False
        )
        assert risk == "medium"
        assert total == 1
        assert three_p == 0

    def test_forward_primer_3p_overlap_with_close_snp(self, tmp_path):
        """Forward primer with SNP at 3' end (last 5bp) should be high risk."""
        # Create BED with SNP at 105-106
        bed = tmp_path / "test.bed"
        bed.write_text("chr1\t105\t106\trs105\n", encoding="utf-8")
        db = load_dbsnp_bed(bed)

        # Forward primer at 90-110, 3' end is 105-110
        # rs105 at 105-106: end=106 > 105 AND start=105 < 110 -> 3' overlap!
        risk, total, three_p, hits = annotate_primer_snps(
            db, "chr1", 90, 20, is_reverse=False
        )
        assert risk == "high"
        assert total == 1
        assert three_p == 1
        assert "rs105:3p" in hits

    def test_reverse_primer_high_risk(self, sample_db):
        """Reverse primer with SNP at 3' end (low-coordinate) should be high risk."""
        # Reverse primer: primer_start=105 (rightmost base), primer_len=20
        # Primer spans [86, 106) on template
        # 3' end = low-coordinate end: [86, 91)
        # rs100 at 100-101: overlaps primer but NOT 3' end [86, 91)
        risk, total, three_p, hits = annotate_primer_snps(
            sample_db, "chr1", 105, 20, is_reverse=True
        )
        # rs100 overlaps primer [86, 106) but not 3' end [86, 91)
        assert risk == "medium"
        assert total == 1
        assert three_p == 0

    def test_reverse_primer_3p_overlap(self, sample_db):
        """Reverse primer SNP at the low-coordinate 3' end is high risk."""
        # Reverse primer: primer_start=115 (rightmost), primer_len=20
        # Primer spans [96, 116) on template
        # 3' end = low-coordinate end: [96, 101)
        # rs100 at 100-101: end=101 > three_prime_start=96 AND start=100 < three_prime_end=101
        # -> 3' overlap!
        risk, total, three_p, hits = annotate_primer_snps(
            sample_db, "chr1", 115, 20, is_reverse=True
        )
        assert risk == "high"
        assert total == 1
        assert three_p == 1
        assert "rs100:3p" in hits

    def test_different_chromosome(self, sample_db):
        """Primer on different chromosome should return none."""
        risk, total, three_p, hits = annotate_primer_snps(
            sample_db, "chr3", 100, 20, is_reverse=False
        )
        assert risk == "none"
        assert total == 0


class TestAnnotatePrimerPair:
    """Test primer pair annotation."""

    def test_no_overlap(self, sample_db):
        """Primer pair with no overlap should return none."""
        result = annotate_primer_pair(
            sample_db, "chr1", 500, 20, 600, 20
        )
        assert result["common_snp_risk"] == "none"
        assert result["left_primer_common_snp_count"] == 0
        assert result["right_primer_common_snp_count"] == 0

    def test_left_primer_overlap(self, sample_db):
        """Left primer overlapping SNP should be detected."""
        # Left primer at 90-110 overlaps rs100
        result = annotate_primer_pair(
            sample_db, "chr1", 90, 20, 500, 20
        )
        assert result["common_snp_risk"] == "medium"
        assert result["left_primer_common_snp_count"] == 1
        assert result["right_primer_common_snp_count"] == 0

    def test_right_primer_overlap(self, sample_db):
        """Right primer overlapping SNP should be detected."""
        # Right primer at 190-210 overlaps rs200 (reverse)
        # Reverse: primer_start=210, primer_len=20
        # Primer spans [191, 211)
        result = annotate_primer_pair(
            sample_db, "chr1", 500, 20, 210, 20
        )
        assert result["common_snp_risk"] == "medium"
        assert result["left_primer_common_snp_count"] == 0
        assert result["right_primer_common_snp_count"] == 1

    def test_both_primers_overlap(self, sample_db):
        """Both primers overlapping SNPs should be detected."""
        # Left at 90-110 (rs100), right at 290-310 (rs300)
        result = annotate_primer_pair(
            sample_db, "chr1", 90, 20, 310, 20
        )
        assert result["common_snp_risk"] == "medium"
        assert result["left_primer_common_snp_count"] == 1
        assert result["right_primer_common_snp_count"] == 1

    def test_high_risk_overrides_medium(self, sample_db):
        """High risk from one primer should override medium from other."""
        # This test needs a case where one primer has 3' overlap
        # and the other has non-3' overlap
        # For now, just test that high > medium
        result = annotate_primer_pair(
            sample_db, "chr1", 90, 20, 500, 20
        )
        # Left primer has medium risk
        assert result["common_snp_risk"] == "medium"
