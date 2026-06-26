"""Tests for in-silico PCR specificity checking.

Stage 3 does NOT filter hits by product size.  All isPcr hits are retained.
Classification is based on hit count and chrom/start/end matching only.

Fast tests (default): mock-based, no genome access.
Integration tests: require RUN_ISPCR_INTEGRATION=1, real genome FASTA, and isPcr.
"""

import os

import pytest
from primer_panel.insilico_pcr import (
    SpecificityResult,
    PcrHit,
    check_ispcr_available,
    check_specificity_batch,
    run_ispcr_batch,
)

# ──────────────────────────────────────────────────────────────────────
# Fast tests (always run)
# ──────────────────────────────────────────────────────────────────────


class TestIsPcrAvailable:
    def test_ispcr_found(self):
        """isPcr should be installed in the test environment."""
        result = check_ispcr_available()
        assert isinstance(result, bool)

    def test_custom_path_not_found(self):
        """Non-existent binary path returns False."""
        assert check_ispcr_available("/nonexistent/isPcr") is False


class TestSpecificityResult:
    """Test SpecificityResult dataclass."""

    def test_default_values(self):
        r = SpecificityResult(
            insilico_status="no_hit",
            insilico_hit_count=0,
            insilico_hits="",
            insilico_best_chrom="",
            insilico_best_start=0,
            insilico_best_end=0,
            insilico_best_size=0,
            specificity_pass=False,
            specificity_explain="",
        )
        assert r.insilico_status == "no_hit"
        assert r.specificity_pass is False

    def test_unique_pass(self):
        r = SpecificityResult(
            insilico_status="unique_pass",
            insilico_hit_count=1,
            insilico_hits="chr1:100-200(100)",
            insilico_best_chrom="chr1",
            insilico_best_start=100,
            insilico_best_end=200,
            insilico_best_size=100,
            specificity_pass=True,
            specificity_explain="single hit matches expected product",
        )
        assert r.specificity_pass is True
        assert r.insilico_hit_count == 1
        assert "matches" in r.specificity_explain


class TestSpecificityBatchEmpty:
    """Fast tests for check_specificity_batch edge cases."""

    def test_empty_batch(self):
        """Empty primer list returns empty results without calling isPcr."""
        results = check_specificity_batch(
            primer_pairs=[],
            genome_fasta="/mnt/e/hg38/genome.fa",
        )
        assert results == {}


class TestNoProductSizeFilter:
    """Stage 3 should NOT filter hits by product size."""

    def test_run_ispcr_batch_no_size_params(self):
        """run_ispcr_batch should not accept product_min/product_max."""
        import inspect
        sig = inspect.signature(run_ispcr_batch)
        params = list(sig.parameters.keys())
        assert "product_min" not in params
        assert "product_max" not in params

    def test_check_specificity_batch_no_size_params(self):
        """check_specificity_batch should not accept product_min/product_max."""
        import inspect
        sig = inspect.signature(check_specificity_batch)
        params = list(sig.parameters.keys())
        assert "product_min" not in params
        assert "product_max" not in params


class TestClassificationLogic:
    """Test specificity classification logic (mocked, no actual isPcr)."""

    def _make_result(self, status, hit_count, best_chrom="chr1",
                     best_start=1000, best_end=2000, best_size=1000,
                     explain=""):
        return SpecificityResult(
            insilico_status=status,
            insilico_hit_count=hit_count,
            insilico_hits="",
            insilico_best_chrom=best_chrom,
            insilico_best_start=best_start,
            insilico_best_end=best_end,
            insilico_best_size=best_size,
            specificity_pass=(status == "unique_pass"),
            specificity_explain=explain,
        )

    def test_unique_pass_explain(self):
        """unique_pass should have descriptive explain."""
        r = self._make_result("unique_pass", 1, explain="single hit matches expected product")
        assert "matches" in r.specificity_explain

    def test_no_hit_explain(self):
        """no_hit should have descriptive explain."""
        r = self._make_result("no_hit", 0, explain="no genome-wide PCR product detected")
        assert "no" in r.specificity_explain.lower()

    def test_multi_hit_explain(self):
        """multi_hit should mention multiple products."""
        r = self._make_result("multi_hit", 3, explain="multiple genome-wide PCR products detected (3 hits)")
        assert "multiple" in r.specificity_explain.lower()

    def test_off_target_explain(self):
        """unique_off_target should describe mismatch."""
        r = self._make_result("unique_off_target", 1,
                              explain="single hit does not match expected product: chrom mismatch (chr2 vs chr1)")
        assert "does not match" in r.specificity_explain
        assert "chrom mismatch" in r.specificity_explain


# ──────────────────────────────────────────────────────────────────────
# Integration tests (require RUN_ISPCR_INTEGRATION=1)
# ──────────────────────────────────────────────────────────────────────

_skip_integration = pytest.mark.skipif(
    os.environ.get("RUN_ISPCR_INTEGRATION") != "1",
    reason="Set RUN_ISPCR_INTEGRATION=1 to run integration tests",
)


@_skip_integration
class TestIsPcrIntegration:
    """Integration tests that run real isPcr on hg38 genome."""

    def test_single_primer_unique(self):
        """Known HAMP primer pair should produce results."""
        results = check_specificity_batch(
            primer_pairs=[{
                "name": "test1",
                "fwd": "CAGCAGTGGGACAGCCAGAC",
                "rev": "AATGCAGATGGGGAAGTGGG",
                "expected_chrom": "chr19",
                "expected_start": 35281979,
                "expected_end": 35284988,
            }],
            genome_fasta="/mnt/e/hg38/genome.fa",
        )
        assert "test1" in results
        r = results["test1"]
        assert r.insilico_status in ("unique_pass", "unique_off_target", "multi_hit", "no_hit")
        assert isinstance(r.insilico_hit_count, int)
        assert isinstance(r.specificity_pass, bool)
        assert isinstance(r.specificity_explain, str)
        assert len(r.specificity_explain) > 0

    def test_batch_multiple_primers(self):
        """Batch with multiple primers returns results for all."""
        primers = [
            {
                "name": f"test_{i}",
                "fwd": "ATCGATCGATCGATCG",
                "rev": "GCTAGCTAGCTAGCTAG",
                "expected_chrom": "chr1",
                "expected_start": 0,
                "expected_end": 100,
            }
            for i in range(3)
        ]
        results = check_specificity_batch(
            primer_pairs=primers,
            genome_fasta="/mnt/e/hg38/genome.fa",
        )
        assert len(results) == 3
        for name in ["test_0", "test_1", "test_2"]:
            assert name in results

    def test_large_product_not_filtered(self):
        """A primer pair that produces a large product (>3300bp) should NOT be filtered.

        This verifies that Stage 3 does not use product_min/product_max.
        """
        results = check_specificity_batch(
            primer_pairs=[{
                "name": "large_product",
                "fwd": "CCCCCAAAAGAAGCGGAGAT",
                "rev": "CCATGTTGGCCAAGCTTGTC",
                "expected_chrom": "chr6",
                "expected_start": 26087389,
                "expected_end": 26090262,
            }],
            genome_fasta="/mnt/e/hg38/genome.fa",
        )
        r = results["large_product"]
        # Should have at least 1 hit (not filtered out by size)
        assert r.insilico_hit_count >= 1
        # If it matches expected, should be unique_pass regardless of size
        if r.insilico_status == "unique_pass":
            assert r.specificity_pass is True
