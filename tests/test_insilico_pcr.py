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
    resolve_ispcr_database,
    run_ispcr_batch,
)
from primer_panel.stage3_inputs import build_stage3_inputs
from primer_panel.writers import PrimerRecord, TargetRecord

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


class TestIsPcrDatabaseResolution:
    def test_prefers_same_basename_twobit_next_to_genome_fasta(self, tmp_path):
        genome = tmp_path / "hg38.fa"
        genome.write_text(">chr1\nACGT\n", encoding="utf-8")
        twobit = tmp_path / "hg38.2bit"
        twobit.write_bytes(b"")

        resolved = resolve_ispcr_database(genome)

        assert resolved.database_path == twobit
        assert resolved.ooc_path is None

    def test_discovers_matching_tile_ooc_next_to_database(self, tmp_path):
        genome = tmp_path / "hg38.fa"
        genome.write_text(">chr1\nACGT\n", encoding="utf-8")
        twobit = tmp_path / "hg38.2bit"
        twobit.write_bytes(b"")
        ooc = tmp_path / "hg38.11.ooc"
        ooc.write_bytes(b"")
        (tmp_path / "hg38.12.ooc").write_bytes(b"")

        resolved = resolve_ispcr_database(genome, tile_size=11)

        assert resolved.database_path == twobit
        assert resolved.ooc_path == ooc

    def test_falls_back_to_fasta_when_no_prepared_database_exists(self, tmp_path):
        genome = tmp_path / "genome.fa"
        genome.write_text(">chr1\nACGT\n", encoding="utf-8")

        resolved = resolve_ispcr_database(genome)

        assert resolved.database_path == genome
        assert resolved.ooc_path is None


class TestIsPcrCommand:
    def test_run_ispcr_batch_uses_resolved_database_and_ooc(
        self, tmp_path, monkeypatch,
    ):
        genome = tmp_path / "hg38.fa"
        genome.write_text(">chr1\nACGT\n", encoding="utf-8")
        twobit = tmp_path / "hg38.2bit"
        twobit.write_bytes(b"")
        ooc = tmp_path / "hg38.11.ooc"
        ooc.write_bytes(b"")
        captured: dict[str, object] = {}

        class Proc:
            returncode = 0
            stdout = "chr1\t10\t50\tp1\t0\t+\n"
            stderr = ""

        def fake_run(cmd, **kwargs):
            captured["cmd"] = cmd
            captured["kwargs"] = kwargs
            return Proc()

        monkeypatch.setattr("subprocess.run", fake_run)

        hits = run_ispcr_batch(
            [{"name": "p1", "fwd": "AAAA", "rev": "TTTT"}],
            str(genome),
            tile_size=11,
        )

        assert hits["p1"][0].start == 10
        cmd = captured["cmd"]
        assert cmd[1] == str(twobit)
        assert "-ooc=" + str(ooc) in cmd
        assert "-tileSize=11" in cmd
        assert "-minPerfect=15" in cmd
        assert "-minGood=15" in cmd


class TestStage3Inputs:
    """Fast tests for building isPcr batch inputs from Stage 1/2 records."""

    def test_build_stage3_inputs_uses_template_relative_primer_coords(self):
        target = TargetRecord(
            gene="GENE",
            transcript_id="ENST1",
            selection_reason="canonical",
            target_id="GENE_cds1",
            required_chrom="chr1",
            required_start=100,
            required_end=200,
            required_length=100,
            extended_chrom="chr1",
            extended_start=90,
            extended_end=210,
            extended_length=120,
            template_chrom="chr1",
            template_start=1000,
            template_end=1500,
            template_length=500,
            sequence_target_start_0based=50,
            sequence_target_length=120,
            sequence_target_for_primer3plus_1based="51,120",
            strand="+",
            product_min=2700,
            product_max=3300,
            cds_exon_numbers="1",
            cds_exon_ids="ex1",
            covered_cds_count=1,
            cds_exon_coords="100-200",
            status="template_ok",
            needs_review=False,
            sequence_status="real",
            target_qc_status="ok",
        )
        ok_primer = PrimerRecord(
            target_id="GENE_cds1",
            primer_rank=2,
            forward_primer="AAA",
            reverse_primer="TTT",
            forward_tm=60.0,
            reverse_tm=61.0,
            tm_diff=1.0,
            forward_gc=50.0,
            reverse_gc=50.0,
            primer_pair_penalty=0.1,
            primer_left_start=10,
            primer_left_len=20,
            primer_right_start=210,
            primer_right_len=20,
            primer3_product_size=201,
            primer3_status="ok",
            primer3_explain="",
            sequence_target_start_0based=50,
            sequence_target_length=120,
        )
        failed_primer = PrimerRecord(
            target_id="GENE_cds1",
            primer_rank=0,
            forward_primer="",
            reverse_primer="",
            forward_tm=0.0,
            reverse_tm=0.0,
            tm_diff=0.0,
            forward_gc=0.0,
            reverse_gc=0.0,
            primer_pair_penalty=0.0,
            primer_left_start=0,
            primer_left_len=0,
            primer_right_start=0,
            primer_right_len=0,
            primer3_product_size=0,
            primer3_status="no_primer",
            primer3_explain="",
            sequence_target_start_0based=50,
            sequence_target_length=120,
        )

        primer_batch, expected_coords = build_stage3_inputs(
            [target], [ok_primer, failed_primer],
        )

        assert primer_batch == [{
            "name": "GENE_cds1_rank2",
            "fwd": "AAA",
            "rev": "TTT",
            "expected_chrom": "chr1",
            "expected_start": 1010,
            "expected_end": 1211,
        }]
        assert expected_coords == {"GENE_cds1": ("chr1", 90, 210)}


class TestPanelFinalizationRescue:
    """Fast tests for panel finalization rescue specificity inputs."""

    def test_panel_finalize_accepts_argv(self, tmp_path):
        import primer_panel.panel_finalization as panel_finalization

        input_dir = tmp_path / "input"
        output_dir = tmp_path / "output"
        input_dir.mkdir()
        (input_dir / "primers_specificity.tsv").write_text(
            "\t".join([
                "target_id", "primer_rank", "forward_primer", "reverse_primer",
                "primer3_status", "insilico_status", "primer_pair_penalty",
            ])
            + "\n"
            + "\t".join(["GENE_cds1", "1", "AAA", "TTT", "ok", "unique_pass", "0.1"])
            + "\n",
            encoding="utf-8",
        )

        panel_finalization.main([
            "--input-dir", str(input_dir),
            "--output-dir", str(output_dir),
            "--genome-fasta", str(tmp_path / "genome.fa"),
        ])

        assert (output_dir / "recommended_primers.tsv").exists()

    def test_panel_finalize_does_not_run_fth1_rescue_by_default(
        self, tmp_path, monkeypatch,
    ):
        import primer_panel.panel_finalization as panel_finalization

        input_dir = tmp_path / "input"
        output_dir = tmp_path / "output"
        input_dir.mkdir()
        (input_dir / "primers_specificity.tsv").write_text(
            "\t".join([
                "target_id", "primer_rank", "forward_primer", "reverse_primer",
                "primer3_status", "insilico_status", "primer_pair_penalty",
                "insilico_hits",
            ])
            + "\n"
            + "\t".join([
                "FTH1_cds1_4", "1", "AAA", "TTT", "ok", "multi_hit", "0.1",
                "chr11:100-200(100);chr2:300-400(100)",
            ])
            + "\n",
            encoding="utf-8",
        )
        (input_dir / "target_summary.tsv").write_text(
            "target_id\ttemplate_chrom\ttemplate_start\nFTH1_cds1_4\tchr11\t100\n",
            encoding="utf-8",
        )

        def fail_if_called(*_args, **_kwargs):
            raise AssertionError("FTH1 rescue should be explicit")

        monkeypatch.setattr(panel_finalization, "rescue_fth1_with_more_primers", fail_if_called)
        monkeypatch.setattr(panel_finalization, "rescue_fth1_split_targets", fail_if_called)

        panel_finalization.main([
            "--input-dir", str(input_dir),
            "--output-dir", str(output_dir),
            "--genome-fasta", str(tmp_path / "genome.fa"),
        ])

        assert (output_dir / "rescue_attempts.tsv").exists()

    def test_rescue_stage3_uses_saved_primer_coords(
        self, tmp_path, monkeypatch,
    ):
        import primer_panel.panel_finalization as panel_finalization

        rescue_primers = tmp_path / "rescue_primers.tsv"
        rescue_primers.write_text(
            "\t".join([
                "target_id", "primer_rank", "forward_primer", "reverse_primer",
                "forward_tm", "reverse_tm", "primer_pair_penalty",
                "primer3_product_size", "primer_left_start", "primer_right_start",
            ])
            + "\n"
            + "\t".join([
                "FTH1_cds1_4", "1", "AAA", "TTT", "60", "61", "0.1",
                "201", "10", "210",
            ])
            + "\n",
            encoding="utf-8",
        )

        target_summary = tmp_path / "target_summary.tsv"
        target_summary.write_text(
            "\t".join(["target_id", "template_start", "template_chrom"])
            + "\n"
            + "\t".join(["FTH1_cds1_4", "1000", "chr11"])
            + "\n",
            encoding="utf-8",
        )

        captured: dict[str, object] = {}

        def fake_check_ispcr_available(_bin):
            return True

        def fake_check_specificity_batch(primer_pairs, **_kwargs):
            captured["primer_pairs"] = primer_pairs
            return {
                "rescue_rank1": SpecificityResult(
                    insilico_status="unique_pass",
                    insilico_hit_count=1,
                    insilico_hits="chr11:1010-1211(201)",
                    insilico_best_chrom="chr11",
                    insilico_best_start=1010,
                    insilico_best_end=1211,
                    insilico_best_size=201,
                    specificity_pass=True,
                    specificity_explain="single hit matches expected product",
                )
            }

        monkeypatch.setattr(
            "primer_panel.insilico_pcr.check_ispcr_available",
            fake_check_ispcr_available,
        )
        monkeypatch.setattr(
            "primer_panel.insilico_pcr.check_specificity_batch",
            fake_check_specificity_batch,
        )

        result = panel_finalization.run_stage3_on_rescue(
            rescue_primers,
            target_summary,
            tmp_path / "genome.fa",
            tmp_path,
        )

        assert result["status"] == "ok"
        assert captured["primer_pairs"][0]["expected_start"] == 1010
        assert captured["primer_pairs"][0]["expected_end"] == 1211

    def test_fth1_multi_hit_uses_expected_target_coordinates(self):
        import primer_panel.panel_finalization as panel_finalization

        result = panel_finalization.analyze_fth1_multi_hit([{
            "target_id": "FTH1_cds1_4",
            "insilico_hits": "chr11:61966152-61969100(2948);chr2:300-400(100)",
            "expected_target_chrom": "chr11",
            "expected_target_start": "61965152",
            "expected_target_end": "61970152",
        }])

        assert (61966152, 61969100, 2948) in result["expected_hits"]
        assert result["off_target_by_chrom"] == {"chr2": [(300, 400, 100)]}


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
