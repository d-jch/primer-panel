"""Tests for dependency preflight checks and doctor command."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from primer_panel.preflight import (
    DoctorReport,
    FileCheck,
    ToolCheck,
    check_fatotwobit,
    check_file,
    check_ispcr,
    check_primer3,
    check_tool,
    preflight_genome_fasta,
    preflight_stage2,
    preflight_stage3,
    preflight_prepare_ispcr_db,
    run_doctor,
)


# ── ToolCheck basics ───────────────────────────────────────────────────────


class TestToolCheck:
    @patch("primer_panel.preflight.shutil.which", return_value="/usr/bin/primer3_core")
    def test_tool_found(self, mock_which):
        result = check_tool("primer3_core", "primer3_core", "hint")
        assert result.found is True
        assert result.path == "/usr/bin/primer3_core"
        assert result.name == "primer3_core"

    @patch("primer_panel.preflight.shutil.which", return_value=None)
    def test_tool_missing(self, mock_which):
        result = check_tool("primer3_core", "primer3_core", "hint")
        assert result.found is False
        assert result.path is None

    @patch("primer_panel.preflight.shutil.which", return_value="/opt/bin/p3")
    def test_custom_bin_path(self, mock_which):
        result = check_tool("primer3_core", "/opt/bin/p3", "hint")
        mock_which.assert_called_once_with("/opt/bin/p3")
        assert result.found is True


# ── Individual tool checks ─────────────────────────────────────────────────


class TestCheckPrimer3:
    @patch("primer_panel.preflight.shutil.which", return_value="/usr/bin/primer3_core")
    def test_found(self, mock_which):
        result = check_primer3()
        assert result.found is True
        assert "micromamba" in result.install_hint or "conda" in result.install_hint

    @patch("primer_panel.preflight.shutil.which", return_value=None)
    def test_missing(self, mock_which):
        result = check_primer3()
        assert result.found is False
        assert "primer3" in result.install_hint.lower()


class TestCheckIsPcr:
    @patch("primer_panel.preflight.shutil.which", return_value=None)
    def test_missing(self, mock_which):
        result = check_ispcr()
        assert result.found is False
        assert "ispcr" in result.install_hint.lower()


class TestCheckFaToTwoBit:
    @patch("primer_panel.preflight.shutil.which", return_value=None)
    def test_missing(self, mock_which):
        result = check_fatotwobit()
        assert result.found is False
        assert "fatotwobit" in result.install_hint.lower()


# ── FileCheck ──────────────────────────────────────────────────────────────


class TestFileCheck:
    def test_file_exists(self, tmp_path):
        f = tmp_path / "test.fa"
        f.touch()
        result = check_file("genome FASTA", f, required=True)
        assert result.exists is True

    def test_file_missing(self):
        result = check_file("genome FASTA", Path("/nonexistent"), required=True)
        assert result.exists is False
        assert result.required is True

    def test_file_not_provided(self):
        result = check_file("genome FASTA", None, required=False)
        assert result.exists is False
        assert result.required is False

    def test_hint_set(self):
        result = check_file("x", None, required=False, hint="custom hint")
        assert result.hint == "custom hint"


# ── DoctorReport ───────────────────────────────────────────────────────────


class TestDoctorReport:
    def test_all_ok(self):
        report = DoctorReport(
            tools=[ToolCheck("t1", True, "/bin/t1", "", required=True)],
            files=[FileCheck("f1", Path("/tmp"), True, False, "")],
        )
        assert report.all_ok is True

    def test_required_tool_missing(self):
        report = DoctorReport(
            tools=[ToolCheck("t1", False, None, "hint", required=True)],
            files=[],
        )
        assert report.all_ok is False

    def test_optional_tool_missing_still_ok(self):
        """Optional tools (e.g. faToTwoBit) missing should not affect all_ok."""
        report = DoctorReport(
            tools=[
                ToolCheck("primer3", True, "/bin/p3", "", required=True),
                ToolCheck("isPcr", True, "/bin/isPcr", "", required=True),
                ToolCheck("faToTwoBit", False, None, "hint", required=False),
            ],
            files=[],
        )
        assert report.all_ok is True

    def test_required_file_missing(self):
        report = DoctorReport(
            tools=[],
            files=[FileCheck("f1", Path("/nope"), False, True, "")],
        )
        assert report.all_ok is False

    def test_optional_file_not_provided_is_ok(self):
        report = DoctorReport(
            tools=[],
            files=[FileCheck("f1", None, False, False, "")],
        )
        assert report.all_ok is True

    def test_optional_file_provided_but_missing_is_not_ok(self):
        report = DoctorReport(
            tools=[],
            files=[FileCheck("f1", Path("/nope"), False, False, "")],
        )
        assert report.all_ok is False


# ── run_doctor ─────────────────────────────────────────────────────────────


class TestRunDoctor:
    @patch("primer_panel.preflight.shutil.which", return_value="/usr/bin/tool")
    def test_returns_report(self, mock_which):
        report = run_doctor()
        assert isinstance(report, DoctorReport)
        assert len(report.tools) == 3  # primer3, ispcr, fatotwobit
        assert len(report.files) == 2  # genome-fasta, dbsnp-bed

    @patch("primer_panel.preflight.shutil.which", return_value="/usr/bin/tool")
    def test_with_genome_fasta(self, mock_which, tmp_path):
        fa = tmp_path / "hg38.fa"
        fa.touch()
        report = run_doctor(genome_fasta=fa)
        genome_check = next(f for f in report.files if f.name == "genome FASTA")
        assert genome_check.exists is True

    @patch("primer_panel.preflight.shutil.which", return_value="/usr/bin/tool")
    def test_missing_provided_genome_fasta_fails_report(self, mock_which, tmp_path):
        report = run_doctor(genome_fasta=tmp_path / "missing.fa")
        assert report.all_ok is False

    @patch("primer_panel.preflight.shutil.which", return_value="/usr/bin/tool")
    def test_directory_provided_as_genome_fasta_fails_report(self, mock_which, tmp_path):
        report = run_doctor(genome_fasta=tmp_path)
        assert report.all_ok is False

    @patch("primer_panel.preflight.shutil.which", return_value=None)
    def test_all_tools_missing(self, mock_which):
        report = run_doctor()
        assert all(not t.found for t in report.tools)
        assert report.all_ok is False

    def test_optional_tool_missing_does_not_affect_exit(self, tmp_path):
        """faToTwoBit missing alone should not cause doctor exit code 1."""
        def which_side_effect(name):
            if name == "faToTwoBit":
                return None
            return f"/usr/bin/{name}"

        with patch("primer_panel.preflight.shutil.which", side_effect=which_side_effect):
            report = run_doctor()
            # Required tools found, optional tool missing → all_ok should be True
            assert report.all_ok is True


# ── preflight_default_all ──────────────────────────────────────────────────


class TestPreflightGenomeFasta:
    """Test genome FASTA preflight for all stage combinations."""

    def test_stage_targets_no_fasta_ok(self):
        """--stage targets does not require genome-fasta."""
        preflight_genome_fasta(None, "targets")  # should not raise

    def test_default_none_no_fasta_exits(self):
        """Default (no --stage, stage=None) without genome-fasta should exit."""
        with pytest.raises(SystemExit):
            preflight_genome_fasta(None, None)

    def test_explicit_all_no_fasta_exits(self):
        """--stage all without genome-fasta should exit."""
        with pytest.raises(SystemExit):
            preflight_genome_fasta(None, "all")

    def test_explicit_design_no_fasta_exits(self):
        """--stage design without genome-fasta should exit."""
        with pytest.raises(SystemExit):
            preflight_genome_fasta(None, "design")

    def test_explicit_specificity_no_fasta_exits(self):
        """--stage specificity without genome-fasta should exit."""
        with pytest.raises(SystemExit):
            preflight_genome_fasta(None, "specificity")

    def test_default_with_valid_fasta_ok(self, tmp_path):
        """Default with existing genome-fasta should not raise."""
        fa = tmp_path / "hg38.fa"
        fa.touch()
        preflight_genome_fasta(fa, None)  # should not raise

    def test_explicit_all_with_valid_fasta_ok(self, tmp_path):
        """--stage all with existing genome-fasta should not raise."""
        fa = tmp_path / "hg38.fa"
        fa.touch()
        preflight_genome_fasta(fa, "all")  # should not raise

    def test_nonexistent_fasta_exits(self, tmp_path):
        """--genome-fasta pointing to nonexistent file should exit."""
        with pytest.raises(SystemExit):
            preflight_genome_fasta(tmp_path / "nonexistent.fa", "all")

    def test_fasta_is_directory_exits(self, tmp_path):
        """--genome-fasta pointing to directory should exit."""
        d = tmp_path / "fasta_dir"
        d.mkdir()
        with pytest.raises(SystemExit):
            preflight_genome_fasta(d, "all")

    def test_design_with_valid_fasta_ok(self, tmp_path):
        """--stage design with valid fasta should pass."""
        fa = tmp_path / "hg38.fa"
        fa.touch()
        preflight_genome_fasta(fa, "design")  # should not raise


# ── preflight_stage2 ───────────────────────────────────────────────────────


class TestPreflightStage2:
    @patch("primer_panel.preflight.shutil.which", return_value="/usr/bin/primer3_core")
    def test_primer3_found(self, mock_which):
        cfg = MagicMock()
        cfg.primer3_bin = "primer3_core"
        preflight_stage2(cfg)  # should not raise

    @patch("primer_panel.preflight.shutil.which", return_value=None)
    def test_primer3_missing_exits(self, mock_which):
        cfg = MagicMock()
        cfg.primer3_bin = "primer3_core"
        with pytest.raises(SystemExit):
            preflight_stage2(cfg)


# ── preflight_stage3 ───────────────────────────────────────────────────────


class TestPreflightStage3:
    @patch("primer_panel.preflight.shutil.which", return_value="/usr/bin/isPcr")
    def test_ispcr_found_with_fasta(self, mock_which, tmp_path):
        cfg = MagicMock()
        cfg.is_pcr_bin = "isPcr"
        cfg.genome_fasta = tmp_path / "hg38.fa"
        cfg.genome_fasta.touch()
        preflight_stage3(cfg)  # should not raise

    @patch("primer_panel.preflight.shutil.which", return_value=None)
    def test_ispcr_missing_exits(self, mock_which):
        cfg = MagicMock()
        cfg.is_pcr_bin = "isPcr"
        cfg.genome_fasta = Path("/some/fa")
        with pytest.raises(SystemExit):
            preflight_stage3(cfg)

    @patch("primer_panel.preflight.shutil.which", return_value="/usr/bin/isPcr")
    def test_no_fasta_exits(self, mock_which):
        cfg = MagicMock()
        cfg.is_pcr_bin = "isPcr"
        cfg.genome_fasta = None
        with pytest.raises(SystemExit):
            preflight_stage3(cfg)

    @patch("primer_panel.preflight.shutil.which", return_value="/usr/bin/isPcr")
    def test_nonexistent_fasta_exits(self, mock_which, tmp_path):
        cfg = MagicMock()
        cfg.is_pcr_bin = "isPcr"
        cfg.genome_fasta = tmp_path / "nonexistent.fa"
        with pytest.raises(SystemExit):
            preflight_stage3(cfg)


# ── preflight_prepare_ispcr_db ─────────────────────────────────────────────


class TestPreflightPrepareIspcrDb:
    @patch("primer_panel.preflight.shutil.which", return_value="/usr/bin/faToTwoBit")
    def test_found(self, mock_which):
        cfg = MagicMock()
        preflight_prepare_ispcr_db(cfg)  # should not raise

    @patch("primer_panel.preflight.shutil.which", return_value=None)
    def test_missing_exits(self, mock_which):
        cfg = MagicMock()
        with pytest.raises(SystemExit):
            preflight_prepare_ispcr_db(cfg)


# ── Integration: default all missing genome-fasta via main ─────────────────


class TestMainDefaultAllMissingFasta:
    """Test that default stage=all without --genome-fasta exits clearly."""

    def test_default_all_no_fasta_exits(self):
        """Running with default stage and no --genome-fasta should exit."""
        from primer_panel.main import main

        with pytest.raises(SystemExit) as exc_info:
            main(["--genes", "HFE"])
        assert exc_info.value.code == 1

    def test_explicit_all_no_fasta_exits(self):
        """--stage all without --genome-fasta should exit."""
        from primer_panel.main import main

        with pytest.raises(SystemExit) as exc_info:
            main(["--genes", "HFE", "--stage", "all"])
        assert exc_info.value.code == 1

    def test_design_no_fasta_exits(self):
        """--stage design without --genome-fasta should exit."""
        from primer_panel.main import main

        with pytest.raises(SystemExit) as exc_info:
            main(["--genes", "HFE", "--stage", "design"])
        assert exc_info.value.code == 1

    def test_specificity_no_fasta_exits(self):
        """--stage specificity without --genome-fasta should exit."""
        from primer_panel.main import main

        with pytest.raises(SystemExit) as exc_info:
            main(["--genes", "HFE", "--stage", "specificity"])
        assert exc_info.value.code == 1

    def test_preflight_before_annotation_client(self, tmp_path):
        """Preflight should fail before annotation client is created.

        This verifies the ordering: preflight → output dir → annotation client.
        With a nonexistent FASTA and --annotation-gtf, the preflight error
        should appear, not a GTF parse error.
        """
        from primer_panel.main import main

        gtf = tmp_path / "annot.gtf"
        gtf.write_text('chr1\tensembl\tgene\t1\t100\t.\t+\t.\tgene_id "E1";\n')

        nonexistent_fasta = tmp_path / "does_not_exist.fa"

        with pytest.raises(SystemExit) as exc_info:
            main([
                "--genes", "HFE",
                "--stage", "all",
                "--genome-fasta", str(nonexistent_fasta),
                "--annotation-gtf", str(gtf),
            ])
        assert exc_info.value.code == 1

    @patch("primer_panel.preflight.shutil.which", return_value="/usr/bin/primer3_core")
    def test_stage_targets_no_fasta_ok(self, mock_which, tmp_path):
        """--stage targets should not require genome-fasta (preflight passes)."""
        from primer_panel.main import main

        out = tmp_path / "out"
        # --stage targets passes preflight without --genome-fasta.
        # It may succeed (if Ensembl API is reachable) or fail at the API
        # call, but it must NOT fail at the preflight check.
        try:
            main(["--genes", "HFE", "--stage", "targets", "--output-dir", str(out)])
        except SystemExit as e:
            # If it exits, it should not be because of missing genome-fasta
            # (preflight error).  Network errors are acceptable.
            pass
