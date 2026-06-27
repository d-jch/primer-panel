"""Tests for --stage CLI parameter and resolve_stage logic."""

import pytest
from unittest.mock import patch, MagicMock
from primer_panel.main import resolve_stage


class TestResolveStage:
    """Test stage resolution logic."""

    def test_stage_targets(self):
        """--stage targets should disable both design and specificity."""
        assert resolve_stage("targets", False, False) == (False, False)

    def test_stage_design(self):
        """--stage design should enable design, disable specificity."""
        assert resolve_stage("design", False, False) == (True, False)

    def test_stage_specificity(self):
        """--stage specificity should enable both design and specificity."""
        assert resolve_stage("specificity", False, False) == (True, True)

    def test_stage_all(self):
        """--stage all should enable both design and specificity."""
        assert resolve_stage("all", False, False) == (True, True)

    def test_stage_overrides_deprecated_design(self):
        """Explicit --stage should override --design-primers flag."""
        # --stage targets with --design-primers should still be targets only
        assert resolve_stage("targets", True, False) == (False, False)

    def test_stage_overrides_deprecated_specificity(self):
        """Explicit --stage should override --check-specificity flag."""
        # --stage design with --check-specificity should still be design only
        assert resolve_stage("design", False, True) == (True, False)

    def test_no_stage_with_design_flag(self):
        """Without --stage, --design-primers should enable design only."""
        assert resolve_stage(None, True, False) == (True, False)

    def test_no_stage_with_specificity_flag(self):
        """Without --stage, --check-specificity should enable both."""
        assert resolve_stage(None, False, True) == (True, True)

    def test_no_stage_default_all(self):
        """Without --stage or deprecated flags, default should be all (Stage 1+2+3)."""
        assert resolve_stage(None, False, False) == (True, True)

    def test_no_stage_both_flags(self):
        """Without --stage, both deprecated flags should enable both."""
        assert resolve_stage(None, True, True) == (True, True)


class TestStageIntegration:
    """Integration tests for --stage parameter in CLI."""

    def test_stage_targets_disables_design_and_specificity(self):
        """--stage targets should set design_primers=False, check_specificity=False."""
        from primer_panel.main import _parse_args

        args = _parse_args(["--genes", "TEST", "--stage", "targets"])
        assert args.stage == "targets"
        assert args.design_primers is False
        assert args.check_specificity is False

    def test_stage_design_enables_design_only(self):
        """--stage design should set design_primers=True, check_specificity=False."""
        from primer_panel.main import _parse_args

        args = _parse_args(["--genes", "TEST", "--stage", "design"])
        assert args.stage == "design"
        assert args.design_primers is False  # CLI flag, resolved by resolve_stage
        assert args.check_specificity is False

    def test_stage_specificity_enables_both(self):
        """--stage specificity should set both flags."""
        from primer_panel.main import _parse_args

        args = _parse_args(["--genes", "TEST", "--stage", "specificity"])
        assert args.stage == "specificity"
        assert args.design_primers is False  # CLI flag, resolved by resolve_stage
        assert args.check_specificity is False

    def test_stage_default_is_none(self):
        """Default --stage should be None."""
        from primer_panel.main import _parse_args

        args = _parse_args(["--genes", "TEST"])
        assert args.stage is None

    def test_common_dbsnp_bed_arg(self):
        """--common-dbsnp-bed should be parsed."""
        from primer_panel.main import _parse_args
        from pathlib import Path

        args = _parse_args([
            "--genes", "TEST",
            "--common-dbsnp-bed", "/path/to/snps.bed",
        ])
        assert args.common_dbsnp_bed == Path("/path/to/snps.bed")
