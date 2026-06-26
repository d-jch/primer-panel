"""Tests for design_template extension logic, SEQUENCE_TARGET, and Primer3 QC."""

import pytest
from pathlib import Path

from primer_panel.config import PipelineConfig
from primer_panel.writers import (
    Target,
    TargetRecord,
    _directional_extend,
    build_records,
)
from primer_panel.primer3_runner import (
    check_target_coverage,
    is_all_n,
    has_n,
    _build_boulder_input,
    _extract_primers,
    _parse_boulder_output,
)
from primer_panel.ensembl_client import TranscriptInfo, CdsExon


# ──────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────

def _make_target(chrom, start, end, strand=1, cds_numbers=None, cds_ids=None, cds_coords=None):
    """Create a minimal Target for testing."""
    return Target(
        chrom=chrom,
        start=start,
        end=end,
        strand=strand,
        cds_exon_numbers=cds_numbers or [],
        cds_exon_ids=cds_ids or [],
        cds_exon_coords=cds_coords or [],
    )


def _make_ti(transcript_id="ENST00000000001"):
    """Create a minimal TranscriptInfo for testing."""
    return TranscriptInfo(
        transcript_id=transcript_id,
        biotype="protein_coding",
        is_mane_select=False,
        is_mane_plus_clinical=False,
        is_canonical=True,
        selection_reason="canonical",
        exons=[],
        cds_exons=[],
    )


def _default_cfg(**kwargs):
    """Create a PipelineConfig with test defaults."""
    defaults = dict(
        product_min=2700,
        product_max=3300,
        primer_flank=500,
    )
    defaults.update(kwargs)
    return PipelineConfig(**defaults)


# ──────────────────────────────────────────────────────────────────────
# A. Short target at gene start edge → extend rightward
# ──────────────────────────────────────────────────────────────────────

class TestEdgeExtension:
    """Short targets at gene edges should extend toward gene interior."""

    def test_gene_start_edge_extends_right(self):
        """HFE_cds1-like: short target at gene start → extend rightward."""
        cfg = _default_cfg()
        # Gene span: 1000..50000
        # Target: 1000..1136 (136bp, at gene start edge)
        target = _make_target("chr1", 1000, 1136, cds_numbers=[1])
        ti = _make_ti()

        records = build_records(
            "TEST", ti, [target], cfg, "placeholder",
            gene_required_start=1000,
            gene_required_end=50000,
        )
        r = records[0]

        # required unchanged
        assert r.required_start == 1000
        assert r.required_end == 1136
        assert r.required_length == 136

        # extended: rightward to product_min (2700)
        assert r.extended_start == 1000  # no left extension
        assert r.extended_length == 2700
        assert r.extended_end == 1000 + 2700

        # template: extended + flank
        assert r.template_start == 1000 - 500  # but clamped to 0? no, 500 >= 0
        assert r.template_length == 2700 + 2 * 500

    def test_gene_end_edge_extends_left(self):
        """Short target at gene end → extend leftward."""
        cfg = _default_cfg()
        # Gene span: 1000..50000
        # Target: 49864..50000 (136bp, at gene end edge)
        target = _make_target("chr1", 49864, 50000, cds_numbers=[5])
        ti = _make_ti()

        records = build_records(
            "TEST", ti, [target], cfg, "placeholder",
            gene_required_start=1000,
            gene_required_end=50000,
        )
        r = records[0]

        assert r.required_length == 136
        # extended: leftward to product_min
        assert r.extended_end == 50000  # no right extension
        assert r.extended_length == 2700
        assert r.extended_start == 50000 - 2700

    def test_middle_target_extends_symmetric(self):
        """Middle target → symmetric extension."""
        cfg = _default_cfg()
        # Gene span: 1000..50000
        # Target: 20000..20136 (136bp, in the middle)
        target = _make_target("chr1", 20000, 20136, cds_numbers=[3])
        ti = _make_ti()

        records = build_records(
            "TEST", ti, [target], cfg, "placeholder",
            gene_required_start=1000,
            gene_required_end=50000,
        )
        r = records[0]

        assert r.required_length == 136
        assert r.extended_length == 2700
        # symmetric: left = extra//2, right = extra - left
        extra = 2700 - 136  # 2564
        left = extra // 2
        right = extra - left
        assert r.extended_start == 20000 - left
        assert r.extended_end == 20136 + right

    def test_gene_start_at_chr0_target_at_start(self):
        """Gene starts at chr0, target at gene start → extend rightward."""
        cfg = _default_cfg()
        # Gene span: 0..500
        # Target: 0..100 (at gene start edge)
        target = _make_target("chr1", 0, 100, cds_numbers=[1])
        ti = _make_ti()

        records = build_records(
            "TEST", ti, [target], cfg, "placeholder",
            gene_required_start=0,
            gene_required_end=500,
        )
        r = records[0]
        assert r.required_length == 100
        assert r.extended_length == 2700
        assert r.extended_start == 0
        assert r.extended_end == 2700

    def test_gene_start_at_chr0_target_at_end(self):
        """Gene starts at chr0, target at gene end → extend left, compensate right."""
        cfg = _default_cfg()
        # Gene span: 0..500
        # Target: 400..500 (at gene end edge)
        target = _make_target("chr1", 400, 500, cds_numbers=[2])
        ti = _make_ti()

        records = build_records(
            "TEST", ti, [target], cfg, "placeholder",
            gene_required_start=0,
            gene_required_end=500,
        )
        r = records[0]
        assert r.required_length == 100
        # Should extend to product_min despite chr0 clamp
        assert r.extended_length == 2700
        assert r.extended_start == 0
        assert r.extended_end == 2700

    def test_middle_near_chr0_left_clamp_compensates(self):
        """Middle target near chr0 with left clamp → compensate right."""
        cfg = _default_cfg()
        # Gene span: 0..10000
        # Target: 100..200 (100bp, in middle but close to 0)
        target = _make_target("chr1", 100, 200, cds_numbers=[1])
        ti = _make_ti()

        records = build_records(
            "TEST", ti, [target], cfg, "placeholder",
            gene_required_start=0,
            gene_required_end=10000,
        )
        r = records[0]
        assert r.required_length == 100
        assert r.extended_length == 2700
        assert r.extended_start >= 0
        # Full length must be achieved
        assert r.extended_end - r.extended_start == 2700

    def test_template_length_always_sufficient(self):
        """template_length >= extended_length + 2*primer_flank, even near chr0."""
        cfg = _default_cfg(product_min=2700, primer_flank=500)
        target = _make_target("chr1", 0, 100, cds_numbers=[1])
        ti = _make_ti()

        records = build_records(
            "TEST", ti, [target], cfg, "placeholder",
            gene_required_start=0,
            gene_required_end=500,
        )
        r = records[0]
        assert r.template_length >= r.extended_length + 2 * cfg.primer_flank


# ──────────────────────────────────────────────────────────────────────
# B. required_length not modified by extension
# ──────────────────────────────────────────────────────────────────────

class TestRequiredLengthPreserved:
    """required_length must always reflect original CDS+buffer."""

    def test_short_target_preserves_required_length(self):
        cfg = _default_cfg()
        target = _make_target("chr1", 1000, 1136, cds_numbers=[1])
        ti = _make_ti()

        records = build_records(
            "TEST", ti, [target], cfg, "placeholder",
            gene_required_start=1000,
            gene_required_end=50000,
        )
        r = records[0]
        assert r.required_length == 136  # never changed

    def test_long_target_preserves_required_length(self):
        cfg = _default_cfg()
        # required already > product_min
        target = _make_target("chr1", 1000, 4000, cds_numbers=[1, 2, 3])
        ti = _make_ti()

        records = build_records(
            "TEST", ti, [target], cfg, "placeholder",
            gene_required_start=1000,
            gene_required_end=50000,
        )
        r = records[0]
        assert r.required_length == 3000  # unchanged
        # extended should equal required (no extension needed)
        assert r.extended_length == 3000


# ──────────────────────────────────────────────────────────────────────
# C. template_length = product_min + 2*primer_flank for short targets
# ──────────────────────────────────────────────────────────────────────

class TestTemplateLength:
    """template_length should be extended_length + 2*primer_flank."""

    def test_short_target_template_length(self):
        cfg = _default_cfg(product_min=2700, primer_flank=500)
        target = _make_target("chr1", 1000, 1136, cds_numbers=[1])
        ti = _make_ti()

        records = build_records(
            "TEST", ti, [target], cfg, "placeholder",
            gene_required_start=1000,
            gene_required_end=50000,
        )
        r = records[0]
        assert r.template_length == 2700 + 2 * 500  # 3700

    def test_long_target_template_length(self):
        cfg = _default_cfg(product_min=2700, primer_flank=500)
        # required = 3000, no extension needed
        target = _make_target("chr1", 1000, 4000, cds_numbers=[1, 2])
        ti = _make_ti()

        records = build_records(
            "TEST", ti, [target], cfg, "placeholder",
            gene_required_start=1000,
            gene_required_end=50000,
        )
        r = records[0]
        # extended = required (3000), template = 3000 + 1000
        assert r.template_length == 3000 + 2 * 500

    def test_near_chrom_zero_clamp(self):
        """Target near chr0 should clamp start to 0 and compensate right."""
        cfg = _default_cfg(product_min=2700, primer_flank=500)
        target = _make_target("chr1", 0, 136, cds_numbers=[1])
        ti = _make_ti()

        records = build_records(
            "TEST", ti, [target], cfg, "placeholder",
            gene_required_start=0,
            gene_required_end=50000,
        )
        r = records[0]
        assert r.template_start >= 0
        # template_length should still be >= product_min + 2*flank
        assert r.template_length >= 2700 + 2 * 500


# ──────────────────────────────────────────────────────────────────────
# D. SEQUENCE_TARGET computation
# ──────────────────────────────────────────────────────────────────────

class TestSequenceTarget:
    """SEQUENCE_TARGET must correctly represent extended_target within design_template."""

    def test_hfe_cds1_sequence_target(self):
        """HFE_cds1: template_length=3700, extended_length=2700, flank=500.

        sequence_target_start_0based = extended_start - template_start
        If gene starts at 26087410, extended_start = 26087410 (gene start edge, rightward),
        template_start = 26087410 - 500 = 26086910.
        seq_target_start_0based = 26087410 - 26086910 = 500.
        seq_target_length = 2700.
        """
        cfg = _default_cfg(product_min=2700, primer_flank=500)
        target = _make_target("chr6", 26087410, 26087546, cds_numbers=[1])
        ti = _make_ti()

        records = build_records(
            "HFE", ti, [target], cfg, "placeholder",
            gene_required_start=26087410,
            gene_required_end=35284988,  # approximate gene end
        )
        r = records[0]

        assert r.template_length == 3700
        assert r.extended_length == 2700
        assert r.sequence_target_start_0based == 500
        assert r.sequence_target_length == 2700
        assert r.sequence_target_for_primer3plus_1based == "501,2700"

    def test_sequence_target_within_template(self):
        """SEQUENCE_TARGET must be fully contained within SEQUENCE_TEMPLATE."""
        cfg = _default_cfg()
        target = _make_target("chr1", 5000, 5136, cds_numbers=[1])
        ti = _make_ti()

        records = build_records(
            "TEST", ti, [target], cfg, "placeholder",
            gene_required_start=5000,
            gene_required_end=50000,
        )
        r = records[0]

        assert r.sequence_target_start_0based >= 0
        assert r.sequence_target_start_0based + r.sequence_target_length <= r.template_length

    def test_sequence_target_covers_required(self):
        """extended_target (SEQUENCE_TARGET region) must cover required_region."""
        cfg = _default_cfg()
        target = _make_target("chr1", 5000, 5136, cds_numbers=[1])
        ti = _make_ti()

        records = build_records(
            "TEST", ti, [target], cfg, "placeholder",
            gene_required_start=5000,
            gene_required_end=50000,
        )
        r = records[0]

        # extended_target = [extended_start, extended_end)
        # It must cover [required_start, required_end)
        assert r.extended_start <= r.required_start
        assert r.extended_end >= r.required_end

    def test_sequence_target_1based_format(self):
        """1-based format for Primer3Plus display."""
        cfg = _default_cfg()
        target = _make_target("chr1", 5000, 5136, cds_numbers=[1])
        ti = _make_ti()

        records = build_records(
            "TEST", ti, [target], cfg, "placeholder",
            gene_required_start=5000,
            gene_required_end=50000,
        )
        r = records[0]

        start_1based = r.sequence_target_start_0based + 1
        expected = f"{start_1based},{r.sequence_target_length}"
        assert r.sequence_target_for_primer3plus_1based == expected


# ──────────────────────────────────────────────────────────────────────
# E. Stage 1 QC
# ──────────────────────────────────────────────────────────────────────

class TestStage1QC:
    """Stage 1 target_qc_status checks."""

    def test_qc_ok_for_normal_target(self):
        """Normal target should have target_qc_status = 'ok'."""
        cfg = _default_cfg()
        target = _make_target("chr1", 5000, 5136, cds_numbers=[1])
        ti = _make_ti()

        records = build_records(
            "TEST", ti, [target], cfg, "real",
            gene_required_start=5000,
            gene_required_end=50000,
        )
        r = records[0]
        assert r.target_qc_status == "ok"

    def test_qc_placeholder_sequence(self):
        """Placeholder sequence should be flagged."""
        cfg = _default_cfg()
        target = _make_target("chr1", 5000, 5136, cds_numbers=[1])
        ti = _make_ti()

        records = build_records(
            "TEST", ti, [target], cfg, "placeholder",
            gene_required_start=5000,
            gene_required_end=50000,
        )
        r = records[0]
        assert "placeholder_sequence" in r.target_qc_status

    def test_qc_tiled_target(self):
        """Tiled target should be flagged."""
        cfg = _default_cfg(product_max=3300)
        # Create a target wider than product_max to force tiling
        target = _make_target("chr1", 5000, 9000, cds_numbers=[1, 2])
        target.tiled = True
        target.status = "tiled"
        ti = _make_ti()

        records = build_records(
            "TEST", ti, [target], cfg, "real",
            gene_required_start=5000,
            gene_required_end=50000,
        )
        r = records[0]
        assert "tiled" in r.target_qc_status


# ──────────────────────────────────────────────────────────────────────
# F. Primer3 coverage QC (check_target_coverage)
# ──────────────────────────────────────────────────────────────────────

class TestTargetCoverageQC:
    """check_target_coverage: QC check for SEQUENCE_TARGET coverage."""

    def test_valid_product_covers_target(self):
        """Product covers SEQUENCE_TARGET → pass."""
        # template-relative: product 0..3200, target 500..3200
        # product_end = right_start + 1 = 3199 + 1 = 3200
        # target_end = 500 + 2700 = 3200
        ok, detail = check_target_coverage(
            primer_left_start=0, primer_left_len=20,
            primer_right_start=3199, primer_right_len=20,
            seq_target_start_0based=500, seq_target_length=2700,
            template_start=0,
        )
        assert ok is True
        assert detail == ""

    def test_product_misses_target_start(self):
        """Product start > target start → QC warning."""
        ok, detail = check_target_coverage(
            primer_left_start=600, primer_left_len=20,  # product starts at 600
            primer_right_start=2999, primer_right_len=20,
            seq_target_start_0based=500, seq_target_length=2700,
            template_start=0,
        )
        assert ok is False
        assert "product_start" in detail

    def test_product_misses_target_end(self):
        """Product end < target end → QC warning."""
        # target_end = 500 + 2700 = 3200
        # product_end = 2000 + 1 = 2001 < 3200
        ok, detail = check_target_coverage(
            primer_left_start=0, primer_left_len=20,
            primer_right_start=2000, primer_right_len=20,
            seq_target_start_0based=500, seq_target_length=2700,
            template_start=0,
        )
        assert ok is False
        assert "product_end" in detail


# ──────────────────────────────────────────────────────────────────────
# G. Primer3 boulder input (Primer3Plus-like)
# ──────────────────────────────────────────────────────────────────────

class TestBoulderInput:
    """Primer3 boulder input should match Primer3Plus defaults."""

    def _get_boulder(self, cfg, seq_target_start=500, seq_target_len=2700):
        """Helper: build boulder input and return (boulder_str, first_base_index)."""
        return _build_boulder_input("TEST", "ATCG" * 1000, cfg, seq_target_start, seq_target_len)

    def test_sequence_target_is_1based(self):
        """SEQUENCE_TARGET must be 1-based (Primer3Plus default PRIMER_FIRST_BASE_INDEX=1)."""
        cfg = _default_cfg()
        boulder, fbi = self._get_boulder(cfg, 500, 2700)
        assert fbi == 1
        assert "SEQUENCE_TARGET=501,2700" in boulder
        assert "SEQUENCE_TARGET=500,2700" not in boulder

    def test_contains_first_base_index(self):
        """Boulder input must include PRIMER_FIRST_BASE_INDEX=1."""
        cfg = _default_cfg()
        boulder, _ = self._get_boulder(cfg)
        assert "PRIMER_FIRST_BASE_INDEX=1" in boulder

    def test_contains_primer_task_generic(self):
        """Boulder input must include PRIMER_TASK=generic (Primer3Plus default)."""
        cfg = _default_cfg()
        boulder, _ = self._get_boulder(cfg)
        assert "PRIMER_TASK=generic" in boulder

    def test_contains_outside_penalty_zero(self):
        """Boulder input must include PRIMER_OUTSIDE_PENALTY=0 (Primer3Plus default)."""
        cfg = _default_cfg()
        boulder, _ = self._get_boulder(cfg)
        assert "PRIMER_OUTSIDE_PENALTY=0" in boulder

    def test_contains_inside_penalty(self):
        """Boulder input must include PRIMER_INSIDE_PENALTY=-1.0 (Primer3Plus default)."""
        cfg = _default_cfg()
        boulder, _ = self._get_boulder(cfg)
        assert "PRIMER_INSIDE_PENALTY=-1.0" in boulder

    def test_contains_pick_internal_oligo_0(self):
        """Boulder input must include PRIMER_PICK_INTERNAL_OLIGO=0."""
        cfg = _default_cfg()
        boulder, _ = self._get_boulder(cfg)
        assert "PRIMER_PICK_INTERNAL_OLIGO=0" in boulder

    def test_contains_primer3plus_product_size_ranges(self):
        """Boulder input must use Primer3Plus default product size ranges."""
        cfg = _default_cfg()
        boulder, _ = self._get_boulder(cfg)
        assert "501-600" in boulder
        assert "1501-3000" in boulder
        assert "10001-20000" in boulder

    def test_contains_primer_params(self):
        """Boulder input must include primer design parameters from Primer3Plus defaults."""
        cfg = _default_cfg()
        boulder, _ = self._get_boulder(cfg)
        assert "PRIMER_NUM_RETURN=10" in boulder
        assert "PRIMER_OPT_SIZE=20" in boulder
        assert "PRIMER_MIN_SIZE=18" in boulder
        assert "PRIMER_MAX_SIZE=27" in boulder
        assert "PRIMER_OPT_TM=60.0" in boulder
        assert "PRIMER_MIN_TM=57.0" in boulder  # Primer3Plus default
        assert "PRIMER_MAX_TM=63.0" in boulder  # Primer3Plus default

    def test_no_p3p_tags(self):
        """Boulder input must not contain P3P_* tags."""
        cfg = _default_cfg()
        boulder, _ = self._get_boulder(cfg)
        for line in boulder.split("\n"):
            if "=" in line:
                key = line.split("=")[0]
                assert not key.startswith("P3P_"), f"P3P_ tag found: {key}"

    def test_no_thermo_path(self):
        """Boulder input must not contain PRIMER_THERMODYNAMIC_PARAMETERS_PATH."""
        cfg = _default_cfg()
        boulder, _ = self._get_boulder(cfg)
        assert "PRIMER_THERMODYNAMIC_PARAMETERS_PATH" not in boulder


class TestCoordinateConversion:
    """Primer3 output coordinate conversion (1-based → 0-based)."""

    def test_extract_primers_converts_to_0based(self):
        """Primer3 1-based coords must be converted to 0-based."""
        parsed = {
            "PRIMER_PAIR_NUM_RETURNED": "1",
            "PRIMER_LEFT_0": "481,20",
            "PRIMER_RIGHT_0": "3353,20",
            "PRIMER_LEFT_0_SEQUENCE": "A" * 20,
            "PRIMER_RIGHT_0_SEQUENCE": "T" * 20,
            "PRIMER_LEFT_0_TM": "60.0",
            "PRIMER_RIGHT_0_TM": "60.0",
            "PRIMER_LEFT_0_GC_PERCENT": "50.0",
            "PRIMER_RIGHT_0_GC_PERCENT": "50.0",
            "PRIMER_PAIR_0_PENALTY": "0.1",
            "PRIMER_PAIR_0_PRODUCT_SIZE": "2873",
            "PRIMER_PAIR_EXPLAIN": "",
        }
        results = _extract_primers(parsed, "TEST", 1, 500, 2700, first_base_index=1)
        r = results[0]
        assert r.primer_left_start == 480   # 481 - 1
        assert r.primer_right_start == 3352  # 3353 - 1
        assert r.primer3_product_size == 2873

    def test_extract_primers_first_base_0(self):
        """With first_base_index=0, coords should not be shifted."""
        parsed = {
            "PRIMER_PAIR_NUM_RETURNED": "1",
            "PRIMER_LEFT_0": "480,20",
            "PRIMER_RIGHT_0": "3352,20",
            "PRIMER_LEFT_0_SEQUENCE": "A" * 20,
            "PRIMER_RIGHT_0_SEQUENCE": "T" * 20,
            "PRIMER_LEFT_0_TM": "60.0",
            "PRIMER_RIGHT_0_TM": "60.0",
            "PRIMER_LEFT_0_GC_PERCENT": "50.0",
            "PRIMER_RIGHT_0_GC_PERCENT": "50.0",
            "PRIMER_PAIR_0_PENALTY": "0.1",
            "PRIMER_PAIR_0_PRODUCT_SIZE": "2873",
            "PRIMER_PAIR_EXPLAIN": "",
        }
        results = _extract_primers(parsed, "TEST", 1, 500, 2700, first_base_index=0)
        r = results[0]
        assert r.primer_left_start == 480
        assert r.primer_right_start == 3352

    def test_stage3_expected_coords_use_0based(self):
        """Stage 3 expected coords must use 0-based PrimerResult coords."""
        # template_start=26086910, left_start=480 (0-based), right_start=3352 (0-based)
        template_start = 26086910
        left_start = 480
        right_start = 3352
        expected_start = template_start + left_start
        expected_end = template_start + right_start + 1
        assert expected_start == 26087390
        assert expected_end == 26090263


class TestTxtSettingsParser:
    """/txt settings file parser must skip P3_FILE_* tags."""

    def test_txt_settings_no_p3_file_tags(self):
        """Loading .txt settings must not include P3_FILE_TYPE, P3_FILE_ID, etc."""
        from primer_panel.primer3plus_core_adapter import load_settings_dict
        from pathlib import Path
        settings_path = Path(
            "/home/djch/micromamba/envs/primer_panel/lib/python3.14/"
            "site-packages/primer3plus_core/settings_files/"
            "primer3plus_2_4_2_default_settings.txt"
        )
        if not settings_path.exists():
            pytest.skip("primer3plus_2_4_2_default_settings.txt not found")
        settings = load_settings_dict(settings_path)
        for key in settings:
            assert not key.startswith("P3_FILE_"), f"P3_FILE_ tag in settings: {key}"
            assert key != "P3_COMMENT", "P3_COMMENT in settings"

    def test_boulder_from_txt_settings_clean(self):
        """Boulder input from .txt settings must not contain P3_FILE_* or P3P_*."""
        from primer_panel.primer3plus_core_adapter import build_primer3plus_input
        from pathlib import Path
        settings_path = Path(
            "/home/djch/micromamba/envs/primer_panel/lib/python3.14/"
            "site-packages/primer3plus_core/settings_files/"
            "primer3plus_2_4_2_default_settings.txt"
        )
        if not settings_path.exists():
            pytest.skip("primer3plus_2_4_2_default_settings.txt not found")
        boulder, _ = build_primer3plus_input(
            sequence_id="TEST",
            sequence_template="ACGT" * 250,
            sequence_target="500,500",
            settings_file=settings_path,
        )
        for line in boulder.split("\n"):
            if "=" in line:
                key = line.split("=")[0]
                assert not key.startswith("P3_FILE_"), f"P3_FILE_ tag: {key}"
                assert not key.startswith("P3P_"), f"P3P_ tag: {key}"


# ──────────────────────────────────────────────────────────────────────
# H. Primer3 placeholder detection
# ──────────────────────────────────────────────────────────────────────

class TestPlaceholderDetection:
    """All-N sequences should be detected."""

    def test_all_n_detected(self):
        assert is_all_n("NNNNNNNN") is True
        assert is_all_n("NNNN NNNN") is True

    def test_real_seq_not_all_n(self):
        assert is_all_n("ATCGATCG") is False

    def test_mixed_has_n(self):
        assert has_n("ATCGNNNN") is True
        assert has_n("ATCGATCG") is False

    def test_all_n_has_n(self):
        assert has_n("NNNN") is True
