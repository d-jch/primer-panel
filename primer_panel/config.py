"""Pipeline configuration with sensible defaults."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class PipelineConfig:
    """All tuneable parameters for the primer panel pipeline.

    Conceptual layers
    -----------------
    required_region  : raw CDS exons — must be covered by the amplicon.
    extended_target  : required_region extended toward gene interior to product_min
                       (only when required_region is shorter than product_min).
    design_template  : extended_target + primer_flank on each side — given to Primer3.
    primer_product   : the actual PCR amplicon — must cover required_region,
                       length constrained by product_min .. product_max.
    """

    # --- product sizing (PCR product length constraints) ---
    product_min: int = 2700           # minimum PCR product length (bp)
    product_max: int = 3300           # maximum PCR product length (bp)

    # --- CDS buffer (deprecated compatibility field; no longer used) ---
    cds_buffer: int = 0

    # --- primer flank (design_template extension beyond extended_target) ---
    primer_flank: int = 300           # bp added to each side of extended_target for Primer3 search

    # --- tiling overlap for oversized required intervals ---
    tile_overlap: int = 200           # bp overlap between tiled targets for a single oversized interval

    # --- Ensembl API ---
    ensembl_base: str = "https://rest.ensembl.org"
    api_timeout: int = 60            # seconds per request

    # --- sequence source ---
    genome_fasta: Path | None = None # path to bgzipped+indexed hg38 FASTA (optional)

    # --- output ---
    output_dir: Path = field(default_factory=lambda: Path("outputs"))

    # --- Primer3 (Stage 2) ---
    design_primers: bool = False       # enable Primer3 design
    primer3_bin: str = "primer3_core"  # path to primer3_core binary

    # Primer3Plus settings file (None = use primer3plus-core bundled default)
    primer3plus_settings_file: Path | None = None

    # Save per-target Primer3 Boulder input to output_dir/primer3_inputs/
    write_primer3_inputs: bool = False

    # User-specified Primer3 overrides (None = use Primer3Plus defaults)
    # These are only set when user explicitly passes CLI flags.
    primer_num_return: int | None = None
    primer_opt_size: int | None = None
    primer_min_size: int | None = None
    primer_max_size: int | None = None
    primer_opt_tm: float | None = None
    primer_min_tm: float | None = None
    primer_max_tm: float | None = None
    primer_max_tm_diff: float | None = None
    primer_min_gc: float | None = None
    primer_max_gc: float | None = None

    # --- In-silico PCR (Stage 3) ---
    check_specificity: bool = False     # enable in-silico PCR specificity check
    is_pcr_bin: str = "isPcr"          # path to isPcr binary
    pcr_tolerance: int = 10             # bp tolerance for coordinate matching

    # isPcr database acceleration (all optional; default = auto-discover)
    ispcr_db: Path | None = None        # explicit .2bit/.nib database path
    ispcr_ooc: Path | None = None       # explicit overused-tile file path
    ispcr_tile_size: int = 11           # tileSize for isPcr (default 11)
    prepare_ispcr_db: bool = False      # create .2bit from FASTA (explicit only)
    make_ispcr_ooc: bool = False        # create .ooc file (explicit only)

    # Common dbSNP annotation (optional)
    common_dbsnp_bed: Path | None = None  # path to common dbSNP BED file

    # --- Local annotation (Stage 1) ---
    annotation_gtf: Path | None = None      # local Ensembl GTF for offline annotation
    annotation_source: str = "auto"         # "auto", "ensembl-api", or "gtf"

    # --- Rescue ---
    rescue_flank: int = 300                 # template extension for rescue (default: same as primer_flank)
    rescue_num_return: int = 20             # PRIMER_NUM_RETURN for rescue (default: more candidates)

    def build_primer3_overrides(self) -> dict[str, str]:
        """Build dict of Primer3 tag overrides from user-specified CLI args.

        Only includes tags where the user explicitly set a value (not None).
        Maps CLI arg names to correct Primer3 tag names.
        """
        overrides: dict[str, str] = {}
        mapping = {
            "primer_num_return": "PRIMER_NUM_RETURN",
            "primer_opt_size": "PRIMER_OPT_SIZE",
            "primer_min_size": "PRIMER_MIN_SIZE",
            "primer_max_size": "PRIMER_MAX_SIZE",
            "primer_opt_tm": "PRIMER_OPT_TM",
            "primer_min_tm": "PRIMER_MIN_TM",
            "primer_max_tm": "PRIMER_MAX_TM",
            "primer_max_tm_diff": "PRIMER_PAIR_MAX_DIFF_TM",
            "primer_min_gc": "PRIMER_MIN_GC",
            "primer_max_gc": "PRIMER_MAX_GC",
        }
        for attr, tag in mapping.items():
            val = getattr(self, attr, None)
            if val is not None:
                overrides[tag] = str(val)
        return overrides
