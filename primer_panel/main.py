"""CLI entry point for the primer panel pipeline."""

from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path

from .config import PipelineConfig
from .ensembl_client import EnsemblClient
from .cds_handler import build_required_intervals
from .preflight import (
    preflight_genome_fasta,
    preflight_prepare_ispcr_db,
    preflight_stage2,
    preflight_stage3,
    print_doctor_report,
    run_doctor,
)
from .target_planner_adapter import plan_targets_with_external_planner
from .stage3_inputs import build_stage3_inputs
from .variant_annotation import load_dbsnp_db, annotate_primer_pair
from .writers import (
    FailedTarget,
    PrimerRecord,
    TargetRecord,
    build_records,
    build_primer_records,
    build_specificity_records,
    write_bed,
    write_failed_targets,
    write_fasta,
    write_primers_tsv,
    write_primers_xlsx,
    write_qc_summary,
    write_required_bed,
    write_specificity_summary,
    write_specificity_tsv,
    write_summary_tsv,
    write_summary_xlsx,
    write_unique_primers,
)

logger = logging.getLogger("primer_panel")


def _run_stage2(all_records: list[TargetRecord], cfg: PipelineConfig) -> list:
    """Run Stage 2: Primer3Plus-like primer design on all targets.

    Primer3 receives SEQUENCE_TEMPLATE (design_template) and SEQUENCE_TARGET
    (extended_target relative coords) from Stage 1.  Product size ranges
    use Primer3Plus defaults.  Coverage check is QC-only (warning, not filter).
    """
    from .primer3_runner import (
        check_target_coverage,
        check_primer3_available,
        is_all_n,
        has_n,
        parse_fasta_sequences,
        run_primer3_for_target,
    )

    # Check Primer3 availability (preflight)
    preflight_stage2(cfg)

    # Load FASTA sequences (design_template sequences)
    fa_path = cfg.output_dir / "targets.fa"
    if not fa_path.exists():
        logger.error("targets.fa not found — run Stage 1 first")
        sys.exit(1)

    sequences = parse_fasta_sequences(fa_path)
    logger.info("Loaded %d sequences from %s", len(sequences), fa_path)

    # Check for placeholder sequences
    placeholder_targets = [tid for tid, seq in sequences.items() if is_all_n(seq)]
    if placeholder_targets:
        logger.error(
            "All sequences are N placeholders (%d targets). "
            "Primer3 requires real sequences. "
            "Please provide --genome-fasta with a bgzipped+indexed hg38 FASTA.",
            len(placeholder_targets),
        )
        sys.exit(1)

    # Build sequence-warning map
    sequence_warnings: dict[str, str] = {}
    for tid, seq in sequences.items():
        if has_n(seq) and not is_all_n(seq):
            sequence_warnings[tid] = "contains_N"
        elif is_all_n(seq):
            sequence_warnings[tid] = "all_N_placeholder"

    # Run Primer3 for each target
    primer_results: dict[str, list] = {}
    success_count = 0
    fail_count = 0

    for rec in all_records:
        seq = sequences.get(rec.target_id)
        if seq is None:
            logger.warning("No sequence found for %s — skipping Primer3", rec.target_id)
            primer_results[rec.target_id] = []
            sequence_warnings.setdefault(rec.target_id, "missing_sequence")
            fail_count += 1
            continue

        try:
            results, boulder_input = run_primer3_for_target(
                target_id=rec.target_id,
                sequence=seq,
                cfg=cfg,
                seq_target_start=rec.sequence_target_start_0based,
                seq_target_len=rec.sequence_target_length,
                return_input=True,
            )

            # Save Primer3 input if requested
            if cfg.write_primer3_inputs and boulder_input:
                input_dir = cfg.output_dir / "primer3_inputs"
                input_dir.mkdir(parents=True, exist_ok=True)
                input_path = input_dir / f"{rec.target_id}.primer3.input"
                input_path.write_text(boulder_input, encoding="utf-8")

            # Post-Primer3 QC (warning only, not filter):
            # Check that primer product covers SEQUENCE_TARGET
            for r in results:
                if r.primer3_status != "ok":
                    continue
                covers, detail = check_target_coverage(
                    r.primer_left_start, r.primer_left_len,
                    r.primer_right_start, r.primer_right_len,
                    rec.sequence_target_start_0based, rec.sequence_target_length,
                    rec.template_start,
                )
                if not covers:
                    r.primer3_explain = f"QC_warning: {detail}"
                    logger.warning("%s rank %d coverage QC: %s", rec.target_id, r.primer_rank, detail)

            primer_results[rec.target_id] = results

            ok_count = sum(1 for r in results if r.primer3_status == "ok")
            if ok_count > 0:
                logger.info(
                    "%s: %d primer pairs (best penalty=%.4f)",
                    rec.target_id, ok_count, results[0].primer_pair_penalty,
                )
                success_count += 1
            else:
                logger.warning(
                    "%s: no valid primers — %s",
                    rec.target_id,
                    results[0].primer3_explain if results[0].primer3_explain else "unknown",
                )
                fail_count += 1

        except Exception as exc:
            logger.error("Primer3 failed for %s: %s", rec.target_id, exc)
            primer_results[rec.target_id] = []
            sequence_warnings.setdefault(rec.target_id, f"primer3_error: {exc}")
            fail_count += 1

    # Build primer records and write outputs
    primer_records = build_primer_records(all_records, primer_results, sequence_warnings)

    primers_tsv_path = cfg.output_dir / "primers.tsv"
    primers_xlsx_path = cfg.output_dir / "primers.xlsx"

    write_primers_tsv(primer_records, primers_tsv_path)
    logger.info("Primers TSV → %s (%d records)", primers_tsv_path, len(primer_records))

    if write_primers_xlsx(all_records, primer_records, primers_xlsx_path):
        logger.info("Primers XLSX → %s", primers_xlsx_path)
    else:
        logger.info("Install openpyxl for primers XLSX output: pip install openpyxl")

    logger.info(
        "Stage 2 done. %d targets with primers, %d failed → %s",
        success_count, fail_count, cfg.output_dir,
    )
    return primer_records


def _run_stage3(
    all_records: list[TargetRecord],
    primer_records: list,
    cfg: PipelineConfig,
) -> None:
    """Run Stage 3: In-silico PCR specificity check using UCSC isPcr.

    Checks genome-wide specificity for all ok primer pairs in one batch.
    Stage 3 is the ONLY stage that produces genomic product coordinates.
    """
    from .insilico_pcr import (
        check_ispcr_available,
        check_specificity_batch,
        make_ispcr_ooc,
        prepare_ispcr_twobit,
    )

    # Check isPcr and genome-fasta availability (preflight)
    preflight_stage3(cfg)

    # Filter to ok primers only
    ok_primers = [pr for pr in primer_records if pr.primer3_status == "ok"]
    logger.info("Checking specificity for %d ok primer pairs …", len(ok_primers))

    if not ok_primers:
        logger.info("No ok primers to check — skipping Stage 3")
        return

    primer_batch, expected_coords = build_stage3_inputs(all_records, primer_records)

    # Optionally prepare isPcr database files (explicit opt-in only)
    ispcr_db_path = str(cfg.ispcr_db) if cfg.ispcr_db else None
    ispcr_ooc_path = str(cfg.ispcr_ooc) if cfg.ispcr_ooc else None

    if cfg.prepare_ispcr_db and cfg.genome_fasta:
        preflight_prepare_ispcr_db(cfg)
        logger.info("Preparing .2bit database from %s …", cfg.genome_fasta)
        twobit = prepare_ispcr_twobit(cfg.genome_fasta)
        ispcr_db_path = str(twobit)
        logger.info(".2bit database → %s", twobit)

    if cfg.make_ispcr_ooc and ispcr_db_path:
        logger.info("Creating .ooc file (tileSize=%d) …", cfg.ispcr_tile_size)
        ooc = make_ispcr_ooc(ispcr_db_path, tile_size=cfg.ispcr_tile_size)
        ispcr_ooc_path = str(ooc)
        logger.info(".ooc file → %s", ooc)

    # Run batch isPcr (no product-size filtering)
    t0 = time.time()
    specificity_results = check_specificity_batch(
        primer_pairs=primer_batch,
        genome_fasta=str(cfg.genome_fasta),
        ispcr_bin=cfg.is_pcr_bin,
        tolerance=cfg.pcr_tolerance,
        ispcr_db=ispcr_db_path,
        ispcr_ooc=ispcr_ooc_path,
        tile_size=cfg.ispcr_tile_size,
    )
    elapsed = time.time() - t0

    pass_count = sum(1 for r in specificity_results.values() if r.specificity_pass)
    logger.info(
        "isPcr done in %.0fs: %d unique_pass / %d total ok",
        elapsed, pass_count, len(ok_primers),
    )

    # Load common dbSNP database if provided
    snp_db = None
    if cfg.common_dbsnp_bed:
        if not cfg.common_dbsnp_bed.exists():
            logger.error("--common-dbsnp-bed path not found: %s", cfg.common_dbsnp_bed)
            sys.exit(1)
        logger.info("Loading common dbSNP from %s …", cfg.common_dbsnp_bed)
        snp_db = load_dbsnp_db(cfg.common_dbsnp_bed)

    # Build SNP annotations for all primer records
    snp_annotations: dict[str, dict] = {}
    if snp_db:
        for pr in primer_records:
            if pr.primer3_status != "ok":
                continue
            # Find the target record to get template coordinates
            target_rec = next(
                (r for r in all_records if r.target_id == pr.target_id), None
            )
            if target_rec is None:
                continue

            # Calculate genomic coordinates for primers
            left_start = target_rec.template_start + pr.primer_left_start
            right_start = target_rec.template_start + pr.primer_right_start

            annotation = annotate_primer_pair(
                snp_db,
                target_rec.template_chrom,
                left_start,
                pr.primer_left_len,
                right_start,
                pr.primer_right_len,
            )
            primer_name = f"{pr.target_id}_rank{pr.primer_rank}"
            snp_annotations[primer_name] = annotation

        logger.info("Annotated %d primer pairs with SNP data", len(snp_annotations))

    # Build specificity records (includes all primer records, not just ok)
    spec_records = build_specificity_records(
        primer_records, specificity_results, expected_coords,
        snp_annotations=snp_annotations if snp_annotations else None,
    )

    # Write outputs
    spec_tsv_path = cfg.output_dir / "primers_specificity.tsv"
    unique_path = cfg.output_dir / "primers_unique.tsv"
    spec_summary_path = cfg.output_dir / "stage3_summary.txt"

    write_specificity_tsv(spec_records, spec_tsv_path)
    write_unique_primers(spec_records, unique_path)
    write_specificity_summary(spec_records, spec_summary_path)

    logger.info("Primers specificity TSV → %s", spec_tsv_path)
    logger.info("Unique primers TSV → %s", unique_path)
    logger.info("Stage 3 summary → %s", spec_summary_path)
    logger.info(
        "Stage 3 done. %d unique_pass / %d total ok → %s",
        pass_count, len(ok_primers), cfg.output_dir,
    )


def _parse_product_size(value: str) -> tuple[int, int]:
    """Parse --product-size MIN-MAX string."""
    parts = value.split("-")
    if len(parts) != 2:
        raise argparse.ArgumentTypeError(
            f"Invalid product-size format '{value}'; expected MIN-MAX (e.g. 2700-3300)"
        )
    try:
        lo, hi = int(parts[0]), int(parts[1])
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"Invalid product-size format '{value}'; both values must be integers"
        )
    if lo >= hi:
        raise argparse.ArgumentTypeError(
            f"product-size MIN ({lo}) must be < MAX ({hi})"
        )
    return lo, hi


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="primer_panel",
        description="Generate PCR primer panel targets covering CDS regions from human gene symbols (hg38).",
    )
    p.add_argument(
        "--genes", nargs="+", default=None,
        help="Gene symbols to process (e.g. HFE HJV TFR2).",
    )
    p.add_argument("--output-dir", type=Path, default=Path("outputs"),
                    help="Directory for output files (default: outputs).")
    p.add_argument("--target-size", "--product-size", type=_parse_product_size,
                    default=(2700, 3300), metavar="MIN-MAX", dest="target_size",
                    help="Target size range for Stage 1 CDS merging/extension (default: 2700-3300). "
                         "Controls required_region grouping and extended_target scaling. "
                         "Does NOT constrain Primer3 product size (Stage 2 uses Primer3Plus defaults).")
    p.add_argument("--cds-buffer", type=int, default=0,
                    help="[Deprecated, ignored] CDS buffer is no longer used. Raw CDS exon coordinates are used directly.")
    p.add_argument("--primer-flank", type=int, default=300,
                    help="Bp added to each side of extended_target for Primer3 search (default: 300).")
    p.add_argument("--genome-fasta", type=Path, default=None,
                    help="Path to bgzipped+indexed hg38 FASTA for real sequence extraction.")

    # --- Stage 2: Primer3 design ---
    p.add_argument("--design-primers", action="store_true",
                    help="Enable Primer3 primer design (Stage 2). Requires real sequences.")
    p.add_argument("--primer3plus-settings", type=Path, default=None, dest="primer3plus_settings",
                    help="Path to Primer3Plus settings JSON/txt for Stage 2. "
                         "Default: primer3plus-core bundled default_settings.json.")
    p.add_argument("--write-primer3-inputs", action="store_true", dest="write_primer3_inputs",
                    help="Save per-target Primer3 Boulder input to output_dir/primer3_inputs/.")
    p.add_argument("--primer-num-return", type=int, default=None,
                    help="Override PRIMER_NUM_RETURN (default: Primer3Plus setting).")
    p.add_argument("--primer-opt-size", type=int, default=None,
                    help="Override PRIMER_OPT_SIZE (default: Primer3Plus setting).")
    p.add_argument("--primer-min-size", type=int, default=None,
                    help="Override PRIMER_MIN_SIZE (default: Primer3Plus setting).")
    p.add_argument("--primer-max-size", type=int, default=None,
                    help="Override PRIMER_MAX_SIZE (default: Primer3Plus setting).")
    p.add_argument("--primer-opt-tm", type=float, default=None,
                    help="Override PRIMER_OPT_TM (default: Primer3Plus setting).")
    p.add_argument("--primer-min-tm", type=float, default=None,
                    help="Override PRIMER_MIN_TM (default: Primer3Plus setting).")
    p.add_argument("--primer-max-tm", type=float, default=None,
                    help="Override PRIMER_MAX_TM (default: Primer3Plus setting).")
    p.add_argument("--primer-max-tm-diff", type=float, default=None,
                    help="Override PRIMER_PAIR_MAX_DIFF_TM (default: Primer3Plus setting).")
    p.add_argument("--primer-min-gc", type=float, default=None,
                    help="Override PRIMER_MIN_GC (default: Primer3Plus setting).")
    p.add_argument("--primer-max-gc", type=float, default=None,
                    help="Override PRIMER_MAX_GC (default: Primer3Plus setting).")
    p.add_argument("--primer3-bin", type=str, default="primer3_core",
                    help="Path to primer3_core binary (default: primer3_core).")

    # --- Stage 3: In-silico PCR specificity ---
    p.add_argument("--check-specificity", action="store_true",
                    help="Enable in-silico PCR specificity check (Stage 3).")
    p.add_argument("--is-pcr-bin", type=str, default="isPcr",
                    help="Path to isPcr binary (default: isPcr).")
    p.add_argument("--pcr-tolerance", type=int, default=10,
                    help="Bp tolerance for coordinate matching in specificity check (default: 10).")
    p.add_argument("--ispcr-db", type=Path, default=None,
                    help="Explicit .2bit/.nib database for isPcr (default: auto-discover alongside genome-fasta).")
    p.add_argument("--ispcr-ooc", type=Path, default=None,
                    help="Explicit overused-tile (.ooc) file for isPcr (default: auto-discover).")
    p.add_argument("--ispcr-tile-size", type=int, default=11,
                    help="Tile size for isPcr (default: 11).")
    p.add_argument("--prepare-ispcr-db", action="store_true",
                    help="Create a .2bit database from the genome FASTA before running isPcr.")
    p.add_argument("--make-ispcr-ooc", action="store_true",
                    help="Create an overused-tile (.ooc) file before running isPcr.")

    # --- Stage control ---
    p.add_argument("--stage", choices=["targets", "design", "specificity", "all"],
                   default=None,
                   help="Pipeline stage to run (default: all). "
                        "targets=Stage1 only, design=Stage1+2, specificity/all=Stage1+2+3. "
                        "Default (no --stage) is 'all' which requires --genome-fasta.")

    # --- Common dbSNP annotation ---
    p.add_argument("--common-dbsnp-bed", type=Path, default=None,
                   help="Path to common dbSNP file (.bed or .bb) for primer risk annotation. "
                        ".bed from N: wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp151Common.txt.gz "
                        "&& zcat snp151Common.txt.gz | cut -f2,3,4,5. "
                        ".bb from UCSC gbdb: wget https://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp155Common.bb "
                        "(requires bigBedToBed: micromamba install -c bioconda ucsc-bigbedtobed). "
                        "See README for details.")

    # --- Local annotation (Stage 1) ---
    p.add_argument("--annotation-gtf", type=Path, default=None,
                   help="Path to local Ensembl GTF (plain or .gtf.gz) for offline annotation. "
                        "Replaces Ensembl REST API lookups for transcript/CDS/exon data.")
    p.add_argument("--annotation-source", choices=["auto", "ensembl-api", "gtf"],
                   default="auto",
                   help="Annotation source: 'auto' uses GTF if --annotation-gtf is provided, "
                        "otherwise Ensembl API.  'gtf' and 'ensembl-api' force that source.")

    # --- Diagnostics ---
    p.add_argument("--doctor", action="store_true",
                   help="Run dependency and environment checks, then exit.")

    p.add_argument("-v", "--verbose", action="store_true",
                    help="Enable debug logging.")
    return p.parse_args(argv)


def resolve_stage(
    stage: str | None,
    deprecated_design_primers: bool,
    deprecated_check_specificity: bool,
) -> tuple[bool, bool]:
    """Resolve which stages to run based on --stage and deprecated flags.

    Returns:
        (run_design, run_specificity)
    """
    if stage is not None:
        # Explicit --stage takes precedence
        if deprecated_design_primers or deprecated_check_specificity:
            logger.warning(
                "--design-primers/--check-specificity ignored when --stage is specified"
            )
        if stage == "targets":
            return False, False
        elif stage == "design":
            return True, False
        else:  # specificity or all
            return True, True

    # No explicit --stage: use deprecated flags for compatibility
    if deprecated_check_specificity:
        return True, True
    elif deprecated_design_primers:
        return True, False
    else:
        # Default: full pipeline (Stage 1+2+3)
        return True, True


def _recompute_sequence_qc(
    record,
    seq_status: str,
    fa_sequences: dict[str, str],
) -> str:
    """Recompute target_qc_status from actual FASTA sequence content.

    Strips any previous sequence-related flags (``placeholder_sequence``,
    ``contains_N``, ``all_N``, ``ok``) and re-adds the correct one based on
    the extracted FASTA content.
    """
    from .primer3_runner import has_n, is_all_n

    stale_flags = {"placeholder_sequence", "contains_N", "all_N", "ok"}
    qc_parts = [p for p in record.target_qc_status.split(";") if p and p not in stale_flags]

    if seq_status == "placeholder":
        qc_parts.append("placeholder_sequence")
    elif seq_status == "real":
        seq = fa_sequences.get(record.target_id, "")
        if is_all_n(seq):
            qc_parts.append("all_N")
        elif has_n(seq):
            qc_parts.append("contains_N")

    return ";".join(qc_parts) if qc_parts else "ok"


def _resolve_annotation_source(args) -> str:
    """Determine which annotation source to use.

    Returns ``"gtf"`` or ``"ensembl-api"``.
    """
    source = args.annotation_source
    if source == "auto":
        return "gtf" if args.annotation_gtf else "ensembl-api"
    if source == "gtf" and not args.annotation_gtf:
        logger.error("--annotation-source gtf requires --annotation-gtf")
        sys.exit(1)
    return source


def _create_annotation_client(args, cfg):
    """Create the appropriate annotation client (GTF or Ensembl API).

    Returns an object with ``select_transcript(symbol)`` and
    ``lookup_gene(symbol)`` methods.
    """
    source = _resolve_annotation_source(args)
    if source == "gtf":
        from .gtf_annotation import GtfAnnotationClient

        gtf_path = args.annotation_gtf
        if not gtf_path.exists():
            logger.error("--annotation-gtf path not found: %s", gtf_path)
            sys.exit(1)
        logger.info("Using local GTF annotation: %s", gtf_path)
        return GtfAnnotationClient(gtf_path)
    else:
        return EnsemblClient(cfg)


def main(argv: list[str] | None = None) -> None:
    args = _parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s %(name)s: %(message)s",
    )

    # --- doctor mode: check environment and exit ---
    if args.doctor:
        report = run_doctor(
            genome_fasta=args.genome_fasta,
            common_dbsnp_bed=args.common_dbsnp_bed,
            primer3_bin=args.primer3_bin,
            ispcr_bin=args.is_pcr_bin,
        )
        print_doctor_report(report)
        sys.exit(0 if report.all_ok else 1)

    # --- Validate --genes is required for pipeline runs ---
    if not args.genes:
        logger.error("--genes is required (e.g. --genes HFE HJV TFR2)")
        sys.exit(1)

    # Resolve stages
    run_design, run_specificity = resolve_stage(
        args.stage, args.design_primers, args.check_specificity
    )

    product_min, product_max = args.target_size

    cfg = PipelineConfig(
        product_min=product_min,
        product_max=product_max,
        cds_buffer=args.cds_buffer,
        primer_flank=args.primer_flank,
        genome_fasta=args.genome_fasta,
        output_dir=args.output_dir,
        # Primer3 config
        design_primers=run_design,
        primer3_bin=args.primer3_bin,
        primer3plus_settings_file=args.primer3plus_settings,
        write_primer3_inputs=args.write_primer3_inputs,
        # User-specified overrides (None = use Primer3Plus defaults)
        primer_num_return=args.primer_num_return,
        primer_opt_size=args.primer_opt_size,
        primer_min_size=args.primer_min_size,
        primer_max_size=args.primer_max_size,
        primer_opt_tm=args.primer_opt_tm,
        primer_min_tm=args.primer_min_tm,
        primer_max_tm=args.primer_max_tm,
        primer_max_tm_diff=args.primer_max_tm_diff,
        primer_min_gc=args.primer_min_gc,
        primer_max_gc=args.primer_max_gc,
        # Stage 3 config
        check_specificity=run_specificity,
        is_pcr_bin=args.is_pcr_bin,
        pcr_tolerance=args.pcr_tolerance,
        ispcr_db=args.ispcr_db,
        ispcr_ooc=args.ispcr_ooc,
        ispcr_tile_size=args.ispcr_tile_size,
        prepare_ispcr_db=args.prepare_ispcr_db,
        make_ispcr_ooc=args.make_ispcr_ooc,
        # Common dbSNP annotation
        common_dbsnp_bed=args.common_dbsnp_bed,
        # Local annotation
        annotation_gtf=args.annotation_gtf,
        annotation_source=args.annotation_source,
    )

    # --- Preflight: validate genome-fasta BEFORE creating output dir / clients ---
    preflight_genome_fasta(cfg.genome_fasta, args.stage)

    cfg.output_dir.mkdir(parents=True, exist_ok=True)

    # Create annotation client (GTF or Ensembl API)
    client = _create_annotation_client(args, cfg)
    all_records: list[TargetRecord] = []
    all_failed: list[FailedTarget] = []

    for gene in args.genes:
        logger.info("Processing %s …", gene)
        try:
            ti, reason = client.select_transcript(gene)
        except Exception as exc:
            logger.error("Failed to process %s: %s", gene, exc)
            all_failed.append(FailedTarget(
                gene=gene,
                reason="api_or_lookup_error",
                detail=str(exc),
            ))
            continue

        # Check for CDS
        if not ti.cds_exons:
            logger.error("%s: no CDS/Translation found — skipping", gene)
            all_failed.append(FailedTarget(
                gene=gene,
                reason="no_cds",
                detail=f"Transcript {ti.transcript_id} has no CDS/Translation data",
            ))
            continue

        # Build required intervals (merge adjacent CDS exons, no buffer)
        required_intervals = build_required_intervals(ti.cds_exons)
        logger.info(
            "%s: %d CDS exons → %d required intervals",
            gene, len(ti.cds_exons), len(required_intervals),
        )

        # Compute gene-level bounds
        gene_req_start = min(ri.start for ri in required_intervals)
        gene_req_end = max(ri.end for ri in required_intervals)

        # Generate targets using primer-target-planner
        gene_data = client.lookup_gene(gene)
        gene_bounds_start = gene_data.get("start", gene_req_start + 1) - 1  # 1-based → 0-based
        gene_bounds_end = gene_data.get("end", gene_req_end)
        targets = plan_targets_with_external_planner(
            required_intervals, cfg,
            gene_start=gene_bounds_start, gene_end=gene_bounds_end,
        )
        logger.info("%s: %d targets generated", gene, len(targets))

        # Build output records (computes design_template with dynamic extension)
        seq_status = "placeholder"  # updated by write_fasta
        all_records.extend(build_records(
            gene, ti, targets, cfg, seq_status,
            gene_required_start=gene_req_start,
            gene_required_end=gene_req_end,
        ))

    if not all_records and not all_failed:
        logger.error("No targets generated — check gene symbols and network connectivity.")
        sys.exit(1)

    # --- Write outputs ---
    bed_path = cfg.output_dir / "targets.bed"
    required_bed_path = cfg.output_dir / "required_regions.bed"
    tsv_path = cfg.output_dir / "target_summary.tsv"
    xlsx_path = cfg.output_dir / "target_summary.xlsx"
    fa_path = cfg.output_dir / "targets.fa"
    failed_path = cfg.output_dir / "failed_targets.tsv"

    write_bed(all_records, bed_path)
    write_required_bed(all_records, required_bed_path)
    write_summary_tsv(all_records, tsv_path)

    if write_summary_xlsx(all_records, xlsx_path):
        logger.info("XLSX written to %s", xlsx_path)
    else:
        logger.info("Install openpyxl for XLSX output: pip install openpyxl")

    seq_status = write_fasta(all_records, cfg, fa_path)
    # Update records with actual sequence status and recompute QC
    from .primer3_runner import parse_fasta_sequences
    fa_sequences: dict[str, str] = {}
    if fa_path.exists():
        fa_sequences = parse_fasta_sequences(fa_path)

    for r in all_records:
        r.sequence_status = seq_status
        r.target_qc_status = _recompute_sequence_qc(r, seq_status, fa_sequences)
    # Re-write TSV with updated status
    write_summary_tsv(all_records, tsv_path)
    write_summary_xlsx(all_records, xlsx_path)  # already logged above if successful

    # Write failed targets (or remove stale file)
    if all_failed:
        write_failed_targets(all_failed, failed_path)
        logger.info("%d failed targets → %s", len(all_failed), failed_path)
    elif failed_path.exists():
        failed_path.unlink()
        logger.debug("Removed stale failed_targets.tsv")

    total_genes = len(args.genes)
    success_genes = len(set(r.gene for r in all_records))
    logger.info(
        "Stage 1 done. %d targets for %d/%d genes → %s",
        len(all_records), success_genes, total_genes, cfg.output_dir,
    )

    # ================================================================
    # Stage 2: Primer3 design (optional)
    # ================================================================
    primer_records = None
    if cfg.design_primers:
        primer_records = _run_stage2(all_records, cfg)

    # ================================================================
    # Stage 3: In-silico PCR specificity (optional)
    # ================================================================
    if cfg.check_specificity:
        if primer_records is None:
            logger.error("--check-specificity requires --design-primers")
            sys.exit(1)
        _run_stage3(all_records, primer_records, cfg)

    # --- QC summary ---
    qc_path = cfg.output_dir / "run_summary.txt"
    write_qc_summary(all_records, all_failed, primer_records, cfg, qc_path)

    logger.info("QC summary → %s", qc_path)


if __name__ == "__main__":
    main()
