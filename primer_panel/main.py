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
    preflight_bigbed_dbsnp,
    preflight_genome_fasta,
    preflight_prepare_ispcr_db,
    preflight_stage2,
    preflight_stage3,
    print_doctor_report,
    run_doctor,
)
from .target_planner_adapter import plan_targets_with_external_planner
from .stage3_inputs import build_stage3_inputs
from .variant_annotation import load_dbsnp_db, annotate_primer_pair, annotate_primer_snps
from .writers import (
    FailedTarget,
    PrimerRecord,
    SpecificityRecord,
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
    write_specificity_clean_tsv,
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
) -> list:
    """Run Stage 3: In-silico PCR specificity check using UCSC isPcr.

    Checks genome-wide specificity for all ok primer pairs in one batch.
    Stage 3 is the ONLY stage that produces genomic product coordinates.

    Returns:
        list[SpecificityRecord] for potential rescue reuse.
    """
    from .insilico_pcr import (
        check_ispcr_available,
        check_specificity_batch,
        ensure_twobit,
        make_ispcr_ooc,
        prepare_ispcr_twobit,
    )

    # Check isPcr and genome-fasta availability (preflight)
    preflight_stage3(cfg)
    # Check bigBedToBed BEFORE Stage 3 isPcr run (avoid wasting computation)
    preflight_bigbed_dbsnp(cfg)

    # Filter to ok primers only
    ok_primers = [pr for pr in primer_records if pr.primer3_status == "ok"]
    logger.info("Checking specificity for %d ok primer pairs …", len(ok_primers))

    if not ok_primers:
        logger.info("No ok primers to check — skipping Stage 3")
        return []

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

    # Auto-discover or auto-create .2bit when no explicit db was specified
    if ispcr_db_path is None and cfg.genome_fasta is not None:
        auto_twobit = ensure_twobit(cfg.genome_fasta)
        if auto_twobit is not None:
            ispcr_db_path = auto_twobit

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

            if target_rec.strand == "-":
                # Minus-strand gene: forward primer (now PRIMER_RIGHT sequence)
                # is at right_start and extends LEFT on the + strand;
                # reverse primer (now PRIMER_LEFT sequence) is at left_start
                # and extends RIGHT.  Use per-primer annotation to get the
                # 3' end direction correct.
                fwd_risk, fwd_total, fwd_3p, fwd_hits = annotate_primer_snps(
                    snp_db, target_rec.template_chrom,
                    right_start, pr.primer_right_len, is_reverse=True,
                )
                rev_risk, rev_total, rev_3p, rev_hits = annotate_primer_snps(
                    snp_db, target_rec.template_chrom,
                    left_start, pr.primer_left_len, is_reverse=False,
                )
                risk_order = {"none": 0, "medium": 1, "high": 2}
                overall_risk = max(fwd_risk, rev_risk, key=lambda r: risk_order.get(r, 0))
                all_hits = []
                if fwd_hits:
                    all_hits.append(f"left:{fwd_hits}")
                if rev_hits:
                    all_hits.append(f"right:{rev_hits}")
                annotation = {
                    "common_snp_risk": overall_risk,
                    "left_primer_common_snp_count": fwd_total,
                    "right_primer_common_snp_count": rev_total,
                    "left_primer_3p_common_snp_count": fwd_3p,
                    "right_primer_3p_common_snp_count": rev_3p,
                    "common_snp_hits": "|".join(all_hits),
                }
            else:
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

    clean_path = cfg.output_dir / "primers_specificity_clean.tsv"
    write_specificity_clean_tsv(spec_records, clean_path)

    logger.info("Primers specificity TSV → %s", spec_tsv_path)
    logger.info("Unique primers TSV → %s", unique_path)
    logger.info("Clean specificity TSV → %s", clean_path)
    logger.info("Stage 3 summary → %s", spec_summary_path)
    logger.info(
        "Stage 3 done. %d unique_pass / %d total ok → %s",
        pass_count, len(ok_primers), cfg.output_dir,
    )
    return spec_records


def _write_rescue_summary(
    primer_records: list,
    spec_records: list,
    rescue_targets: list[str],
    path: Path,
    rescue_flank: int,
    rescue_num_return: int,
) -> None:
    """Write a rescue summary with per-target outcome."""
    from collections import Counter
    rescue_set = set(rescue_targets)
    lines = ["# Rescue Summary", ""]

    for tid in sorted(rescue_set):
        tid_pr = [pr for pr in primer_records
                  if pr.target_id == tid and pr.primer_rank >= 100]
        tid_sp = [sr for sr in spec_records
                  if sr.target_id == tid and sr.primer_rank >= 100]
        unique_pass = [sr for sr in tid_sp if sr.insilico_status == "unique_pass"]
        clean_pass = [sr for sr in unique_pass if sr.common_snp_risk == "none"]

        lines.append(f"## {tid}")
        lines.append(f"  rescue_primers: {len(tid_pr)}")
        lines.append(f"  unique_pass: {len(unique_pass)}")
        lines.append(f"  clean_pass (no SNP): {len(clean_pass)}")
        statuses = Counter(sr.insilico_status for sr in tid_sp)
        for s, c in sorted(statuses.items()):
            lines.append(f"  status_{s}: {c}")
        snp_risks = Counter(sr.common_snp_risk for sr in tid_sp)
        for r, c in sorted(snp_risks.items()):
            lines.append(f"  snp_risk_{r}: {c}")
        lines.append("")

    lines.append("# Config")
    lines.append(f"  flank: {rescue_flank}")
    lines.append(f"  num_return: {rescue_num_return}")
    lines.append("")

    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    logger.info("Rescue summary → %s", path)


def _run_rescue_standalone(
    args: argparse.Namespace,
    cfg: PipelineConfig,
) -> None:
    """Standalone rescue: load existing outputs, re-run Stage 2+3, merge.

    Called when --rescue-target is specified without --genes.
    """
    import csv

    out = cfg.output_dir

    # 1. Validate prerequisites (preflight checks first)
    preflight_stage3(cfg)
    preflight_bigbed_dbsnp(cfg)
    if cfg.genome_fasta is None:
        logger.error("Rescue requires --genome-fasta")
        sys.exit(1)
    if cfg.common_dbsnp_bed and not cfg.common_dbsnp_bed.exists():
        logger.error("--common-dbsnp-bed path not found: %s", cfg.common_dbsnp_bed)
        sys.exit(1)

    tsv_path = out / "target_summary.tsv"
    if not tsv_path.exists():
        logger.error("target_summary.tsv not found in %s — run Stage 1 first", out)
        sys.exit(1)

    all_records: list[TargetRecord] = []
    with open(tsv_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            all_records.append(TargetRecord(
                gene=row["gene"],
                transcript_id=row["transcript_id"],
                selection_reason=row.get("selection_reason", ""),
                target_id=row["target_id"],
                required_chrom=row["required_chrom"],
                required_start=int(row["required_start"]),
                required_end=int(row["required_end"]),
                required_length=int(row["required_length"]),
                extended_chrom=row["extended_chrom"],
                extended_start=int(row["extended_start"]),
                extended_end=int(row["extended_end"]),
                extended_length=int(row["extended_length"]),
                template_chrom=row["template_chrom"],
                template_start=int(row["template_start"]),
                template_end=int(row["template_end"]),
                template_length=int(row["template_length"]),
                sequence_target_start_0based=int(row["sequence_target_start_0based"]),
                sequence_target_length=int(row["sequence_target_length"]),
                sequence_target_for_primer3plus_1based=row["sequence_target_for_primer3plus_1based"],
                strand=row["strand"],
                product_min=int(row["product_min"]),
                product_max=int(row["product_max"]),
                cds_exon_numbers=row["cds_exon_numbers"],
                cds_exon_ids=row["cds_exon_ids"],
                covered_cds_count=int(row["covered_cds_count"]),
                cds_exon_coords=row["cds_exon_coords"],
                status=row["status"],
                needs_review=row.get("needs_review", "False") == "True",
                sequence_status=row.get("sequence_status", ""),
                target_qc_status=row.get("target_qc_status", ""),
            ))
    logger.info("Loaded %d targets from %s", len(all_records), tsv_path)

    # 2. Load primer_records from primers.tsv
    primers_path = out / "primers.tsv"
    primer_records: list = []
    if primers_path.exists():
        with open(primers_path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                primer_records.append(PrimerRecord(
                    target_id=row["target_id"],
                    primer_rank=int(row["primer_rank"]),
                    forward_primer=row["forward_primer"],
                    reverse_primer=row["reverse_primer"],
                    forward_tm=float(row["forward_tm"]),
                    reverse_tm=float(row["reverse_tm"]),
                    tm_diff=float(row["tm_diff"]),
                    forward_gc=float(row["forward_gc"]),
                    reverse_gc=float(row["reverse_gc"]),
                    primer_pair_penalty=float(row["primer_pair_penalty"]),
                    primer_left_start=int(row["primer_left_start"]),
                    primer_left_len=int(row["primer_left_len"]),
                    primer_right_start=int(row["primer_right_start"]),
                    primer_right_len=int(row["primer_right_len"]),
                    primer3_product_size=int(row["primer3_product_size"]),
                    primer3_status=row["primer3_status"],
                    primer3_explain=row.get("primer3_explain", ""),
                    sequence_target_start_0based=int(row["sequence_target_start_0based"]),
                    sequence_target_length=int(row["sequence_target_length"]),
                ))
    logger.info("Loaded %d primer records from %s", len(primer_records), primers_path)

    # 3. Load spec_records from primers_specificity.tsv (may not exist).
    #    Coordinate columns (primer_left_start etc.) are not in the
    #    specificity TSV; pull them from the primer_records lookup.
    primer_by_name: dict[tuple[str, int], PrimerRecord] = {
        (pr.target_id, pr.primer_rank): pr for pr in primer_records
    }
    spec_path = out / "primers_specificity.tsv"
    if not spec_path.exists():
        logger.error(
            "primers_specificity.tsv not found in %s — run the full pipeline "
            "(Stage 1-3) before rescue", out,
        )
        sys.exit(1)
    spec_records: list = []
    with open(spec_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            key = (row["target_id"], int(row["primer_rank"]))
            pr = primer_by_name.get(key)
            spec_records.append(SpecificityRecord(
                    target_id=row["target_id"],
                    primer_rank=int(row["primer_rank"]),
                    forward_primer=row["forward_primer"],
                    reverse_primer=row["reverse_primer"],
                    forward_tm=float(row["forward_tm"]),
                    reverse_tm=float(row["reverse_tm"]),
                    tm_diff=float(row["tm_diff"]),
                    forward_gc=float(row["forward_gc"]),
                    reverse_gc=float(row["reverse_gc"]),
                    primer_pair_penalty=float(row["primer_pair_penalty"]),
                    primer_left_start=pr.primer_left_start if pr else 0,
                    primer_left_len=pr.primer_left_len if pr else 0,
                    primer_right_start=pr.primer_right_start if pr else 0,
                    primer_right_len=pr.primer_right_len if pr else 0,
                    primer3_product_size=int(row["primer3_product_size"]),
                    primer3_status=row["primer3_status"],
                    primer3_explain=row.get("primer3_explain", ""),
                    sequence_target_start_0based=int(row["sequence_target_start_0based"]),
                    sequence_target_length=int(row["sequence_target_length"]),
                    insilico_status=row["insilico_status"],
                    insilico_hit_count=int(row["insilico_hit_count"]),
                    insilico_hits=row["insilico_hits"],
                    insilico_best_chrom=row["insilico_best_chrom"],
                    insilico_best_start=int(row["insilico_best_start"]) if row["insilico_best_start"] else 0,
                    insilico_best_end=int(row["insilico_best_end"]) if row["insilico_best_end"] else 0,
                    insilico_best_size=int(row["insilico_best_size"]) if row["insilico_best_size"] else 0,
                    specificity_pass=row.get("specificity_pass", "False") == "True",
                    expected_target_chrom=row["expected_target_chrom"],
                    expected_target_start=int(row["expected_target_start"]) if row["expected_target_start"] else 0,
                    expected_target_end=int(row["expected_target_end"]) if row["expected_target_end"] else 0,
                    specificity_explain=row.get("specificity_explain", ""),
                    common_snp_risk=row.get("common_snp_risk", "none"),
                    left_primer_common_snp_count=int(row.get("left_primer_common_snp_count", "0")),
                    right_primer_common_snp_count=int(row.get("right_primer_common_snp_count", "0")),
                    left_primer_3p_common_snp_count=int(row.get("left_primer_3p_common_snp_count", "0")),
                    right_primer_3p_common_snp_count=int(row.get("right_primer_3p_common_snp_count", "0")),
                    common_snp_hits=row.get("common_snp_hits", ""),
            ))
    logger.info("Loaded %d specificity records from %s", len(spec_records), spec_path)

    # 4. Run rescue
    rescue_targets = args.rescue_target
    primer_records, spec_records = _run_rescue(
        all_records, primer_records, spec_records, rescue_targets, cfg,
    )

    # 5. Re-write outputs with merged rescue results
    write_primers_tsv(primer_records, out / "primers.tsv")
    write_primers_xlsx(all_records, primer_records, out / "primers.xlsx")
    write_specificity_tsv(spec_records, out / "primers_specificity.tsv")
    write_unique_primers(spec_records, out / "primers_unique.tsv")
    write_specificity_clean_tsv(spec_records, out / "primers_specificity_clean.tsv")
    write_specificity_summary(spec_records, out / "stage3_summary.txt")
    _write_rescue_summary(
        primer_records, spec_records, rescue_targets,
        out / "rescue_summary.txt",
        cfg.rescue_flank, cfg.rescue_num_return,
    )

    logger.info("Rescue done. Outputs rewritten in %s", out)


def _run_rescue(
    all_records: list[TargetRecord],
    primer_records: list,
    spec_records: list,
    rescue_targets: list[str],
    cfg: PipelineConfig,
) -> tuple[list, list]:
    """Rescue specific targets with relaxed Primer3 parameters.

    1. Validate rescue_targets against all_records.
    2. Rebuild template with rescue_flank, re-extract FASTA.
    3. Run Primer3 with rescue_num_return, isPcr, dbSNP.
    4. Merge rescue primers into primer_records and spec_records.

    Returns:
        (updated_primer_records, updated_spec_records).
    """
    import dataclasses

    from .primer3_runner import run_primer3_for_target
    from .insilico_pcr import check_specificity_batch

    # 1. Validate rescue targets (deduplicate, preserve order)
    valid_tids = {r.target_id for r in all_records}
    seen: set[str] = set()
    deficient: list[str] = []
    for tid in rescue_targets:
        if tid not in seen and tid in valid_tids:
            seen.add(tid)
            deficient.append(tid)
    missing = sorted(set(tid for tid in rescue_targets if tid not in valid_tids))
    if missing:
        logger.warning("Rescue: target IDs not found in records — skipping: %s", ", ".join(missing))
    if not deficient:
        logger.warning("Rescue: none of the specified targets exist — nothing to do")
        return primer_records, spec_records

    logger.info("Rescue: re-running Stage 2+3 for %d targets: %s", len(deficient), ", ".join(deficient))

    if cfg.genome_fasta is None:
        logger.error("Rescue requires --genome-fasta")
        sys.exit(1)

    # 2. Map target_id → original TargetRecord
    original_records: dict[str, TargetRecord] = {}
    for r in all_records:
        original_records[r.target_id] = r

    # Build rescue config: apply rescue-specific Primer3 overrides
    rescue_overrides = {}
    if cfg.rescue_min_tm is not None:
        rescue_overrides["primer_min_tm"] = cfg.rescue_min_tm
    if cfg.rescue_max_tm is not None:
        rescue_overrides["primer_max_tm"] = cfg.rescue_max_tm
    if cfg.rescue_min_gc is not None:
        rescue_overrides["primer_min_gc"] = cfg.rescue_min_gc
    if cfg.rescue_max_gc is not None:
        rescue_overrides["primer_max_gc"] = cfg.rescue_max_gc
    if cfg.rescue_min_size is not None:
        rescue_overrides["primer_min_size"] = cfg.rescue_min_size
    if cfg.rescue_opt_size is not None:
        rescue_overrides["primer_opt_size"] = cfg.rescue_opt_size
    if cfg.rescue_max_size is not None:
        rescue_overrides["primer_max_size"] = cfg.rescue_max_size
    rescue_cfg = dataclasses.replace(
        cfg,
        primer_flank=cfg.rescue_flank,
        primer_num_return=cfg.rescue_num_return,
        **rescue_overrides,
    )

    # Recompute template boundaries with rescue_flank
    rescue_target_records: list[TargetRecord] = []
    for tid in deficient:
        orig = original_records[tid]
        new_flank = cfg.rescue_flank
        ext_start, ext_end = orig.extended_start, orig.extended_end
        desired = (ext_end - ext_start) + 2 * new_flank
        new_template_start = max(0, ext_start - new_flank)
        new_template_end = ext_end + new_flank
        if (new_template_end - new_template_start) < desired and new_template_start == 0:
            new_template_end += desired - (new_template_end - new_template_start)
        new_template_len = new_template_end - new_template_start
        new_seq_target_start = ext_start - new_template_start
        new_seq_target_len = ext_end - ext_start
        new_seq_target_1based = f"{new_seq_target_start + 1},{new_seq_target_len}"

        rescue_target_records.append(dataclasses.replace(
            orig,
            template_start=new_template_start,
            template_end=new_template_end,
            template_length=new_template_len,
            sequence_target_start_0based=new_seq_target_start,
            sequence_target_length=new_seq_target_len,
            sequence_target_for_primer3plus_1based=new_seq_target_1based,
            product_min=cfg.product_min,
            product_max=cfg.product_max,
        ))

    # 3. Extract real FASTA
    try:
        from pyfaidx import Fasta
    except ImportError:
        logger.error("Rescue requires pyfaidx for FASTA extraction")
        return primer_records, spec_records

    genome = Fasta(str(cfg.genome_fasta))
    rescue_sequences: dict[str, str] = {}
    for rec in rescue_target_records:
        try:
            rescue_sequences[rec.target_id] = str(
                genome[rec.template_chrom][rec.template_start:rec.template_end]
            )
        except Exception:
            logger.warning("Rescue: cannot extract sequence for %s — skipping", rec.target_id)

    if not rescue_sequences:
        logger.warning("Rescue: no sequences extracted — aborting")
        return primer_records, spec_records

    logger.info("Rescue: extracted %d FASTA sequences with flank=%d",
                 len(rescue_sequences), cfg.rescue_flank)

    # 4. Run Primer3 on rescue targets
    rescue_primer_results: dict[str, list] = {}
    for rec in rescue_target_records:
        seq = rescue_sequences.get(rec.target_id)
        if seq is None:
            continue
        try:
            results = run_primer3_for_target(
                target_id=rec.target_id,
                sequence=seq,
                cfg=rescue_cfg,
                seq_target_start=rec.sequence_target_start_0based,
                seq_target_len=rec.sequence_target_length,
            )
            rescue_primer_results[rec.target_id] = results
            ok_count = sum(1 for r in results if r.primer3_status == "ok")
            logger.info("Rescue %s: %d primer pairs", rec.target_id, ok_count)
        except Exception as exc:
            logger.warning("Rescue Primer3 failed for %s: %s", rec.target_id, exc)

    if not rescue_primer_results:
        logger.warning("Rescue: no Primer3 results — aborting")
        return primer_records, spec_records

    # Remove any previous rescue primers for targets that were successfully rescued
    rescued_tid_set = set(rescue_primer_results.keys())
    primer_records = [pr for pr in primer_records
                      if not (pr.target_id in rescued_tid_set and pr.primer_rank >= 100)]
    spec_records = [sr for sr in spec_records
                    if not (sr.target_id in rescued_tid_set and sr.primer_rank >= 100)]

    # Determine rank offset per target (safe against repeated rescues)
    rank_offset: dict[str, int] = {}
    for tid in deficient:
        max_rank = max([0] + [pr.primer_rank for pr in primer_records if pr.target_id == tid])
        rank_offset[tid] = max(100, max_rank + 1)

    # Build rescue PrimerRecords with strand-aware swap (rescue-template-relative coords).
    # Coordinates are translated to original-template-relative AFTER isPcr/SNP.
    rescue_primer_records: list[PrimerRecord] = []
    for tr in rescue_target_records:
        results = rescue_primer_results.get(tr.target_id, [])
        offset = rank_offset[tr.target_id]
        for pr in results:
            if pr.primer3_status != "ok":
                continue

            if tr.strand == "-":
                # Minus-strand: swap so forward follows transcription direction
                rec = PrimerRecord(
                    target_id=tr.target_id,
                    primer_rank=pr.primer_rank + offset,
                    forward_primer=pr.reverse_primer,
                    reverse_primer=pr.forward_primer,
                    forward_tm=pr.reverse_tm,
                    reverse_tm=pr.forward_tm,
                    tm_diff=pr.tm_diff,
                    forward_gc=pr.reverse_gc,
                    reverse_gc=pr.forward_gc,
                    primer_pair_penalty=pr.primer_pair_penalty,
                    primer_left_start=pr.primer_left_start,
                    primer_left_len=pr.primer_left_len,
                    primer_right_start=pr.primer_right_start,
                    primer_right_len=pr.primer_right_len,
                    primer3_product_size=pr.primer3_product_size,
                    primer3_status=pr.primer3_status,
                    primer3_explain=pr.primer3_explain,
                    sequence_target_start_0based=pr.sequence_target_start_0based,
                    sequence_target_length=pr.sequence_target_length,
                )
            else:
                rec = PrimerRecord(
                    target_id=tr.target_id,
                    primer_rank=pr.primer_rank + offset,
                    forward_primer=pr.forward_primer,
                    reverse_primer=pr.reverse_primer,
                    forward_tm=pr.forward_tm,
                    reverse_tm=pr.reverse_tm,
                    tm_diff=pr.tm_diff,
                    forward_gc=pr.forward_gc,
                    reverse_gc=pr.reverse_gc,
                    primer_pair_penalty=pr.primer_pair_penalty,
                    primer_left_start=pr.primer_left_start,
                    primer_left_len=pr.primer_left_len,
                    primer_right_start=pr.primer_right_start,
                    primer_right_len=pr.primer_right_len,
                    primer3_product_size=pr.primer3_product_size,
                    primer3_status=pr.primer3_status,
                    primer3_explain=pr.primer3_explain,
                    sequence_target_start_0based=pr.sequence_target_start_0based,
                    sequence_target_length=pr.sequence_target_length,
                )
            rescue_primer_records.append(rec)

    logger.info("Rescue: %d rescue primer pairs built", len(rescue_primer_records))

    # Skip expensive isPcr if no rescue primers were produced
    if not rescue_primer_records:
        logger.warning("Rescue: no ok primer pairs produced — aborting")
        return primer_records, spec_records

    # 5. Prepare isPcr database (same logic as _run_stage3)
    from .insilico_pcr import ensure_twobit, make_ispcr_ooc, prepare_ispcr_twobit
    ispcr_db_path = str(cfg.ispcr_db) if cfg.ispcr_db else None
    ispcr_ooc_path = str(cfg.ispcr_ooc) if cfg.ispcr_ooc else None

    if cfg.prepare_ispcr_db and cfg.genome_fasta:
        preflight_prepare_ispcr_db(cfg)
        logger.info("Preparing .2bit database from %s …", cfg.genome_fasta)
        twobit = prepare_ispcr_twobit(cfg.genome_fasta)
        ispcr_db_path = str(twobit)
        logger.info(".2bit database → %s", twobit)

    if ispcr_db_path is None and cfg.genome_fasta is not None:
        auto_twobit = ensure_twobit(cfg.genome_fasta)
        if auto_twobit is not None:
            ispcr_db_path = auto_twobit

    if cfg.make_ispcr_ooc and ispcr_db_path:
        logger.info("Creating .ooc file (tileSize=%d) …", cfg.ispcr_tile_size)
        ooc = make_ispcr_ooc(ispcr_db_path, tile_size=cfg.ispcr_tile_size)
        ispcr_ooc_path = str(ooc)
        logger.info(".ooc file → %s", ooc)

    # Run isPcr with rescue-template-relative coords
    rescue_batch, rescue_expected_coords = build_stage3_inputs(
        rescue_target_records, rescue_primer_records,
    )

    rescue_specificity = check_specificity_batch(
        primer_pairs=rescue_batch,
        genome_fasta=str(cfg.genome_fasta),
        ispcr_bin=cfg.is_pcr_bin,
        tolerance=cfg.pcr_tolerance,
        ispcr_db=ispcr_db_path,
        ispcr_ooc=ispcr_ooc_path,
        tile_size=cfg.ispcr_tile_size,
    )

    pass_count = sum(1 for r in rescue_specificity.values() if r.specificity_pass)
    logger.info("Rescue isPcr: %d unique_pass / %d rescue primers", pass_count, len(rescue_primer_records))

    # 6. dbSNP annotation — uses rescue-template-relative coords + rescue TargetRecords
    rescue_by_tid_for_snp = {r.target_id: r for r in rescue_target_records}
    snp_annotations: dict[str, dict] = {}
    if cfg.common_dbsnp_bed and cfg.common_dbsnp_bed.exists():
        snp_db = load_dbsnp_db(cfg.common_dbsnp_bed)
        for pr in rescue_primer_records:
            if pr.primer3_status != "ok":
                continue
            target_rec = rescue_by_tid_for_snp.get(pr.target_id)
            if target_rec is None:
                continue

            left_start = target_rec.template_start + pr.primer_left_start
            right_start = target_rec.template_start + pr.primer_right_start

            if target_rec.strand == "-":
                fwd_risk, fwd_total, fwd_3p, fwd_hits = annotate_primer_snps(
                    snp_db, target_rec.template_chrom,
                    right_start, pr.primer_right_len, is_reverse=True,
                )
                rev_risk, rev_total, rev_3p, rev_hits = annotate_primer_snps(
                    snp_db, target_rec.template_chrom,
                    left_start, pr.primer_left_len, is_reverse=False,
                )
                risk_order = {"none": 0, "medium": 1, "high": 2}
                overall_risk = max(fwd_risk, rev_risk, key=lambda r: risk_order.get(r, 0))
                all_hits = []
                if fwd_hits:
                    all_hits.append(f"left:{fwd_hits}")
                if rev_hits:
                    all_hits.append(f"right:{rev_hits}")
                annotation = {
                    "common_snp_risk": overall_risk,
                    "left_primer_common_snp_count": fwd_total,
                    "right_primer_common_snp_count": rev_total,
                    "left_primer_3p_common_snp_count": fwd_3p,
                    "right_primer_3p_common_snp_count": rev_3p,
                    "common_snp_hits": "|".join(all_hits),
                }
            else:
                annotation = annotate_primer_pair(
                    snp_db,
                    target_rec.template_chrom,
                    left_start, pr.primer_left_len,
                    right_start, pr.primer_right_len,
                )
            primer_name = f"{pr.target_id}_rank{pr.primer_rank}"
            snp_annotations[primer_name] = annotation
        logger.info("Rescue: annotated %d primer pairs with SNP data", len(snp_annotations))

    # 7. Translate rescue primer coordinates to original-template-relative
    #    so they are consistent with the original target_summary.tsv template.
    #    genomic = rescue_start + rescue_offset = orig_start + orig_offset
    #    → orig_offset = resc_offset + (rescue_start - orig_start)
    rescue_by_tid = {r.target_id: r for r in rescue_target_records}
    for pr in rescue_primer_records:
        orig = original_records[pr.target_id]
        tr = rescue_by_tid[pr.target_id]
        coord_delta = tr.template_start - orig.template_start
        pr.primer_left_start += coord_delta
        pr.primer_right_start += coord_delta
        pr.sequence_target_start_0based += coord_delta

    # 7. Build rescue specificity records
    rescue_spec_records = build_specificity_records(
        rescue_primer_records, rescue_specificity,
        expected_coords=rescue_expected_coords,
        snp_annotations=snp_annotations if snp_annotations else None,
    )

    rescue_clean = sum(
        1 for r in rescue_spec_records
        if r.insilico_status == "unique_pass" and r.common_snp_risk == "none"
    )
    logger.info("Rescue: %d clean primers added for %s", rescue_clean, ", ".join(deficient))

    updated_primer_records = list(primer_records) + rescue_primer_records
    updated_spec_records = list(spec_records) + rescue_spec_records

    return updated_primer_records, updated_spec_records


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

    # --- Rescue ---
    p.add_argument("--rescue-target", nargs="+", default=None,
                   help="Target IDs to rescue (e.g. DENND3_cds18). Re-runs Stage 2+3 "
                        "with relaxed parameters for the specified targets.")
    p.add_argument("--rescue-flank", type=int, default=None,
                   help="Template extension for rescue (default: same as --primer-flank).")
    p.add_argument("--rescue-num-return", type=int, default=20,
                   help="PRIMER_NUM_RETURN for rescue runs (default: 20).")
    p.add_argument("--rescue-min-tm", type=float, default=None,
                   help="Tighter PRIMER_MIN_TM for rescue (default: same as --primer-min-tm).")
    p.add_argument("--rescue-max-tm", type=float, default=None,
                   help="Tighter PRIMER_MAX_TM for rescue (default: same as --primer-max-tm).")
    p.add_argument("--rescue-min-gc", type=float, default=None,
                   help="Tighter PRIMER_MIN_GC for rescue (default: same as --primer-min-gc).")
    p.add_argument("--rescue-max-gc", type=float, default=None,
                   help="Tighter PRIMER_MAX_GC for rescue (default: same as --primer-max-gc).")
    p.add_argument("--rescue-min-size", type=int, default=None,
                   help="Tighter PRIMER_MIN_SIZE for rescue (default: same as --primer-min-size).")
    p.add_argument("--rescue-opt-size", type=int, default=None,
                   help="Tighter PRIMER_OPT_SIZE for rescue (default: same as --primer-opt-size).")
    p.add_argument("--rescue-max-size", type=int, default=None,
                   help="Tighter PRIMER_MAX_SIZE for rescue (default: same as --primer-max-size).")

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

    # --- Mutual exclusion: --rescue-target vs --genes ---
    if args.rescue_target and args.genes:
        logger.error("--rescue-target and --genes are mutually exclusive")
        sys.exit(1)

    # --- Resolve stages (shared by both modes) ---
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
        # Rescue
        rescue_flank=args.rescue_flank if args.rescue_flank is not None else args.primer_flank,
        rescue_num_return=args.rescue_num_return,
        rescue_min_tm=args.rescue_min_tm,
        rescue_max_tm=args.rescue_max_tm,
        rescue_min_gc=args.rescue_min_gc,
        rescue_max_gc=args.rescue_max_gc,
        rescue_min_size=args.rescue_min_size,
        rescue_opt_size=args.rescue_opt_size,
        rescue_max_size=args.rescue_max_size,
    )

    # --- Standalone rescue mode ---
    if args.rescue_target:
        _run_rescue_standalone(args, cfg)
        return

    # --- Validate --genes is required for pipeline runs ---
    if not args.genes:
        logger.error("--genes is required (e.g. --genes HFE HJV TFR2)")
        sys.exit(1)

    # Deduplicate genes (preserve order)
    seen: set[str] = set()
    deduped_genes: list[str] = []
    for gene in args.genes:
        if gene not in seen:
            seen.add(gene)
            deduped_genes.append(gene)
    if len(deduped_genes) < len(args.genes):
        duplicates = [g for g in args.genes if args.genes.count(g) > 1]
        dup_set = sorted(set(duplicates))
        logger.warning(
            "Duplicate genes detected and removed: %s. "
            "%d unique genes will be processed.",
            ", ".join(dup_set), len(deduped_genes),
        )
    args.genes = deduped_genes

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
    spec_records = None
    if cfg.check_specificity:
        if primer_records is None:
            logger.error("--check-specificity requires --design-primers")
            sys.exit(1)
        spec_records = _run_stage3(all_records, primer_records, cfg)

    # --- QC summary ---
    qc_path = cfg.output_dir / "run_summary.txt"
    write_qc_summary(all_records, all_failed, primer_records, cfg, qc_path)

    logger.info("QC summary → %s", qc_path)


if __name__ == "__main__":
    main()
