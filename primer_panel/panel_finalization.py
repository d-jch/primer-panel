"""Panel Finalization: select recommended primers and rescue failed targets.

Reads Stage 2+3 outputs and selects the best unique primer per target.
Experimental rescue attempts only run when explicitly requested.

Usage:
    primer-panel-finalize \\
        --input-dir outputs/hcc6_primers \\
        --output-dir outputs/panel_final \\
        --genome-fasta /path/to/hg38.fa
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
import time
from collections import defaultdict
from pathlib import Path

logger = logging.getLogger("panel_finalization")


# ----------------------------------------------------------------------
# Data loading
# ----------------------------------------------------------------------

def load_specificity_tsv(path: Path) -> list[dict]:
    """Load primers_specificity.tsv."""
    rows = []
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            rows.append(row)
    return rows


def load_target_summary(path: Path) -> dict[str, dict]:
    """Load target_summary.tsv keyed by target_id."""
    targets = {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            targets[row["target_id"]] = row
    return targets


# ----------------------------------------------------------------------
# Recommended primer selection
# ----------------------------------------------------------------------

def select_recommended(rows: list[dict]) -> tuple[list[dict], list[str]]:
    """Select best unique_pass primer per target.

    Returns (recommended_list, failed_target_ids).
    """
    # Group by target_id, filter to unique_pass
    by_target: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        if row.get("insilico_status") == "unique_pass" and row.get("primer3_status") == "ok":
            by_target[row["target_id"]].append(row)

    all_targets = sorted(set(row["target_id"] for row in rows))
    recommended = []
    failed = []

    for tid in all_targets:
        candidates = by_target.get(tid, [])
        if candidates:
            # Sort by primer_pair_penalty (lowest = best)
            candidates.sort(key=lambda r: float(r.get("primer_pair_penalty", "999")))
            best = candidates[0]
            recommended.append(best)
        else:
            failed.append(tid)

    return recommended, failed


# ----------------------------------------------------------------------
# FTH1 rescue analysis
# ----------------------------------------------------------------------

def _parse_optional_int(value: object) -> int | None:
    if value in (None, ""):
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def analyze_fth1_multi_hit(rows: list[dict], expected_tolerance: int = 1000) -> dict:
    """Analyze FTH1_cds1_4 multi-hit pattern."""
    fth1_rows = [r for r in rows if r["target_id"] == "FTH1_cds1_4"]

    expected_intervals = []
    for row in fth1_rows:
        chrom = row.get("expected_target_chrom", "")
        start = _parse_optional_int(row.get("expected_target_start"))
        end = _parse_optional_int(row.get("expected_target_end"))
        if chrom and start is not None and end is not None:
            expected_intervals.append((chrom, start, end))

    # Collect all unique hit locations
    all_locations: dict[str, set] = defaultdict(set)  # chrom -> set of (start, end, size)
    for row in fth1_rows:
        hits_str = row.get("insilico_hits", "")
        if not hits_str:
            continue
        for hit in hits_str.split(";"):
            # Parse 'chr11:61964652-61967597(2945)'
            try:
                loc_part, size_part = hit.split("(")
                chrom, coords = loc_part.split(":")
                start, end = coords.split("-")
                size = int(size_part.rstrip(")"))
                all_locations[chrom].add((int(start), int(end), size))
            except ValueError:
                continue

    def is_expected_hit(chrom: str, start: int, end: int) -> bool:
        if expected_intervals:
            return any(
                hit_chrom == chrom
                and start <= interval_end + expected_tolerance
                and end >= interval_start - expected_tolerance
                for hit_chrom, interval_start, interval_end in expected_intervals
            )

        return chrom == "chr11" and abs(start - 61964652) < expected_tolerance

    # Separate expected locus hits from off-targets.
    expected_hits = set()
    off_target_by_chrom: dict[str, list] = defaultdict(list)

    for chrom, locs in all_locations.items():
        for start, end, size in sorted(locs):
            if is_expected_hit(chrom, start, end):
                expected_hits.add((start, end, size))
            else:
                off_target_by_chrom[chrom].append((start, end, size))

    return {
        "total_primers": len(fth1_rows),
        "expected_hits": expected_hits,
        "off_target_by_chrom": dict(off_target_by_chrom),
        "off_target_chroms": sorted(off_target_by_chrom.keys()),
    }


# ----------------------------------------------------------------------
# Rescue: re-run Primer3 with different parameters
# ----------------------------------------------------------------------

def rescue_fth1_with_more_primers(
    target_summary: dict[str, dict],
    output_dir: Path,
    genome_fasta: Path,
    primer_num_return: int = 50,
    primer_min_size: int = 20,
    primer_opt_size: int = 22,
    primer_max_size: int = 28,
    primer_opt_tm: float = 61.0,
    primer_min_tm: float = 59.0,
    primer_max_tm: float = 65.0,
    primer_min_gc: float = 40.0,
    primer_max_gc: float = 60.0,
) -> dict:
    """Re-design FTH1_cds1_4 with stricter primer parameters and more candidates.

    Returns dict with rescue results.
    """
    from primer_panel.config import PipelineConfig
    from primer_panel.primer3_runner import (
        check_primer3_available,
        run_primer3_for_target,
        is_all_n,
    )
    from primer_panel.writers import Target, build_records
    from primer_panel.ensembl_client import TranscriptInfo

    rescue_output = output_dir / "rescue_fth1"
    rescue_output.mkdir(parents=True, exist_ok=True)

    # Get FTH1 target info
    fth1_info = target_summary.get("FTH1_cds1_4")
    if not fth1_info:
        return {"status": "error", "detail": "FTH1_cds1_4 not found in target_summary"}

    # Build Target object
    cds_coords = []
    for part in fth1_info["cds_exon_coords"].split(","):
        s, e = part.split("-")
        cds_coords.append((int(s), int(e)))

    cds_numbers = [int(n) for n in fth1_info["cds_exon_numbers"].split(",")]

    t = Target(
        chrom=fth1_info["required_chrom"],
        start=int(fth1_info["required_start"]),
        end=int(fth1_info["required_end"]),
        strand=1 if fth1_info["strand"] == "+" else -1,
        cds_exon_numbers=cds_numbers,
        cds_exon_ids=fth1_info.get("cds_exon_ids", "").split(","),
        cds_exon_coords=cds_coords,
    )

    ti = TranscriptInfo(
        transcript_id=fth1_info["transcript_id"],
        biotype="protein_coding",
        is_mane_select=False,
        is_mane_plus_clinical=False,
        is_canonical=True,
        selection_reason=fth1_info.get("selection_reason", ""),
        exons=[],
        cds_exons=[],
    )

    # Build config with stricter parameters
    cfg = PipelineConfig(
        product_min=2700,
        product_max=3300,
        primer_flank=300,
        genome_fasta=genome_fasta,
        output_dir=rescue_output,
        design_primers=True,
        primer_num_return=primer_num_return,
        primer_opt_size=primer_opt_size,
        primer_min_size=primer_min_size,
        primer_max_size=primer_max_size,
        primer_opt_tm=primer_opt_tm,
        primer_min_tm=primer_min_tm,
        primer_max_tm=primer_max_tm,
        primer_max_tm_diff=2.0,
        primer_min_gc=primer_min_gc,
        primer_max_gc=primer_max_gc,
    )

    # Build records to get sequence_target coordinates
    records = build_records(
        "FTH1", ti, [t], cfg, "real",
        gene_required_start=int(fth1_info["required_start"]),
        gene_required_end=int(fth1_info["required_end"]),
    )
    rec = records[0]

    # Extract real sequence from genome FASTA
    try:
        from pyfaidx import Fasta
        genome = Fasta(str(genome_fasta))
        seq = str(genome[rec.template_chrom][rec.template_start:rec.template_end])
    except Exception as exc:
        return {"status": "error", "detail": f"Cannot extract sequence: {exc}"}

    if is_all_n(seq):
        return {"status": "error", "detail": "Sequence is all N"}

    # Run Primer3
    if not check_primer3_available(cfg.primer3_bin):
        return {"status": "error", "detail": "primer3_core not found"}

    logger.info("Rescue FTH1_cds1_4: running Primer3 with %d return, stricter params", primer_num_return)

    results = run_primer3_for_target(
        target_id="FTH1_cds1_4",
        sequence=seq,
        cfg=cfg,
        seq_target_start=rec.sequence_target_start_0based,
        seq_target_len=rec.sequence_target_length,
    )

    ok_results = [r for r in results if r.primer3_status == "ok"]
    logger.info("  Primer3 returned %d ok pairs", len(ok_results))

    # Write rescue primers TSV
    rescue_primers_path = rescue_output / "rescue_primers.tsv"
    with open(rescue_primers_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "target_id", "primer_rank", "forward_primer", "reverse_primer",
            "forward_tm", "reverse_tm", "tm_diff", "forward_gc", "reverse_gc",
            "primer_pair_penalty", "primer3_product_size",
            "primer_left_start", "primer_left_len", "primer_right_start", "primer_right_len",
        ])
        for r in ok_results:
            writer.writerow([
                "FTH1_cds1_4", r.primer_rank,
                r.forward_primer, r.reverse_primer,
                r.forward_tm, r.reverse_tm, r.tm_diff,
                r.forward_gc, r.reverse_gc,
                r.primer_pair_penalty, r.primer3_product_size,
                r.primer_left_start, r.primer_left_len,
                r.primer_right_start, r.primer_right_len,
            ])

    return {
        "status": "ok",
        "primer_count": len(ok_results),
        "rescue_primers_path": str(rescue_primers_path),
        "results": ok_results,
        "config": {
            "primer_num_return": primer_num_return,
            "primer_min_size": primer_min_size,
            "primer_opt_size": primer_opt_size,
            "primer_max_size": primer_max_size,
            "primer_opt_tm": primer_opt_tm,
            "primer_min_tm": primer_min_tm,
            "primer_max_tm": primer_max_tm,
            "primer_min_gc": primer_min_gc,
            "primer_max_gc": primer_max_gc,
        },
    }


def run_stage3_on_rescue(
    rescue_primers_path: Path,
    target_summary_path: Path,
    genome_fasta: Path,
    output_dir: Path,
    is_pcr_bin: str = "isPcr",
    tolerance: int = 10,
) -> dict:
    """Run Stage 3 on rescue primers."""
    from primer_panel.insilico_pcr import check_ispcr_available, check_specificity_batch

    if not check_ispcr_available(is_pcr_bin):
        return {"status": "error", "detail": "isPcr not found"}

    # Load rescue primers
    primers = []
    with open(rescue_primers_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            primers.append(row)

    if not primers:
        return {"status": "error", "detail": "No rescue primers found"}

    # Load target info for expected coords
    target_info = {}
    with open(target_summary_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            target_info[row["target_id"]] = row

    # Build batch
    primer_batch = []
    for p in primers:
        target_id = p.get("target_id", "")
        tgt = target_info.get(target_id)
        if tgt is None:
            return {"status": "error", "detail": f"{target_id} not found in target_summary"}
        template_chrom = tgt.get("template_chrom", "")
        if not template_chrom:
            return {"status": "error", "detail": f"{target_id} missing template_chrom"}
        template_start = int(tgt.get("template_start", 0))
        left_start = int(p.get("primer_left_start", 0))
        right_start = int(p.get("primer_right_start", 0))
        primer_batch.append({
            "name": f"rescue_rank{p['primer_rank']}",
            "fwd": p["forward_primer"],
            "rev": p["reverse_primer"],
            "expected_chrom": template_chrom,
            "expected_start": template_start + left_start,
            "expected_end": template_start + right_start + 1,
        })

    logger.info("Running isPcr on %d rescue primers ...", len(primer_batch))

    t0 = time.time()
    results = check_specificity_batch(
        primer_pairs=primer_batch,
        genome_fasta=str(genome_fasta),
        ispcr_bin=is_pcr_bin,
        tolerance=tolerance,
    )
    elapsed = time.time() - t0
    logger.info("isPcr done in %.0fs", elapsed)

    # Analyze results
    unique_pass = []
    multi_hit = []
    no_hit = []
    primers_by_rank = {p["primer_rank"]: p for p in primers}
    for name, r in results.items():
        rank = int(name.replace("rescue_rank", ""))
        primer = primers_by_rank.get(str(rank))
        if primer is None:
            continue

        entry = {
            "target_id": primer.get("target_id", ""),
            "rank": rank,
            "fwd": primer["forward_primer"],
            "rev": primer["reverse_primer"],
            "fwd_tm": primer["forward_tm"],
            "rev_tm": primer["reverse_tm"],
            "penalty": primer["primer_pair_penalty"],
            "product_size": primer["primer3_product_size"],
            "insilico_status": r.insilico_status,
            "hit_count": r.insilico_hit_count,
            "specificity_explain": r.specificity_explain,
            "best_chrom": r.insilico_best_chrom,
            "best_start": r.insilico_best_start,
            "best_end": r.insilico_best_end,
            "best_size": r.insilico_best_size,
        }

        if r.insilico_status == "unique_pass":
            unique_pass.append(entry)
        elif r.insilico_status == "multi_hit":
            multi_hit.append(entry)
        elif r.insilico_status == "no_hit":
            no_hit.append(entry)

    return {
        "status": "ok",
        "total": len(primers),
        "unique_pass": unique_pass,
        "multi_hit": multi_hit,
        "no_hit": no_hit,
        "elapsed": elapsed,
    }


# ----------------------------------------------------------------------
# Split target rescue: FTH1_cds1_2 + FTH1_cds3_4
# ----------------------------------------------------------------------

def rescue_fth1_split_targets(
    genome_fasta: Path,
    output_dir: Path,
    is_pcr_bin: str = "isPcr",
) -> dict:
    """Try splitting FTH1_cds1_4 into FTH1_cds1_2 and FTH1_cds3_4."""
    from primer_panel.config import PipelineConfig
    from primer_panel.ensembl_client import EnsemblClient
    from primer_panel.cds_handler import build_required_intervals
    from primer_panel.target_planner_adapter import plan_targets_with_external_planner
    from primer_panel.writers import build_records, write_fasta
    from primer_panel.primer3_runner import (
        check_primer3_available,
        run_primer3_for_target,
        parse_fasta_sequences,
        is_all_n,
    )
    from primer_panel.insilico_pcr import check_ispcr_available, check_specificity_batch

    rescue_dir = output_dir / "rescue_fth1_split"
    rescue_dir.mkdir(parents=True, exist_ok=True)

    # Use pipeline to get FTH1 CDS exons
    cfg = PipelineConfig(
        product_min=2700,
        product_max=3300,
        primer_flank=300,
        genome_fasta=genome_fasta,
        output_dir=rescue_dir,
        primer_num_return=50,
        primer_min_size=20,
        primer_opt_size=22,
        primer_max_size=28,
        primer_opt_tm=61.0,
        primer_min_tm=59.0,
        primer_max_tm=65.0,
        primer_min_gc=40.0,
        primer_max_gc=60.0,
    )

    client = EnsemblClient(cfg)
    ti, reason = client.select_transcript("FTH1")

    if not ti.cds_exons:
        return {"status": "error", "detail": "No CDS exons found for FTH1"}

    logger.info("FTH1: %d CDS exons", len(ti.cds_exons))

    # Build required intervals
    required_intervals = build_required_intervals(ti.cds_exons)

    # Use a smaller target window to force splitting into smaller targets
    cfg_split = PipelineConfig(
        product_min=1200,
        product_max=1500,
        primer_flank=300,
        genome_fasta=genome_fasta,
        output_dir=rescue_dir,
        primer_num_return=50,
        primer_min_size=20,
        primer_opt_size=22,
        primer_max_size=28,
        primer_opt_tm=61.0,
        primer_min_tm=59.0,
        primer_max_tm=65.0,
        primer_min_gc=40.0,
        primer_max_gc=60.0,
    )

    gene_req_start = min(ri.start for ri in required_intervals)
    gene_req_end = max(ri.end for ri in required_intervals)
    gene_data = client.lookup_gene("FTH1")
    gene_bounds_start = gene_data.get("start", gene_req_start + 1) - 1
    gene_bounds_end = gene_data.get("end", gene_req_end)
    targets = plan_targets_with_external_planner(
        required_intervals,
        cfg_split,
        gene_start=gene_bounds_start,
        gene_end=gene_bounds_end,
    )
    logger.info("FTH1 with split: %d targets generated", len(targets))

    if len(targets) < 2:
        return {"status": "error", "detail": f"Split produced only {len(targets)} targets, need >= 2"}

    # Build records
    records = build_records(
        "FTH1", ti, targets, cfg_split, "real",
        gene_required_start=gene_req_start,
        gene_required_end=gene_req_end,
    )

    logger.info("Split targets:")
    for rec in records:
        logger.info("  %s: required=%s:%d-%d (%dbp), extended=%d-%d (%dbp), template=%d-%d (%dbp)",
                    rec.target_id, rec.required_chrom, rec.required_start, rec.required_end,
                    rec.required_length, rec.extended_start, rec.extended_end,
                    rec.extended_length, rec.template_start, rec.template_end,
                    rec.template_length)

    # Write FASTA and extract sequences
    fa_path = rescue_dir / "split_targets.fa"
    write_fasta(records, cfg_split, fa_path)
    sequences = parse_fasta_sequences(fa_path)

    # Check primer3 and isPcr availability
    if not check_primer3_available(cfg_split.primer3_bin):
        return {"status": "error", "detail": "primer3_core not found"}
    if not check_ispcr_available(is_pcr_bin):
        return {"status": "error", "detail": "isPcr not found"}

    # Run Primer3 on each split target
    all_results = {}
    for rec in records:
        seq = sequences.get(rec.target_id)
        if seq is None or is_all_n(seq):
            logger.warning("No valid sequence for %s", rec.target_id)
            continue

        logger.info("Running Primer3 for %s ...", rec.target_id)
        p3_results = run_primer3_for_target(
            target_id=rec.target_id,
            sequence=seq,
            cfg=cfg_split,
            seq_target_start=rec.sequence_target_start_0based,
            seq_target_len=rec.sequence_target_length,
        )
        ok_results = [r for r in p3_results if r.primer3_status == "ok"]
        logger.info("  %s: %d ok primer pairs", rec.target_id, len(ok_results))
        all_results[rec.target_id] = (rec, ok_results)

    # Run isPcr on all primers
    primer_batch = []
    for tid, (rec, results) in all_results.items():
        for r in results:
            exp_start = rec.template_start + r.primer_left_start
            exp_end = rec.template_start + r.primer_right_start + 1
            primer_batch.append({
                "name": f"{tid}_rank{r.primer_rank}",
                "fwd": r.forward_primer,
                "rev": r.reverse_primer,
                "expected_chrom": rec.template_chrom,
                "expected_start": exp_start,
                "expected_end": exp_end,
            })

    logger.info("Running isPcr on %d split primers ...", len(primer_batch))
    t0 = time.time()
    specificity_results = check_specificity_batch(
        primer_pairs=primer_batch,
        genome_fasta=str(genome_fasta),
        ispcr_bin=is_pcr_bin,
        tolerance=10,
    )
    elapsed = time.time() - t0
    logger.info("isPcr done in %.0fs", elapsed)

    # Analyze per-target
    split_summary = {}
    for tid, (rec, p3_results) in all_results.items():
        unique = []
        multi = []
        no_hit = []
        for r in p3_results:
            name = f"{tid}_rank{r.primer_rank}"
            sp = specificity_results.get(name)
            if sp is None:
                continue
            entry = {
                "rank": r.primer_rank,
                "fwd": r.forward_primer,
                "rev": r.reverse_primer,
                "penalty": r.primer_pair_penalty,
                "product_size": r.primer3_product_size,
                "status": sp.insilico_status,
                "hit_count": sp.insilico_hit_count,
                "explain": sp.specificity_explain,
            }
            if sp.insilico_status == "unique_pass":
                unique.append(entry)
            elif sp.insilico_status == "multi_hit":
                multi.append(entry)
            else:
                no_hit.append(entry)

        split_summary[tid] = {
            "total": len(p3_results),
            "unique_pass": unique,
            "multi_hit": multi,
            "no_hit": no_hit,
        }

    return {
        "status": "ok",
        "targets": split_summary,
        "elapsed": elapsed,
    }


# ----------------------------------------------------------------------
# Output writers
# ----------------------------------------------------------------------

def write_recommended_primers(recommended: list[dict], output_dir: Path) -> None:
    """Write recommended_primers.tsv and .xlsx."""
    cols = [
        "target_id", "primer_rank", "forward_primer", "reverse_primer",
        "forward_tm", "reverse_tm", "tm_diff", "forward_gc", "reverse_gc",
        "primer_pair_penalty", "primer3_product_size",
        "insilico_best_chrom", "insilico_best_start", "insilico_best_end",
        "insilico_best_size", "insilico_hit_count", "specificity_explain",
        # Common dbSNP annotation columns
        "common_snp_risk", "left_primer_common_snp_count",
        "right_primer_common_snp_count", "left_primer_3p_common_snp_count",
        "right_primer_3p_common_snp_count", "common_snp_hits",
    ]

    tsv_path = output_dir / "recommended_primers.tsv"
    with open(tsv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for r in recommended:
            writer.writerow(r)
    logger.info("Recommended primers TSV -> %s (%d rows)", tsv_path, len(recommended))

    # XLSX
    try:
        from openpyxl import Workbook
        wb = Workbook()
        ws = wb.active
        ws.title = "Recommended Primers"
        ws.append(cols)
        for r in recommended:
            ws.append([r.get(c, "") for c in cols])
        xlsx_path = output_dir / "recommended_primers.xlsx"
        wb.save(xlsx_path)
        logger.info("Recommended primers XLSX -> %s", xlsx_path)
    except ImportError:
        logger.info("openpyxl not installed -- skipping XLSX")


def write_failed_targets(failed_tids: list[str], rows: list[dict], output_dir: Path) -> None:
    """Write failed_or_needs_review_targets.tsv."""
    path = output_dir / "failed_or_needs_review_targets.tsv"
    with open(path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["target_id", "reason", "primer_count", "unique_pass_count",
                         "multi_hit_count", "no_hit_count", "best_explain"])
        for tid in failed_tids:
            tid_rows = [r for r in rows if r["target_id"] == tid]
            unique = sum(1 for r in tid_rows if r.get("insilico_status") == "unique_pass")
            multi = sum(1 for r in tid_rows if r.get("insilico_status") == "multi_hit")
            no_hit = sum(1 for r in tid_rows if r.get("insilico_status") == "no_hit")
            best_explain = tid_rows[0].get("specificity_explain", "") if tid_rows else ""
            reason = "no_unique_primer"
            writer.writerow([tid, reason, len(tid_rows), unique, multi, no_hit, best_explain])
    logger.info("Failed targets -> %s (%d targets)", path, len(failed_tids))


def write_panel_summary(
    recommended: list[dict],
    failed_tids: list[str],
    rescue_results: dict,
    output_dir: Path,
) -> None:
    """Write panel_summary.txt."""
    lines = []
    lines.append("# Panel Finalization Summary")
    lines.append("")

    # Overall stats
    total_targets = len(recommended) + len(failed_tids)
    lines.append(f"total_targets\t{total_targets}")
    lines.append(f"targets_with_unique_primer\t{len(recommended)}")
    lines.append(f"targets_without_unique_primer\t{len(failed_tids)}")
    lines.append("")

    # Recommended primers
    lines.append("# Recommended Primers")
    lines.append("target_id\tprimer_rank\tfwd\trev\tpenalty\tproduct_size\tbest_hit")
    for r in recommended:
        lines.append(
            f"{r['target_id']}\t{r['primer_rank']}\t{r['forward_primer']}\t{r['reverse_primer']}"
            f"\t{r.get('primer_pair_penalty', '')}\t{r.get('primer3_product_size', '')}"
            f"\t{r.get('insilico_best_chrom', '')}:{r.get('insilico_best_start', '')}"
            f"-{r.get('insilico_best_end', '')}"
        )
    lines.append("")

    # Failed targets
    if failed_tids:
        lines.append("# Targets Without Unique Primer")
        for tid in failed_tids:
            lines.append(f"  {tid}")
        lines.append("")

    # Rescue results
    if rescue_results:
        lines.append("# Rescue Attempts")
        for attempt_name, result in rescue_results.items():
            lines.append(f"\n## {attempt_name}")
            lines.append(f"  status: {result.get('status', 'unknown')}")

            if result.get("status") == "ok":
                if "unique_pass" in result:
                    # Single-target rescue
                    lines.append(f"  total_primers: {result.get('total', 0)}")
                    lines.append(f"  unique_pass: {len(result.get('unique_pass', []))}")
                    lines.append(f"  multi_hit: {len(result.get('multi_hit', []))}")
                    lines.append(f"  no_hit: {len(result.get('no_hit', []))}")
                    if result.get("unique_pass"):
                        lines.append("  RESCUED -- unique primers found:")
                        for p in result["unique_pass"][:3]:
                            lines.append(f"    rank {p['rank']}: {p['fwd']} / {p['rev']} (penalty={p['penalty']})")
                elif "targets" in result:
                    # Split-target rescue
                    for tid, stats in result["targets"].items():
                        lines.append(f"\n  {tid}:")
                        lines.append(f"    total: {stats['total']}")
                        lines.append(f"    unique_pass: {len(stats['unique_pass'])}")
                        lines.append(f"    multi_hit: {len(stats['multi_hit'])}")
                        if stats["unique_pass"]:
                            lines.append("    RESCUED:")
                            for p in stats["unique_pass"][:3]:
                                lines.append(f"      rank {p['rank']}: {p['fwd']} / {p['rev']} (penalty={p['penalty']})")
                        else:
                            lines.append("    NOT rescued -- no unique primers")
            else:
                lines.append(f"  error: {result.get('detail', 'unknown')}")

    path = output_dir / "panel_summary.txt"
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    logger.info("Panel summary -> %s", path)


def write_rescue_attempts(rescue_results: dict, output_dir: Path) -> None:
    """Write rescue_attempts.tsv with all rescue primer results."""
    path = output_dir / "rescue_attempts.tsv"
    with open(path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "attempt", "target_id", "rank", "fwd", "rev", "penalty",
            "product_size", "insilico_status", "hit_count", "specificity_explain",
        ])
        for attempt_name, result in rescue_results.items():
            if result.get("status") != "ok":
                continue

            if "unique_pass" in result:
                # Single-target rescue
                for status_key in ["unique_pass", "multi_hit", "no_hit"]:
                    for p in result.get(status_key, []):
                        writer.writerow([
                            attempt_name, p.get("target_id", "FTH1_cds1_4"), p["rank"],
                            p["fwd"], p["rev"], p["penalty"], p["product_size"],
                            p["insilico_status"], p["hit_count"], p["specificity_explain"],
                        ])
            elif "targets" in result:
                # Split-target rescue
                for tid, stats in result["targets"].items():
                    for status_key in ["unique_pass", "multi_hit", "no_hit"]:
                        for p in stats.get(status_key, []):
                            writer.writerow([
                                attempt_name, tid, p["rank"],
                                p["fwd"], p["rev"], p["penalty"], p["product_size"],
                                p["status"], p["hit_count"], p["explain"],
                            ])

    logger.info("Rescue attempts -> %s", path)


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Panel finalization")
    parser.add_argument("--input-dir", type=Path, required=True,
                        help="Directory with Stage 2+3 outputs")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="Output directory for finalization results")
    parser.add_argument("--genome-fasta", type=Path, required=True,
                        help="Path to hg38 FASTA")
    parser.add_argument("--is-pcr-bin", type=str, default="isPcr",
                        help="Path to isPcr binary")
    parser.add_argument("--rescue-target", action="append", default=[],
                        help="Explicit failed target to rescue. Currently only FTH1_cds1_4 has an implemented rescue strategy.")
    parser.add_argument("--rescue-all", action="store_true",
                        help="Run implemented rescue strategies for all supported failed targets.")
    parser.add_argument("--experimental-split-rescue", action="store_true",
                        help="Also try the experimental FTH1 split-target rescue.")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Enable debug logging.")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s %(name)s: %(message)s",
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    spec_path = args.input_dir / "primers_specificity.tsv"
    target_path = args.input_dir / "target_summary.tsv"

    if not spec_path.exists():
        logger.error("primers_specificity.tsv not found in %s", args.input_dir)
        sys.exit(1)

    logger.info("Loading specificity results from %s", spec_path)
    rows = load_specificity_tsv(spec_path)
    logger.info("Loaded %d primer records", len(rows))

    target_summary = {}
    if target_path.exists():
        target_summary = load_target_summary(target_path)
        logger.info("Loaded %d targets from %s", len(target_summary), target_path)

    # 1. Select recommended primers
    recommended, failed_tids = select_recommended(rows)
    logger.info("Recommended: %d targets, Failed: %d targets", len(recommended), len(failed_tids))

    requested_rescue = set(failed_tids if args.rescue_all else args.rescue_target)
    requested_rescue &= set(failed_tids)

    # 2. Analyze FTH1 multi-hit only when FTH1 rescue was requested.
    if "FTH1_cds1_4" in requested_rescue:
        logger.info("\n=== FTH1_cds1_4 Multi-hit Analysis ===")
        fth1_analysis = analyze_fth1_multi_hit(rows)
        logger.info("  Off-target chromosomes: %s", fth1_analysis["off_target_chroms"])
        logger.info("  Expected hits: %d locations", len(fth1_analysis["expected_hits"]))
        for chrom, locs in fth1_analysis["off_target_by_chrom"].items():
            logger.info("  %s: %d off-target locations", chrom, len(locs))
    elif "FTH1_cds1_4" in failed_tids:
        logger.info(
            "FTH1_cds1_4 failed; skipping rescue by default. "
            "Pass --rescue-target FTH1_cds1_4 to enable it."
        )

    # 3. Rescue attempts
    rescue_results = {}

    for tid in sorted(tid for tid in requested_rescue if tid != "FTH1_cds1_4"):
        rescue_results[f"unsupported_{tid}"] = {
            "status": "error",
            "detail": f"No implemented rescue strategy for {tid}",
        }

    if "FTH1_cds1_4" in requested_rescue:
        # Attempt 1: More primers with stricter params
        logger.info("\n=== Rescue Attempt 1: More primers (50), stricter params ===")
        rescue1 = rescue_fth1_with_more_primers(
            target_summary, args.output_dir, args.genome_fasta,
            primer_num_return=50,
            primer_min_size=20, primer_opt_size=22, primer_max_size=28,
            primer_opt_tm=61.0, primer_min_tm=59.0, primer_max_tm=65.0,
            primer_min_gc=40.0, primer_max_gc=60.0,
        )

        if rescue1["status"] == "ok" and rescue1["primer_count"] > 0:
            logger.info("  Primer3 returned %d ok pairs, running Stage 3 ...", rescue1["primer_count"])
            s3_result = run_stage3_on_rescue(
                Path(rescue1["rescue_primers_path"]),
                target_path, args.genome_fasta, args.output_dir,
                is_pcr_bin=args.is_pcr_bin,
            )
            rescue_results["attempt1_more_primers_strict"] = s3_result
            if s3_result.get("unique_pass"):
                logger.info("  [OK] RESCUED: %d unique_pass primers found", len(s3_result["unique_pass"]))
            else:
                logger.info("  [FAIL] Not rescued: %d multi_hit, %d no_hit",
                           len(s3_result.get("multi_hit", [])), len(s3_result.get("no_hit", [])))

        if args.experimental_split_rescue:
            # Attempt 2: Split targets
            logger.info("\n=== Rescue Attempt 2: Split FTH1_cds1_4 into smaller targets ===")
            rescue2 = rescue_fth1_split_targets(
                args.genome_fasta, args.output_dir, is_pcr_bin=args.is_pcr_bin,
            )
            rescue_results["attempt2_split_targets"] = rescue2

            if rescue2.get("status") == "ok":
                for tid, stats in rescue2.get("targets", {}).items():
                    if stats["unique_pass"]:
                        logger.info("  rescued %s: %d unique_pass", tid, len(stats["unique_pass"]))
                    else:
                        logger.info(
                            "  %s not rescued: %d multi_hit, %d no_hit",
                            tid, len(stats["multi_hit"]), len(stats["no_hit"]),
                        )

    # 4. Write outputs
    logger.info("\n=== Writing outputs ===")
    write_recommended_primers(recommended, args.output_dir)
    write_failed_targets(failed_tids, rows, args.output_dir)
    write_rescue_attempts(rescue_results, args.output_dir)
    write_panel_summary(recommended, failed_tids, rescue_results, args.output_dir)

    # 5. Print final summary
    print(f"\n{'='*60}")
    print(f"Panel Finalization Results")
    print(f"{'='*60}")
    print(f"Total targets:           {len(recommended) + len(failed_tids)}")
    print(f"With unique primer:      {len(recommended)}")
    print(f"Without unique primer:   {len(failed_tids)}")
    if failed_tids:
        print(f"Failed targets:          {', '.join(failed_tids)}")

    if rescue_results:
        print(f"\nRescue attempts: {len(rescue_results)}")
        for name, result in rescue_results.items():
            status = result.get("status", "unknown")
            if status == "ok":
                if "unique_pass" in result:
                    n = len(result["unique_pass"])
                    print(f"  {name}: {'[OK] RESCUED' if n > 0 else '[FAIL] failed'} ({n} unique)")
                elif "targets" in result:
                    for tid, stats in result["targets"].items():
                        n = len(stats["unique_pass"])
                        print(f"  {name}/{tid}: {'[OK] RESCUED' if n > 0 else '[FAIL] failed'} ({n} unique)")
            else:
                print(f"  {name}: error -- {result.get('detail', '')}")

    print(f"\nOutputs -> {args.output_dir}")


if __name__ == "__main__":
    main()
