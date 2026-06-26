"""Standalone Stage 3 CLI: run in-silico PCR from an existing primers.tsv.

Usage:
    python -m primer_panel.insilico_pcr_cli \
        --primers-tsv outputs/stage2_real_hcc6/primers.tsv \
        --targets-tsv outputs/stage2_real_hcc6/target_summary.tsv \
        --genome-fasta /mnt/e/hg38/genome.fa \
        --output-dir outputs/stage3_from_existing

This reads an existing primers.tsv (from Stage 2) and target_summary.tsv
(from Stage 1), runs isPcr on all ok primer pairs, and writes specificity
results without re-running Stage 1 or Stage 2.  No Ensembl API access required.

Stage 3 does NOT filter hits by product size.  All isPcr hits are retained.
Stage 3 is the ONLY stage that produces genomic product coordinates.
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
import time
from pathlib import Path

logger = logging.getLogger("primer_panel.insilico_pcr_cli")


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="primer_panel.insilico_pcr_cli",
        description="Run Stage 3 in-silico PCR from an existing primers.tsv.",
    )
    p.add_argument("--primers-tsv", type=Path, required=True,
                    help="Path to primers.tsv from Stage 2.")
    p.add_argument("--targets-tsv", type=Path, default=None,
                    help="Path to target_summary.tsv from Stage 1 (for expected coords).")
    p.add_argument("--genome-fasta", type=Path, required=True,
                    help="Path to hg38 FASTA (bgzipped or plain).")
    p.add_argument("--output-dir", type=Path, required=True,
                    help="Output directory for specificity results.")
    p.add_argument("--is-pcr-bin", type=str, default="isPcr",
                    help="Path to isPcr binary (default: isPcr).")
    p.add_argument("--pcr-tolerance", type=int, default=10,
                    help="Bp tolerance for coordinate matching (default: 10).")
    p.add_argument("--ispcr-db", type=Path, default=None,
                    help="Explicit .2bit/.nib database for isPcr (default: auto-discover).")
    p.add_argument("--ispcr-ooc", type=Path, default=None,
                    help="Explicit overused-tile (.ooc) file for isPcr (default: auto-discover).")
    p.add_argument("--ispcr-tile-size", type=int, default=11,
                    help="Tile size for isPcr (default: 11).")
    p.add_argument("--prepare-ispcr-db", action="store_true",
                    help="Create a .2bit database from the genome FASTA before running isPcr.")
    p.add_argument("--make-ispcr-ooc", action="store_true",
                    help="Create an overused-tile (.ooc) file before running isPcr.")
    p.add_argument("--product-size", type=str, default=None, metavar="MIN-MAX",
                    help="[DEPRECATED, ignored] Stage 3 does not filter by product size. "
                         "Kept for backward compatibility only.")
    p.add_argument("-v", "--verbose", action="store_true",
                    help="Enable debug logging.")
    return p.parse_args(argv)


def _load_target_coords(targets_tsv: Path) -> dict[str, dict]:
    """Load target_summary.tsv to get template_start for genomic coord computation."""
    coords: dict[str, dict] = {}
    with open(targets_tsv) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            coords[row["target_id"]] = {
                "template_chrom": row.get("template_chrom", ""),
                "template_start": int(row.get("template_start", 0)),
                "extended_chrom": row.get("extended_chrom", ""),
                "extended_start": int(row.get("extended_start", 0)),
                "extended_end": int(row.get("extended_end", 0)),
            }
    return coords


def main(argv: list[str] | None = None) -> None:
    args = _parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s %(name)s: %(message)s",
    )

    if args.product_size is not None:
        logger.warning("--product-size is deprecated and ignored. "
                        "Stage 3 does not filter hits by product size.")

    from .insilico_pcr import (
        check_ispcr_available,
        check_specificity_batch,
        make_ispcr_ooc,
        prepare_ispcr_twobit,
    )
    from .stage3_inputs import build_stage3_inputs_from_target_coords
    from .writers import (
        PrimerRecord,
        build_specificity_records,
        write_specificity_tsv,
        write_unique_primers,
        write_specificity_summary,
    )

    # Check isPcr
    if not check_ispcr_available(args.is_pcr_bin):
        logger.error("%s not found; pass --is-pcr-bin or install isPcr", args.is_pcr_bin)
        sys.exit(1)

    # Check inputs
    if not args.primers_tsv.exists():
        logger.error("primers.tsv not found: %s", args.primers_tsv)
        sys.exit(1)
    if not args.genome_fasta.exists():
        logger.error("genome FASTA not found: %s", args.genome_fasta)
        sys.exit(1)

    # Load target coords if available (for expected region and template_start)
    target_coords: dict[str, dict] = {}
    if args.targets_tsv and args.targets_tsv.exists():
        target_coords = _load_target_coords(args.targets_tsv)
        logger.info("Loaded %d target coords from %s", len(target_coords), args.targets_tsv)

    # Load primer records from TSV
    logger.info("Loading primers from %s", args.primers_tsv)
    primer_records: list[PrimerRecord] = []
    with open(args.primers_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            # Support both new and legacy TSV formats
            try:
                rec = PrimerRecord(
                    target_id=row.get("target_id", ""),
                    primer_rank=int(row.get("primer_rank", 0)),
                    forward_primer=row.get("forward_primer", ""),
                    reverse_primer=row.get("reverse_primer", ""),
                    forward_tm=float(row.get("forward_tm", 0)),
                    reverse_tm=float(row.get("reverse_tm", 0)),
                    tm_diff=float(row.get("tm_diff", 0)),
                    forward_gc=float(row.get("forward_gc", 0)),
                    reverse_gc=float(row.get("reverse_gc", 0)),
                    primer_pair_penalty=float(row.get("primer_pair_penalty", 0)),
                    primer_left_start=int(row.get("primer_left_start", 0)),
                    primer_left_len=int(row.get("primer_left_len", 0)),
                    primer_right_start=int(row.get("primer_right_start", 0)),
                    primer_right_len=int(row.get("primer_right_len", 0)),
                    primer3_product_size=int(row.get(
                        "primer3_product_size",
                        row.get("primer_product_size", 0),  # legacy compat
                    )),
                    primer3_status=row.get("primer3_status", ""),
                    primer3_explain=row.get("primer3_explain", ""),
                    sequence_target_start_0based=int(row.get("sequence_target_start_0based", 0)),
                    sequence_target_length=int(row.get("sequence_target_length", 0)),
                )
            except (ValueError, KeyError) as exc:
                logger.warning("Skipping row: %s", exc)
                continue
            primer_records.append(rec)

    ok_primers = [pr for pr in primer_records if pr.primer3_status == "ok"]
    logger.info("Loaded %d records, %d ok primers", len(primer_records), len(ok_primers))

    if not ok_primers:
        logger.warning("No ok primers found — nothing to check")
        return

    primer_batch, expected_coords = build_stage3_inputs_from_target_coords(
        target_coords,
        primer_records,
    )

    # Optionally prepare isPcr database files (explicit opt-in only)
    ispcr_db_path = str(args.ispcr_db) if args.ispcr_db else None
    ispcr_ooc_path = str(args.ispcr_ooc) if args.ispcr_ooc else None

    if args.prepare_ispcr_db:
        logger.info("Preparing .2bit database from %s …", args.genome_fasta)
        twobit = prepare_ispcr_twobit(args.genome_fasta)
        ispcr_db_path = str(twobit)
        logger.info(".2bit database → %s", twobit)

    if args.make_ispcr_ooc and ispcr_db_path:
        logger.info("Creating .ooc file (tileSize=%d) …", args.ispcr_tile_size)
        ooc = make_ispcr_ooc(ispcr_db_path, tile_size=args.ispcr_tile_size)
        ispcr_ooc_path = str(ooc)
        logger.info(".ooc file → %s", ooc)

    # Run isPcr (no product-size filtering)
    logger.info("Running isPcr on %d primers …", len(primer_batch))
    t0 = time.time()
    specificity_results = check_specificity_batch(
        primer_pairs=primer_batch,
        genome_fasta=str(args.genome_fasta),
        ispcr_bin=args.is_pcr_bin,
        tolerance=args.pcr_tolerance,
        ispcr_db=ispcr_db_path,
        ispcr_ooc=ispcr_ooc_path,
        tile_size=args.ispcr_tile_size,
    )
    elapsed = time.time() - t0

    pass_count = sum(1 for r in specificity_results.values() if r.specificity_pass)
    logger.info("isPcr done in %.0fs: %d unique_pass / %d total ok", elapsed, pass_count, len(ok_primers))

    # Write outputs
    args.output_dir.mkdir(parents=True, exist_ok=True)

    spec_records = build_specificity_records(
        primer_records, specificity_results, expected_coords,
    )

    spec_tsv_path = args.output_dir / "primers_specificity.tsv"
    unique_path = args.output_dir / "primers_unique.tsv"
    spec_summary_path = args.output_dir / "stage3_summary.txt"

    write_specificity_tsv(spec_records, spec_tsv_path)
    write_unique_primers(spec_records, unique_path)
    write_specificity_summary(spec_records, spec_summary_path)

    logger.info("Primers specificity TSV → %s", spec_tsv_path)
    logger.info("Unique primers TSV → %s", unique_path)
    logger.info("Stage 3 summary → %s", spec_summary_path)

    # Print summary
    from collections import Counter
    tid_unique: Counter[str] = Counter()
    for key, r in specificity_results.items():
        if r.specificity_pass:
            tid = key.rsplit("_rank", 1)[0]
            tid_unique[tid] += 1

    all_tids = sorted(set(pr.target_id for pr in ok_primers))
    print(f"\n{'='*60}")
    print(f"Stage 3 Results")
    print(f"{'='*60}")
    print(f"Total ok primers:    {len(ok_primers)}")
    print(f"Unique pass:         {pass_count}")
    print(f"Multi-hit:           {sum(1 for r in specificity_results.values() if r.insilico_status == 'multi_hit')}")
    print(f"No hit:              {sum(1 for r in specificity_results.values() if r.insilico_status == 'no_hit')}")
    print(f"Off-target:          {sum(1 for r in specificity_results.values() if r.insilico_status == 'unique_off_target')}")
    print(f"\n{'Target':25s} {'Unique':>6s}")
    print("-" * 35)
    for tid in all_tids:
        cnt = tid_unique.get(tid, 0)
        mark = "✓" if cnt > 0 else "✗"
        print(f"{mark} {tid:23s} {cnt:6d}")

    logger.info("Done. %d unique_pass / %d total ok → %s", pass_count, len(ok_primers), args.output_dir)


if __name__ == "__main__":
    main()
