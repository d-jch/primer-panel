"""Stage 3: In-silico PCR specificity check using UCSC isPcr.

Uses the UCSC Kent tools isPcr binary for genome-wide PCR simulation.
  1. Write a query file with all primer pairs.
  2. Run isPcr once on the genome FASTA.
  3. Parse BED output to get PCR products.
  4. Classify each primer pair as unique_pass / multi_hit / no_hit /
     unique_off_target / pcr_error.

Stage 3 does NOT filter hits by product size.  All isPcr hits are retained.
Product size is reported as information only, not used for pass/fail.
"""

from __future__ import annotations

import csv
import logging
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

logger = logging.getLogger(__name__)


# ──────────────────────────────────────────────────────────────────────
# Data classes
# ──────────────────────────────────────────────────────────────────────

@dataclass
class PcrHit:
    """One in-silico PCR product."""
    chrom: str
    start: int          # 0-based, inclusive
    end: int            # 0-based, exclusive
    size: int
    strand: str         # "+" or "-"


@dataclass
class SpecificityResult:
    """Specificity result for one primer pair."""
    insilico_status: str        # unique_pass / multi_hit / no_hit / unique_off_target / pcr_error
    insilico_hit_count: int
    insilico_hits: str          # semicolon-separated "chrom:start-end(size)" strings
    insilico_best_chrom: str
    insilico_best_start: int
    insilico_best_end: int
    insilico_best_size: int
    specificity_pass: bool
    specificity_explain: str    # human-readable explanation of the classification


# ──────────────────────────────────────────────────────────────────────
# isPcr availability check
# ──────────────────────────────────────────────────────────────────────

def check_ispcr_available(bin_path: str = "isPcr") -> bool:
    """Check whether isPcr is callable."""
    return shutil.which(bin_path) is not None


# ──────────────────────────────────────────────────────────────────────
# Batch isPcr runner
# ──────────────────────────────────────────────────────────────────────

def run_ispcr_batch(
    primers: list[dict],  # list of {name, fwd, rev}
    genome_fasta: str,
    ispcr_bin: str = "isPcr",
    min_perfect: int = 15,
    min_good: int = 15,
) -> dict[str, list[PcrHit]]:
    """Run isPcr on a batch of primer pairs and return ALL hits.

    No product-size filtering — all isPcr hits are retained.

    Each primer dict should have: name (str), fwd (str), rev (str).
    Returns {primer_name: [PcrHit, ...]}.
    """
    if not primers:
        return {}

    # Write query file
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".txt", delete=False, prefix="ispcr_query_"
    ) as qf:
        query_path = qf.name
        for p in primers:
            qf.write(f"{p['name']}\t{p['fwd']}\t{p['rev']}\n")

    # Run isPcr
    try:
        cmd = [
            ispcr_bin,
            genome_fasta,
            query_path,
            "stdout",
            "-out=bed",
        ]
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=600,  # 10 min timeout for large genome
        )
    finally:
        Path(query_path).unlink(missing_ok=True)

    if proc.returncode != 0:
        logger.error("isPcr failed (exit %d): %s", proc.returncode, proc.stderr[:500])
        raise RuntimeError(f"isPcr failed: {proc.stderr[:200]}")

    # Parse BED output — retain ALL hits, no size filter
    hits_by_name: dict[str, list[PcrHit]] = {p["name"]: [] for p in primers}

    for line in proc.stdout.strip().splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 6:
            continue

        chrom = parts[0]
        start = int(parts[1])   # BED is 0-based start
        end = int(parts[2])     # BED is 0-based exclusive end
        name = parts[3]
        strand = parts[5]
        size = end - start

        if name in hits_by_name:
            hits_by_name[name].append(PcrHit(
                chrom=chrom, start=start, end=end, size=size, strand=strand,
            ))

    return hits_by_name


# ──────────────────────────────────────────────────────────────────────
# Main API
# ──────────────────────────────────────────────────────────────────────

def check_specificity_batch(
    primer_pairs: list[dict],  # list of {name, fwd, rev, expected_chrom, expected_start, expected_end}
    genome_fasta: str,
    ispcr_bin: str = "isPcr",
    tolerance: int = 10,
    min_perfect: int = 15,
    min_good: int = 15,
) -> dict[str, SpecificityResult]:
    """Check genome-wide PCR specificity for a batch of primer pairs.

    No product-size filtering.  Classification is based on:
    - hit count
    - best hit chrom/start/end vs expected chrom/start/end (within tolerance)

    Returns {primer_name: SpecificityResult}.
    """
    # Run isPcr batch
    try:
        hits_by_name = run_ispcr_batch(
            primers=[{"name": p["name"], "fwd": p["fwd"], "rev": p["rev"]} for p in primer_pairs],
            genome_fasta=genome_fasta,
            ispcr_bin=ispcr_bin,
            min_perfect=min_perfect,
            min_good=min_good,
        )
    except Exception as exc:
        logger.error("isPcr batch failed: %s", exc)
        # Return pcr_error for all
        return {
            p["name"]: SpecificityResult(
                insilico_status="pcr_error",
                insilico_hit_count=0,
                insilico_hits="",
                insilico_best_chrom="",
                insilico_best_start=0,
                insilico_best_end=0,
                insilico_best_size=0,
                specificity_pass=False,
                specificity_explain=f"isPcr error: {exc}",
            )
            for p in primer_pairs
        }

    # Classify each primer pair
    results: dict[str, SpecificityResult] = {}

    for p in primer_pairs:
        name = p["name"]
        hits = hits_by_name.get(name, [])
        expected_chrom = p["expected_chrom"]
        expected_start = p["expected_start"]
        expected_end = p["expected_end"]

        hit_count = len(hits)
        hits_str = ";".join(f"{h.chrom}:{h.start}-{h.end}({h.size})" for h in hits)

        if hit_count == 0:
            results[name] = SpecificityResult(
                insilico_status="no_hit",
                insilico_hit_count=0,
                insilico_hits="",
                insilico_best_chrom="",
                insilico_best_start=0,
                insilico_best_end=0,
                insilico_best_size=0,
                specificity_pass=False,
                specificity_explain="no genome-wide PCR product detected",
            )
            continue

        # Find best hit (closest to expected coordinates)
        def _hit_distance(h: PcrHit) -> int:
            return abs(h.start - expected_start) + abs(h.end - expected_end)

        best = min(hits, key=_hit_distance)

        # Check if best hit matches expected location
        chrom_match = best.chrom == expected_chrom
        start_delta = abs(best.start - expected_start)
        end_delta = abs(best.end - expected_end)
        coord_match = start_delta <= tolerance and end_delta <= tolerance
        expected_match = chrom_match and coord_match

        if hit_count == 1 and expected_match:
            results[name] = SpecificityResult(
                insilico_status="unique_pass",
                insilico_hit_count=1,
                insilico_hits=hits_str,
                insilico_best_chrom=best.chrom,
                insilico_best_start=best.start,
                insilico_best_end=best.end,
                insilico_best_size=best.size,
                specificity_pass=True,
                specificity_explain="single hit matches expected product",
            )
        elif hit_count == 1 and not expected_match:
            # Explain what didn't match
            reasons: list[str] = []
            if not chrom_match:
                reasons.append(f"chrom mismatch ({best.chrom} vs {expected_chrom})")
            if start_delta > tolerance:
                reasons.append(f"start delta={start_delta}bp (tolerance={tolerance})")
            if end_delta > tolerance:
                reasons.append(f"end delta={end_delta}bp (tolerance={tolerance})")
            explain = "single hit does not match expected product: " + "; ".join(reasons)

            results[name] = SpecificityResult(
                insilico_status="unique_off_target",
                insilico_hit_count=1,
                insilico_hits=hits_str,
                insilico_best_chrom=best.chrom,
                insilico_best_start=best.start,
                insilico_best_end=best.end,
                insilico_best_size=best.size,
                specificity_pass=False,
                specificity_explain=explain,
            )
        else:
            # hit_count >= 2 — multi_hit regardless of whether one matches
            results[name] = SpecificityResult(
                insilico_status="multi_hit",
                insilico_hit_count=hit_count,
                insilico_hits=hits_str,
                insilico_best_chrom=best.chrom,
                insilico_best_start=best.start,
                insilico_best_end=best.end,
                insilico_best_size=best.size,
                specificity_pass=False,
                specificity_explain=f"multiple genome-wide PCR products detected ({hit_count} hits)",
            )

    return results
