"""Common dbSNP variant annotation for primer binding sites.

Provides efficient overlap queries between primer coordinates and
common dbSNP variants loaded from a BED file.

BED format (0-based, half-open):
    chrom   start   end   rsid
"""

from __future__ import annotations

import csv
import logging
from bisect import bisect_left
from dataclasses import dataclass, field
from pathlib import Path

logger = logging.getLogger(__name__)


@dataclass
class SnpEntry:
    """A single SNP entry from BED file."""

    chrom: str
    start: int  # 0-based
    end: int    # half-open
    rsid: str


@dataclass
class ChromIndex:
    """Pre-computed index for one chromosome."""

    entries: list[SnpEntry]
    starts: list[int]  # sorted start positions
    max_end_prefix: list[int]  # max(entries[0..i].end) for early termination


@dataclass
class SnpDatabase:
    """Indexed SNP database for efficient overlap queries."""

    chroms: dict[str, ChromIndex] = field(default_factory=dict)


def load_dbsnp_bed(bed_path: Path) -> SnpDatabase:
    """Load common dbSNP BED file and build per-chromosome index.

    BED format: chrom start end name (0-based, half-open)
    Lines with fewer than 4 columns are skipped.
    """
    chrom_entries: dict[str, list[SnpEntry]] = {}

    with open(bed_path) as f:
        for row in csv.reader(f, delimiter="\t"):
            if len(row) < 4:
                continue
            # Skip comment/header lines
            if row[0].startswith("#"):
                continue
            chrom, start_str, end_str, rsid = row[0], row[1], row[2], row[3]
            try:
                start, end = int(start_str), int(end_str)
            except ValueError:
                continue
            if chrom not in chrom_entries:
                chrom_entries[chrom] = []
            chrom_entries[chrom].append(SnpEntry(chrom, start, end, rsid))

    # Build sorted index per chromosome
    db = SnpDatabase()
    for chrom, entries in chrom_entries.items():
        entries.sort(key=lambda e: e.start)
        starts = [e.start for e in entries]
        max_end_prefix: list[int] = []
        running_max = 0
        for e in entries:
            running_max = max(running_max, e.end)
            max_end_prefix.append(running_max)
        db.chroms[chrom] = ChromIndex(
            entries=entries, starts=starts, max_end_prefix=max_end_prefix,
        )

    total_snps = sum(len(v.entries) for v in db.chroms.values())
    logger.info(
        "Loaded %d SNPs across %d chromosomes from %s",
        total_snps,
        len(db.chroms),
        bed_path,
    )
    return db


def _find_overlapping_snps(
    db: SnpDatabase,
    chrom: str,
    start: int,
    end: int,
) -> list[SnpEntry]:
    """Find SNPs overlapping [start, end) interval.

    Uses binary search on sorted starts and a prefix-maximum of ends
    for efficient early termination.  Correctly handles intervals of any
    length, not just short SNPs.
    """
    if chrom not in db.chroms:
        return []

    idx = db.chroms[chrom]
    if not idx.entries:
        return []

    # Upper bound: entries with start >= end cannot overlap
    upper = bisect_left(idx.starts, end)

    result = []
    # Scan backwards from upper-1; use max_end_prefix for early termination
    for i in range(upper - 1, -1, -1):
        entry = idx.entries[i]
        if entry.end > start:
            result.append(entry)
        # If the maximum end among entries[0..i] doesn't reach start,
        # no earlier entry can overlap either.
        if idx.max_end_prefix[i] <= start:
            break

    return result


def annotate_primer_snps(
    db: SnpDatabase,
    chrom: str,
    primer_start: int,
    primer_len: int,
    is_reverse: bool = False,
) -> tuple[str, int, int, str]:
    """Annotate a single primer with overlapping SNPs.

    Args:
        db: SNP database
        chrom: chromosome name
        primer_start: primer start position (template-relative)
        primer_len: primer length
        is_reverse: True for reverse primer, False for forward

    Returns:
        (risk_level, total_count, three_prime_count, hits_str)
        - risk_level: "high" if 3' end overlap, "medium" if other overlap, "none"
        - total_count: total overlapping SNPs
        - three_prime_count: SNPs overlapping 3' end (last 5bp)
        - hits_str: semicolon-separated hit descriptions
    """
    if is_reverse:
        # Reverse primer: primer_start is the rightmost (high-coordinate) base
        # on the template (Primer3 PRIMER_RIGHT product-end coordinate).
        # Template region: [primer_start - primer_len + 1, primer_start + 1)
        # Synthesis goes right-to-left, so the 3' end is the LOW-coordinate end.
        # High-risk 3' region: first 5 bases of the primer on the template.
        primer_end = primer_start + 1
        primer_start_genomic = primer_start - primer_len + 1
        three_prime_start = primer_start_genomic
        three_prime_end = primer_start_genomic + 5
    else:
        # Forward primer: 3' end is at primer_start + primer_len - 1
        # Primer extends right: [primer_start, primer_start + primer_len)
        primer_start_genomic = primer_start
        primer_end = primer_start + primer_len
        three_prime_start = primer_start + primer_len - 5
        three_prime_end = primer_start + primer_len

    # Find all overlapping SNPs
    overlapping = _find_overlapping_snps(db, chrom, primer_start_genomic, primer_end)

    if not overlapping:
        return "none", 0, 0, ""

    # Classify hits by 3' overlap
    three_prime_hits = []
    other_hits = []

    for snp in overlapping:
        if snp.end > three_prime_start and snp.start < three_prime_end:
            three_prime_hits.append(snp)
        else:
            other_hits.append(snp)

    # Determine risk level
    if three_prime_hits:
        risk = "high"
    elif other_hits:
        risk = "medium"
    else:
        risk = "none"

    # Build hits string
    hits = []
    for snp in three_prime_hits:
        hits.append(f"{snp.rsid}:3p")
    for snp in other_hits:
        hits.append(f"{snp.rsid}:other")

    return risk, len(overlapping), len(three_prime_hits), ";".join(hits)


def annotate_primer_pair(
    db: SnpDatabase,
    chrom: str,
    left_start: int,
    left_len: int,
    right_start: int,
    right_len: int,
) -> dict[str, str | int]:
    """Annotate a primer pair with SNP overlap information.

    Returns dict with keys matching SpecificityRecord fields:
        common_snp_risk, left_primer_common_snp_count,
        right_primer_common_snp_count, left_primer_3p_common_snp_count,
        right_primer_3p_common_snp_count, common_snp_hits
    """
    left_risk, left_total, left_3p, left_hits = annotate_primer_snps(
        db, chrom, left_start, left_len, is_reverse=False
    )
    right_risk, right_total, right_3p, right_hits = annotate_primer_snps(
        db, chrom, right_start, right_len, is_reverse=True
    )

    # Overall risk is highest of left/right
    risk_order = {"none": 0, "medium": 1, "high": 2}
    overall_risk = "none"
    if risk_order.get(left_risk, 0) > risk_order.get(overall_risk, 0):
        overall_risk = left_risk
    if risk_order.get(right_risk, 0) > risk_order.get(overall_risk, 0):
        overall_risk = right_risk

    # Combine hits
    all_hits = []
    if left_hits:
        all_hits.append(f"left:{left_hits}")
    if right_hits:
        all_hits.append(f"right:{right_hits}")

    return {
        "common_snp_risk": overall_risk,
        "left_primer_common_snp_count": left_total,
        "right_primer_common_snp_count": right_total,
        "left_primer_3p_common_snp_count": left_3p,
        "right_primer_3p_common_snp_count": right_3p,
        "common_snp_hits": "|".join(all_hits),
    }
