"""Common dbSNP variant annotation for primer binding sites.

Provides efficient overlap queries between primer coordinates and
common dbSNP variants, supporting both plain BED and bigBed (.bb) files.

BED format (0-based, half-open):
    chrom   start   end   rsid

bigBed files are queried directly via ``bigBedToBed`` region queries —
no intermediate conversion or full-text dump is needed.
"""

from __future__ import annotations

import csv
import logging
import shutil
import subprocess
from bisect import bisect_left
from dataclasses import dataclass, field
from pathlib import Path
from typing import Protocol

logger = logging.getLogger(__name__)


# ── data structures ──────────────────────────────────────────────────────────


@dataclass
class SnpEntry:
    """A single SNP entry."""

    chrom: str
    start: int  # 0-based
    end: int    # half-open
    rsid: str


class SnpDatabase(Protocol):
    """Interface for SNP overlap queries.

    Implementations:
        BedSnpDatabase  — in-memory index from a plain BED file
        BigBedSnpDatabase — region queries against a .bb file via bigBedToBed
    """

    def query_region(self, chrom: str, start: int, end: int) -> list[SnpEntry]:
        """Return SNPs overlapping [start, end)."""
        ...


@dataclass
class ChromIndex:
    """Pre-computed index for one chromosome."""

    entries: list[SnpEntry]
    starts: list[int]  # sorted start positions
    max_end_prefix: list[int]  # max(entries[0..i].end) for early termination


@dataclass
class BedSnpDatabase:
    """In-memory SNP index loaded from a plain BED file."""

    chroms: dict[str, ChromIndex] = field(default_factory=dict)
    _total_snps: int = 0

    def query_region(self, chrom: str, start: int, end: int) -> list[SnpEntry]:
        """Find SNPs overlapping [start, end) using in-memory B-tree."""
        if chrom not in self.chroms:
            return []

        idx = self.chroms[chrom]
        if not idx.entries:
            return []

        upper = bisect_left(idx.starts, end)
        result = []
        for i in range(upper - 1, -1, -1):
            entry = idx.entries[i]
            if entry.end > start:
                result.append(entry)
            if idx.max_end_prefix[i] <= start:
                break

        return result


@dataclass
class BigBedSnpDatabase:
    """bigBed-backed SNP database queried via ``bigBedToBed``.

    Each ``query_region`` call spawns a ``bigBedToBed`` subprocess with
    ``-chrom -start -end`` flags.  bigBed's internal B+tree index makes
    each region query fast — only the relevant data blocks are read.
    """

    path: Path

    def query_region(self, chrom: str, start: int, end: int) -> list[SnpEntry]:
        """Query bigBed for SNPs in [start, end)."""
        cmd = [
            "bigBedToBed",
            "-chrom=" + chrom,
            "-start=" + str(start),
            "-end=" + str(end),
            str(self.path),
            "stdout",
        ]
        try:
            proc = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=30,
            )
        except FileNotFoundError:
            raise RuntimeError(
                "bigBedToBed is required to query .bb files.\n"
                "  Install via micromamba/conda:\n"
                "    micromamba install -c bioconda ucsc-bigbedtobed\n"
                "  Or convert the .bb to .bed first:\n"
                "    bigBedToBed input.bb output.bed"
            ) from None

        if proc.returncode != 0:
            stderr = proc.stderr.strip()
            raise RuntimeError(
                f"bigBedToBed failed (exit {proc.returncode}): {stderr}"
            )

        results: list[SnpEntry] = []
        for line in proc.stdout.splitlines():
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                continue
            try:
                results.append(SnpEntry(
                    chrom=fields[0],
                    start=int(fields[1]),
                    end=int(fields[2]),
                    rsid=fields[3],
                ))
            except (ValueError, IndexError):
                continue

        return results


# ── factory ──────────────────────────────────────────────────────────────────


def load_dbsnp_db(path: Path) -> SnpDatabase:
    """Load a SNP database from a BED or bigBed file.

    - ``.bed`` files are loaded into an in-memory index.
    - ``.bb``  files are opened for on-demand region queries via
      ``bigBedToBed`` (requires ``ucsc-bigbedtobed``).
    """
    if path.suffix == ".bb":
        _check_bigbedtobed()
        logger.info("Opening bigBed SNP database: %s", path)
        return BigBedSnpDatabase(path=path)

    logger.info("Loading BED SNP database from %s …", path)
    return _load_bed(path)


# ── BED parser (internal) ────────────────────────────────────────────────────


def _load_bed(bed_path: Path) -> BedSnpDatabase:
    """Load a plain BED file and build per-chromosome index."""
    chrom_entries: dict[str, list[SnpEntry]] = {}

    with open(bed_path, encoding="utf-8") as f:
        for row in csv.reader(f, delimiter="\t"):
            if len(row) < 4:
                continue
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

    db = BedSnpDatabase()
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
    db._total_snps = sum(len(v.entries) for v in db.chroms.values())

    logger.info(
        "Loaded %d SNPs across %d chromosomes from %s",
        db._total_snps,
        len(db.chroms),
        bed_path,
    )
    return db


# ── helpers ──────────────────────────────────────────────────────────────────


_BIGBEDTOBED_HINT = (
    "bigBedToBed is required to query .bb SNP databases.\n"
    "  Install via micromamba/conda:\n"
    "    micromamba install -c bioconda ucsc-bigbedtobed"
)


def check_bigbedtobed() -> str | None:
    """Return install hint if bigBedToBed is missing, else None."""
    if shutil.which("bigBedToBed") is None:
        return _BIGBEDTOBED_HINT
    return None


def _check_bigbedtobed() -> None:
    """Raise RuntimeError if bigBedToBed is not available."""
    hint = check_bigbedtobed()
    if hint:
        raise RuntimeError(hint)


# ── overlap query ────────────────────────────────────────────────────────────


def _find_overlapping_snps(
    db: SnpDatabase,
    chrom: str,
    start: int,
    end: int,
) -> list[SnpEntry]:
    """Find SNPs overlapping [start, end) interval.

    Delegates to ``db.query_region`` so the same code works for both
    in-memory (BED) and on-demand (bigBed) backends.
    """
    return db.query_region(chrom, start, end)


# ── primer annotation ────────────────────────────────────────────────────────


def annotate_primer_snps(
    db: SnpDatabase,
    chrom: str,
    primer_start: int,
    primer_len: int,
    is_reverse: bool = False,
) -> tuple[str, int, int, str]:
    """Annotate a single primer with overlapping SNPs.

    Args:
        db: SNP database (BED or bigBed backend)
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
    overall_risk = max(left_risk, right_risk, key=lambda r: risk_order.get(r, 0))

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


# ── deprecated alias ─────────────────────────────────────────────────────────


def load_dbsnp_bed(bed_path: Path) -> SnpDatabase:
    """Deprecated: use ``load_dbsnp_db`` instead.

    Kept for backward compatibility with existing callers.
    """
    logger.warning(
        "load_dbsnp_bed is deprecated — use load_dbsnp_db which also "
        "supports .bb (bigBed) files."
    )
    return load_dbsnp_db(bed_path)
