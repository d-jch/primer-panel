"""CDS interval processing: merge adjacent CDS exons into required intervals."""

from __future__ import annotations

from dataclasses import dataclass, field

from .ensembl_client import CdsExon
from .utils import chrom_sort_key


@dataclass
class CdsRequiredInterval:
    """A required interval covering one or more adjacent CDS exons.

    These intervals represent the regions that *must* be covered by PCR targets.
    Adjacent CDS exons on the same chromosome are merged into a single
    CdsRequiredInterval.
    """

    chrom: str
    start: int          # 0-based, CDS exon start
    end: int            # exclusive, CDS exon end
    strand: int
    cds_exons: list[CdsExon] = field(default_factory=list)

    @property
    def length(self) -> int:
        return self.end - self.start


def build_required_intervals(
    cds_exons: list[CdsExon],
) -> list[CdsRequiredInterval]:
    """Build required intervals from CDS exons by merging adjacent ones.

    Adjacent/overlapping CDS exons on the same chromosome are merged.
    No buffer is added — the raw CDS exon coordinates are used directly.

    Returns intervals sorted by (chrom, start).
    """
    if not cds_exons:
        return []

    # Wrap each CDS exon as a CdsRequiredInterval
    intervals: list[CdsRequiredInterval] = []
    for cex in cds_exons:
        intervals.append(CdsRequiredInterval(
            chrom=cex.chrom,
            start=cex.start,
            end=cex.end,
            strand=cex.strand,
            cds_exons=[cex],
        ))

    # Sort by (chrom, start)
    intervals.sort(key=lambda iv: (chrom_sort_key(iv.chrom), iv.start, iv.end))

    # Merge overlapping/adjacent on same chrom
    merged: list[CdsRequiredInterval] = [intervals[0]]
    for iv in intervals[1:]:
        prev = merged[-1]
        if iv.chrom == prev.chrom and iv.start <= prev.end:
            # overlap or adjacent — merge
            prev.end = max(prev.end, iv.end)
            prev.cds_exons.extend(iv.cds_exons)
        else:
            merged.append(iv)

    return merged
