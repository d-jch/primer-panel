"""Exon interval processing: buffer, sort, merge."""

from __future__ import annotations

from dataclasses import dataclass

from .ensembl_client import ExonCoord
from .utils import chrom_sort_key


@dataclass
class GenomicInterval:
    """A merged genomic interval covering one or more buffered exons."""

    chrom: str
    start: int          # 0-based
    end: int            # exclusive
    strand: int         # 1 or -1
    exon_ids: list[str] # constituent exon IDs

    @property
    def length(self) -> int:
        return self.end - self.start


def buffer_and_merge(
    exons: list[ExonCoord],
    splice_buffer: int = 30,
) -> list[GenomicInterval]:
    """Add splice buffers to each exon, then merge overlapping intervals.

    Input exons may be unsorted; output is sorted by (chrom, start).
    """
    if not exons:
        return []

    # Apply buffer
    buffered = [
        GenomicInterval(
            chrom=e.chrom,
            start=max(0, e.start - splice_buffer),
            end=e.end + splice_buffer,
            strand=e.strand,
            exon_ids=[e.exon_id],
        )
        for e in exons
    ]

    # Sort
    buffered.sort(key=lambda iv: (chrom_sort_key(iv.chrom), iv.start, iv.end))

    # Merge overlapping on same chrom
    merged: list[GenomicInterval] = [buffered[0]]
    for iv in buffered[1:]:
        prev = merged[-1]
        if iv.chrom == prev.chrom and iv.start <= prev.end:
            # overlap or adjacent — merge
            prev.end = max(prev.end, iv.end)
            prev.exon_ids.extend(iv.exon_ids)
        else:
            merged.append(iv)

    return merged
