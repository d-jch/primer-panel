"""Small helpers shared across modules."""

from __future__ import annotations

import re


_CHROM_ORDER: dict[str, int] = {
    f"chr{i}": i for i in range(1, 23)
}
_CHROM_ORDER.update({"chrX": 23, "chrY": 24, "chrM": 25})


def chrom_sort_key(chrom: str) -> tuple[int, str]:
    """Return a sortable key for chromosome names (chr1..chr22, chrX, chrY, chrM, then alpha)."""
    if chrom in _CHROM_ORDER:
        return (0, "", _CHROM_ORDER[chrom])
    # fallback: try to extract numeric part
    m = re.match(r"chr(\d+)", chrom)
    if m:
        return (0, "", int(m.group(1)))
    return (1, chrom, 0)
