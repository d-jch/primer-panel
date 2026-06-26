"""Ensembl REST API client — gene lookup, transcript selection, CDS extraction."""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass, field

import requests

from .config import PipelineConfig

logger = logging.getLogger(__name__)


@dataclass
class ExonCoord:
    """A single exon coordinate from Ensembl."""

    exon_id: str
    chrom: str          # e.g. "chr6"
    start: int          # 0-based
    end: int            # exclusive
    strand: int         # 1 or -1


@dataclass
class CdsExon:
    """A CDS segment within one exon.

    ``cds_exon_number`` is assigned in transcription order:
      - positive strand: genome ascending → 1, 2, 3 …
      - negative strand: genome descending → 1, 2, 3 …
    """

    cds_exon_number: int
    cds_exon_id: str     # typically the Ensembl exon ID
    chrom: str
    start: int           # 0-based, CDS start within this exon
    end: int             # exclusive, CDS end within this exon
    strand: int


@dataclass
class TranscriptInfo:
    """Metadata for the selected transcript."""

    transcript_id: str
    biotype: str
    is_mane_select: bool
    is_mane_plus_clinical: bool
    is_canonical: bool
    selection_reason: str
    exons: list[ExonCoord]
    cds_exons: list[CdsExon] = field(default_factory=list)
    translation_id: str | None = None


def _strip_version(ensembl_id: str) -> str:
    """ENST00000340047.8 → ENST00000340047"""
    return ensembl_id.split(".")[0]


def _chrom_from_region(region: str) -> str:
    """Convert Ensembl region '6:...' → 'chr6'."""
    chrom = region.split(":")[0]
    if not chrom.startswith("chr"):
        chrom = f"chr{chrom}"
    return chrom


class EnsemblClient:
    """Thin wrapper around Ensembl REST /lookup endpoints."""

    def __init__(self, cfg: PipelineConfig) -> None:
        self.cfg = cfg
        self.session = requests.Session()
        self.session.headers.update({"Content-Type": "application/json"})
        self._last_request_ts: float = 0.0

    # ------------------------------------------------------------------
    # Rate limiting (15 req/s for unauthenticated)
    # ------------------------------------------------------------------

    def _throttle(self) -> None:
        elapsed = time.time() - self._last_request_ts
        if elapsed < 0.075:  # ~13 req/s, safe margin
            time.sleep(0.075 - elapsed)
        self._last_request_ts = time.time()

    # ------------------------------------------------------------------
    # API calls
    # ------------------------------------------------------------------

    def _get_json(self, endpoint: str) -> dict:
        url = f"{self.cfg.ensembl_base}{endpoint}"
        self._throttle()
        resp = self.session.get(url, timeout=self.cfg.api_timeout)
        if resp.status_code == 429:
            retry_after = float(resp.headers.get("Retry-After", 1))
            logger.warning("Rate limited; sleeping %.1fs", retry_after)
            time.sleep(retry_after)
            resp = self.session.get(url, timeout=self.cfg.api_timeout)
        resp.raise_for_status()
        return resp.json()

    def lookup_gene(self, symbol: str) -> dict:
        """Full gene lookup with expanded transcripts."""
        return self._get_json(
            f"/lookup/symbol/homo_sapiens/{symbol}?expand=1;MANE=1"
        )

    # ------------------------------------------------------------------
    # Transcript selection
    # ------------------------------------------------------------------

    def select_transcript(self, symbol: str) -> tuple[TranscriptInfo, str]:
        """Return the best transcript for *symbol* and the reason string."""
        data = self.lookup_gene(symbol)

        if "error" in data:
            raise ValueError(f"Ensembl lookup failed for {symbol}: {data['error']}")

        transcripts: list[dict] = data.get("Transcript", [])
        if not transcripts:
            raise ValueError(f"No transcripts found for {symbol}")

        # Filter to protein-coding
        pc = [t for t in transcripts if t.get("biotype") == "protein_coding"]
        if not pc:
            # Fallback: use all transcripts (still pick longest)
            logger.warning("%s: no protein_coding transcripts; using all", symbol)
            pc = transcripts

        # Build candidates
        candidates: list[tuple[int, str, TranscriptInfo, int]] = []
        for t in pc:
            tid = t["id"]

            # MANE annotation may come as a list or dict depending on API version
            is_mane_select = False
            is_mane_plus = False
            mane_list = t.get("MANE", [])
            if isinstance(mane_list, list):
                for m in mane_list:
                    if isinstance(m, dict):
                        if m.get("type") == "MANE_Select":
                            is_mane_select = True
                        elif m.get("type") == "MANE_Plus_Clinical":
                            is_mane_plus = True
            elif isinstance(mane_list, dict):
                if mane_list.get("type") == "MANE_Select":
                    is_mane_select = True
                elif mane_list.get("type") == "MANE_Plus_Clinical":
                    is_mane_plus = True

            is_canonical = t.get("is_canonical", 0) == 1

            exons = self._extract_exons(t)
            cds_exons = self._extract_cds_exons(t, exons)
            total_bp = sum(e.end - e.start for e in exons)

            ti = TranscriptInfo(
                transcript_id=tid,
                biotype=t.get("biotype", "unknown"),
                is_mane_select=is_mane_select,
                is_mane_plus_clinical=is_mane_plus,
                is_canonical=is_canonical,
                selection_reason="",  # filled below
                exons=exons,
                cds_exons=cds_exons,
                translation_id=t.get("Translation", {}).get("id") if t.get("Translation") else None,
            )

            if is_mane_select:
                candidates.append((0, "MANE_Select", ti, total_bp))
            elif is_mane_plus:
                candidates.append((1, "MANE_Plus_Clinical", ti, total_bp))
            elif is_canonical:
                candidates.append((2, "canonical_protein_coding", ti, total_bp))
            else:
                candidates.append((3, "longest_protein_coding", ti, total_bp))

        if not candidates:
            raise ValueError(f"No suitable transcripts for {symbol}")

        # Sort by priority, then longest
        candidates.sort(key=lambda c: (c[0], -c[3]))
        _, reason, best, _ = candidates[0]
        best.selection_reason = reason

        logger.info(
            "%s: selected %s (%s, %d exons, %d CDS exons)",
            symbol, best.transcript_id, reason, len(best.exons), len(best.cds_exons),
        )
        return best, reason

    # ------------------------------------------------------------------
    # Exon extraction
    # ------------------------------------------------------------------

    def _extract_exons(self, transcript: dict) -> list[ExonCoord]:
        """Parse exon list from a transcript dict, converting to 0-based coords."""
        exons: list[ExonCoord] = []
        strand_val = transcript.get("strand", 1)

        for exon_data in transcript.get("Exon", []):
            start = exon_data["start"] - 1   # Ensembl is 1-based → 0-based
            end = exon_data["end"]
            chrom = _chrom_from_region(exon_data.get("region", transcript.get("seq_region_name", "unknown")))
            exons.append(ExonCoord(
                exon_id=exon_data.get("id", "unknown"),
                chrom=chrom,
                start=start,
                end=end,
                strand=strand_val,
            ))

        return exons

    # ------------------------------------------------------------------
    # CDS extraction
    # ------------------------------------------------------------------

    def _extract_cds_exons(self, transcript: dict, exons: list[ExonCoord]) -> list[CdsExon]:
        """Extract CDS segments from transcript Translation data.

        For each exon, compute the intersection with the CDS (translation) region.
        Returns CDS exon objects numbered in transcription order.
        Returns empty list if no Translation data (non-coding transcript).
        """
        translation = transcript.get("Translation")
        if not translation:
            return []

        # Translation start/end in Ensembl 1-based inclusive coords → 0-based
        cds_start = translation["start"] - 1
        cds_end = translation["end"]  # end is inclusive in Ensembl API, becomes exclusive in 0-based

        if cds_start >= cds_end:
            logger.warning("Invalid CDS range: start=%d >= end=%d", cds_start, cds_end)
            return []

        strand_val = transcript.get("strand", 1)
        chrom = exons[0].chrom if exons else "unknown"

        # Build CDS exon segments by intersecting each exon with CDS region
        raw_segments: list[tuple[int, int, str]] = []  # (start, end, exon_id)
        for exon in exons:
            # Intersection of [exon.start, exon.end) and [cds_start, cds_end)
            isect_start = max(exon.start, cds_start)
            isect_end = min(exon.end, cds_end)
            if isect_start < isect_end:
                raw_segments.append((isect_start, isect_end, exon.exon_id))

        if not raw_segments:
            logger.warning("Translation exists but no exon overlaps with CDS region")
            return []

        # Number in transcription order
        # Positive strand: genome ascending → 1, 2, 3 …
        # Negative strand: genome descending → 1, 2, 3 …
        if strand_val == 1:
            raw_segments.sort(key=lambda s: s[0])
        else:
            raw_segments.sort(key=lambda s: -s[0])

        cds_exons = []
        for i, (start, end, exon_id) in enumerate(raw_segments, 1):
            cds_exons.append(CdsExon(
                cds_exon_number=i,
                cds_exon_id=exon_id,
                chrom=chrom,
                start=start,
                end=end,
                strand=strand_val,
            ))

        return cds_exons
