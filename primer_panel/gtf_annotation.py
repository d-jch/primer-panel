"""Local Ensembl GTF annotation client.

Parses an Ensembl GTF file (plain or gzip) and provides gene/transcript
lookup without network access.  Coordinate convention: GTF uses 1-based
inclusive; this module converts to 0-based half-open to match
``EnsemblClient`` output structures.
"""

from __future__ import annotations

import gzip
import logging
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from .ensembl_client import CdsExon, ExonCoord, TranscriptInfo

logger = logging.getLogger(__name__)

# ── GTF attribute parser ───────────────────────────────────────────────────

_ATTR_RE = re.compile(r'(\w+)\s+"([^"]*)"')


def _parse_attributes(attr_str: str) -> dict[str, str]:
    """Parse GTF attribute column into a dict.

    Handles both comma-separated and repeated ``tag`` attributes:

    - ``tag "basic"; tag "MANE_Select";``  (Ensembl real format)
    - ``tag "MANE_Select,Ensembl_canonical";``  (single-value comma format)

    The ``tag`` key is returned as a comma-joined string of all values found.
    All other keys use last-wins semantics (they are unique in valid GTF).
    """
    result: dict[str, str] = {}
    for key, value in _ATTR_RE.findall(attr_str):
        if key == "tag":
            # Accumulate: merge with existing value if present
            existing = result.get("tag", "")
            result["tag"] = f"{existing},{value}" if existing else value
        else:
            result[key] = value
    return result


# ── internal data model (GTF-native, before coordinate conversion) ─────────


@dataclass
class _GtfExon:
    """Raw exon record from GTF (0-based half-open after conversion)."""

    exon_id: str
    chrom: str
    start: int  # 0-based
    end: int  # exclusive
    strand: int
    exon_number: int | None = None


@dataclass
class _GtfTranscript:
    """Aggregated data for one transcript from the GTF."""

    transcript_id: str
    gene_id: str
    gene_name: str
    biotype: str
    chrom: str
    strand: int
    tags: set[str] = field(default_factory=set)
    exons: list[_GtfExon] = field(default_factory=list)
    cds_starts: list[int] = field(default_factory=list)
    cds_ends: list[int] = field(default_factory=list)
    # CDS records keyed by (start, end) for exon matching
    cds_records: list[tuple[int, int]] = field(default_factory=list)


@dataclass
class _GtfGene:
    """Aggregated data for one gene."""

    gene_id: str
    gene_name: str
    chrom: str
    strand: int
    start: int  # min exon start (0-based)
    end: int  # max exon end (exclusive)
    transcripts: dict[str, _GtfTranscript] = field(default_factory=dict)


# ── main parser ────────────────────────────────────────────────────────────


class GtfAnnotationClient:
    """Parse an Ensembl GTF and provide ``EnsemblClient``-compatible lookups."""

    def __init__(self, gtf_path: str | Path) -> None:
        self._path = Path(gtf_path)
        self._genes_by_name: dict[str, _GtfGene] = {}
        self._genes_by_id: dict[str, _GtfGene] = {}
        self._parse()

    # ── public API (matches EnsemblClient shape) ───────────────────────

    def lookup_gene(self, symbol: str) -> dict:
        """Return a dict compatible with ``EnsemblClient.lookup_gene``.

        Provides ``start``, ``end``, and ``Transcript`` keys used by the
        pipeline's Stage 1 target planning.
        """
        gene = self._resolve_gene(symbol)
        transcripts_list = []
        for t in gene.transcripts.values():
            transcripts_list.append({
                "id": t.transcript_id,
                "biotype": t.biotype,
                "is_canonical": 1 if "Ensembl_canonical" in t.tags else 0,
                "strand": t.strand,
                "Exon": [
                    {
                        "id": e.exon_id,
                        "start": e.start + 1,  # back to 1-based for compat
                        "end": e.end,
                        "region": e.chrom,
                    }
                    for e in t.exons
                ],
                "Translation": self._build_translation_dict(t),
                "MANE": self._build_mane_list(t),
            })
        return {
            "start": gene.start + 1,  # 1-based for pipeline compat
            "end": gene.end,
            "Transcript": transcripts_list,
        }

    def select_transcript(self, symbol: str) -> tuple[TranscriptInfo, str]:
        """Select the best transcript for *symbol*.

        Priority (mirrors ``EnsemblClient.select_transcript``):
          1. MANE_Select (protein_coding)
          2. MANE_Plus_Clinical (protein_coding)
          3. Ensembl_canonical (protein_coding)
          4. Longest CDS (protein_coding)
          5. Fallback: longest transcript (any biotype, only if no protein_coding)

        Returns ``(TranscriptInfo, reason_string)``.
        """
        gene = self._resolve_gene(symbol)
        if not gene.transcripts:
            raise ValueError(f"No transcripts found for {symbol} in GTF")

        # Filter to protein_coding transcripts first
        pc_transcripts = [
            t for t in gene.transcripts.values() if t.biotype == "protein_coding"
        ]
        # Fallback: use all transcripts if no protein_coding found
        if not pc_transcripts:
            logger.warning("%s: no protein_coding transcripts in GTF; using all", symbol)
            pc_transcripts = list(gene.transcripts.values())

        candidates: list[tuple[int, str, TranscriptInfo, int]] = []

        for t in pc_transcripts:
            ti = self._build_transcript_info(t)
            if not ti.cds_exons:
                continue  # skip transcripts without CDS

            cds_total = sum(e.end - e.start for e in ti.cds_exons)

            if "MANE_Select" in t.tags:
                candidates.append((0, "MANE_Select", ti, cds_total))
            elif "MANE_Plus_Clinical" in t.tags:
                candidates.append((1, "MANE_Plus_Clinical", ti, cds_total))
            elif ti.is_canonical:
                candidates.append((2, "canonical_protein_coding", ti, cds_total))
            else:
                candidates.append((3, "longest_protein_coding", ti, cds_total))

        if not candidates:
            # Fallback: try all transcripts (including non-protein_coding with exons)
            for t in gene.transcripts.values():
                ti = self._build_transcript_info(t)
                exon_total = sum(e.end - e.start for e in ti.exons)
                candidates.append((4, "longest_transcript", ti, exon_total))

        if not candidates:
            raise ValueError(f"No suitable transcripts for {symbol} in GTF")

        candidates.sort(key=lambda c: (c[0], -c[3]))
        _, reason, best, _ = candidates[0]
        best.selection_reason = reason

        logger.info(
            "%s: selected %s from GTF (%s, %d exons, %d CDS exons)",
            symbol,
            best.transcript_id,
            reason,
            len(best.exons),
            len(best.cds_exons),
        )
        return best, reason

    # ── internal: parse GTF ────────────────────────────────────────────

    def _parse(self) -> None:
        """Read the GTF and build in-memory indexes."""
        opener = gzip.open if self._path.suffix == ".gz" else open
        logger.info("Parsing GTF: %s", self._path)

        with opener(self._path, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                seqname, source, feature, start, end, score, strand, frame, attrs = parts
                if feature not in ("gene", "transcript", "exon", "CDS"):
                    continue

                attr = _parse_attributes(attrs)
                chrom = seqname if seqname.startswith("chr") else f"chr{seqname}"
                strand_int = 1 if strand == "+" else -1

                # GTF: 1-based inclusive → 0-based half-open
                gtf_start = int(start) - 1
                gtf_end = int(end)

                if feature == "gene":
                    self._handle_gene(attr, chrom, strand_int, gtf_start, gtf_end)
                elif feature == "transcript":
                    self._handle_transcript(attr, chrom, strand_int)
                elif feature == "exon":
                    self._handle_exon(attr, chrom, strand_int, gtf_start, gtf_end)
                elif feature == "CDS":
                    self._handle_cds(attr, gtf_start, gtf_end)

        # Finalize gene bounds
        for gene in self._genes_by_id.values():
            if gene.transcripts:
                all_starts = []
                all_ends = []
                for t in gene.transcripts.values():
                    for e in t.exons:
                        all_starts.append(e.start)
                        all_ends.append(e.end)
                if all_starts:
                    gene.start = min(all_starts)
                    gene.end = max(all_ends)

        logger.info(
            "GTF parsed: %d genes (%d by name, %d by id)",
            len(set(id(g) for g in self._genes_by_id.values())),
            len(self._genes_by_name),
            len(self._genes_by_id),
        )

    def _handle_gene(
        self, attr: dict, chrom: str, strand: int, start: int, end: int
    ) -> None:
        gene_id = attr.get("gene_id", "")
        gene_name = attr.get("gene_name", gene_id)
        if not gene_id:
            return
        gene = _GtfGene(
            gene_id=gene_id,
            gene_name=gene_name,
            chrom=chrom,
            strand=strand,
            start=start,
            end=end,
        )
        self._genes_by_id[gene_id] = gene
        # Index by name (case-insensitive for lookup convenience)
        key = gene_name.upper()
        if key not in self._genes_by_name:
            self._genes_by_name[key] = gene

    def _handle_transcript(
        self, attr: dict, chrom: str, strand: int
    ) -> None:
        transcript_id = attr.get("transcript_id", "")
        gene_id = attr.get("gene_id", "")
        if not transcript_id or not gene_id:
            return
        gene = self._genes_by_id.get(gene_id)
        if gene is None:
            return
        biotype = attr.get("transcript_biotype", attr.get("gene_biotype", "unknown"))
        tags_raw = attr.get("tag", "")
        # Split on comma; handles both repeated-attribute joins and inline commas
        tags = {t.strip() for t in tags_raw.split(",") if t.strip()} if tags_raw else set()
        t = _GtfTranscript(
            transcript_id=transcript_id,
            gene_id=gene_id,
            gene_name=gene.gene_name,
            biotype=biotype,
            chrom=chrom,
            strand=strand,
            tags=tags,
        )
        gene.transcripts[transcript_id] = t

    def _handle_exon(
        self, attr: dict, chrom: str, strand: int, start: int, end: int
    ) -> None:
        transcript_id = attr.get("transcript_id", "")
        gene_id = attr.get("gene_id", "")
        if not transcript_id:
            return
        gene = self._genes_by_id.get(gene_id)
        if gene is None:
            return
        t = gene.transcripts.get(transcript_id)
        if t is None:
            return
        exon_id = attr.get("exon_id", f"{transcript_id}_exon")
        exon_number = None
        if "exon_number" in attr:
            try:
                exon_number = int(attr["exon_number"])
            except ValueError:
                pass
        t.exons.append(
            _GtfExon(
                exon_id=exon_id,
                chrom=chrom,
                start=start,
                end=end,
                strand=strand,
                exon_number=exon_number,
            )
        )

    def _handle_cds(self, attr: dict, start: int, end: int) -> None:
        transcript_id = attr.get("transcript_id", "")
        gene_id = attr.get("gene_id", "")
        if not transcript_id:
            return
        gene = self._genes_by_id.get(gene_id)
        if gene is None:
            return
        t = gene.transcripts.get(transcript_id)
        if t is None:
            return
        t.cds_records.append((start, end))

    # ── internal: resolve gene ─────────────────────────────────────────

    def _resolve_gene(self, symbol: str) -> _GtfGene:
        """Look up a gene by name or ID."""
        gene = self._genes_by_name.get(symbol.upper())
        if gene is not None:
            return gene
        gene = self._genes_by_id.get(symbol)
        if gene is not None:
            return gene
        raise ValueError(
            f"Gene '{symbol}' not found in GTF.  "
            f"Check that the gene name/ID matches the GTF annotation."
        )

    # ── internal: build TranscriptInfo ─────────────────────────────────

    def _build_transcript_info(self, t: _GtfTranscript) -> TranscriptInfo:
        """Convert internal _GtfTranscript to pipeline TranscriptInfo."""
        exons = [
            ExonCoord(
                exon_id=e.exon_id,
                chrom=e.chrom,
                start=e.start,
                end=e.end,
                strand=e.strand,
            )
            for e in sorted(t.exons, key=lambda e: e.start)
        ]

        cds_exons = self._build_cds_exons(t)

        is_mane_select = "MANE_Select" in t.tags
        is_mane_plus = "MANE_Plus_Clinical" in t.tags
        is_canonical = "Ensembl_canonical" in t.tags

        return TranscriptInfo(
            transcript_id=t.transcript_id,
            biotype=t.biotype,
            is_mane_select=is_mane_select,
            is_mane_plus_clinical=is_mane_plus,
            is_canonical=is_canonical,
            selection_reason="",
            exons=exons,
            cds_exons=cds_exons,
            translation_id=None,
        )

    def _build_cds_exons(self, t: _GtfTranscript) -> list[CdsExon]:
        """Build CDS exons by intersecting CDS records with exons.

        For each exon, find overlapping CDS segments.  Number in
        transcription order (positive: ascending; negative: descending).
        """
        if not t.cds_records:
            return []

        # Build CDS lookup: list of (start, end) sorted by start
        cds_list = sorted(t.cds_records, key=lambda c: c[0])

        raw_segments: list[tuple[int, int, str]] = []
        for exon in sorted(t.exons, key=lambda e: e.start):
            for cds_start, cds_end in cds_list:
                # Intersection of [exon.start, exon.end) and [cds_start, cds_end)
                isect_start = max(exon.start, cds_start)
                isect_end = min(exon.end, cds_end)
                if isect_start < isect_end:
                    raw_segments.append((isect_start, isect_end, exon.exon_id))

        if not raw_segments:
            return []

        # Number in transcription order
        if t.strand == 1:
            raw_segments.sort(key=lambda s: s[0])
        else:
            raw_segments.sort(key=lambda s: -s[0])

        cds_exons = []
        for i, (start, end, exon_id) in enumerate(raw_segments, 1):
            cds_exons.append(
                CdsExon(
                    cds_exon_number=i,
                    cds_exon_id=exon_id,
                    chrom=t.chrom,
                    start=start,
                    end=end,
                    strand=t.strand,
                )
            )

        return cds_exons

    # ── internal: compat helpers ───────────────────────────────────────

    @staticmethod
    def _build_translation_dict(t: _GtfTranscript) -> dict | None:
        """Build a Translation-like dict from CDS records."""
        if not t.cds_records:
            return None
        cds_start = min(s for s, _ in t.cds_records)
        cds_end = max(e for _, e in t.cds_records)
        # Return in 1-based inclusive (Ensembl API convention)
        return {"start": cds_start + 1, "end": cds_end}

    @staticmethod
    def _build_mane_list(t: _GtfTranscript) -> list[dict]:
        """Build a MANE-like list from transcript tags."""
        mane = []
        if "MANE_Select" in t.tags:
            mane.append({"type": "MANE_Select"})
        if "MANE_Plus_Clinical" in t.tags:
            mane.append({"type": "MANE_Plus_Clinical"})
        return mane
