"""Tests for local Ensembl GTF annotation client."""

from __future__ import annotations

import gzip
import textwrap
from pathlib import Path
from unittest.mock import patch

import pytest

from primer_panel.gtf_annotation import GtfAnnotationClient, _parse_attributes
from primer_panel.ensembl_client import CdsExon, ExonCoord, TranscriptInfo


# ── Fixture GTF data ───────────────────────────────────────────────────────

# Positive strand gene with 3 CDS exons
POSITIVE_GTF = textwrap.dedent("""\
chr6\tensembl\tgene\t28000001\t28004000\t.\t+\t.\tgene_id "ENSG0000001"; gene_name "TESTPOS"; gene_biotype "protein_coding";
chr6\tensembl\ttranscript\t28000001\t28004000\t.\t+\t.\tgene_id "ENSG0000001"; transcript_id "ENST0000001"; gene_name "TESTPOS"; transcript_biotype "protein_coding"; tag "Ensembl_canonical";
chr6\tensembl\texon\t28000001\t28000500\t.\t+\t.\tgene_id "ENSG0000001"; transcript_id "ENST0000001"; exon_id "ENSE0000001"; exon_number "1";
chr6\tensembl\texon\t28001001\t28001500\t.\t+\t.\tgene_id "ENSG0000001"; transcript_id "ENST0000001"; exon_id "ENSE0000002"; exon_number "2";
chr6\tensembl\texon\t28002001\t28004000\t.\t+\t.\tgene_id "ENSG0000001"; transcript_id "ENST0000001"; exon_id "ENSE0000003"; exon_number "3";
chr6\tensembl\tCDS\t28000201\t28000500\t.\t+\t0\tgene_id "ENSG0000001"; transcript_id "ENST0000001";
chr6\tensembl\tCDS\t28001001\t28001300\t.\t+\t0\tgene_id "ENSG0000001"; transcript_id "ENST0000001";
chr6\tensembl\tCDS\t28002001\t28002500\t.\t+\t0\tgene_id "ENSG0000001"; transcript_id "ENST0000001";
""")

# Negative strand gene with 3 CDS exons
NEGATIVE_GTF = textwrap.dedent("""\
chr17\tensembl\tgene\t41200001\t41204000\t.\t-\t.\tgene_id "ENSG0000002"; gene_name "TESTNEG"; gene_biotype "protein_coding";
chr17\tensembl\ttranscript\t41200001\t41204000\t.\t-\t.\tgene_id "ENSG0000002"; transcript_id "ENST0000002"; gene_name "TESTNEG"; transcript_biotype "protein_coding"; tag "Ensembl_canonical";
chr17\tensembl\texon\t41200001\t41200500\t.\t-\t.\tgene_id "ENSG0000002"; transcript_id "ENST0000002"; exon_id "ENSE0000010"; exon_number "1";
chr17\tensembl\texon\t41201001\t41201500\t.\t-\t.\tgene_id "ENSG0000002"; transcript_id "ENST0000002"; exon_id "ENSE0000011"; exon_number "2";
chr17\tensembl\texon\t41202001\t41204000\t.\t-\t.\tgene_id "ENSG0000002"; transcript_id "ENST0000002"; exon_id "ENSE0000012"; exon_number "3";
chr17\tensembl\tCDS\t41202001\t41202500\t.\t-\t0\tgene_id "ENSG0000002"; transcript_id "ENST0000002";
chr17\tensembl\tCDS\t41201001\t41201300\t.\t-\t0\tgene_id "ENSG0000002"; transcript_id "ENST0000002";
chr17\tensembl\tCDS\t41200201\t41200500\t.\t-\t0\tgene_id "ENSG0000002"; transcript_id "ENST0000002";
""")

# Multi-transcript gene: canonical vs longer non-canonical
MULTI_TX_GTF = textwrap.dedent("""\
chr1\tensembl\tgene\t1000001\t1003000\t.\t+\t.\tgene_id "ENSG0000003"; gene_name "TESTMULTI"; gene_biotype "protein_coding";
chr1\tensembl\ttranscript\t1000001\t1001000\t.\t+\t.\tgene_id "ENSG0000003"; transcript_id "ENST0000010"; gene_name "TESTMULTI"; transcript_biotype "protein_coding"; tag "Ensembl_canonical";
chr1\tensembl\texon\t1000001\t1000500\t.\t+\t.\tgene_id "ENSG0000003"; transcript_id "ENST0000010"; exon_id "ENSE0000020"; exon_number "1";
chr1\tensembl\texon\t1000701\t1001000\t.\t+\t.\tgene_id "ENSG0000003"; transcript_id "ENST0000010"; exon_id "ENSE0000021"; exon_number "2";
chr1\tensembl\tCDS\t1000101\t1000500\t.\t+\t0\tgene_id "ENSG0000003"; transcript_id "ENST0000010";
chr1\tensembl\tCDS\t1000701\t1000900\t.\t+\t0\tgene_id "ENSG0000003"; transcript_id "ENST0000010";
chr1\tensembl\ttranscript\t1000001\t1003000\t.\t+\t.\tgene_id "ENSG0000003"; transcript_id "ENST0000011"; gene_name "TESTMULTI"; transcript_biotype "protein_coding";
chr1\tensembl\texon\t1000001\t1000500\t.\t+\t.\tgene_id "ENSG0000003"; transcript_id "ENST0000011"; exon_id "ENSE0000022"; exon_number "1";
chr1\tensembl\texon\t1001001\t1001500\t.\t+\t.\tgene_id "ENSG0000003"; transcript_id "ENST0000011"; exon_id "ENSE0000023"; exon_number "2";
chr1\tensembl\texon\t1002001\t1003000\t.\t+\t.\tgene_id "ENSG0000003"; transcript_id "ENST0000011"; exon_id "ENSE0000024"; exon_number "3";
chr1\tensembl\tCDS\t1000101\t1000500\t.\t+\t0\tgene_id "ENSG0000003"; transcript_id "ENST0000011";
chr1\tensembl\tCDS\t1001001\t1001500\t.\t+\t0\tgene_id "ENSG0000003"; transcript_id "ENST0000011";
chr1\tensembl\tCDS\t1002001\t1002800\t.\t+\t0\tgene_id "ENSG0000003"; transcript_id "ENST0000011";
""")

# MANE_Select transcript
MANE_GTF = textwrap.dedent("""\
chr7\tensembl\tgene\t55000001\t55002000\t.\t+\t.\tgene_id "ENSG0000004"; gene_name "TESTMANE"; gene_biotype "protein_coding";
chr7\tensembl\ttranscript\t55000001\t55001000\t.\t+\t.\tgene_id "ENSG0000004"; transcript_id "ENST0000020"; gene_name "TESTMANE"; transcript_biotype "protein_coding"; tag "Ensembl_canonical";
chr7\tensembl\texon\t55000001\t55000500\t.\t+\t.\tgene_id "ENSG0000004"; transcript_id "ENST0000020"; exon_id "ENSE0000030"; exon_number "1";
chr7\tensembl\tCDS\t55000101\t55000500\t.\t+\t0\tgene_id "ENSG0000004"; transcript_id "ENST0000020";
chr7\tensembl\ttranscript\t55000001\t55002000\t.\t+\t.\tgene_id "ENSG0000004"; transcript_id "ENST0000021"; gene_name "TESTMANE"; transcript_biotype "protein_coding"; tag "MANE_Select,Ensembl_canonical";
chr7\tensembl\texon\t55000001\t55000500\t.\t+\t.\tgene_id "ENSG0000004"; transcript_id "ENST0000021"; exon_id "ENSE0000031"; exon_number "1";
chr7\tensembl\texon\t55001001\t55002000\t.\t+\t.\tgene_id "ENSG0000004"; transcript_id "ENST0000021"; exon_id "ENSE0000032"; exon_number "2";
chr7\tensembl\tCDS\t55000101\t55000500\t.\t+\t0\tgene_id "ENSG0000004"; transcript_id "ENST0000021";
chr7\tensembl\tCDS\t55001001\t55001500\t.\t+\t0\tgene_id "ENSG0000004"; transcript_id "ENST0000021";
""")

# MANE with real repeated tag format (as found in Ensembl GTF files)
MANE_REAL_TAGS_GTF = textwrap.dedent("""\
chr7\tensembl\tgene\t55000001\t55002000\t.\t+\t.\tgene_id "ENSG0000005"; gene_name "TESTREAL"; gene_biotype "protein_coding";
chr7\tensembl\ttranscript\t55000001\t55001000\t.\t+\t.\tgene_id "ENSG0000005"; transcript_id "ENST0000030"; gene_name "TESTREAL"; transcript_biotype "protein_coding"; tag "basic"; tag "Ensembl_canonical";
chr7\tensembl\texon\t55000001\t55000500\t.\t+\t.\tgene_id "ENSG0000005"; transcript_id "ENST0000030"; exon_id "ENSE0000040"; exon_number "1";
chr7\tensembl\tCDS\t55000101\t55000500\t.\t+\t0\tgene_id "ENSG0000005"; transcript_id "ENST0000030";
chr7\tensembl\ttranscript\t55000001\t55002000\t.\t+\t.\tgene_id "ENSG0000005"; transcript_id "ENST0000031"; gene_name "TESTREAL"; transcript_biotype "protein_coding"; tag "basic"; tag "MANE_Select"; tag "appris_principal_1";
chr7\tensembl\texon\t55000001\t55000500\t.\t+\t.\tgene_id "ENSG0000005"; transcript_id "ENST0000031"; exon_id "ENSE0000041"; exon_number "1";
chr7\tensembl\texon\t55001001\t55002000\t.\t+\t.\tgene_id "ENSG0000005"; transcript_id "ENST0000031"; exon_id "ENSE0000042"; exon_number "2";
chr7\tensembl\tCDS\t55000101\t55000500\t.\t+\t0\tgene_id "ENSG0000005"; transcript_id "ENST0000031";
chr7\tensembl\tCDS\t55001001\t55001500\t.\t+\t0\tgene_id "ENSG0000005"; transcript_id "ENST0000031";
""")

# MANE_Plus_Clinical over canonical
MANE_PLUS_GTF = textwrap.dedent("""\
chr3\tensembl\tgene\t10000001\t10002000\t.\t+\t.\tgene_id "ENSG0000006"; gene_name "TESTMANEPLUS"; gene_biotype "protein_coding";
chr3\tensembl\ttranscript\t10000001\t10001000\t.\t+\t.\tgene_id "ENSG0000006"; transcript_id "ENST0000040"; gene_name "TESTMANEPLUS"; transcript_biotype "protein_coding"; tag "Ensembl_canonical";
chr3\tensembl\texon\t10000001\t10000500\t.\t+\t.\tgene_id "ENSG0000006"; transcript_id "ENST0000040"; exon_id "ENSE0000050"; exon_number "1";
chr3\tensembl\tCDS\t10000101\t10000500\t.\t+\t0\tgene_id "ENSG0000006"; transcript_id "ENST0000040";
chr3\tensembl\ttranscript\t10000001\t10002000\t.\t+\t.\tgene_id "ENSG0000006"; transcript_id "ENST0000041"; gene_name "TESTMANEPLUS"; transcript_biotype "protein_coding"; tag "MANE_Plus_Clinical";
chr3\tensembl\texon\t10000001\t10000500\t.\t+\t.\tgene_id "ENSG0000006"; transcript_id "ENST0000041"; exon_id "ENSE0000051"; exon_number "1";
chr3\tensembl\texon\t10001001\t10002000\t.\t+\t.\tgene_id "ENSG0000006"; transcript_id "ENST0000041"; exon_id "ENSE0000052"; exon_number "2";
chr3\tensembl\tCDS\t10000101\t10000500\t.\t+\t0\tgene_id "ENSG0000006"; transcript_id "ENST0000041";
chr3\tensembl\tCDS\t10001001\t10001500\t.\t+\t0\tgene_id "ENSG0000006"; transcript_id "ENST0000041";
""")

# Gene with non-protein_coding transcript that has CDS, and protein_coding transcript
MIXED_BIOTYPE_GTF = textwrap.dedent("""\
chr5\tensembl\tgene\t50000001\t50002000\t.\t+\t.\tgene_id "ENSG0000007"; gene_name "TESTMIXED"; gene_biotype "protein_coding";
chr5\tensembl\ttranscript\t50000001\t50002000\t.\t+\t.\tgene_id "ENSG0000007"; transcript_id "ENST0000050"; gene_name "TESTMIXED"; transcript_biotype "nonsense_mediated_decay";
chr5\tensembl\texon\t50000001\t50000500\t.\t+\t.\tgene_id "ENSG0000007"; transcript_id "ENST0000050"; exon_id "ENSE0000060"; exon_number "1";
chr5\tensembl\texon\t50001001\t50002000\t.\t+\t.\tgene_id "ENSG0000007"; transcript_id "ENST0000050"; exon_id "ENSE0000061"; exon_number "2";
chr5\tensembl\tCDS\t50000101\t50000500\t.\t+\t0\tgene_id "ENSG0000007"; transcript_id "ENST0000050";
chr5\tensembl\tCDS\t50001001\t50001500\t.\t+\t0\tgene_id "ENSG0000007"; transcript_id "ENST0000050";
chr5\tensembl\ttranscript\t50000001\t50001000\t.\t+\t.\tgene_id "ENSG0000007"; transcript_id "ENST0000051"; gene_name "TESTMIXED"; transcript_biotype "protein_coding";
chr5\tensembl\texon\t50000001\t50000500\t.\t+\t.\tgene_id "ENSG0000007"; transcript_id "ENST0000051"; exon_id "ENSE0000062"; exon_number "1";
chr5\tensembl\tCDS\t50000101\t50000500\t.\t+\t0\tgene_id "ENSG0000007"; transcript_id "ENST0000051";
""")

# Gene with only non-protein_coding transcripts (fallback test)
NON_CODING_GTF = textwrap.dedent("""\
chr9\tensembl\tgene\t30000001\t30001000\t.\t+\t.\tgene_id "ENSG0000008"; gene_name "TESTNONCODING"; gene_biotype "lncRNA";
chr9\tensembl\ttranscript\t30000001\t30001000\t.\t+\t.\tgene_id "ENSG0000008"; transcript_id "ENST0000060"; gene_name "TESTNONCODING"; transcript_biotype "lncRNA";
chr9\tensembl\texon\t30000001\t30000500\t.\t+\t.\tgene_id "ENSG0000008"; transcript_id "ENST0000060"; exon_id "ENSE0000070"; exon_number "1";
chr9\tensembl\texon\t30000701\t30001000\t.\t+\t.\tgene_id "ENSG0000008"; transcript_id "ENST0000060"; exon_id "ENSE0000071"; exon_number "2";
""")


# ── Helpers ────────────────────────────────────────────────────────────────


def _write_gtf(tmp_path: Path, content: str, name: str = "test.gtf") -> Path:
    p = tmp_path / name
    p.write_text(content)
    return p


def _write_gtf_gz(tmp_path: Path, content: str, name: str = "test.gtf.gz") -> Path:
    p = tmp_path / name
    with gzip.open(p, "wt") as f:
        f.write(content)
    return p


# ── Attribute parser ───────────────────────────────────────────────────────


class TestParseAttributes:
    def test_basic(self):
        attrs = _parse_attributes('gene_id "ENSG001"; gene_name "TP53";')
        assert attrs == {"gene_id": "ENSG001", "gene_name": "TP53"}

    def test_tag_field(self):
        attrs = _parse_attributes('tag "MANE_Select,Ensembl_canonical";')
        assert attrs["tag"] == "MANE_Select,Ensembl_canonical"

    def test_empty(self):
        assert _parse_attributes("") == {}

    def test_no_quotes(self):
        assert _parse_attributes("malformed") == {}

    def test_repeated_tag_attributes(self):
        """Real Ensembl GTF uses separate tag lines, not comma-joined."""
        attrs = _parse_attributes(
            'gene_id "ENSG001"; tag "basic"; tag "MANE_Select"; tag "Ensembl_canonical";'
        )
        assert "basic" in attrs["tag"]
        assert "MANE_Select" in attrs["tag"]
        assert "Ensembl_canonical" in attrs["tag"]

    def test_repeated_tag_parsed_as_set(self):
        """Tags from repeated attributes must be individually addressable."""
        attrs = _parse_attributes('tag "basic"; tag "MANE_Select";')
        tags = {t.strip() for t in attrs["tag"].split(",") if t.strip()}
        assert tags == {"basic", "MANE_Select"}


# ── Positive strand gene ───────────────────────────────────────────────────


class TestPositiveStrand:
    def test_select_transcript(self, tmp_path):
        gtf = _write_gtf(tmp_path, POSITIVE_GTF)
        client = GtfAnnotationClient(gtf)
        ti, reason = client.select_transcript("TESTPOS")

        assert isinstance(ti, TranscriptInfo)
        assert ti.transcript_id == "ENST0000001"
        assert ti.biotype == "protein_coding"
        assert ti.is_canonical is True
        assert ti.cds_exons[0].strand == 1
        assert len(ti.exons) == 3
        assert len(ti.cds_exons) == 3

    def test_cds_exon_order_positive(self, tmp_path):
        gtf = _write_gtf(tmp_path, POSITIVE_GTF)
        client = GtfAnnotationClient(gtf)
        ti, _ = client.select_transcript("TESTPOS")

        # Positive strand: CDS exons numbered in ascending genome order
        assert ti.cds_exons[0].cds_exon_number == 1
        assert ti.cds_exons[1].cds_exon_number == 2
        assert ti.cds_exons[2].cds_exon_number == 3
        assert ti.cds_exons[0].start < ti.cds_exons[1].start

    def test_cds_coordinates(self, tmp_path):
        gtf = _write_gtf(tmp_path, POSITIVE_GTF)
        client = GtfAnnotationClient(gtf)
        ti, _ = client.select_transcript("TESTPOS")

        # First CDS: chr6:28000201-28000500 (0-based: 28000200-28000500)
        cds1 = ti.cds_exons[0]
        assert cds1.chrom == "chr6"
        assert cds1.start == 28000200  # 1-based 28000201 → 0-based 28000200
        assert cds1.end == 28000500
        assert cds1.strand == 1

    def test_gene_id_lookup(self, tmp_path):
        gtf = _write_gtf(tmp_path, POSITIVE_GTF)
        client = GtfAnnotationClient(gtf)
        ti, _ = client.select_transcript("ENSG0000001")
        assert ti.transcript_id == "ENST0000001"

    def test_lookup_gene(self, tmp_path):
        gtf = _write_gtf(tmp_path, POSITIVE_GTF)
        client = GtfAnnotationClient(gtf)
        data = client.lookup_gene("TESTPOS")

        assert "start" in data
        assert "end" in data
        assert "Transcript" in data
        assert len(data["Transcript"]) >= 1


# ── Negative strand gene ───────────────────────────────────────────────────


class TestNegativeStrand:
    def test_select_transcript(self, tmp_path):
        gtf = _write_gtf(tmp_path, NEGATIVE_GTF)
        client = GtfAnnotationClient(gtf)
        ti, reason = client.select_transcript("TESTNEG")

        assert ti.cds_exons[0].strand == -1
        assert len(ti.cds_exons) == 3

    def test_cds_exon_order_negative(self, tmp_path):
        gtf = _write_gtf(tmp_path, NEGATIVE_GTF)
        client = GtfAnnotationClient(gtf)
        ti, _ = client.select_transcript("TESTNEG")

        # Negative strand: CDS exons numbered in descending genome order
        # Genome positions: 41200200-41200500, 41201000-41201300, 41202000-41202500
        # Descending: 41202000 (1), 41201000 (2), 41200200 (3)
        assert ti.cds_exons[0].cds_exon_number == 1
        assert ti.cds_exons[0].start > ti.cds_exons[1].start  # descending genome order
        assert ti.cds_exons[1].cds_exon_number == 2
        assert ti.cds_exons[2].cds_exon_number == 3


# ── Multi-transcript selection ─────────────────────────────────────────────


class TestMultiTranscript:
    def test_canonical_selected_over_longer(self, tmp_path):
        """Canonical transcript should be selected even if shorter."""
        gtf = _write_gtf(tmp_path, MULTI_TX_GTF)
        client = GtfAnnotationClient(gtf)
        ti, reason = client.select_transcript("TESTMULTI")

        assert ti.transcript_id == "ENST0000010"
        assert reason == "canonical_protein_coding"
        assert ti.is_canonical is True

    def test_longest_when_no_canonical(self, tmp_path):
        """Without canonical tag, longest CDS should win."""
        # Remove canonical tag from the shorter transcript
        modified = MULTI_TX_GTF.replace(
            'tag "Ensembl_canonical"', 'tag "basic"'
        )
        gtf = _write_gtf(tmp_path, modified)
        client = GtfAnnotationClient(gtf)
        ti, reason = client.select_transcript("TESTMULTI")

        # ENST0000011 has more CDS total
        assert ti.transcript_id == "ENST0000011"
        assert reason == "longest_protein_coding"


# ── MANE_Select priority ───────────────────────────────────────────────────


class TestManeSelect:
    def test_mane_select_wins(self, tmp_path):
        """MANE_Select transcript should be preferred over canonical."""
        gtf = _write_gtf(tmp_path, MANE_GTF)
        client = GtfAnnotationClient(gtf)
        ti, reason = client.select_transcript("TESTMANE")

        assert ti.transcript_id == "ENST0000021"
        assert reason == "MANE_Select"
        assert ti.is_mane_select is True

    def test_mane_select_real_repeated_tags(self, tmp_path):
        """Real Ensembl GTF with repeated tag lines must still detect MANE_Select."""
        gtf = _write_gtf(tmp_path, MANE_REAL_TAGS_GTF)
        client = GtfAnnotationClient(gtf)
        ti, reason = client.select_transcript("TESTREAL")

        assert ti.transcript_id == "ENST0000031"
        assert reason == "MANE_Select"
        assert ti.is_mane_select is True

    def test_mane_plus_clinical_over_canonical(self, tmp_path):
        """MANE_Plus_Clinical should rank above Ensembl_canonical."""
        gtf = _write_gtf(tmp_path, MANE_PLUS_GTF)
        client = GtfAnnotationClient(gtf)
        ti, reason = client.select_transcript("TESTMANEPLUS")

        assert ti.transcript_id == "ENST0000041"
        assert reason == "MANE_Plus_Clinical"
        assert ti.is_mane_plus_clinical is True

    def test_protein_coding_preferred_over_non_coding_with_cds(self, tmp_path):
        """protein_coding transcript should be selected over NMD transcript with CDS."""
        gtf = _write_gtf(tmp_path, MIXED_BIOTYPE_GTF)
        client = GtfAnnotationClient(gtf)
        ti, reason = client.select_transcript("TESTMIXED")

        assert ti.transcript_id == "ENST0000051"
        assert ti.biotype == "protein_coding"

    def test_fallback_non_protein_coding(self, tmp_path):
        """When no protein_coding exists, fallback to longest transcript."""
        gtf = _write_gtf(tmp_path, NON_CODING_GTF)
        client = GtfAnnotationClient(gtf)
        ti, reason = client.select_transcript("TESTNONCODING")

        assert ti.transcript_id == "ENST0000060"
        assert reason == "longest_transcript"


# ── Gzip support ───────────────────────────────────────────────────────────


class TestGzipSupport:
    def test_read_gzip(self, tmp_path):
        gtf = _write_gtf_gz(tmp_path, POSITIVE_GTF)
        client = GtfAnnotationClient(gtf)
        ti, _ = client.select_transcript("TESTPOS")
        assert ti.transcript_id == "ENST0000001"


# ── Error cases ────────────────────────────────────────────────────────────


class TestErrors:
    def test_gene_not_found(self, tmp_path):
        gtf = _write_gtf(tmp_path, POSITIVE_GTF)
        client = GtfAnnotationClient(gtf)
        with pytest.raises(ValueError, match="not found in GTF"):
            client.select_transcript("NONEXISTENT")

    def test_gene_not_found_by_name(self, tmp_path):
        gtf = _write_gtf(tmp_path, POSITIVE_GTF)
        client = GtfAnnotationClient(gtf)
        with pytest.raises(ValueError, match="not found in GTF"):
            client.select_transcript("FAKEGENE")


# ── TranscriptInfo structure matches EnsemblClient output ──────────────────


class TestTranscriptInfoCompat:
    """Ensure GTF output is compatible with pipeline expectations."""

    def test_has_exons(self, tmp_path):
        gtf = _write_gtf(tmp_path, POSITIVE_GTF)
        client = GtfAnnotationClient(gtf)
        ti, _ = client.select_transcript("TESTPOS")
        assert all(isinstance(e, ExonCoord) for e in ti.exons)

    def test_has_cds_exons(self, tmp_path):
        gtf = _write_gtf(tmp_path, POSITIVE_GTF)
        client = GtfAnnotationClient(gtf)
        ti, _ = client.select_transcript("TESTPOS")
        assert all(isinstance(c, CdsExon) for c in ti.cds_exons)

    def test_exon_coords_zero_based(self, tmp_path):
        gtf = _write_gtf(tmp_path, POSITIVE_GTF)
        client = GtfAnnotationClient(gtf)
        ti, _ = client.select_transcript("TESTPOS")
        # First exon: GTF 28000001-28000500 → 0-based 28000000-28000500
        assert ti.exons[0].start == 28000000
        assert ti.exons[0].end == 28000500


# ── Integration: --annotation-gtf skips EnsemblClient ──────────────────────


class TestMainGtfIntegration:
    """Test that --annotation-gtf uses GTF client instead of Ensembl API."""

    @patch("primer_panel.preflight.shutil.which", return_value="/usr/bin/primer3_core")
    def test_gtf_mode_no_api_calls(self, mock_which, tmp_path):
        """With --annotation-gtf, EnsemblClient should not be instantiated."""
        gtf = _write_gtf(tmp_path, POSITIVE_GTF, "annot.gtf")
        out = tmp_path / "out"

        # We can't run the full pipeline without network + Primer3, but we
        # can verify that main() creates a GtfAnnotationClient, not an
        # EnsemblClient, by checking the import path.
        from primer_panel.main import _create_annotation_client, _resolve_annotation_source

        # Simulate args
        args = type("Args", (), {
            "annotation_gtf": gtf,
            "annotation_source": "auto",
        })()
        cfg = type("Cfg", (), {})()

        source = _resolve_annotation_source(args)
        assert source == "gtf"

        client = _create_annotation_client(args, cfg)
        assert isinstance(client, GtfAnnotationClient)

    def test_annotation_source_auto_no_gtf(self):
        """auto without --annotation-gtf should default to ensembl-api."""
        from primer_panel.main import _resolve_annotation_source

        args = type("Args", (), {
            "annotation_gtf": None,
            "annotation_source": "auto",
        })()

        source = _resolve_annotation_source(args)
        assert source == "ensembl-api"

    def test_annotation_source_gtf_without_path_exits(self, tmp_path):
        """--annotation-source gtf without --annotation-gtf should exit."""
        from primer_panel.main import _resolve_annotation_source

        args = type("Args", (), {
            "annotation_gtf": None,
            "annotation_source": "gtf",
        })()

        with pytest.raises(SystemExit):
            _resolve_annotation_source(args)
