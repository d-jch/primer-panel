# Primer Panel

Primer Panel is a Python CLI for building human hg38 PCR primer panels from
gene symbols. It separates target planning, Primer3 primer design, and
genome-wide in-silico PCR specificity checks into clear pipeline stages.

## What It Does

| Stage | Purpose | Main output |
| --- | --- | --- |
| 1. Target generation | Select coding transcripts from Ensembl and build CDS target windows | `target_summary.tsv`, `targets.bed`, `targets.fa` |
| 2. Primer design | Run Primer3 with Primer3Plus-like defaults | `primers.tsv`, `primers.xlsx` |
| 3. Specificity check | Run UCSC `isPcr` and classify genome-wide hits | `primers_specificity.tsv`, `primers_unique.tsv` |

The pipeline targets CDS regions, not full exons or UTRs. Stage 1 creates the
template and `SEQUENCE_TARGET`; Stage 2 records Primer3 design metrics; Stage 3
is the only stage that reports genomic PCR product coordinates.

## Installation

**Primer Panel does not bundle `primer3_core` or UCSC `isPcr` binaries.**

```bash
# Create environment with all tools
micromamba create -n primer_panel -c conda-forge -c bioconda \
  python=3.11 primer3 ispcr ucsc-fatotwobit \
  requests openpyxl pyfaidx -y
micromamba activate primer_panel

# Install primer-panel itself
pip install primer-panel
```

Or for development: `pip install -e .`

### Verify

```bash
primer-panel --doctor
```

### Tool summary

| Tool | Stage | Install | Required? |
| --- | --- | --- | --- |
| `primer3_core` | 2 | `micromamba install -c bioconda primer3` | Required |
| `isPcr` | 3 | `micromamba install -c bioconda ispcr` | Required |
| `faToTwoBit` | 3 | `micromamba install -c bioconda ucsc-fatotwobit` | Optional (auto-`.2bit`) |
| `bigBedToBed` | dbSNP | `micromamba install -c bioconda ucsc-bigbedtobed` | Optional (.bb support) |
| `pyfaidx` | 1 | `pip install pyfaidx` | Required |
| `openpyxl` | any | `pip install openpyxl` | Optional (XLSX) |

### Docker / Apptainer

```dockerfile
FROM mambaorg/micromamba:latest
RUN micromamba install -y -c conda-forge -c bioconda \
    python=3.11 primer3 ispcr ucsc-fatotwobit \
    requests openpyxl pyfaidx && \
    micromamba clean -afy
RUN pip install primer-panel
```

## Quick Start

```bash
# Stage 1 only: target coordinates, no genome needed
primer-panel --genes HFE HJV TFR2 SLC40A1 HAMP FTH1 \
  --stage targets --output-dir outputs/targets

# Full pipeline
primer-panel --genes HFE HJV TFR2 SLC40A1 HAMP FTH1 \
  --target-size 2700-3300 --primer-flank 500 \
  --genome-fasta /path/to/hg38.fa \
  --common-dbsnp-bed /path/to/dbSnp155Common.bb
```

## Common Options

| Option | Default | Description |
| --- | --- | --- |
| `--genes` | (required) | Gene symbols to process. |
| `--stage` | `all` | `targets`, `design`, `specificity`, or `all`. |
| `--target-size MIN-MAX` | `2700-3300` | CDS target extension size. |
| `--primer-flank N` | `300` | Flanking bases for primer search. |
| `--genome-fasta PATH` | — | hg38 FASTA (required for Stage 2+). |
| `--annotation-gtf PATH` | — | Local Ensembl GTF for offline annotation. |
| `--common-dbsnp-bed PATH` | — | dbSNP file (`.bed` or `.bb`). See [dbSNP](#dbsnp-annotation). |
| `--ispcr-db PATH` | auto | Database for isPcr (auto-`.2bit` if `faToTwoBit` installed). |
| `--ispcr-ooc PATH` | auto | Overused-tile file. |
| `--prepare-ispcr-db` | off | Force (re)create `.2bit`. |
| `--make-ispcr-ooc` | off | Create `.ooc` file. |
| `--doctor` | off | Run dependency checks and exit. |

## isPcr Acceleration

isPcr over a full genome is slow.  A `.2bit` binary database makes it fast.
The pipeline **auto-generates** it on first run if `faToTwoBit` is installed
(one-time, ~1-3 min).  If `faToTwoBit` is missing, it falls back to FASTA with
a one-time hint.

Also supported: `.ooc` overused-tile files (auto-discovered), and `.nib`
databases (legacy).  No flags needed — everything is automatic.

To force rebuild or pre-create: `--prepare-ispcr-db`
(see `primer-panel --help`).

## dbSNP Annotation

`--common-dbsnp-bed` annotates primers with common variants (MAF ≥ 1%) from
dbSNP, flagging primers that overlap high-frequency SNPs — especially at the
3' end where mismatches impair extension.  Accepts `.bed` and `.bb` (bigBed) files.

### Method 1: BED via UCSC table dump (dbSNP 151, easy)

```bash
# Download and convert (~9 MB compressed)
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp151Common.txt.gz
zcat snp151Common.txt.gz | cut -f2,3,4,5 > snp151Common_hg38.bed

primer-panel --genes HFE HJV TFR2 --genome-fasta /path/to/hg38.fa \
  --common-dbsnp-bed snp151Common_hg38.bed --output-dir outputs/run
```

### Method 2: bigBed directly (dbSNP 155+, requires `bigBedToBed`)

```bash
# Install once
micromamba install -c bioconda ucsc-bigbedtobed

# Download and use directly — no conversion (~2 GB)
wget https://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp155Common.bb

primer-panel --genes HFE HJV TFR2 --genome-fasta /path/to/hg38.fa \
  --common-dbsnp-bed dbSnp155Common.bb --output-dir outputs/run
```

Available builds: `dbSnp153Common.bb`, `dbSnp155Common.bb`.
Check `https://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/` for newer releases.

> snp151's ~6M common variants are generally sufficient — common variant calls
> (MAF ≥ 1%) are stable between builds.

### Output annotation

`primers_specificity.tsv` includes: `common_snp_risk` (none/medium/high),
`*_common_snp_count`, `*_3p_common_snp_count`, `common_snp_hits`.

## Key Outputs

| File | Stage | Description |
| --- | --- | --- |
| `target_summary.tsv` / `.xlsx` | 1 | Coordinates, transcript, CDS coverage, QC status. |
| `targets.bed` | 1 | Design-template BED6. |
| `required_regions.bed` | 1 | Raw CDS regions that must be covered. |
| `targets.fa` | 1 | Template sequences (N-placeholder unless `--genome-fasta`). |
| `failed_targets.tsv` | 1 | Genes that could not be processed (if any). |
| `primers.tsv` / `.xlsx` | 2 | Primer3 pairs and design metrics. |
| `primers_specificity.tsv` | 3 | Primers + isPcr hit classification + dbSNP annotation. |
| `primers_unique.tsv` | 3 | Subset: `unique_pass` only. |
| `stage3_summary.txt` | 3 | Per-target specificity summary. |
| `run_summary.txt` | all | QC summary for the full run. |

## Rescue

Re-run Stage 2+3 for specific targets with relaxed parameters.  Useful when
a target has too few clean primers after the main run.

```bash
# Run rescue against an existing output directory (no --genes)
primer-panel --output-dir outputs/full \
  --genome-fasta /path/to/hg38.fa \
  --common-dbsnp-bed /path/to/dbSnp155Common.bb \
  --rescue-target DENND3_cds18 --rescue-flank 500 --rescue-num-return 20
```

| Flag | Default | Effect |
| --- | --- | --- |
| `--rescue-target ID [ID ...]` | (required) | Target IDs to rescue, e.g. `DENND3_cds18`. |
| `--rescue-flank` | same as `--primer-flank` | Wider template for Primer3 search. |
| `--rescue-num-return` | 20 | More Primer3 candidates. |

Rescue loads existing outputs from `--output-dir`, re-runs Stage 2+3 for the
specified targets, and merges results back into all output files.  Primers from
the rescue pass have ranks offset by +100.

## Specificity Status

| Status | Meaning |
| --- | --- |
| `unique_pass` | One hit, matching expected chromosome and coordinates. |
| `unique_off_target` | One hit, but not at the expected location. |
| `multi_hit` | Multiple genome-wide products. |
| `no_hit` | No product detected. |
| `pcr_error` | `isPcr` failed for the batch. |

## Dependency Doctor

```bash
primer-panel --doctor
```

Reports `primer3_core`, `isPcr`, `faToTwoBit`, `bigBedToBed`, genome FASTA,
and common dbSNP status.  Does not require `--genes`.

## Offline Annotation (Ensembl GTF)

By default, Stage 1 queries the Ensembl REST API.  For offline use (or to pin
a specific release), provide a local Ensembl GTF:

```bash
primer-panel --genes HFE HJV TFR2 --stage targets \
  --annotation-gtf /path/to/Homo_sapiens.GRCh38.110.gtf.gz \
  --output-dir outputs/offline
```

GTF provides annotation only — for the full pipeline, also pass
`--genome-fasta`.  Use GTF and FASTA from the same Ensembl release.

Transcript selection priority: `MANE_Select` → `Ensembl_canonical` → longest
CDS among `protein_coding` transcripts.

## Development

```bash
python -m pytest
python -m compileall -q primer_panel tests
python -m build --no-isolation --sdist --wheel --outdir /tmp/primer-panel-build
```

## License

MIT. See [LICENSE](LICENSE).
