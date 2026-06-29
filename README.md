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

## Recommended Installation

**Primer Panel does not bundle `primer3_core` or UCSC `isPcr` binaries.**
These tools have separate licenses and platform-specific build requirements.
Install them via micromamba/conda or a container runtime.

### Create the environment (recommended)

```bash
# Create environment with all tools
micromamba create -n primer_panel -c conda-forge -c bioconda \
  python=3.11 primer3 ispcr ucsc-fatotwobit \
  requests openpyxl pyfaidx -y

micromamba activate primer_panel

# Install primer-panel itself
pip install primer-panel
```

Or for development:

```bash
pip install -e .
```

### Verify installation

```bash
# Check everything is in order
primer-panel --doctor

# Or check individual tools
primer3_core --version
isPcr
faToTwoBit
```

### Tool summary

| Tool | Stage | Install | Required? |
| --- | --- | --- | --- |
| `primer3_core` | 2 (primer design) | `micromamba install -c bioconda primer3` | For Stage 2+ |
| `isPcr` | 3 (specificity) | `micromamba install -c bioconda ispcr` | For Stage 3 |
| `faToTwoBit` | 3 (db prep) | `micromamba install -c bioconda ucsc-fatotwobit` | Optional |
| `pyfaidx` | 1 (sequences) | `pip install pyfaidx` | For real sequences |
| `openpyxl` | any (XLSX) | `pip install openpyxl` | Optional |

### Alternative: Docker / Apptainer

For containerized environments, you can build an image that includes all
bioinformatics tools.  This avoids platform-specific binary issues:

```dockerfile
FROM mambaorg/micromamba:latest
RUN micromamba install -y -c conda-forge -c bioconda \
    python=3.11 primer3 ispcr ucsc-fatotwobit \
    requests openpyxl pyfaidx && \
    micromamba clean -afy
RUN pip install primer-panel
```

## Quick Start

Run only Stage 1 (target coordinates, no primer design):

```bash
primer-panel \
  --genes HFE HJV TFR2 SLC40A1 HAMP FTH1 \
  --stage targets \
  --output-dir outputs/hcc6_targets
```

Run the full pipeline (Stage 1+2+3, this is the default when no `--stage` is given):

```bash
primer-panel \
  --genes HFE HJV TFR2 SLC40A1 HAMP FTH1 \
  --target-size 2700-3300 \
  --primer-flank 500 \
  --genome-fasta /path/to/hg38.fa \
  --output-dir outputs/hcc6_primers
```

> **Note:** The default stage is `all` (Stage 1+2+3).  `--genome-fasta` is
> required because Stage 2 needs real sequences for Primer3 and Stage 3 needs
> the genome for specificity checks.  If you only want target coordinates, use
> `--stage targets`.

Run Stage 1+2 (targets + primer design, no specificity check):

```bash
primer-panel \
  --genes HFE HJV TFR2 SLC40A1 HAMP FTH1 \
  --stage design \
  --genome-fasta /path/to/hg38.fa \
  --output-dir outputs/hcc6_primers
```

Annotate primers with common dbSNP variants:

```bash
primer-panel \
  --genes HFE HJV TFR2 SLC40A1 HAMP FTH1 \
  --genome-fasta /path/to/hg38.fa \
  --common-dbsnp-bed /path/to/common_snps.bed \
  --output-dir outputs/hcc6_primers
```

Finalize a panel from Stage 2 and Stage 3 outputs:

```bash
primer-panel-finalize \
  --input-dir outputs/hcc6_primers \
  --output-dir outputs/panel_final \
  --genome-fasta /path/to/hg38.fa
```

## Common Options

| Option | Default | Description |
| --- | --- | --- |
| `--genes` | required | Gene symbols to process. |
| `--stage STAGE` | `all` | Pipeline stage: `targets`, `design`, `specificity`, or `all`. Default runs full pipeline (requires `--genome-fasta`). |
| `--target-size MIN-MAX` | `2700-3300` | Stage 1 CDS grouping and target extension size. This does not constrain Primer3 product size. |
| `--primer-flank N` | `300` | Bases added on both sides of each target for primer search. |
| `--genome-fasta PATH` | required for full pipeline | hg38 FASTA used for real sequence extraction and specificity checks. |
| `--annotation-gtf PATH` | none | Local Ensembl GTF for offline annotation (replaces Ensembl REST API). |
| `--annotation-source` | `auto` | Annotation source: `auto`, `ensembl-api`, or `gtf`. |
| `--common-dbsnp-bed PATH` | none | Common dbSNP BED file for primer risk annotation. |
| `--doctor` | off | Run dependency checks and exit. |
| `--design-primers` | off | [Deprecated] Use `--stage` instead. |
| `--check-specificity` | off | [Deprecated] Use `--stage` instead. |
| `--primer3-bin PATH` | `primer3_core` | Primer3 executable. |
| `--is-pcr-bin PATH` | `isPcr` | UCSC `isPcr` executable. |
| `--ispcr-db PATH` | auto | Explicit `.2bit`/`.nib` database for isPcr. |
| `--ispcr-ooc PATH` | auto | Explicit overused-tile (`.ooc`) file for isPcr. |
| `--ispcr-tile-size N` | `11` | Tile size for isPcr. |
| `--prepare-ispcr-db` | off | Create a `.2bit` database from the genome FASTA. |
| `--make-ispcr-ooc` | off | Create an overused-tile (`.ooc`) file. |

## isPcr Acceleration

Running isPcr on a full hg38 genome can be slow. Two pre-computed files speed
it up significantly:

- **`.2bit` database** — a compact binary version of the genome FASTA. isPcr
  reads this much faster than a gzipped FASTA.
- **`.ooc` (overused-tile) file** — pre-computed mask for repetitive k-mers.
  Skips alignment against high-frequency tiles during the genome scan.

### Auto-discovery (default)

When you pass `--genome-fasta /path/to/hg38.fa`, the pipeline automatically
looks for these files in the same directory:

| File searched | Example |
| --- | --- |
| Same-basename `.2bit` | `hg38.2bit` |
| Same-basename `.nib` | `hg38.nib` |
| `*.<tileSize>.ooc` | `hg38.11.ooc` |

If found, they are used automatically. No extra flags needed.

### Explicit creation

To create these files yourself (one-time setup):

```bash
# Create .2bit from FASTA (requires UCSC faToTwoBit)
primer-panel \
  --genes HFE --genome-fasta /path/to/hg38.fa \
  --prepare-ispcr-db --check-specificity --design-primers \
  --output-dir outputs/hcc6

# Create .ooc overused-tile file
primer-panel \
  --genes HFE --genome-fasta /path/to/hg38.fa \
  --make-ispcr-ooc --check-specificity --design-primers \
  --output-dir outputs/hcc6
```

You can also pass explicit paths if the files are elsewhere:

```bash
primer-panel \
  --genes HFE --genome-fasta /path/to/hg38.fa \
  --ispcr-db /data/hg38.2bit --ispcr-ooc /data/hg38.11.ooc \
  --check-specificity --design-primers \
  --output-dir outputs/hcc6
```

## Key Outputs

### Stage 1 — Target Generation

| File | Description |
| --- | --- |
| `targets.bed` | Design-template coordinates in BED6 format. |
| `required_regions.bed` | Raw CDS regions that must be covered. |
| `target_summary.tsv` | Target coordinates, selected transcript, CDS coverage, and QC status. |
| `targets.fa` | Design-template FASTA. Uses placeholder `N` sequence unless `--genome-fasta` is provided. |
| `target_summary.xlsx` | Same as TSV, Excel format (requires openpyxl). |
| `failed_targets.tsv` | Genes that could not be processed (only if failures exist). |

### Stage 2 — Primer Design

| File | Description |
| --- | --- |
| `primers.tsv` | Primer3 primer pairs and design metrics. No genomic product coordinates. |
| `primers.xlsx` | Same as TSV with Targets and Primers sheets (requires openpyxl). |

### Stage 3 — Specificity Check

| File | Description |
| --- | --- |
| `primers_specificity.tsv` | Primer records plus in-silico PCR hit classification and genomic coordinates. |
| `primers_unique.tsv` | Primer pairs classified as `unique_pass`. |
| `stage3_summary.txt` | Per-target specificity summary. |

### Pipeline Summary

| File | Description |
| --- | --- |
| `run_summary.txt` | QC summary for the full pipeline run. |

## primer-panel-finalize

The `primer-panel-finalize` command selects the best unique primer per target
from Stage 3 outputs and writes a clean panel recommendation.

**Default behavior** (no rescue):

```bash
primer-panel-finalize \
  --input-dir outputs/hcc6_primers \
  --output-dir outputs/panel_final \
  --genome-fasta /path/to/hg38.fa
```

This produces:

| File | Description |
| --- | --- |
| `recommended_primers.tsv` | Best unique_pass primer per target. |
| `recommended_primers.xlsx` | Same, Excel format (requires openpyxl). |
| `failed_or_needs_review_targets.tsv` | Targets without any unique_pass primer. |
| `rescue_attempts.tsv` | Rescue results (empty if no rescue requested). |
| `panel_summary.txt` | Human-readable summary. |

### Rescue (experimental)

Rescue attempts re-run Primer3 with different parameters for targets that have
no unique primer. This is an experimental feature and must be explicitly
requested:

```bash
# Rescue a specific target
primer-panel-finalize \
  --input-dir outputs/hcc6_primers \
  --output-dir outputs/panel_final \
  --genome-fasta /path/to/hg38.fa \
  --rescue-target FTH1_cds1_4

# Rescue all supported targets
primer-panel-finalize \
  --input-dir outputs/hcc6_primers \
  --output-dir outputs/panel_final \
  --genome-fasta /path/to/hg38.fa \
  --rescue-all

# Enable experimental split-target rescue
primer-panel-finalize \
  --input-dir outputs/hcc6_primers \
  --output-dir outputs/panel_final \
  --genome-fasta /path/to/hg38.fa \
  --rescue-target FTH1_cds1_4 \
  --experimental-split-rescue
```

Currently only `FTH1_cds1_4` has an implemented rescue strategy.

## Specificity Status

Stage 3 reports one of these statuses for each primer pair:

| Status | Meaning |
| --- | --- |
| `unique_pass` | One hit, matching the expected chromosome and coordinates within tolerance. |
| `unique_off_target` | One hit, but not at the expected location. |
| `multi_hit` | Multiple genome-wide products. |
| `no_hit` | No product detected. |
| `pcr_error` | `isPcr` failed for the batch. |

## Dependency Doctor

Check your environment for required tools and files:

```bash
primer-panel --doctor
```

This reports:

- `primer3_core` — found / missing
- `isPcr` — found / missing
- `faToTwoBit` — found / missing (optional)
- genome FASTA — provided and exists / not provided
- common dbSNP BED — provided and exists / not provided

The `--doctor` flag does **not** require `--genes` and does not run any
pipeline stage.

## Offline / Cached Annotation with Ensembl GTF

By default, Stage 1 queries the Ensembl REST API for transcript and CDS
annotation.  To avoid network access (or to pin a specific annotation
release), you can provide a local Ensembl GTF file:

```bash
primer-panel \
  --genes HFE HJV TFR2 \
  --stage targets \
  --annotation-gtf /path/to/Homo_sapiens.GRCh38.110.gtf.gz \
  --output-dir outputs/offline_run
```

The GTF must be an Ensembl-format GTF (plain or gzip-compressed).  Both
`gene_name` and `gene_id` lookups are supported.

### Combining GTF with genome FASTA

GTF provides annotation (coordinates), not sequence.  For the full pipeline
you still need `--genome-fasta`:

```bash
primer-panel \
  --genes HFE HJV TFR2 \
  --annotation-gtf /path/to/Homo_sapiens.GRCh38.110.gtf.gz \
  --genome-fasta /path/to/hg38.fa \
  --output-dir outputs/full_offline
```

### Assembly / release consistency

Use GTF and FASTA from the same Ensembl release (e.g. both GRCh38.110).
Mixing releases may cause coordinate mismatches.

### Transcript selection priority

When using a GTF, transcript selection follows the same priority as the
Ensembl API mode:

1. `MANE_Select` tag
2. `Ensembl_canonical` tag
3. Longest CDS among `protein_coding` transcripts

## Development

Run the test suite:

```bash
python -m pytest
```

Run basic packaging and syntax checks:

```bash
python -m compileall -q primer_panel tests
python -m build --no-isolation --sdist --wheel --outdir /tmp/primer-panel-build
```

## License

MIT. See [LICENSE](LICENSE).
