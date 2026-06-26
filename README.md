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

## Install

Create an environment with the required Python and bioinformatics tools:

```bash
micromamba create -n primer_panel -c conda-forge -c bioconda \
  python=3.11 requests openpyxl pyfaidx primer3 -y

micromamba activate primer_panel
pip install -e .
```

Check the CLI entry points:

```bash
primer-panel --help
primer-panel-finalize --help
primer3_core --version
```

Optional tools:

- `pyfaidx`: required when extracting real sequence from a genome FASTA.
- `openpyxl`: enables XLSX output.
- UCSC `isPcr`: required for Stage 3 specificity checks.
- UCSC `faToTwoBit`: only needed if you want to create `.2bit` databases.

## Quick Start

Generate Stage 1 target files without designing primers:

```bash
primer-panel \
  --genes HFE HJV TFR2 SLC40A1 HAMP FTH1 \
  --target-size 2700-3300 \
  --primer-flank 500 \
  --output-dir outputs/hcc6_targets
```

Run the full pipeline with primer design and specificity checks:

```bash
primer-panel \
  --genes HFE HJV TFR2 SLC40A1 HAMP FTH1 \
  --target-size 2700-3300 \
  --primer-flank 500 \
  --genome-fasta /path/to/hg38.fa \
  --design-primers \
  --check-specificity \
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
| `--target-size MIN-MAX` | `2700-3300` | Stage 1 CDS grouping and target extension size. This does not constrain Primer3 product size. |
| `--primer-flank N` | `300` | Bases added on both sides of each target for primer search. |
| `--genome-fasta PATH` | none | hg38 FASTA used for real sequence extraction and specificity checks. |
| `--design-primers` | off | Enable Primer3 primer design. |
| `--check-specificity` | off | Enable UCSC `isPcr` specificity checks. Requires `--design-primers`. |
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
