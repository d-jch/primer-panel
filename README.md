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

## Key Outputs

| File | Created by | Description |
| --- | --- | --- |
| `targets.bed` | Stage 1 | Design-template coordinates in BED6 format. |
| `required_regions.bed` | Stage 1 | Raw CDS regions that must be covered. |
| `target_summary.tsv` | Stage 1 | Target coordinates, selected transcript, CDS coverage, and QC status. |
| `targets.fa` | Stage 1 | Design-template FASTA. Uses placeholder `N` sequence unless `--genome-fasta` is provided. |
| `primers.tsv` | Stage 2 | Primer3 primer pairs and design metrics. No genomic product coordinates. |
| `primers_specificity.tsv` | Stage 3 | Primer records plus in-silico PCR hit classification and genomic coordinates. |
| `primers_unique.tsv` | Stage 3 | Primer pairs classified as `unique_pass`. |
| `run_summary.txt` | Pipeline | Short QC summary for the run. |

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
