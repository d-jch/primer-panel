# Changelog

## 0.4.0 (2026-06-24)

### CLI
- Add `primer-panel` console entry point (install with `pip install -e .`)
- Add `primer-panel-finalize` console entry point for panel finalization
- Add `__main__.py` guard for `python -m primer_panel`

### Stage 1 — Target Generation
- Integrate `primer-target-planner` as default grouper (replaces internal v1/v2)
- Strand-aware grouping: +strand sweeps right, -strand sweeps left
- Remove `cds_buffer` (raw CDS exon coordinates used directly)
- Pass gene-level bounds from Ensembl for correct `terminal_reverse` logic
- Add `SEQUENCE_TARGET` fields (0-based and 1-based) to target output
- Add `target_qc_status` for Stage 1 quality checks

### Stage 2 — Primer Design
- Use `primer3plus-core` package for Primer3Plus-like defaults
- `PRIMER_TASK=generic`, `PRIMER_OUTSIDE_PENALTY=0`, `PRIMER_INSIDE_PENALTY=-1.0`
- Primer3Plus default product size ranges (501-600 ... 10001-20000)
- Convert 0-based `SEQUENCE_TARGET` to 1-based for Primer3 input
- Convert Primer3 1-based output coordinates to 0-based internally
- Remove `product_min`/`product_max` hard filtering of primer products
- Stage 2 output: only Primer3 design fields, no genomic product coordinates

### Stage 3 — Specificity
- Remove `product_min`/`product_max` filtering from isPcr hits
- All isPcr hits retained; classification based on chrom/start/end matching
- Add `specificity_explain` to all specificity results
- Stage 3 is the only stage producing genomic product coordinates

### Configuration
- `--target-size` replaces `--product-size` (backward compatible)
- `--primer3plus-settings` for custom Primer3 settings file
- `--write-primer3-inputs` for debugging
- Primer params default to `None` (only override when explicitly specified)

### Dependencies
- Add `primer3plus_core` package
- Add `primer-target-planner` package
