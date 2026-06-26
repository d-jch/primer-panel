#!/usr/bin/env bash
# Smoke test: run HAMP (small gene, few exons) through the pipeline
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

echo "=== Primer Panel Smoke Test ==="
echo ""

# Clean previous run
rm -rf outputs/smoke_test

# Run pipeline with HAMP only
echo "Running: python3 -m primer_panel --genes HAMP --output-dir outputs/smoke_test"
python3 -m primer_panel --genes HAMP --output-dir outputs/smoke_test

echo ""
echo "=== Checking outputs ==="

PASS=true

# Check targets.bed
if [[ -f outputs/smoke_test/targets.bed ]]; then
    LINES=$(wc -l < outputs/smoke_test/targets.bed)
    echo "[OK] targets.bed exists ($LINES lines)"
else
    echo "[FAIL] targets.bed missing"
    PASS=false
fi

# Check target_summary.tsv
if [[ -f outputs/smoke_test/target_summary.tsv ]]; then
    LINES=$(wc -l < outputs/smoke_test/target_summary.tsv)
    echo "[OK] target_summary.tsv exists ($LINES lines)"
else
    echo "[FAIL] target_summary.tsv missing"
    PASS=false
fi

# Check targets.fa
if [[ -f outputs/smoke_test/targets.fa ]]; then
    echo "[OK] targets.fa exists"
else
    echo "[FAIL] targets.fa missing"
    PASS=false
fi

echo ""
if $PASS; then
    echo "=== SMOKE TEST PASSED ==="
else
    echo "=== SMOKE TEST FAILED ==="
    exit 1
fi
