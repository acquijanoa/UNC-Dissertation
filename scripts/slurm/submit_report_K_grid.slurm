#!/bin/bash
#SBATCH --job-name=report_K_grid
#SBATCH --output=logs/report_K_grid_%j.out
#SBATCH --error=logs/report_K_grid_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=general

# ─── Generate K-grid coverage report ────────────────────────────────────────
# Merge all K-grid results and generate a report showing coverage by K value.
#
# Usage (after all K-grid jobs complete):
#   sbatch scripts/slurm/submit_report_K_grid.slurm
# ───────────────────────────────────────────────────────────────────────────

set -euo pipefail

PROJ_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJ_ROOT}"

module purge
module load r/4.5.0

RESULTS_DIR="${RESULTS_DIR:-${PROJ_ROOT}/data/friedman_K_grid}"

echo "=== Generating K-grid report | $(date) ==="
echo "Results directory: ${RESULTS_DIR}"

if [ ! -d "${RESULTS_DIR}" ]; then
  echo "ERROR: ${RESULTS_DIR} not found"
  exit 1
fi

# Count result files
NUM_FILES=$(find "${RESULTS_DIR}" -name "result_*.csv" | wc -l)
echo "Found ${NUM_FILES} result files"

if [ "${NUM_FILES}" -lt 1000 ]; then
  echo "WARNING: Few result files found (${NUM_FILES})"
fi

# Generate report
echo ""
echo "Generating report..."
Rscript scripts/report_friedman_recovery.R \
  "${RESULTS_DIR}" \
  "${RESULTS_DIR}/report_table.csv" \
  "${RESULTS_DIR}/report.md"

echo ""
echo "=== Report complete | $(date) ==="
echo "Results:"
echo "  - CSV: ${RESULTS_DIR}/report_table.csv"
echo "  - MD:  ${RESULTS_DIR}/report.md"
