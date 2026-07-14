#!/bin/bash
#SBATCH --job-name=report_A2_K
#SBATCH --output=logs/report_A2_K_%j.out
#SBATCH --error=logs/report_A2_K_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=general

# Report for Method A K-grid with weight_norm + no per-draw renorm

set -euo pipefail

PROJ_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJ_ROOT}"

module purge
module load r/4.5.0

RESULTS_DIR="${PROJ_ROOT}/data/friedman_K_grid_A2"

echo "=== Generating Method A2 K-grid report | $(date) ==="
echo "Results directory: ${RESULTS_DIR}"

NUM_FILES=$(find "${RESULTS_DIR}" -name "result_*.csv" | wc -l)
echo "Found ${NUM_FILES} result files"

Rscript scripts/report_friedman_recovery.R \
  "${RESULTS_DIR}" \
  "${RESULTS_DIR}/report_table.csv" \
  "${RESULTS_DIR}/report.md"

echo "=== Done | $(date) ==="
cat "${RESULTS_DIR}/report.md"
