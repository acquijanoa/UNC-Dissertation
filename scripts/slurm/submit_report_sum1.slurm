#!/bin/bash
#SBATCH --job-name=report_sum1
#SBATCH --output=logs/report_sum1_%j.out
#SBATCH --error=logs/report_sum1_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=general

set -euo pipefail
PROJ_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJ_ROOT}"
module purge
module load r/4.5.0

RESULTS_DIR="${PROJ_ROOT}/data/friedman_results_sum1"
echo "=== sum1 report | files=$(find "${RESULTS_DIR}" -name 'result_*.csv' | wc -l) | $(date) ==="
Rscript scripts/report_friedman_recovery.R \
  "${RESULTS_DIR}" \
  "${RESULTS_DIR}/report_table.csv" \
  "${RESULTS_DIR}/report.md"
cat "${RESULTS_DIR}/report.md"
