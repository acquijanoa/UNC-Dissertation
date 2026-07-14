#!/bin/bash
#SBATCH --job-name=report_friedman
#SBATCH --output=logs/report_friedman_%j.out
#SBATCH --error=logs/report_friedman_%j.err
#SBATCH --time=00:20:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=general

# Merge result CSVs and write MC report (bias, rbias, empSE, estSE, coverage)
#   sbatch scripts/slurm/submit_report_friedman.slurm

set -euo pipefail

PROJ_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJ_ROOT}"
mkdir -p logs data/friedman_results

module purge
module load r/4.5.0

echo "=== Merge + Friedman report | $(date) ==="

Rscript scripts/merge_friedman_results.R \
  data/friedman_results \
  data/friedman_results/summary.csv

Rscript scripts/report_friedman_recovery.R \
  data/friedman_results \
  data/friedman_results/report_table.csv \
  data/friedman_results/report.md

echo "=== Done | $(date) ==="
echo "Report: data/friedman_results/report.md"
