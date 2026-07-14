#!/bin/bash
#SBATCH --job-name=merge_friedman
#SBATCH --output=logs/merge_friedman_%j.out
#SBATCH --error=logs/merge_friedman_%j.err
#SBATCH --time=00:15:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=general

# Merge Friedman recovery CSVs into data/friedman_results/summary.csv
#   sbatch scripts/slurm/submit_merge_friedman.slurm

set -euo pipefail

PROJ_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJ_ROOT}"
mkdir -p logs

module purge
module load r/4.5.0

echo "=== Merging Friedman results | $(date) ==="
echo "Project root: ${PROJ_ROOT}"

Rscript scripts/merge_friedman_results.R \
  data/friedman_results \
  data/friedman_results/summary.csv

echo "=== Done | $(date) ==="
