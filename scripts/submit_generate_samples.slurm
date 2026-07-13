#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2g
#SBATCH -t 00:05:00
#SBATCH -o logs/sample_%A_%a.out
#SBATCH -e logs/sample_%A_%a.err
#SBATCH -p general
#SBATCH --array=1-1000

# ─── Generate Samples — Longleaf SLURM Submit Script ────────────────────────

set -euo pipefail

# Purge modules and load R
module purge
module load r/4.5.0

# Print status info
echo "=== Starting sample generation task ${SLURM_ARRAY_TASK_ID} | $(date) ==="

# Execute sample R script with the task ID as sample index
Rscript scripts/sample.R "${SLURM_ARRAY_TASK_ID}"

echo "=== Done | $(date) ==="
