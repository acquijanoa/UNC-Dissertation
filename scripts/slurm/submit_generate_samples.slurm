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
# Submit from anywhere:
#   sbatch scripts/slurm/submit_generate_samples.slurm

set -euo pipefail

# SLURM copies the batch script to spool; use the directory where sbatch was run
PROJ_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJ_ROOT}"
mkdir -p logs data/samples

module purge
module load r/4.5.0

echo "=== Starting sample generation task ${SLURM_ARRAY_TASK_ID} | $(date) ==="
echo "Project root: ${PROJ_ROOT}"

Rscript scripts/sample.R "${SLURM_ARRAY_TASK_ID}"

echo "=== Done | $(date) ==="
