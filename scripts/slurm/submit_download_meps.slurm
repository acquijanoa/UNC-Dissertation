#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH -t 00:30:00
#SBATCH -o logs/download_meps_%j.out
#SBATCH -e logs/download_meps_%j.err
#SBATCH -p general

# ─── Download MEPS Data — Longleaf SLURM Submit Script ──────────────────────
# Submit from anywhere:
#   sbatch scripts/slurm/submit_download_meps.slurm

set -euo pipefail

# SLURM copies the batch script to spool; use the directory where sbatch was run
PROJ_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJ_ROOT}"
mkdir -p logs

module purge
module load r/4.4.0

echo "=== Starting download_meps.R | $(date) ==="
echo "Project root: ${PROJ_ROOT}"

Rscript scripts/download_meps.R

echo "=== Done | $(date) ==="
