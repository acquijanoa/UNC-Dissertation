#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH -t 00:30:00
#SBATCH -o logs/download_meps_%j.out
#SBATCH -e logs/download_meps_%j.err
#SBATCH -p general

# ─── Download MEPS Data — Longleaf SLURM Submit Script ──────────────────────

set -euo pipefail

# Purge modules and load R 4.4.0
module purge
module load r/4.4.0

# Print status info
echo "=== Starting download_meps.R | $(date) ==="

# Execute R script
Rscript scripts/download_meps.R

echo "=== Done | $(date) ==="
