#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH -t 00:10:00
#SBATCH -o logs/gen_pop_%j.out
#SBATCH -e logs/gen_pop_%j.err
#SBATCH -p general

# ─── Generate Population — Longleaf SLURM Submit Script ──────────────────────

set -euo pipefail

# Purge modules and load R
module purge
module load r/4.5.0

# Print status info
echo "=== Starting gen_pop.R | $(date) ==="

# Execute R script
# Note: R CMD BATCH generates gen_pop.Rout in the current directory or wherever specified
Rscript scripts/gen_pop.R

echo "=== Done | $(date) ==="
