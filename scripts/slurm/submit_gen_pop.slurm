#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH -t 00:10:00
#SBATCH -o logs/gen_pop_%j.out
#SBATCH -e logs/gen_pop_%j.err
#SBATCH -p general

# ─── Generate Population — Longleaf SLURM Submit Script ──────────────────────
# Submit from anywhere:
#   sbatch scripts/slurm/submit_gen_pop.slurm

set -euo pipefail

# SLURM copies the batch script to spool; use the directory where sbatch was run
PROJ_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJ_ROOT}"
mkdir -p logs

module purge
module load r/4.5.0

echo "=== Starting gen_pop.R | $(date) ==="
echo "Project root: ${PROJ_ROOT}"

Rscript scripts/gen_pop.R

echo "=== Done | $(date) ==="
