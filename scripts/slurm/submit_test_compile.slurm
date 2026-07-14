#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5g
#SBATCH -t 00:10:00
#SBATCH -o logs/job_%j.out
#SBATCH -e logs/job_%j.err
#SBATCH -p general

# ─── Test Compile — Longleaf SLURM Submit Script ────────────────────────────
# Submit from anywhere:
#   sbatch scripts/slurm/submit_test_compile.slurm

set -euo pipefail

# SLURM copies the batch script to spool; use the directory where sbatch was run
PROJ_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJ_ROOT}"
mkdir -p logs

module purge
module load r/4.5.0

echo "=== Starting test_compile.R | $(date) ==="
echo "Project root: ${PROJ_ROOT}"

Rscript scripts/test_compile.R

echo "=== Done | $(date) ==="
