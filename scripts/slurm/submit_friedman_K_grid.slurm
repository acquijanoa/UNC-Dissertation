#!/bin/bash
#SBATCH --job-name=friedman_K
#SBATCH --output=logs/friedman_K_%A_%a.out
#SBATCH --error=logs/friedman_K_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-1000
#SBATCH --partition=general

# ─── Friedman recovery K-grid experiment ────────────────────────────────────
# Method A: weight_norm + Kim-Rao, NO per-draw Σw*=n (K can inflate Term 2)
# Method B: raw weight + Kim-Rao, no Σw*=n
#
# Usage:
#   sbatch --export=ALL,METHOD=A,K=1.0 scripts/slurm/submit_friedman_K_grid.slurm
#   sbatch --export=ALL,METHOD=A,K=2.0 scripts/slurm/submit_friedman_K_grid.slurm
#   sbatch --export=ALL,METHOD=B,K=10.0 scripts/slurm/submit_friedman_K_grid.slurm
# ───────────────────────────────────────────────────────────────────────────

set -euo pipefail

# SLURM copies the batch script to spool; use the directory where sbatch was run
PROJ_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJ_ROOT}"

module purge
module load r/4.5.0

METHOD="${METHOD:-A}"
K="${K:-1.0}"
# Default to A2 dir so new Method A (no per-draw renorm) does not overwrite
# the earlier K-grid where A used normalize_weights=TRUE.
RESULTS_DIR="${RESULTS_DIR:-${PROJ_ROOT}/data/friedman_K_grid_A2}"
mkdir -p logs "${RESULTS_DIR}"

SAMPLE_IDX=$(printf "%04d" "${SLURM_ARRAY_TASK_ID}")
SAMPLE_FILE="${PROJ_ROOT}/data/samples/sample_${SAMPLE_IDX}.csv"

# Format K for filename: replace . with p (e.g., 1.0 -> 1p0, 0.5 -> 0p5)
K_STR=$(echo "${K}" | sed 's/\./_/g')
OUT_FILE="${RESULTS_DIR}/result_${SAMPLE_IDX}_${METHOD}_K${K_STR}.csv"

echo "=== Friedman K-grid | method=${METHOD} | K=${K} | sample=${SAMPLE_IDX} | $(date) ==="
echo "Project root: ${PROJ_ROOT}"

Rscript scripts/simulation_friedman_recovery.R \
    --single_sample "${SAMPLE_FILE}" \
    --out_file      "${OUT_FILE}" \
    --method        "${METHOD}" \
    --K             "${K}" \
    --n_mcmc        2000 \
    --burn_in       500 \
    --n_trees       200 \
    --n_test        2000 \
    --seed          "${SLURM_ARRAY_TASK_ID}"

echo "=== Done | $(date) ==="
