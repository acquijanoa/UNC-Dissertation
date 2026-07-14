#!/bin/bash
#SBATCH --job-name=friedman_rec
#SBATCH --output=logs/friedman_%A_%a.out
#SBATCH --error=logs/friedman_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-1000
#SBATCH --partition=general

# ─── Friedman recovery (y_f ~ x_f*, held-out vs mu_f) ───────────────────────
# Methods: A (weight_norm Σw=n + Kim-Rao, no per-draw renorm),
#          B (raw weight + Kim-Rao, no Σw*=n),
#          sum1 (w/Σw so Σw=1 + Kim-Rao, no per-draw renorm),
#          chipman (BART::wbart, x_f only / no design IDs)
#
#   sbatch --export=ALL,METHOD=A scripts/slurm/submit_friedman_recovery.slurm
#   sbatch --export=ALL,METHOD=B scripts/slurm/submit_friedman_recovery.slurm
#   sbatch --export=ALL,METHOD=sum1 scripts/slurm/submit_friedman_recovery.slurm
#   sbatch --export=ALL,METHOD=chipman scripts/slurm/submit_friedman_recovery.slurm
# ───────────────────────────────────────────────────────────────────────────

set -euo pipefail

# SLURM copies the batch script to spool; use the directory where sbatch was run
PROJ_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJ_ROOT}"

module purge
module load r/4.5.0

METHOD="${METHOD:-A}"
K="${K:-1.0}"
RESULTS_DIR="${RESULTS_DIR:-${PROJ_ROOT}/data/friedman_results}"
mkdir -p logs "${RESULTS_DIR}"

# Ensure BART is available for the Chipman arm only
if [ "${METHOD}" = "chipman" ]; then
  Rscript -e 'if (!requireNamespace("BART", quietly=TRUE)) install.packages("BART", repos="https://cloud.r-project.org")'
fi

SAMPLE_IDX=$(printf "%04d" "${SLURM_ARRAY_TASK_ID}")
SAMPLE_FILE="${PROJ_ROOT}/data/samples/sample_${SAMPLE_IDX}.csv"
OUT_FILE="${RESULTS_DIR}/result_${SAMPLE_IDX}_${METHOD}.csv"

echo "=== Friedman recovery | method=${METHOD} | K=${K} | sample=${SAMPLE_IDX} | $(date) ==="
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
