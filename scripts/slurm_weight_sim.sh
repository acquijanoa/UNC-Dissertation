#!/bin/bash
#SBATCH --job-name=wt_sim
#SBATCH --output=logs/wt_sim_%A_%a.out
#SBATCH --error=logs/wt_sim_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-1000
#SBATCH --partition=general

# ─── Weight Condition Simulation — Longleaf Array Job ──────────────────────
# Each array task handles one replication (1 sample × 3 conditions).
# Results are written as separate CSV files and merged afterward (see scripts/archive/).
#
# Submit from the project root on Longleaf:
#   cd /proj/sta790/alvaro/dissertation   (adjust path)
#   mkdir -p logs
#   sbatch scripts/slurm_weight_sim.sh
# ───────────────────────────────────────────────────────────────────────────

set -euo pipefail

# ── R environment ───────────────────────────────────────────────
module purge
module load r/4.5.0

# ── Paths ────────────────────────────────────────────────────────────────────
PROJ_ROOT="$(pwd)"
SCRIPT="${PROJ_ROOT}/scripts/simulation_weight_conditions.R"
RESULTS_DIR="${PROJ_ROOT}/data/sim_results"
mkdir -p "${RESULTS_DIR}"

# ── Per-task arguments ───────────────────────────────────────────────────────
# Default to Condition A if not set via --export
COND="${COND:-A}"

# SLURM_ARRAY_TASK_ID is 1-indexed (1..1000) → maps to sample file index
SAMPLE_IDX=$(printf "%04d" "${SLURM_ARRAY_TASK_ID}")
SAMPLE_FILE="${PROJ_ROOT}/data/samples/sample_${SAMPLE_IDX}.csv"
OUT_FILE="${RESULTS_DIR}/result_${SAMPLE_IDX}_${COND}.csv"

echo "=== Task ${SLURM_ARRAY_TASK_ID} | Sample ${SAMPLE_IDX} | Condition ${COND} | $(date) ==="

Rscript "${SCRIPT}" \
    --single_sample "${SAMPLE_FILE}" \
    --out_file      "${OUT_FILE}"    \
    --condition     "${COND}"        \
    --n_mcmc        2000             \
    --burn_in       500              \
    --n_trees       200              \
    --seed          "${SLURM_ARRAY_TASK_ID}"

echo "=== Done | $(date) ==="
