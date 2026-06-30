#!/bin/bash
#SBATCH --job-name=wt_sim
#SBATCH --output=logs/wt_sim_%A_%a.out
#SBATCH --error=logs/wt_sim_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-50
#SBATCH --partition=general

# ─── Weight Condition Simulation — Longleaf Array Job ──────────────────────
# Each array task handles one replication (1 sample × 3 conditions).
# Results are written as separate CSV files and merged by merge_results.py.
#
# Submit from the project root on Longleaf:
#   cd /proj/sta790/alvaro/dissertation   (adjust path)
#   mkdir -p logs
#   sbatch scripts/slurm_weight_sim.sh
# ───────────────────────────────────────────────────────────────────────────

set -euo pipefail

# ── Conda / Python environment ───────────────────────────────────────────────
module purge
module load anaconda/2023.09

# Activate your conda environment (adjust name if different)
conda activate dissertation 2>/dev/null || \
    source activate dissertation 2>/dev/null || \
    echo "Using base conda environment"

# ── Paths ────────────────────────────────────────────────────────────────────
PROJ_ROOT="$(pwd)"
SCRIPT="${PROJ_ROOT}/scripts/simulation_weight_conditions.py"
RESULTS_DIR="${PROJ_ROOT}/data/sim_results"
mkdir -p "${RESULTS_DIR}"

# ── Per-task arguments ───────────────────────────────────────────────────────
# SLURM_ARRAY_TASK_ID is 1-indexed (1..50) → maps to sample file index
SAMPLE_IDX=$(printf "%04d" "${SLURM_ARRAY_TASK_ID}")
SAMPLE_FILE="${PROJ_ROOT}/data/samples/sample_${SAMPLE_IDX}.csv"
OUT_FILE="${RESULTS_DIR}/result_${SAMPLE_IDX}.csv"

echo "=== Task ${SLURM_ARRAY_TASK_ID} | Sample ${SAMPLE_IDX} | $(date) ==="

python3 "${SCRIPT}" \
    --single_sample "${SAMPLE_FILE}" \
    --out_file      "${OUT_FILE}"    \
    --n_mcmc        2000             \
    --burn_in       500              \
    --n_trees       200              \
    --seed          "${SLURM_ARRAY_TASK_ID}"

echo "=== Done | $(date) ==="
