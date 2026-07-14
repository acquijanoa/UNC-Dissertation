#!/bin/bash
#SBATCH --job-name=orchestrate_K
#SBATCH --output=logs/orchestrate_K_%j.out
#SBATCH --error=logs/orchestrate_K_%j.err
#SBATCH --time=00:05:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --partition=general

# ─── Orchestrate K-grid experiment ─────────────────────────────────────────
# Submit array jobs for multiple K values for methods A and B.
# Tests K ∈ {0.5, 1, 2, 5, 10, 20} to find optimal coverage.
# ───────────────────────────────────────────────────────────────────────────

set -euo pipefail

PROJ_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJ_ROOT}"

# Check that samples exist
if [ ! -d "data/samples" ]; then
  echo "ERROR: data/samples/ not found. Run gen_pop + samples first."
  exit 1
fi

# K values to test
K_VALUES=(0.5 1 2 5 10 20)

echo "=== Submitting K-grid experiment ==="
echo "Methods: A, B"
echo "K values: ${K_VALUES[*]}"
echo ""

# Method A with different K values
for K in "${K_VALUES[@]}"; do
  echo "Submitting METHOD=A K=${K}..."
  sbatch --export=ALL,METHOD=A,K=${K} scripts/slurm/submit_friedman_K_grid.slurm
done

# Method B with different K values
for K in "${K_VALUES[@]}"; do
  echo "Submitting METHOD=B K=${K}..."
  sbatch --export=ALL,METHOD=B,K=${K} scripts/slurm/submit_friedman_K_grid.slurm
done

echo ""
echo "=== All K-grid jobs submitted ==="
echo "Total: 12 array jobs (6 K values × 2 methods)"
echo "Each array job has 1000 replications"
echo ""
echo "Monitor with: squeue -u \$USER"
echo "Results will be in: data/friedman_K_grid/"
