#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16g
#SBATCH -t 05:00:00
#SBATCH --array=1-8
#SBATCH -o logs/estimate_causal_meps_%A_%a.out
#SBATCH -e logs/estimate_causal_meps_%A_%a.err
#SBATCH -p general

# ─── Estimate Causal Effect of USC on MEPS — Longleaf SLURM Submit Script ────

set -euo pipefail

# Purge modules and load R
module purge
module load r/4.5.0

if [ $SLURM_ARRAY_TASK_ID -le 4 ]; then
    OUTCOME="totslf21"
    CHAIN_ID=$SLURM_ARRAY_TASK_ID
else
    OUTCOME="totexp21"
    CHAIN_ID=$((SLURM_ARRAY_TASK_ID - 4))
fi

# Print status info
echo "=== Starting estimate_causal_meps.R ($OUTCOME, Chain $CHAIN_ID) | $(date) ==="
echo "Node: $SLURMD_NODENAME"

# Execute R script
Rscript scripts/estimate_causal_meps.R \
    --n_mcmc 1000 \
    --burn_in 500 \
    --n_trees 200 \
    --outcome "$OUTCOME" \
    --chain_id "$CHAIN_ID"

echo "=== Done | $(date) ==="
