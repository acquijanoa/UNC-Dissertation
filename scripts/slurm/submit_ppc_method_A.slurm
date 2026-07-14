#!/bin/bash
#SBATCH --job-name=ppc_method_A
#SBATCH --output=logs/ppc_method_A_%a.out
#SBATCH --error=logs/ppc_method_A_%a.err
#SBATCH --time=02:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=general

# Method A OOB posterior predictive check (95% coverage of y)
# Usage:
#   sbatch --array=1-1000 --export=MODE=calibrate scripts/slurm/submit_ppc_method_A.slurm
#   sbatch --array=1-1000 --export=MODE=check,K_VAL=1.0 scripts/slurm/submit_ppc_method_A.slurm

PROJ_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$PROJ_ROOT" || exit 1

module load r/4.5.0

MODE="${MODE:-calibrate}"
K_VAL="${K_VAL:-1.0}"
RESULTS_DIR="${RESULTS_DIR:-data/friedman_ppc_A}"
SAMPLE_ID=$(printf "%04d" "${SLURM_ARRAY_TASK_ID}")
SAMPLE_FILE="data/samples/sample_${SAMPLE_ID}.csv"
OUT_FILE="${RESULTS_DIR}/ppc_${SAMPLE_ID}_${MODE}.csv"

mkdir -p logs "${RESULTS_DIR}"

echo "Method A OOB PPC | mode=${MODE} | sample=${SAMPLE_ID} | K=${K_VAL}"
echo "Output: ${OUT_FILE}"

Rscript scripts/ppc_method_A.R \
  --single_sample "${SAMPLE_FILE}" \
  --out_file "${OUT_FILE}" \
  --mode "${MODE}" \
  --K "${K_VAL}" \
  --K_grid "0.01,0.05,0.1,0.2,0.5,1,2" \
  --target 0.95 \
  --level 0.95 \
  --response y_f \
  --n_mcmc 2000 \
  --burn_in 500 \
  --n_trees 200 \
  --seed $((42 + SLURM_ARRAY_TASK_ID))

exit 0
