#!/bin/bash
#SBATCH --job-name=friedman_oob_K
#SBATCH --output=logs/friedman_oob_K_%a.out
#SBATCH --error=logs/friedman_oob_K_%a.err
#SBATCH --time=02:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=general

# Submit Friedman OOB K-calibration for a single method
# Usage: sbatch --array=1-100 --export=METHOD=A scripts/slurm/submit_friedman_oob_K.slurm

PROJ_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$PROJ_ROOT" || exit 1

module load r/4.3.1

METHOD="${METHOD:-A}"
SAMPLE_ID=$(printf "%04d" ${SLURM_ARRAY_TASK_ID})
SAMPLE_FILE="data/samples/sample_${SAMPLE_ID}.csv"
RESULTS_DIR="${RESULTS_DIR:-data/friedman_results_oob}"
OUT_FILE="${RESULTS_DIR}/result_${SAMPLE_ID}_${METHOD}_OOB.csv"

# OOB K-grid: focus on region around optimal from grid search
# For method A: optimal was K≈0.1, so search 0.01-1
# For method B: original K=1 was good, so search 0.5-5
if [[ "$METHOD" == "A" ]]; then
  K_GRID="0.01,0.05,0.1,0.2,0.5,1.0"
elif [[ "$METHOD" == "B" ]]; then
  K_GRID="0.5,1.0,2.0,5.0,10.0"
elif [[ "$METHOD" == "sum1" ]]; then
  K_GRID="0.1,0.5,1.0,2.0,5.0"
else
  echo "Unknown METHOD: $METHOD"
  exit 1
fi

echo "Running OOB K-calibration for sample ${SAMPLE_ID}, method ${METHOD}"
echo "K grid: ${K_GRID}"
echo "Output: ${OUT_FILE}"

Rscript scripts/simulation_friedman_oob_K.R \
  --single_sample "${SAMPLE_FILE}" \
  --out_file "${OUT_FILE}" \
  --method "${METHOD}" \
  --K_grid "${K_GRID}" \
  --calib_trees 50 \
  --calib_iter 500 \
  --calib_burn 150 \
  --n_test 2000 \
  --n_mcmc 2000 \
  --burn_in 500 \
  --n_trees 200 \
  --seed $((42 + SLURM_ARRAY_TASK_ID))

exit 0
