#!/bin/bash
#SBATCH --job-name=report_oob_K
#SBATCH --output=logs/report_oob_K.out
#SBATCH --error=logs/report_oob_K.err
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=general

PROJ_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$PROJ_ROOT" || exit 1

module load r/4.3.1

RESULTS_DIR="${RESULTS_DIR:-data/friedman_results_oob}"
REPORT_CSV="${RESULTS_DIR}/report_table.csv"
REPORT_MD="${RESULTS_DIR}/report.md"

echo "Generating OOB K-calibration report from ${RESULTS_DIR}"
echo "CSV output: ${REPORT_CSV}"
echo "MD output: ${REPORT_MD}"

Rscript scripts/report_friedman_oob_K.R \
  --results_dir "${RESULTS_DIR}" \
  --out_csv "${REPORT_CSV}" \
  --out_md "${REPORT_MD}"

echo "Report generation complete."
cat "${REPORT_MD}"

exit 0
