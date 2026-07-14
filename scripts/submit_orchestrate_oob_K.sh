#!/bin/bash
# submit_orchestrate_oob_K.sh
# Submit OOB K-calibration jobs for methods A and B (100 samples each)

set -e

PROJ_ROOT="$(pwd)"
cd "$PROJ_ROOT" || exit 1

mkdir -p logs data/friedman_results_oob

echo "=== Submitting OOB K-calibration jobs ==="
echo "Project root: $PROJ_ROOT"
echo ""

# Method A: 1000 samples (full MC)
echo "Submitting Method A (OOB) ..."
JOB_A=$(sbatch --parsable --array=1-1000 \
  --export=ALL,METHOD=A,RESULTS_DIR=data/friedman_results_oob \
  scripts/slurm/submit_friedman_oob_K.slurm)
echo "  Job ID: $JOB_A"

# Method B: 1000 samples (full MC)
echo "Submitting Method B (OOB) ..."
JOB_B=$(sbatch --parsable --array=1-1000 \
  --export=ALL,METHOD=B,RESULTS_DIR=data/friedman_results_oob \
  scripts/slurm/submit_friedman_oob_K.slurm)
echo "  Job ID: $JOB_B"

echo ""
echo "=== Submitting report job (depends on A and B completion) ==="
JOB_REPORT=$(sbatch --parsable --dependency=afterok:${JOB_A}:${JOB_B} \
  scripts/slurm/submit_report_oob_K.slurm)
echo "  Report Job ID: $JOB_REPORT"

echo ""
echo "All jobs submitted. Monitor with: squeue -u $USER"
echo "  Method A (OOB): $JOB_A"
echo "  Method B (OOB): $JOB_B"
echo "  Report: $JOB_REPORT"
