#!/bin/bash
# submit_K_grid_with_report.sh
# Submit K-grid jobs and automatically generate report when they complete

set -euo pipefail

echo "=== Submitting K-grid experiment with auto-report ==="

# Submit all K-grid jobs and collect their job IDs
K_VALUES=(0.5 1 2 5 10 20)
JOB_IDS=()

echo "Submitting Method A jobs..."
for K in "${K_VALUES[@]}"; do
  JOB_ID=$(sbatch --parsable --export=ALL,METHOD=A,K=${K} scripts/slurm/submit_friedman_K_grid.slurm)
  JOB_IDS+=("${JOB_ID}")
  echo "  K=${K}: ${JOB_ID}"
done

echo "Submitting Method B jobs..."
for K in "${K_VALUES[@]}"; do
  JOB_ID=$(sbatch --parsable --export=ALL,METHOD=B,K=${K} scripts/slurm/submit_friedman_K_grid.slurm)
  JOB_IDS+=("${JOB_ID}")
  echo "  K=${K}: ${JOB_ID}"
done

# Build dependency string
DEPENDENCY=$(IFS=:; echo "${JOB_IDS[*]}")

echo ""
echo "Submitting report job with dependency on all K-grid jobs..."
REPORT_JOB=$(sbatch --parsable --dependency=afterok:${DEPENDENCY} scripts/slurm/submit_report_K_grid.slurm)
echo "Report job: ${REPORT_JOB}"

echo ""
echo "=== Setup complete ==="
echo "K-grid jobs: ${#JOB_IDS[@]} array jobs"
echo "Report will run automatically after all jobs complete"
echo ""
echo "Monitor with: squeue -u \$USER"
