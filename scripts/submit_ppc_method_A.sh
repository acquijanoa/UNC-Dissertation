#!/bin/bash
# Submit Method A OOB PPC: calibrate K to ~95% PP coverage of y, then check
set -e
PROJ_ROOT="$(pwd)"
cd "$PROJ_ROOT" || exit 1
mkdir -p logs data/friedman_ppc_A

echo "=== Method A OOB Posterior Predictive Check (target 95% cover of y) ==="
JOB=$(sbatch --parsable --array=1-1000 \
  --export=ALL,MODE=calibrate,RESULTS_DIR=data/friedman_ppc_A \
  scripts/slurm/submit_ppc_method_A.slurm)
echo "Submitted job: $JOB"
echo "Monitor: squeue -u \$USER"
echo "Results: data/friedman_ppc_A/"
