#!/bin/bash
#SBATCH --job-name=pipeline_friedman
#SBATCH --output=logs/pipeline_friedman_%j.out
#SBATCH --error=logs/pipeline_friedman_%j.err
#SBATCH --time=00:05:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --partition=general

# Orchestrator: after gen_pop completes, submit samples then Friedman + report.
# Usage: started automatically with --dependency after gen_pop, OR:
#   POP_JOB=$(sbatch --parsable scripts/slurm/submit_gen_pop.slurm)
#   sbatch --dependency=afterok:$POP_JOB scripts/slurm/submit_pipeline_friedman_clean.slurm

set -euo pipefail
PROJ_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJ_ROOT}"
mkdir -p logs

echo "=== Submitting sample generation + Friedman recovery | $(date) ==="

SAMP_JOB=$(sbatch --parsable scripts/slurm/submit_generate_samples.slurm)
echo "samples: ${SAMP_JOB}"

JOB_A=$(sbatch --parsable --dependency=afterok:${SAMP_JOB} --export=ALL,METHOD=A scripts/slurm/submit_friedman_recovery.slurm)
JOB_B=$(sbatch --parsable --dependency=afterok:${SAMP_JOB} --export=ALL,METHOD=B scripts/slurm/submit_friedman_recovery.slurm)
JOB_C=$(sbatch --parsable --dependency=afterok:${SAMP_JOB} --export=ALL,METHOD=chipman scripts/slurm/submit_friedman_recovery.slurm)
JOB_R=$(sbatch --parsable --dependency=afterok:${JOB_A}:${JOB_B}:${JOB_C} scripts/slurm/submit_report_friedman.slurm)

echo "A=${JOB_A} B=${JOB_B} chipman=${JOB_C} report=${JOB_R}"
echo "=== Orchestrator done | $(date) ==="
