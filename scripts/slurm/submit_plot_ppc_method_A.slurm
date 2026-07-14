#!/bin/bash
#SBATCH --job-name=plot_ppc_A
#SBATCH --output=logs/plot_ppc_method_A.out
#SBATCH --error=logs/plot_ppc_method_A.err
#SBATCH --time=01:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=general

PROJ_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$PROJ_ROOT" || exit 1

module load r/4.5.0
mkdir -p logs figures

Rscript scripts/plot_ppc_method_A.R \
  --single_sample data/samples/sample_0001.csv \
  --out_plot figures/ppc_method_A_sample0001.png \
  --mode calibrate \
  --K_grid 0.01,0.05,0.1,0.2,0.5,1,2 \
  --response y_f \
  --n_mcmc 2000 \
  --burn_in 500 \
  --n_trees 200 \
  --seed 42

ls -la figures/ppc_method_A_sample0001.*
exit 0
