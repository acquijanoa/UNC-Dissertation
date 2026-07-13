#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5g
#SBATCH -t 00:10:00
#SBATCH -o logs/job_%j.out
#SBATCH -e logs/job_%j.err

mkdir -p logs

module purge
module load r/4.5.0

R CMD BATCH --no-save test_compile.R
