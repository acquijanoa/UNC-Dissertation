#!/usr/bin/env python3
"""
merge_slurm_results.py
======================
Combines individual result CSVs from the SLURM array job on Longleaf
into a single final results file, and prints the summary table.

Usage:
    python3 scripts/merge_slurm_results.py
"""

import sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = ROOT / "data" / "sim_results"
FINAL_OUT = ROOT / "data" / "simulation_weight_conditions.csv"

def main():
    if not RESULTS_DIR.exists():
        print(f"Error: Results directory {RESULTS_DIR} does not exist.")
        sys.exit(1)

    result_files = sorted(RESULTS_DIR.glob("result_*.csv"))
    if not result_files:
        print(f"No result CSVs found in {RESULTS_DIR}")
        sys.exit(1)

    print(f"Found {len(result_files)} result files. Merging...")
    
    dfs = []
    for f in result_files:
        try:
            df = pd.read_csv(f)
            dfs.append(df)
        except Exception as e:
            print(f"Warning: Could not read {f.name}: {e}")

    if not dfs:
        print("No valid dataframes to merge.")
        sys.exit(1)

    df_all = pd.concat(dfs, ignore_index=True)
    df_all.to_csv(FINAL_OUT, index=False)
    print(f"Successfully saved merged results to: {FINAL_OUT}")

    # Print summary table
    print("\n" + "=" * 78)
    print("SIMULATION SUMMARY")
    print("=" * 78)
    summary = df_all.groupby("condition").agg(
        pp_cov_mean   =("pp_coverage", "mean"),
        pp_cov_std    =("pp_coverage", "std"),
        pp_width_mean =("pp_width",    "mean"),
        f_sd_mean     =("mcmc_sd_f",   "mean"),
        f_sd_std      =("mcmc_sd_f",   "std"),
        f_ci_mean     =("f_ci_width",  "mean"),
        rmse_mean     =("rmse",        "mean"),
        deff_mean     =("deff",        "mean"),
    ).round(4)
    print(summary.to_string())

    sd_A = df_all[df_all.condition == "A"]["mcmc_sd_f"].mean()
    sd_B = df_all[df_all.condition == "B"]["mcmc_sd_f"].mean()
    sd_C = df_all[df_all.condition == "C"]["mcmc_sd_f"].mean()
    deff_mean = df_all["deff"].mean()
    print(f"\nInterpretation:")
    print(f"  MCMC-SD ratio  A/B = {sd_A/sd_B:.3f}  "
          f"(theory: sqrt(1+1/Deff) = {(1+1/deff_mean)**0.5:.3f})")
    print(f"  MCMC-SD ratio  A/C = {sd_A/sd_C:.3f}")

if __name__ == "__main__":
    main()
