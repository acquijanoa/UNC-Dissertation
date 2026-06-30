#!/usr/bin/env python3
"""
Survey BART Imputation Simulation under MAR Missingness.

This script loads a sample, induces Missing at Random (MAR) nonresponse
on the outcome 'y' at specified levels (5%, 10%, 20%, and 50%),
fits Survey BART on the respondents, and evaluates the imputation accuracy
and population mean estimation.
"""

import argparse
import sys
import time
from pathlib import Path
import numpy as np
import pandas as pd
from numpy.random import default_rng

# Add the directory containing this script to sys.path
script_dir = Path(__file__).parent.resolve()
sys.path.append(str(script_dir))

try:
    from survey_bart import SurveyBART, BARTPrior
except ImportError:
    print("ERROR: Could not import SurveyBART. Ensure survey_bart.py is in the same directory.", file=sys.stderr)
    sys.exit(1)


def logit_func(x):
    return 1.0 / (1.0 + np.exp(-x))


def induce_mar_missingness(df, missing_rate, rng, seed=None):
    """
    Induce Missing at Random (MAR) missingness on the outcome variable 'y'.
    The probability of response depends on covariates x_cont_2 and x_bin.

    P(Response = 1) = logit^{-1}(gamma_0 + 0.5*x_cont_2 - 0.5*x_bin)

    We calibrate gamma_0 using binary search to match the target missing_rate.
    """
    target_response_rate = 1.0 - missing_rate
    x_c2 = df["x_cont_2"].values
    x_b = df["x_bin"].values

    # Calibrate gamma_0
    low, high = -15.0, 15.0
    for _ in range(100):
        mid = (low + high) / 2.0
        p_resp = logit_func(mid + 0.5 * x_c2 - 0.5 * x_b)
        mean_p = p_resp.mean()
        if mean_p < target_response_rate:
            low = mid
        else:
            high = mid
    gamma_0 = (low + high) / 2.0
    p_resp = logit_func(gamma_0 + 0.5 * x_c2 - 0.5 * x_b)

    # Sample response indicator R_i (1 = observed, 0 = missing)
    R = rng.binomial(1, p_resp).astype(bool)
    actual_missing_rate = 1.0 - R.mean()
    return R, p_resp, actual_missing_rate


def parse_args():
    parser = argparse.ArgumentParser(
        description="Survey BART Imputation Simulation under MAR"
    )
    parser.add_argument(
        "--sample", type=int, default=1,
        help="Sample index (1-1000). Loads data/samples/sample_XXXX.csv"
    )
    parser.add_argument(
        "--sample_file", type=str, default=None,
        help="Directly specify sample CSV path (overrides --sample)"
    )
    parser.add_argument(
        "--B", type=int, default=50,
        help="Number of Kim-24 bootstrap replicates (default: 50)"
    )
    parser.add_argument(
        "--n_mcmc", type=int, default=2000,
        help="Total MCMC iterations (default: 2000)"
    )
    parser.add_argument(
        "--burn_in", type=int, default=500,
        help="Burn-in iterations (default: 500)"
    )
    parser.add_argument(
        "--m", type=int, default=50,
        help="Number of trees (default: 50)"
    )
    parser.add_argument(
        "--outcome", type=str, default="y",
        help="Outcome column (default: y)"
    )
    parser.add_argument(
        "--features", type=str, nargs="+",
        default=["x_inv", "x_cont_2", "x_bin", "x_cat"],
        help="Covariates to use in BART"
    )
    parser.add_argument(
        "--strata_col", type=str, default="strata"
    )
    parser.add_argument(
        "--psu_col", type=str, default="psuid"
    )
    parser.add_argument(
        "--weight_col", type=str, default="weight"
    )
    parser.add_argument(
        "--seed", type=int, default=20260623,
        help="Random seed"
    )
    parser.add_argument(
        "--output_dir", type=str, default="data/bart_results",
        help="Output directory"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # ---- Load population data to calculate ground truth ----
    pop_path = Path("data/population.csv")
    if pop_path.exists():
        df_pop = pd.read_csv(pop_path)
        true_pop_mean = df_pop[args.outcome].mean()
        print(f"True Population Mean (ground truth): {true_pop_mean:.6f}")
    else:
        true_pop_mean = 22.52220
        print(f"WARNING: population.csv not found. Using default true population mean: {true_pop_mean:.6f}")

    # ---- Resolve sample file ----
    if args.sample_file:
        sample_path = Path(args.sample_file)
    else:
        sample_path = Path(f"data/samples/sample_{args.sample:04d}.csv")

    if not sample_path.exists():
        print(f"ERROR: Sample file not found: {sample_path}", file=sys.stderr)
        sys.exit(1)

    print(f"Loading sample data: {sample_path}")
    df = pd.read_csv(sample_path)
    n = len(df)
    print(f"Sample size: n={n}")

    # Compute fully observed sample survey-weighted mean of y
    y_full = df[args.outcome].values
    weights_full = df[args.weight_col].values
    full_sample_weighted_mean = np.average(y_full, weights=weights_full)
    print(f"Fully Observed Design-Weighted Sample Mean: {full_sample_weighted_mean:.6f}")
    print(f"Fully Observed Unweighted Sample Mean: {y_full.mean():.6f}")

    # One-hot encode the categorical covariates
    available_feats = [f for f in args.features if f in df.columns]
    df_enc = pd.get_dummies(df[available_feats], drop_first=True)
    X = df_enc.values.astype(float)
    y = df[args.outcome].values.astype(float)

    rng = default_rng(args.seed)

    # Missingness rates to test
    missing_rates = [0.05, 0.10, 0.20, 0.50]
    results = []

    print("\nStarting simulation over missingness rates...")
    for m_rate in missing_rates:
        print(f"\n--- Inducing {m_rate*100:.0f}% expected missingness ---")

        # Induce MAR missingness
        R, p_resp, actual_m_rate = induce_mar_missingness(df, m_rate, rng)
        print(f"Actual missingness rate: {actual_m_rate*100:.2f}% (nonrespondents: {np.sum(~R)})")

        # Split sample into respondents and missing
        X_resp, y_resp = X[R], y[R]
        X_miss, y_miss = X[~R], y[~R]
        df_resp = df[R].reset_index(drop=True)

        w_resp = df[R][args.weight_col].values
        w_miss = df[~R][args.weight_col].values
        w_sum = weights_full.sum()

        # Initialize Survey BART model
        prior = BARTPrior(m=args.m)
        model = SurveyBART(
            prior=prior,
            n_mcmc=args.n_mcmc,
            burn_in=args.burn_in,
            seed=args.seed + int(m_rate * 1000)
        )

        # Register missing units as the test set to collect imputed draws
        model.register_test_set(X_miss)

        # Fit model on respondents
        t_start = time.time()
        model.fit(
            X_resp, y_resp, df_resp,
            strata_col=args.strata_col,
            psu_col=args.psu_col,
            weight_col=args.weight_col
        )
        elapsed = time.time() - t_start
        print(f"Fitted Survey BART in {elapsed:.1f}s")

        # Retrieve stochastic predictive draws for missing units: shape (n_post, n_missing)
        test_draws = model._test_draws_

        # 1. Imputation Accuracy (RMSE)
        # We impute using the posterior mean prediction for each missing unit
        imputed_y_mean = test_draws.mean(axis=0)
        impute_rmse = np.sqrt(np.mean((imputed_y_mean - y_miss) ** 2))
        impute_mae = np.mean(np.abs(imputed_y_mean - y_miss))

        # 2. Imputed Population Mean Estimator
        # For each MCMC draw t, the completed-sample design-weighted mean is:
        # y_bar_imp^(t) = (sum_{resp} w_i * y_i + sum_{miss} w_j * y_pred_j^(t)) / sum_S w_k
        w_y_resp_sum = np.dot(w_resp, y_resp)
        w_y_miss_draws = np.dot(test_draws, w_miss)
        y_bar_imp_draws = (w_y_resp_sum + w_y_miss_draws) / w_sum

        # Calculate posterior statistics of the imputed population mean estimator
        est_mean = y_bar_imp_draws.mean()
        est_sd = y_bar_imp_draws.std()
        ci_lower = np.percentile(y_bar_imp_draws, 2.5)
        ci_upper = np.percentile(y_bar_imp_draws, 97.5)
        ci_width = ci_upper - ci_lower

        # Check coverage
        covers_pop = (ci_lower <= true_pop_mean <= ci_upper)
        covers_full = (ci_lower <= full_sample_weighted_mean <= ci_upper)

        bias_vs_pop = est_mean - true_pop_mean
        bias_vs_full = est_mean - full_sample_weighted_mean

        print(f"Imputation RMSE           : {impute_rmse:.4f}")
        print(f"Imputed Pop Mean Estimate : {est_mean:.4f} (95% CI: [{ci_lower:.4f}, {ci_upper:.4f}])")
        print(f"Bias vs Population        : {bias_vs_pop:.4f}")
        print(f"Covers True Pop Mean?     : {covers_pop}")

        results.append({
            "target_missing_rate": m_rate,
            "actual_missing_rate": actual_m_rate,
            "n_respondents": len(X_resp),
            "n_missing": len(X_miss),
            "impute_rmse": impute_rmse,
            "impute_mae": impute_mae,
            "est_mean": est_mean,
            "est_sd": est_sd,
            "ci_lower": ci_lower,
            "ci_upper": ci_upper,
            "ci_width": ci_width,
            "covers_pop": int(covers_pop),
            "covers_full": int(covers_full),
            "bias_vs_pop": bias_vs_pop,
            "bias_vs_full": bias_vs_full,
            "elapsed_s": elapsed
        })

    # ---- Output results ----
    df_res = pd.DataFrame(results)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    stem = Path(args.sample_file).stem if args.sample_file else f"sample_{args.sample:04d}"
    out_csv = out_dir / f"{stem}_missing_summary.csv"
    df_res.to_csv(out_csv, index=False)

    print(f"\n\n{'='*80}")
    print("                      COMPARATIVE SIMULATION RESULTS")
    print(f"{'='*80}")
    print(f"True Population Mean                   : {true_pop_mean:.6f}")
    print(f"Fully Observed Design-Weighted Mean    : {full_sample_weighted_mean:.6f}")
    print(f"{'='*80}")
    print(df_res[[
        "target_missing_rate", "actual_missing_rate", "impute_rmse",
        "est_mean", "ci_lower", "ci_upper", "covers_pop", "elapsed_s"
    ]].to_string(index=False, formatters={
        "target_missing_rate": "{:.0%}".format,
        "actual_missing_rate": "{:.2%}".format,
        "impute_rmse": "{:.4f}".format,
        "est_mean": "{:.4f}".format,
        "ci_lower": "{:.4f}".format,
        "ci_upper": "{:.4f}".format,
        "covers_pop": lambda x: "YES" if x else "NO",
        "elapsed_s": "{:.1f}".format
    }))
    print(f"{'='*80}")
    print(f"Results saved to: {out_csv}\n")


if __name__ == "__main__":
    main()
