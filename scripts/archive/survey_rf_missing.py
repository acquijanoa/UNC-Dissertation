#!/usr/bin/env python3
"""
Survey Random Forest Imputation Simulation under MAR Missingness.

This script loads a sample, induces Missing at Random (MAR) nonresponse
on the outcome 'y' at specified levels (5%, 10%, 20%, and 50%),
fits a Survey Random Forest (Decision Trees grown on Kim-24 survey bootstrap replicates
with random feature selection) to impute the missing outcomes, and evaluates 
imputation accuracy and population mean estimation.

It compares:
1. Imputation under original weights (where respondents' weights are fixed).
2. Imputation under replicate weights (which propagates the design variance of observed units).
"""

import argparse
import sys
import time
from pathlib import Path
import numpy as np
import pandas as pd
from numpy.random import default_rng
from sklearn.tree import DecisionTreeRegressor

# Add the directory containing this script to sys.path
script_dir = Path(__file__).parent.resolve()
sys.path.append(str(script_dir))

try:
    from survey_bart import generate_kim24_weights
except ImportError:
    print("ERROR: Could not import generate_kim24_weights. Ensure survey_bart.py is in the same directory.", file=sys.stderr)
    sys.exit(1)


def logit_func(x):
    return 1.0 / (1.0 + np.exp(-x))


def induce_mar_missingness(df, missing_rate, rng):
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
        description="Survey Random Forest Imputation Simulation under MAR"
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
        "--B", type=int, default=100,
        help="Number of bootstrap replicates/trees (default: 100)"
    )
    parser.add_argument(
        "--max_features", type=str, default="sqrt",
        help="Feature selection parameter for DecisionTreeRegressor (default: 'sqrt')"
    )
    parser.add_argument(
        "--outcome", type=str, default="y",
        help="Outcome column (default: y)"
    )
    parser.add_argument(
        "--features", type=str, nargs="+",
        default=["x_inv", "x_cont_2", "x_bin", "x_cat"],
        help="Covariates to use in the model"
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
        "--output_dir", type=str, default="data/rf_results",
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

        w_resp = df[R][args.weight_col].values
        w_miss = df[~R][args.weight_col].values
        w_sum = weights_full.sum()

        test_draws = np.zeros((args.B, len(X_miss)))
        y_bar_imp_orig = np.zeros(args.B)
        y_bar_imp_rep = np.zeros(args.B)

        # Fit B trees using survey bootstrap replicates (Kim-24)
        t_start = time.time()
        for b in range(args.B):
            # Generate fresh bootstrap weights for the entire dataset
            w_star = generate_kim24_weights(
                df, strata_col=args.strata_col, psu_col=args.psu_col, weight_col=args.weight_col, rng=rng
            )
            w_star_resp = w_star[R]
            w_star_miss = w_star[~R]
            w_star_sum = w_star.sum()

            # Filter training data to units with positive bootstrap weight
            in_bag = w_star_resp > 0
            if in_bag.sum() < 2:
                # Fallback to avoid empty sets (extremely rare)
                X_train, y_train, w_train = X_resp, y_resp, w_resp
            else:
                X_train = X_resp[in_bag]
                y_train = y_resp[in_bag]
                w_train = w_star_resp[in_bag]

            # Fit one tree (a constituent tree of our survey random forest)
            tree_seed = int(rng.integers(0, 2**31 - 1))
            tree = DecisionTreeRegressor(max_features=args.max_features, random_state=tree_seed)
            tree.fit(X_train, y_train, sample_weight=w_train)

            # Predict mean for missing units
            y_pred_miss = tree.predict(X_miss)

            # Calculate weighted residual variance of this tree on respondents
            y_pred_resp = tree.predict(X_resp)
            residuals_resp = y_resp - y_pred_resp
            weighted_rss = np.sum(w_star_resp * (residuals_resp ** 2))
            weighted_w_sum = w_star_resp.sum()
            
            if weighted_w_sum > 0:
                sig2 = weighted_rss / weighted_w_sum
            else:
                sig2 = np.var(residuals_resp)

            # Add random residual noise to make the imputation proper (stochastic imputation)
            sig = np.sqrt(max(sig2, 1e-10))
            epsilon = rng.normal(0, sig, size=len(X_miss))
            imputed_draw = y_pred_miss + epsilon
            test_draws[b] = imputed_draw

            # 1. Original Weights Estimator (Hold respondents fixed, like typical Bayesian MCMC setups)
            # y_bar = (sum_resp w_i * y_i + sum_miss w_j * y_imp_j) / sum_S w_k
            y_bar_orig = (np.dot(w_resp, y_resp) + np.dot(w_miss, imputed_draw)) / w_sum
            y_bar_imp_orig[b] = y_bar_orig

            # 2. Replicate Weights Estimator (Vary weights for both respondents and nonrespondents; design-consistent)
            # y_bar* = (sum_resp w*_i * y_i + sum_miss w*_j * y_imp_j) / sum_S w*_k
            y_bar_rep = (np.dot(w_star_resp, y_resp) + np.dot(w_star_miss, imputed_draw)) / w_star_sum
            y_bar_imp_rep[b] = y_bar_rep

        elapsed = time.time() - t_start
        print(f"Fitted Survey Random Forest ({args.B} trees) and imputed in {elapsed:.1f}s")

        # 1. Imputation Accuracy (RMSE of the posterior mean)
        imputed_y_mean = test_draws.mean(axis=0)
        impute_rmse = np.sqrt(np.mean((imputed_y_mean - y_miss) ** 2))
        impute_mae = np.mean(np.abs(imputed_y_mean - y_miss))

        # 2. Evaluation: Original-weighted completed mean (Fixing observed units)
        est_mean_orig = y_bar_imp_orig.mean()
        est_sd_orig = y_bar_imp_orig.std()
        ci_lower_orig = np.percentile(y_bar_imp_orig, 2.5)
        ci_upper_orig = np.percentile(y_bar_imp_orig, 97.5)
        covers_pop_orig = (ci_lower_orig <= true_pop_mean <= ci_upper_orig)

        # 3. Evaluation: Replicate-weighted completed mean (Varying observed and missing weights; design-consistent)
        est_mean_rep = y_bar_imp_rep.mean()
        est_sd_rep = y_bar_imp_rep.std()
        ci_lower_rep = np.percentile(y_bar_imp_rep, 2.5)
        ci_upper_rep = np.percentile(y_bar_imp_rep, 97.5)
        covers_pop_rep = (ci_lower_rep <= true_pop_mean <= ci_upper_rep)

        bias_vs_pop_rep = est_mean_rep - true_pop_mean

        print(f"Imputation RMSE                     : {impute_rmse:.4f}")
        print(f"Original-Weighted CI (fixed resp)   : [{ci_lower_orig:.4f}, {ci_upper_orig:.4f}] (covers pop? {covers_pop_orig})")
        print(f"Replicate-Weighted CI (varying resp): [{ci_lower_rep:.4f}, {ci_upper_rep:.4f}] (covers pop? {covers_pop_rep})")
        print(f"CI Width Comparison                 : Original={ci_upper_orig-ci_lower_orig:.4f} vs. Replicate={ci_upper_rep-ci_lower_rep:.4f}")

        results.append({
            "target_missing_rate": m_rate,
            "actual_missing_rate": actual_m_rate,
            "impute_rmse": impute_rmse,
            "impute_mae": impute_mae,
            "est_mean_orig": est_mean_orig,
            "est_sd_orig": est_sd_orig,
            "ci_lower_orig": ci_lower_orig,
            "ci_upper_orig": ci_upper_orig,
            "covers_pop_orig": int(covers_pop_orig),
            "est_mean_rep": est_mean_rep,
            "est_sd_rep": est_sd_rep,
            "ci_lower_rep": ci_lower_rep,
            "ci_upper_rep": ci_upper_rep,
            "covers_pop_rep": int(covers_pop_rep),
            "bias_vs_pop": bias_vs_pop_rep,
            "elapsed_s": elapsed
        })

    # ---- Save results ----
    df_res = pd.DataFrame(results)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = Path(args.sample_file).stem if args.sample_file else f"sample_{args.sample:04d}"
    out_csv = out_dir / f"{stem}_rf_missing_summary.csv"
    df_res.to_csv(out_csv, index=False)

    print(f"\n\n{'='*95}")
    print("                    SURVEY RANDOM FOREST COMPARATIVE RESULTS")
    print(f"{'='*95}")
    print(f"True Population Mean                      : {true_pop_mean:.6f}")
    print(f"Fully Observed Design-Weighted Mean       : {full_sample_weighted_mean:.6f}")
    print(f"{'='*95}")
    print(df_res[[
        "target_missing_rate", "impute_rmse",
        "est_mean_orig", "covers_pop_orig", 
        "est_mean_rep", "covers_pop_rep", "elapsed_s"
    ]].to_string(index=False, formatters={
        "target_missing_rate": "{:.0%}".format,
        "impute_rmse": "{:.4f}".format,
        "est_mean_orig": "{:.4f}".format,
        "covers_pop_orig": lambda x: "YES" if x else "NO",
        "est_mean_rep": "{:.4f}".format,
        "covers_pop_rep": lambda x: "YES" if x else "NO",
        "elapsed_s": "{:.1f}".format
    }))
    print(f"{'='*95}")
    print(f"Results saved to: {out_csv}\n")


if __name__ == "__main__":
    main()
