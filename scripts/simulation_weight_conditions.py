#!/usr/bin/env python3
"""
simulation_weight_conditions.py
================================
Empirical validation of the one-step variance theory from weight_calibration_proposal.md.

Three weight conditions are tested on the SAME set of n_samples pre-drawn samples:

  Condition A — Normalized to n (mean=1): w* <- w* * n / sum(w*)
      Theory: Term2 = O(1/n) > 0  =>  Var(theta) ≈ Σ_sandwich + V_model
              => credible intervals TOO WIDE (over-coverage).

  Condition B — Raw Kim-Rao weights (K=1, no rescaling):
      Theory: Term2 = O(1/N_hat) ≈ 0  =>  Var(theta) ≈ Σ_sandwich
              => credible intervals target nominal 95% coverage.

  Condition C — K-scaled weights (K = n / (w_bar * n_min)):
      Theory: same as B but with stronger Term2 collapse.
              K does not affect Term1 (scale-invariant estimator).
              => coverage should also target 95%.

Metrics (across B replications):
  - 95% CI coverage rate  (target: 0.95 for B, C; expected >0.95 for A)
  - Median 95% CI width   (narrower = more efficient; B and C should be sharper)
  - RMSE of posterior mean (should be identical A≈B≈C — scale-invariance)
  - Design effect (Deff)   (measures weight heterogeneity per sample)

Population: data/population.csv (N=65,360).
Test set:   2000 population units NOT in the training sample, drawn fresh each rep.

Usage:
    /Users/alvaroquijano/anaconda3/bin/python3 scripts/simulation_weight_conditions.py \
        --n_samples 50 --n_mcmc 800 --burn_in 250 --n_trees 50
"""

import argparse
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from numpy.random import default_rng

warnings.filterwarnings("ignore")

# ── Path setup ────────────────────────────────────────────────────────────────
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "archive"))

from survey_bart import generate_kim24_weights, BARTPrior, BARTree

# =============================================================================
# Constants
# =============================================================================
FEATURE_COLS = ["x_inv", "x_cont_2", "x_bin", "x_cat"]
TARGET_COL   = "y"
STRATA_COL   = "strata"
PSU_COL      = "psuid"
WEIGHT_COL   = "weight"
N_TEST       = 2000   # held-out population units per replication


def design_effect(weights: np.ndarray) -> float:
    """Kish (1965) design effect."""
    n = len(weights)
    return float((n * (weights ** 2).sum()) / (weights.sum() ** 2))


# =============================================================================
# Fast test set construction — fully vectorized, no Python-level row iteration
# =============================================================================

def build_test_set(pop: pd.DataFrame, df_train: pd.DataFrame,
                   n_test: int, rng: np.random.Generator) -> pd.DataFrame:
    """
    Return up to n_test population units that are NOT in df_train.
    Key matching is vectorized using string concatenation on strata+psuid+subid.
    """
    # Vectorized composite key (O(N) string ops, no Python lambda)
    train_keys = set(
        df_train[STRATA_COL].astype(str) + "_" +
        df_train[PSU_COL].astype(str)    + "_" +
        df_train["subid"].astype(str)
    )
    pop_keys = (
        pop[STRATA_COL].astype(str) + "_" +
        pop[PSU_COL].astype(str)    + "_" +
        pop["subid"].astype(str)
    )
    mask_out = ~pop_keys.isin(train_keys)
    pool = pop[mask_out]
    n_draw = min(n_test, len(pool))
    return pool.sample(n=n_draw, random_state=int(rng.integers(1_000_000)))


# =============================================================================
# Core: fit one BART model under a given weight condition
# =============================================================================

def run_one_condition(
    df_train: pd.DataFrame,
    X_test: np.ndarray,
    y_test: np.ndarray,
    condition: str,
    n_mcmc: int,
    burn_in: int,
    n_trees: int,
    rng: np.random.Generator,
) -> dict:
    """
    Replicate the SurveyBART backfitting loop with condition-specific weight
    scaling, then evaluate 95% CI coverage and width on (X_test, y_test).

    Weight transforms:
      A: w* <- w* * n / sum(w*)          [normalize to mean=1]
      B: w* <- w*                         [raw, no rescaling]
      C: w* <- K * w*   where K = n / (w_bar * n_min)
    """
    n = len(df_train)
    w_orig  = df_train[WEIGHT_COL].values.astype(float)
    w_bar   = w_orig.mean()
    n_min   = 5                                      # minimum leaf floor for K
    K_c     = max(1.0, n / (w_bar * n_min))          # §2.5 formula

    X_train = df_train[FEATURE_COLS].values.astype(float)
    y_train = df_train[TARGET_COL].values.astype(float)

    # ── Prior and initialisation ──────────────────────────────────────────────
    prior = BARTPrior(m=n_trees, alpha=0.95, beta=2.0, k=2.0, nu=3.0, q=0.90)

    y_shift = y_train.mean()
    y_scale = y_train.std() if y_train.std() > 0 else 1.0
    y_norm  = (y_train - y_shift) / y_scale
    prior.calibrate(y_norm)

    trees = [BARTree(prior, rng) for _ in range(n_trees)]
    init_mu = y_norm.mean() / n_trees
    for tree in trees:
        tree.root.mu = init_mu

    sigma2 = max(y_norm.var(), 1e-6)

    # Fallback weights for zero-weight edge case
    w_fallback = w_orig * (n / w_orig.sum())

    # ── Storage ───────────────────────────────────────────────────────────────
    n_post = n_mcmc - burn_in
    f_draws   = np.zeros((n_post, len(y_test)))   # f(x) draws  — no obs noise
    pp_draws  = np.zeros((n_post, len(y_test)))   # posterior predictive — with noise

    # ── MCMC loop ─────────────────────────────────────────────────────────────
    for t in range(n_mcmc):
        # Draw fresh Kim-24 bootstrap weights
        w_raw = generate_kim24_weights(
            df_train, STRATA_COL, PSU_COL, WEIGHT_COL, rng
        )
        w_sum = w_raw.sum()

        if w_sum < 1e-10:
            w_star = w_fallback.copy()
        elif condition == "A":
            w_star = w_raw * (n / w_sum)   # normalize to sum = n
        elif condition == "B":
            w_star = w_raw                  # raw weights
        elif condition == "C":
            w_star = K_c * w_raw            # K-scaled raw weights
        else:
            raise ValueError(f"Unknown condition: {condition}")

        # Bayesian backfitting
        f_all = np.array([tree.predict(X_train) for tree in trees])
        for j, tree in enumerate(trees):
            R_j = y_norm - (f_all.sum(axis=0) - f_all[j])
            if rng.random() < 0.5:
                tree._grow_proposal(X_train, R_j, w_star, sigma2)
            else:
                tree._prune_proposal(X_train, R_j, w_star, sigma2)
            tree.gibbs_sample_leaves(X_train, R_j, w_star, sigma2)
            f_all[j] = tree.predict(X_train)

        # sigma² Gibbs update — always use n-scaled weights for scale stability
        y_hat = f_all.sum(axis=0)
        resid = y_norm - y_hat
        w_sigma = w_star / w_star.sum() * n   # n-scaled regardless of condition
        wrss    = np.dot(w_sigma, resid ** 2)
        nu_post = prior.nu + n
        ss_post = prior.nu * prior.lam + wrss
        sigma2  = ss_post / rng.chisquare(nu_post)

        # Store post-burn-in draws
        if t >= burn_in:
            idx    = t - burn_in
            f_test = sum(tree.predict(X_test) for tree in trees)
            f_test_orig = f_test * y_scale + y_shift

            # f(x) draws — no observation noise (measures Term1 + Term2 directly)
            f_draws[idx] = f_test_orig

            # Posterior predictive — add obs noise (for y-coverage)
            noise = rng.normal(0, np.sqrt(sigma2), size=len(y_test))
            pp_draws[idx] = f_test_orig + noise * y_scale

    # ── Metrics ───────────────────────────────────────────────────────────────
    # 1. Posterior mean RMSE (scale-invariant, should be equal A=B=C)
    post_mean = f_draws.mean(axis=0)
    rmse      = float(np.sqrt(np.mean((post_mean - y_test) ** 2)))

    # 2. MCMC SD of f(x) — directly measures sqrt(Term1 + Term2)
    #    Larger for A (Term2 > 0), should be equal for B and C
    mcmc_sd_f = float(f_draws.std(axis=0).mean())

    # 3. 95% CI for f(x) — width reflects Term1 + Term2 (no noise)
    f_lower = np.percentile(f_draws, 2.5,  axis=0)
    f_upper = np.percentile(f_draws, 97.5, axis=0)
    f_ci_width = float((f_upper - f_lower).mean())

    # 4. Posterior predictive 95% coverage for individual y_i
    #    Should be ~95% for all conditions; any deviation reveals miscalibration
    pp_lower = np.percentile(pp_draws, 2.5,  axis=0)
    pp_upper = np.percentile(pp_draws, 97.5, axis=0)
    pp_coverage = float(((y_test >= pp_lower) & (y_test <= pp_upper)).mean())
    pp_width    = float((pp_upper - pp_lower).mean())

    deff   = design_effect(w_orig)
    K_used = K_c if condition == "C" else (1.0 if condition == "B" else np.nan)

    return {
        "condition":   condition,
        "rmse":        rmse,
        "mcmc_sd_f":   mcmc_sd_f,     # sd of f(x) MCMC draws = sqrt(Term1+Term2)
        "f_ci_width":  f_ci_width,    # CI width for f(x)  — theory comparison
        "pp_coverage": pp_coverage,   # posterior predictive coverage for y
        "pp_width":    pp_width,      # posterior predictive CI width for y
        "deff":        deff,
        "K":           K_used,
        "N_hat":       float(w_orig.sum()),
        "n":           n,
    }



# =============================================================================
# Main simulation loop
# =============================================================================

def main():
    parser = argparse.ArgumentParser()
    # ── Local loop mode ──────────────────────────────────────────────────────
    parser.add_argument("--n_samples",     type=int,  default=50,
                        help="Number of sample files to loop over (local mode)")
    # ── SLURM array mode (single sample) ─────────────────────────────────────
    parser.add_argument("--single_sample", type=str,  default=None,
                        help="Path to a single sample CSV (SLURM array mode)")
    parser.add_argument("--out_file",      type=str,  default=None,
                        help="Output CSV path for this task (SLURM array mode)")
    # ── Shared MCMC settings ─────────────────────────────────────────────────
    parser.add_argument("--n_mcmc",    type=int, default=800)
    parser.add_argument("--burn_in",   type=int, default=250)
    parser.add_argument("--n_trees",   type=int, default=50)
    parser.add_argument("--seed",      type=int, default=42)
    args = parser.parse_args()


    rng        = default_rng(args.seed)
    sample_dir = ROOT / "data" / "samples"
    pop_path   = ROOT / "data" / "population.csv"
    conditions = ["A", "B", "C"]

    # Load population ONCE
    print(f"Loading population from {pop_path} ...", flush=True)
    t0 = time.time()
    pop = pd.read_csv(pop_path)
    print(f"  N = {len(pop):,} rows loaded in {time.time()-t0:.1f}s\n", flush=True)

    sample_files = sorted(sample_dir.glob("sample_*.csv"))[:args.n_samples]
    print(f"Simulation: {len(sample_files)} replications × {len(conditions)} conditions")
    print(f"MCMC: {args.n_mcmc} iters | burn-in: {args.burn_in} | trees: {args.n_trees}\n",
          flush=True)

    results = []
    for rep_idx, sample_path in enumerate(sample_files):
        df_train = pd.read_csv(sample_path)
        n        = len(df_train)

        # Fast test set (vectorized key matching)
        df_test = build_test_set(pop, df_train, N_TEST,
                                 default_rng(int(rng.integers(1_000_000))))
        X_test  = df_test[FEATURE_COLS].values.astype(float)
        y_test  = df_test[TARGET_COL].values.astype(float)

        t_rep = time.time()
        for cond in conditions:
            res = run_one_condition(
                df_train=df_train,
                X_test=X_test,
                y_test=y_test,
                condition=cond,
                n_mcmc=args.n_mcmc,
                burn_in=args.burn_in,
                n_trees=args.n_trees,
                rng=default_rng(int(rng.integers(1_000_000_000))),
            )
            res["rep"] = rep_idx + 1
            results.append(res)
            print(f"  Rep {rep_idx+1:3d} | Cond {cond} | "
                  f"PP-Cov={100*res['pp_coverage']:.1f}%  "
                  f"f-SD={res['mcmc_sd_f']:.3f}  "
                  f"f-CI={res['f_ci_width']:.3f}  "
                  f"RMSE={res['rmse']:.4f}  "
                  f"Deff={res['deff']:.3f}  "
                  f"K={res['K']:.1f}", flush=True)

        elapsed = time.time() - t_rep
        print(f"  Rep {rep_idx+1} done in {elapsed:.1f}s\n", flush=True)

    # ── Save and summarise ────────────────────────────────────────────────────
    df_res   = pd.DataFrame(results)
    out_path = ROOT / "data" / "simulation_weight_conditions.csv"
    df_res.to_csv(out_path, index=False)
    print(f"\nResults saved to: {out_path}")

    print("\n" + "=" * 78)
    print("SIMULATION SUMMARY")
    print("=" * 78)
    summary = df_res.groupby("condition").agg(
        pp_cov_mean   =("pp_coverage", "mean"),
        pp_cov_std    =("pp_coverage", "std"),
        pp_width_mean =("pp_width",    "mean"),
        f_sd_mean     =("mcmc_sd_f",   "mean"),   # sqrt(Term1+Term2) — key theory metric
        f_sd_std      =("mcmc_sd_f",   "std"),
        f_ci_mean     =("f_ci_width",  "mean"),
        rmse_mean     =("rmse",        "mean"),
        deff_mean     =("deff",        "mean"),
    ).round(4)
    print(summary.to_string())

    print("\nInterpretation:")
    # f(x) MCMC-SD ratio A/B — should equal sqrt(1 + 1/Deff) per theory
    sd_A = df_res[df_res.condition=="A"]["mcmc_sd_f"].mean()
    sd_B = df_res[df_res.condition=="B"]["mcmc_sd_f"].mean()
    sd_C = df_res[df_res.condition=="C"]["mcmc_sd_f"].mean()
    deff_mean = df_res["deff"].mean()
    print(f"  MCMC-SD ratio  A/B = {sd_A/sd_B:.3f}  "
          f"(theory: sqrt(1+1/Deff) = {(1+1/deff_mean)**0.5:.3f})")
    print(f"  MCMC-SD ratio  A/C = {sd_A/sd_C:.3f}")
    print("  PP-Coverage  : ~95% for all conditions (obs noise dominates)")
    print("  f-SD / f-CI  : Term1+Term2 — A > B ≈ C validates the proposal")
    print("  RMSE         : equal A≈B≈C (scale-invariance of weighted mean)")


if __name__ == "__main__":
    main()

