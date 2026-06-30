#!/usr/bin/env python3
"""
survey_bart.py
==============
Random Weight Survey BART — Application 2 of the dissertation proposal.

This script implements the Random Weight MCMC framework (proposal §PART III)
applied to Bayesian Additive Regression Trees (BART). At each MCMC iteration,
a fresh set of Kim, Rao & Wang (2024) bootstrap weights w* is drawn and used to:

  1. Evaluate the weighted pseudo-likelihood in MH tree structure proposals.
  2. Compute w*-weighted sufficient statistics in the Gibbs terminal node draws.

The Kim-24 bootstrap (Example 1, stratified multi-stage sampling):
  For each stratum h with n_h sampled PSUs:
    - Draw r*_h ~ Multinomial(n_h - 1, 1/n_h * 1_{n_h})
    - Set k_h = n_h / (n_h - 1)
    - Broadcast w*_i = k_h * r*_{hi} * w_i to all units in PSU i

This satisfies Kim-24 Conditions C6-C7:
    E*{S*_w(θ)} = S_w(θ)  and  V*{S*_w(θ)} = V^{S_w(θ) | F_N}

Usage:
    python survey_bart.py --sample 1 --B 50 --n_mcmc 2000 --burn_in 500
    python survey_bart.py --sample_file data/samples/sample_0001.csv --B 50

Output:
    data/bart_results/sample_XXXX_bart.npz
        - 'y_pred_mean'     : posterior mean predictions (n_test,)
        - 'y_pred_draws'    : post-burn-in posterior draws (n_draws, n_test)
        - 'weights_matrix'  : Kim-24 bootstrap weight matrix (B, n_train)
        - 'oob_mask_matrix' : OOB indicator matrix (B, n_train) — ~37% True
        - 'rmse_train'      : in-sample RMSE per MCMC step
        - 'coverage_95'     : 95% credible interval coverage on test set
"""

import argparse
import os
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from numpy.random import default_rng

warnings.filterwarnings("ignore")


# =============================================================================
# 1.  KIM-24 BOOTSTRAP WEIGHT GENERATOR
# =============================================================================

def generate_kim24_weights(df: pd.DataFrame,
                            strata_col: str = "strata",
                            psu_col: str = "psuid",
                            weight_col: str = "weight",
                            rng: np.random.Generator = None) -> np.ndarray:
    """
    Generate one set of Kim-24 bootstrap weights for a stratified cluster sample.

    Algorithm (Kim, Rao & Wang 2024, Example 1 / Section 3):
      For each stratum h:
        1. Identify the n_h distinct sampled PSUs.
        2. Draw r*_h ~ Multinomial(n_h - 1, (1/n_h, ..., 1/n_h)).
        3. Rescaling factor: k_h = n_h / (n_h - 1).
        4. For every unit belonging to PSU i in stratum h:
               w*_i = k_h * r*_{hi} * w_i

    Properties guaranteed:
      E*[w*_i] = w_i           (unbiasedness)
      V*[S*_w] = V^{S_w | F_N} (correct design variance replication)

    Returns
    -------
    w_star : np.ndarray of shape (n,)
        Bootstrap weights. Approximately 1/e ≈ 37% of PSUs get weight 0
        (the 'OOB' dropout property exploited by the Random Weight MCMC).
    """
    if rng is None:
        rng = default_rng()

    n = len(df)
    w_star = np.zeros(n, dtype=float)

    # Iterate over strata
    for h, stratum_idx in df.groupby(strata_col).groups.items():
        stratum_df = df.loc[stratum_idx]

        # Unique PSUs in this stratum
        psus = stratum_df[psu_col].unique()
        n_h = len(psus)

        if n_h == 1:
            # Edge case: single PSU — cannot compute variance; pass original weight
            w_star[stratum_idx] = stratum_df[weight_col].values
            continue

        # Kim-24: Multinomial(n_h - 1, uniform) rescaling draw
        # r*_h is a vector of length n_h summing to n_h - 1
        r_h = rng.multinomial(n_h - 1, np.ones(n_h) / n_h)  # shape (n_h,)
        k_h = n_h / (n_h - 1)  # rescaling factor

        # Broadcast to individual rows
        psu_to_r = dict(zip(psus, r_h))
        r_expanded = stratum_df[psu_col].map(psu_to_r).values.astype(float)
        original_weights = stratum_df[weight_col].values

        w_star[stratum_idx] = k_h * r_expanded * original_weights

    return w_star


def generate_kim24_weight_matrix(df: pd.DataFrame,
                                  B: int,
                                  strata_col: str = "strata",
                                  psu_col: str = "psuid",
                                  weight_col: str = "weight",
                                  seed: int = None) -> np.ndarray:
    """
    Generate the full B x n matrix of Kim-24 bootstrap weights.

    Returns
    -------
    W : np.ndarray of shape (B, n)
        Row b is one complete bootstrap weight draw w*^(b).
    oob_mask : np.ndarray of shape (B, n), dtype bool
        True where w*^(b)_i == 0 (the OOB 'dropout' units).
    """
    rng = default_rng(seed)
    n = len(df)
    W = np.zeros((B, n), dtype=float)

    for b in range(B):
        W[b] = generate_kim24_weights(df, strata_col, psu_col, weight_col, rng)

    oob_mask = (W == 0)
    return W, oob_mask


# =============================================================================
# 2.  BART PRIOR & LIKELIHOOD UTILITIES
# =============================================================================

class BARTPrior:
    """
    Standard Chipman-George-McCulloch (CGM) BART prior.

    Tree prior parameters:
        alpha : probability a node at depth d is non-terminal = alpha / (1 + d)^beta
        beta  : depth penalisation exponent
        k     : leaf parameter prior = mu_jl ~ N(0, sigma_mu^2)
                where sigma_mu = (y_max - y_min) / (2 * k * sqrt(m))

    Variance prior:
        sigma^2 ~ Inv-Chi^2(nu, lambda)
        Chosen so P(sigma < sigma_hat) = 0.9 (data-driven calibration).
    """

    def __init__(self, m: int = 200, alpha: float = 0.95, beta: float = 2.0,
                 k: float = 2.0, nu: float = 3.0, q: float = 0.90):
        self.m = m          # number of trees
        self.alpha = alpha  # tree prior: base
        self.beta = beta    # tree prior: power
        self.k = k          # leaf width control
        self.nu = nu        # sigma^2 prior df
        self.q = q          # sigma^2 prior quantile calibration

    def calibrate(self, y_train: np.ndarray):
        """Calibrate sigma_mu and lambda from training data."""
        y_range = y_train.max() - y_train.min()
        self.sigma_mu = y_range / (2.0 * self.k * np.sqrt(self.m))
        # Rough OLS residual sigma as sigma_hat
        sigma_hat = np.std(y_train) * 0.9
        from scipy.stats import chi2 as _chi2
        self.lam = (sigma_hat**2 * _chi2.ppf(1 - self.q, self.nu)) / self.nu
        return self


# =============================================================================
# 3.  DECISION TREE NODE (lightweight)
# =============================================================================

class Node:
    """Lightweight binary tree node for BART."""
    __slots__ = ["is_leaf", "feature", "threshold", "left", "right",
                 "mu", "depth", "node_id"]

    def __init__(self, depth: int = 0, node_id: int = 0):
        self.is_leaf = True
        self.feature = None
        self.threshold = None
        self.left = None
        self.right = None
        self.mu = 0.0
        self.depth = depth
        self.node_id = node_id

    def predict(self, X: np.ndarray) -> np.ndarray:
        """Route observations and return leaf mu values."""
        n = X.shape[0]
        out = np.empty(n)
        if self.is_leaf:
            out[:] = self.mu
            return out
        mask = X[:, self.feature] <= self.threshold
        if mask.any():
            out[mask] = self.left.predict(X[mask])
        if (~mask).any():
            out[~mask] = self.right.predict(X[~mask])
        return out

    def get_leaves(self):
        """Return list of all leaf nodes."""
        if self.is_leaf:
            return [self]
        return self.left.get_leaves() + self.right.get_leaves()

    def get_internal_nodes(self):
        """Return list of all non-leaf nodes."""
        if self.is_leaf:
            return []
        return [self] + self.left.get_internal_nodes() + self.right.get_internal_nodes()

    def route(self, X: np.ndarray):
        """Return leaf index for each row of X."""
        n = X.shape[0]
        result = np.empty(n, dtype=int)
        if self.is_leaf:
            result[:] = self.node_id
            return result
        mask = X[:, self.feature] <= self.threshold
        if mask.any():
            result[mask] = self.left.route(X[mask])
        if (~mask).any():
            result[~mask] = self.right.route(X[~mask])
        return result

    def count_nodes(self):
        if self.is_leaf:
            return 1
        return 1 + self.left.count_nodes() + self.right.count_nodes()


# =============================================================================
# 4.  SINGLE BART TREE
# =============================================================================

class BARTree:
    """One regression tree in the BART ensemble."""

    def __init__(self, prior: BARTPrior, rng: np.random.Generator):
        self.prior = prior
        self.rng = rng
        self._node_counter = 0
        self.root = Node(depth=0, node_id=self._next_id())

    def _next_id(self) -> int:
        nid = self._node_counter
        self._node_counter += 1
        return nid

    def predict(self, X: np.ndarray) -> np.ndarray:
        return self.root.predict(X)

    # ------------------------------------------------------------------
    #  MH proposal: GROW — split a random leaf into two children
    # ------------------------------------------------------------------
    def _grow_proposal(self, X: np.ndarray, R: np.ndarray,
                        w_star: np.ndarray, sigma2: float) -> tuple:
        """
        Propose growing a random leaf node into an internal node.

        Uses w*-weighted log-likelihood for the MH acceptance ratio:
            log p_w*(R | T*, M*) - log p_w*(R | T, M)
        where the weighted log-likelihood of a leaf l with unit indices idx is:
            -0.5 * sum_i w*_i * (R_i - mu_l)^2 / sigma2
            (evaluated at the posterior mean mu_l for numerical stability)

        Returns (accepted, modified_root_copy) or (False, None).
        """
        leaves = self.root.get_leaves()
        # Only grow leaves that have enough data (at least 2 units with w*>0)
        eligible = [l for l in leaves
                    if self._leaf_effective_n(X, l, w_star) >= 2]
        if not eligible:
            return False, None

        leaf = self.rng.choice(eligible)
        # Choose a random split feature and threshold
        idx = self._get_leaf_indices(X, leaf)
        X_leaf = X[idx]
        feat = self.rng.integers(0, X.shape[1])
        vals = np.unique(X_leaf[:, feat])
        if len(vals) < 2:
            return False, None
        threshold = self.rng.choice(vals[:-1])

        # MH log-ratio: log-likelihood gain from the split
        w_l = w_star[idx]
        R_l = R[idx]
        mask_left = X_leaf[:, feat] <= threshold
        mask_right = ~mask_left

        if mask_left.sum() == 0 or mask_right.sum() == 0:
            return False, None

        # Weighted log-likelihood of current leaf (no split)
        ll_before = self._weighted_leaf_ll(R_l, w_l, sigma2)

        # Weighted log-likelihood after split
        ll_left  = self._weighted_leaf_ll(R_l[mask_left],  w_l[mask_left],  sigma2)
        ll_right = self._weighted_leaf_ll(R_l[mask_right], w_l[mask_right], sigma2)
        ll_after = ll_left + ll_right

        # Prior ratio: (alpha / (1 + depth)^beta) * (1 - alpha)^2
        depth = leaf.depth
        prior_grow = self.prior.alpha / (1 + depth) ** self.prior.beta
        prior_ratio = np.log(prior_grow) + 2 * np.log(1 - self.prior.alpha /
                      (1 + depth + 1) ** self.prior.beta) - np.log(1 - prior_grow)

        # Transition ratio (simplified: uniform leaf/feature/threshold selection)
        log_accept = ll_after - ll_before + prior_ratio
        if np.log(self.rng.random()) < log_accept:
            # Apply grow
            leaf.is_leaf = False
            leaf.feature = feat
            leaf.threshold = threshold
            leaf.left  = Node(depth=leaf.depth + 1, node_id=self._next_id())
            leaf.right = Node(depth=leaf.depth + 1, node_id=self._next_id())
            return True, None
        return False, None

    # ------------------------------------------------------------------
    #  MH proposal: PRUNE — collapse a random pair of sibling leaves
    # ------------------------------------------------------------------
    def _prune_proposal(self, X: np.ndarray, R: np.ndarray,
                         w_star: np.ndarray, sigma2: float) -> tuple:
        """
        Propose pruning an internal node whose both children are leaves.
        MH ratio is the reverse of GROW.
        """
        internals = self.root.get_internal_nodes()
        prunable = [n for n in internals
                    if n.left.is_leaf and n.right.is_leaf]
        if not prunable:
            return False, None

        node = self.rng.choice(prunable)
        idx  = self._get_leaf_indices(X, node)
        X_n  = X[idx]
        w_n  = w_star[idx]
        R_n  = R[idx]

        mask_left  = X_n[:, node.feature] <= node.threshold
        mask_right = ~mask_left

        ll_after  = self._weighted_leaf_ll(R_n, w_n, sigma2)
        ll_left   = self._weighted_leaf_ll(R_n[mask_left],  w_n[mask_left],  sigma2)
        ll_right  = self._weighted_leaf_ll(R_n[mask_right], w_n[mask_right], sigma2)
        ll_before = ll_left + ll_right

        depth = node.depth
        prior_grow  = self.prior.alpha / (1 + depth) ** self.prior.beta
        prior_ratio = np.log(1 - prior_grow) - np.log(prior_grow) \
                      - 2 * np.log(1 - self.prior.alpha /
                                   (1 + depth + 1) ** self.prior.beta)

        log_accept = ll_after - ll_before + prior_ratio
        if np.log(self.rng.random()) < log_accept:
            node.is_leaf  = True
            node.feature  = None
            node.threshold = None
            node.left     = None
            node.right    = None
            return True, None
        return False, None

    # ------------------------------------------------------------------
    #  Gibbs: sample terminal node parameters w*-weighted
    # ------------------------------------------------------------------
    def gibbs_sample_leaves(self, X: np.ndarray, R: np.ndarray,
                             w_star: np.ndarray, sigma2: float):
        """
        Draw leaf parameters mu_jl from their w*-weighted Gibbs conditionals.

        The weighted posterior for each leaf l:
            mu_l | R, w*, sigma2, T ~ N(mu_l_hat, V_l)

        where (Kim-24 weighted sufficient statistics):
            W_l     = sum_{i in l} w*_i            (effective weight)
            R_bar_l = sum_{i in l} w*_i * R_i / W_l (weighted mean)
            V_l     = 1 / (W_l / sigma2 + 1 / sigma_mu^2)
            mu_l_hat= V_l * W_l * R_bar_l / sigma2

        Because w*_i ~ 0 for ~37% of PSUs (OOB dropout), those units are
        automatically excluded from the sufficient statistics at each step,
        providing native regularization against cluster overfitting.
        """
        leaves = self.root.get_leaves()
        sigma_mu2 = self.prior.sigma_mu ** 2

        for leaf in leaves:
            idx = self._get_leaf_indices(X, leaf)
            w_l = w_star[idx]
            R_l = R[idx]

            W_l = w_l.sum()
            if W_l < 1e-10:
                # All units in this leaf got OOB weight 0 → draw from prior
                leaf.mu = self.rng.normal(0, self.prior.sigma_mu)
                continue

            # w*-weighted mean residual
            R_bar_l = np.dot(w_l, R_l) / W_l

            # Posterior variance and mean
            V_l     = 1.0 / (W_l / sigma2 + 1.0 / sigma_mu2)
            mu_hat  = V_l * (W_l * R_bar_l / sigma2)
            leaf.mu = self.rng.normal(mu_hat, np.sqrt(V_l))

    # ------------------------------------------------------------------
    #  Utilities
    # ------------------------------------------------------------------
    def _get_leaf_indices(self, X: np.ndarray, leaf: Node) -> np.ndarray:
        """Return integer indices of rows routed to `leaf`."""
        leaf_id = leaf.node_id
        routed  = self.root.route(X)
        return np.where(routed == leaf_id)[0]

    def _leaf_effective_n(self, X: np.ndarray, leaf: Node,
                          w_star: np.ndarray) -> int:
        """Count units in leaf with positive bootstrap weight."""
        idx = self._get_leaf_indices(X, leaf)
        return (w_star[idx] > 0).sum()

    @staticmethod
    def _weighted_leaf_ll(R: np.ndarray, w: np.ndarray,
                           sigma2: float) -> float:
        """
        Weighted Gaussian log-likelihood for a leaf, evaluated at the
        w-weighted posterior mean (integrating out mu analytically is
        preferred but computationally expensive; we use the MAP instead).
        """
        W = w.sum()
        if W < 1e-10:
            return 0.0
        mu_hat = np.dot(w, R) / W
        return -0.5 * np.dot(w, (R - mu_hat) ** 2) / sigma2


# =============================================================================
# 5.  SURVEY BART SAMPLER (Random Weight MCMC)
# =============================================================================

class SurveyBART:
    """
    Random Weight Survey BART.

    Implements the Bayesian backfitting algorithm with Kim-24 bootstrap weight
    injection at every MCMC transition step (proposal §PART III, §1.

    Algorithm per iteration t:
      1. Draw fresh w* ~ Kim-24(S, strata, PSUs)       ← new at each step
      2. For each tree j = 1, ..., m:
           a. Compute partial residuals R_j = y - sum_{j'≠j} g(X; T_{j'}, M_{j'})
           b. Propose GROW or PRUNE to T_j using MH with w*-weighted likelihood
           c. Gibbs-draw mu_jl for all leaves l of T_j using w*-weighted moments
      3. Gibbs-draw sigma^2 using w*-weighted residual sum of squares
    """

    def __init__(self, prior: BARTPrior = None, n_mcmc: int = 2000,
                 burn_in: int = 500, seed: int = None, verbose: bool = True):
        self.prior   = prior or BARTPrior()
        self.n_mcmc  = n_mcmc
        self.burn_in = burn_in
        self.rng     = default_rng(seed)
        self.verbose = verbose
        self.trees_  = []
        self.sigma2_ = 1.0
        self.y_scale_: tuple = (0.0, 1.0)  # (shift, scale) for y normalisation
        self._X_test_store = None           # set by register_test_set()
        self._test_draws_  = None           # populated during fit()

    # ------------------------------------------------------------------
    #  Public API
    # ------------------------------------------------------------------
    def fit(self, X_train: np.ndarray, y_train: np.ndarray,
            df_train: pd.DataFrame,
            strata_col: str = "strata",
            psu_col:    str = "psuid",
            weight_col: str = "weight") -> "SurveyBART":
        """
        Fit the Random Weight Survey BART model.

        Parameters
        ----------
        X_train    : (n, p) covariate matrix
        y_train    : (n,)   outcome vector
        df_train   : DataFrame with survey design columns (strata, PSU, weights)
        """
        n, p = X_train.shape
        m    = self.prior.m

        # --- Normalise y to roughly [-0.5, 0.5] for numerical stability ---
        y_shift = y_train.mean()
        y_scale = y_train.std() if y_train.std() > 0 else 1.0
        self.y_scale_ = (y_shift, y_scale)
        y_norm = (y_train - y_shift) / y_scale

        # --- Calibrate prior on normalised scale ---
        self.prior.calibrate(y_norm)

        # --- Initialise trees as single-leaf stumps ---
        self.trees_ = [BARTree(self.prior, self.rng) for _ in range(m)]
        # Initialise all leaf means to y_norm.mean() / m
        init_mu = y_norm.mean() / m
        for tree in self.trees_:
            tree.root.mu = init_mu

        # --- Initialise sigma^2 ---
        self.sigma2_ = y_norm.var() if y_norm.var() > 0 else 1.0

        # --- Storage ---
        n_post  = self.n_mcmc - self.burn_in
        self.posterior_draws_ = np.zeros((n_post, n))  # in-sample
        self.sigma2_draws_    = np.zeros(n_post)
        self.rmse_train_      = []

        # Kim-24 weight matrix (B=n_mcmc draws, one per MCMC step)
        # We store the w* at each iteration for diagnostics
        self.weight_draws_    = np.zeros((n_post, n))
        self.oob_mask_draws_  = np.zeros((n_post, n), dtype=bool)

        # --- Original survey weights (normalised to sum=n) ---
        w_orig = df_train[weight_col].values.copy().astype(float)
        w_orig_norm = w_orig * (n / w_orig.sum())

        # Optional test-set tracking — preserve if register_test_set() was called
        # before fit(). Only reset _test_draws_ storage if X_test is registered.
        if self._X_test_store is not None:
            # Re-initialise storage with correct n_post (may differ from prior call)
            self._test_draws_ = np.zeros((n_post, self._X_test_store.shape[0]))

        t0 = time.time()
        for t in range(self.n_mcmc):
            # ==============================================================
            # Step 1: Draw fresh Kim-24 bootstrap weights for this iteration
            # ==============================================================
            w_star = generate_kim24_weights(
                df_train, strata_col, psu_col, weight_col, self.rng
            )
            # Normalise bootstrap weights to sum = n (preserves scale)
            w_sum = w_star.sum()
            if w_sum > 0:
                w_star_norm = w_star * (n / w_sum)
            else:
                w_star_norm = w_orig_norm.copy()

            oob = (w_star == 0)

            # ==============================================================
            # Step 2: Bayesian backfitting over all m trees
            # ==============================================================
            # Current ensemble predictions
            f_all = np.array([tree.predict(X_train) for tree in self.trees_])
            # f_all shape: (m, n)

            for j, tree in enumerate(self.trees_):
                # ---- Partial residuals (exclude tree j) ----
                R_j = y_norm - (f_all.sum(axis=0) - f_all[j])

                # ---- MH: propose GROW or PRUNE ----
                if self.rng.random() < 0.5:
                    tree._grow_proposal(X_train, R_j, w_star_norm, self.sigma2_)
                else:
                    tree._prune_proposal(X_train, R_j, w_star_norm, self.sigma2_)

                # ---- Gibbs: sample leaf parameters ----
                tree.gibbs_sample_leaves(X_train, R_j, w_star_norm, self.sigma2_)

                # Update this tree's predictions in f_all
                f_all[j] = tree.predict(X_train)

            # ==============================================================
            # Step 3: Gibbs-draw sigma^2 (w*-weighted residuals)
            # ==============================================================
            y_hat   = f_all.sum(axis=0)
            resid   = y_norm - y_hat
            # Weighted RSS
            wrss    = np.dot(w_star_norm, resid ** 2)
            # Inverse-chi-squared: sigma^2 | data ~ Inv-Chi^2(nu + n, ...)
            nu_post = self.prior.nu + n
            ss_post = self.prior.nu * self.prior.lam + wrss
            self.sigma2_ = ss_post / self.rng.chisquare(nu_post)

            # ==============================================================
            # Diagnostics & storage
            # ==============================================================
            rmse = np.sqrt(np.mean(resid ** 2))
            self.rmse_train_.append(rmse)

            if t >= self.burn_in:
                idx = t - self.burn_in
                # Back-transform to original scale
                self.posterior_draws_[idx] = y_hat * y_scale + y_shift
                self.sigma2_draws_[idx]    = self.sigma2_ * y_scale ** 2
                self.weight_draws_[idx]    = w_star
                self.oob_mask_draws_[idx]  = oob
                # Also track test predictions if X_test was registered
                if self._X_test_store is not None:
                    f_test = sum(tree.predict(self._X_test_store)
                                 for tree in self.trees_)
                    # Independent noise per test unit — correct posterior predictive
                    noise = self.rng.normal(0, np.sqrt(self.sigma2_),
                                            size=self._X_test_store.shape[0])
                    self._test_draws_[idx] = f_test * y_scale + y_shift + noise * y_scale

            if self.verbose and (t + 1) % 200 == 0:
                elapsed = time.time() - t0
                pct     = 100 * (t + 1) / self.n_mcmc
                print(f"  MCMC {t+1:>5}/{self.n_mcmc}  ({pct:5.1f}%)  "
                      f"RMSE={rmse:.4f}  sigma={np.sqrt(self.sigma2_):.4f}  "
                      f"OOB%={100*oob.mean():.1f}%  "
                      f"elapsed={elapsed:.1f}s")

        return self

    def register_test_set(self, X_test: np.ndarray):
        """
        Register a test set so that posterior predictive draws are collected
        during fit() rather than only using the final tree state.
        Call this BEFORE fit().
        """
        n_post = self.n_mcmc - self.burn_in
        self._X_test_store = X_test
        self._test_draws_  = np.zeros((n_post, X_test.shape[0]))
        return self

    def predict(self, X_test: np.ndarray = None) -> dict:
        """
        Compute posterior predictive distribution.

        If register_test_set() was called before fit(), uses the full MCMC
        posterior draws collected during sampling (recommended).
        Otherwise, uses the final tree state + sigma2 posterior for uncertainty.

        Returns
        -------
        dict with:
            'mean'     : posterior mean (n_test,)
            'std'      : posterior std (n_test,)
            'lower_95' : 2.5th percentile (n_test,)
            'upper_95' : 97.5th percentile (n_test,)
            'draws'    : all post-burn-in draws (n_post, n_test)
        """
        y_shift, y_scale = self.y_scale_

        # --- Case 1: Proper posterior draws collected during MCMC ---
        if self._test_draws_ is not None and X_test is None:
            pred_draws = self._test_draws_
        elif self._test_draws_ is not None and X_test is not None:
            pred_draws = self._test_draws_
        else:
            # --- Case 2: Use final tree state + sigma2 uncertainty ---
            if X_test is None:
                raise ValueError("Provide X_test or call register_test_set() before fit()")
            y_hat_test = sum(tree.predict(X_test) * y_scale + y_shift
                             for tree in self.trees_)
            sigma2_post = self.sigma2_draws_
            noise_draws = self.rng.normal(
                0, np.sqrt(sigma2_post)[:, None],
                size=(len(sigma2_post), X_test.shape[0])
            )
            pred_draws = y_hat_test[None, :] + noise_draws

        return {
            "mean":     pred_draws.mean(axis=0),
            "std":      pred_draws.std(axis=0),
            "lower_95": np.percentile(pred_draws, 2.5, axis=0),
            "upper_95": np.percentile(pred_draws, 97.5, axis=0),
            "draws":    pred_draws,
        }

    def coverage_95(self, y_test: np.ndarray) -> float:
        """Compute 95% credible interval coverage on test set."""
        preds = self.predict()
        covered = ((y_test >= preds["lower_95"]) &
                   (y_test <= preds["upper_95"]))
        return covered.mean()


# =============================================================================
# 6.  STANDALONE KIM-24 BOOTSTRAP WEIGHT GENERATOR  (B replicates only)
# =============================================================================

def generate_and_save_bootstrap_weights(df: pd.DataFrame,
                                         B: int,
                                         output_path: str,
                                         seed: int = None,
                                         strata_col: str = "strata",
                                         psu_col: str = "psuid",
                                         weight_col: str = "weight"):
    """
    Generate B Kim-24 bootstrap weight replicates and save to CSV and NPZ.

    This is useful as a standalone preprocessing step before running BART
    (or any other method) on the bootstrap replicates.

    Output files:
        {output_path}_weights.csv  — (n x B) weight matrix, columns = boot_1..boot_B
        {output_path}_weights.npz  — numpy archive with arrays W, oob_mask, w_orig
    """
    W, oob_mask = generate_kim24_weight_matrix(
        df, B=B,
        strata_col=strata_col, psu_col=psu_col, weight_col=weight_col,
        seed=seed
    )

    # --- Save NPZ ---
    npz_path = output_path + "_weights.npz"
    np.savez_compressed(
        npz_path,
        W=W,
        oob_mask=oob_mask,
        w_orig=df[weight_col].values,
        strata=df[strata_col].values,
        psu=df[psu_col].values,
    )

    # --- Save human-readable CSV (W transposed: n rows, B cols) ---
    csv_path = output_path + "_weights.csv"
    col_names = [f"boot_{b+1:04d}" for b in range(B)]
    df_w = pd.DataFrame(W.T, columns=col_names)
    df_w.insert(0, strata_col, df[strata_col].values)
    df_w.insert(1, psu_col, df[psu_col].values)
    df_w.insert(2, "w_orig", df[weight_col].values)
    df_w.to_csv(csv_path, index=False)

    # --- Summary statistics ---
    oob_rate = oob_mask.mean()
    w_cv = W[W > 0].std() / W[W > 0].mean()  # coeff of variation of positive weights

    print(f"\n{'='*60}")
    print(f"Kim-24 Bootstrap Weight Summary (B={B})")
    print(f"{'='*60}")
    print(f"  Sample size n         : {len(df)}")
    print(f"  Number of strata      : {df[strata_col].nunique()}")
    print(f"  Mean PSUs per stratum : {df.groupby(strata_col)[psu_col].nunique().mean():.1f}")
    print(f"  OOB rate (mean)       : {100*oob_rate:.1f}%  (theory: ~36.8%)")
    print(f"  CV of positive w*     : {w_cv:.3f}")
    print(f"  Saved NPZ             : {npz_path}")
    print(f"  Saved CSV             : {csv_path}")
    print(f"{'='*60}\n")

    return W, oob_mask, npz_path, csv_path


# =============================================================================
# 7.  MAIN
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Random Weight Survey BART with Kim-24 bootstrap weights",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument(
        "--sample", type=int, default=1,
        help="Sample index (1-1000). Loads data/samples/sample_XXXX.csv"
    )
    parser.add_argument(
        "--sample_file", type=str, default=None,
        help="Override --sample: directly specify a sample CSV path"
    )
    parser.add_argument(
        "--B", type=int, default=50,
        help="Number of Kim-24 bootstrap weight replicates (default: 50)"
    )
    parser.add_argument(
        "--n_mcmc", type=int, default=2000,
        help="Total MCMC iterations (default: 2000)"
    )
    parser.add_argument(
        "--burn_in", type=int, default=500,
        help="Burn-in iterations to discard (default: 500)"
    )
    parser.add_argument(
        "--m", type=int, default=50,
        help="Number of trees in BART ensemble (default: 50; paper uses 200)"
    )
    parser.add_argument(
        "--outcome", type=str, default="y",
        help="Outcome column name (default: y)"
    )
    parser.add_argument(
        "--features", type=str, nargs="+",
        default=["x_inv", "x_cont_2", "x_bin", "x_cat"],
        help="Covariate column names"
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
        "--weights_only", action="store_true",
        help="Only generate and save Kim-24 bootstrap weights, skip BART"
    )
    parser.add_argument(
        "--test_frac", type=float, default=0.2,
        help="Fraction of data for test set (default: 0.2)"
    )
    parser.add_argument(
        "--seed", type=int, default=20260623,
        help="Random seed"
    )
    parser.add_argument(
        "--output_dir", type=str, default="data/bart_results",
        help="Output directory (default: data/bart_results)"
    )
    parser.add_argument(
        "--verbose", action="store_true", default=True
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # ---- Resolve input file ----
    if args.sample_file:
        sample_path = Path(args.sample_file)
    else:
        sample_path = Path(f"data/samples/sample_{args.sample:04d}.csv")

    if not sample_path.exists():
        print(f"ERROR: Sample file not found: {sample_path}", file=sys.stderr)
        sys.exit(1)

    print(f"\n{'='*60}")
    print(f"  Random Weight Survey BART")
    print(f"  Kim, Rao & Wang (2024) bootstrap  (B={args.B})")
    print(f"{'='*60}")
    print(f"  Input  : {sample_path}")
    print(f"  Trees  : m={args.m}")
    print(f"  MCMC   : {args.n_mcmc} iterations ({args.burn_in} burn-in)")
    print(f"  Seed   : {args.seed}")

    # ---- Load data ----
    df = pd.read_csv(sample_path)
    print(f"  n      : {len(df)} observations")
    print(f"  Strata : {df[args.strata_col].nunique()}")
    print(f"  PSUs   : {df[args.psu_col].nunique()}")

    # ---- Create output directory ----
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Stem for output files
    if args.sample_file:
        stem = Path(args.sample_file).stem
    else:
        stem = f"sample_{args.sample:04d}"

    out_stem = str(out_dir / stem)

    # ======================================================================
    # STEP A: Generate Kim-24 bootstrap weights (B replicates)
    # ======================================================================
    print(f"\n--- Generating Kim-24 bootstrap weights (B={args.B}) ---")
    W, oob_mask, npz_wpath, csv_wpath = generate_and_save_bootstrap_weights(
        df, B=args.B,
        output_path=out_stem,
        seed=args.seed,
        strata_col=args.strata_col,
        psu_col=args.psu_col,
        weight_col=args.weight_col
    )

    if args.weights_only:
        print("--weights_only flag set. Exiting after weight generation.")
        return

    # ======================================================================
    # STEP B: Prepare X, y matrices
    # ======================================================================
    # Validate feature columns
    available_feats = [f for f in args.features if f in df.columns]
    missing_feats   = [f for f in args.features if f not in df.columns]
    if missing_feats:
        print(f"  WARNING: Features not found in data: {missing_feats}")
    if not available_feats:
        print("ERROR: No valid feature columns found.", file=sys.stderr)
        sys.exit(1)

    # One-hot encode categorical columns
    df_enc = pd.get_dummies(df[available_feats], drop_first=True)
    X = df_enc.values.astype(float)
    y = df[args.outcome].values.astype(float)

    print(f"\n--- Feature matrix ---")
    print(f"  Features : {list(df_enc.columns)}")
    print(f"  X shape  : {X.shape}")
    print(f"  y range  : [{y.min():.2f}, {y.max():.2f}]")

    # Train/test split (stratified: sample within each stratum to preserve design)
    rng_split = default_rng(args.seed + 1)
    test_mask = np.zeros(len(df), dtype=bool)
    for h, idx in df.groupby(args.strata_col).groups.items():
        idx_arr   = np.array(idx)
        n_test_h  = max(1, int(len(idx_arr) * args.test_frac))
        test_idx  = rng_split.choice(idx_arr, size=n_test_h, replace=False)
        test_mask[test_idx] = True

    X_train, y_train = X[~test_mask], y[~test_mask]
    X_test,  y_test  = X[test_mask],  y[test_mask]
    df_train          = df[~test_mask].reset_index(drop=True)

    print(f"\n--- Train/test split ---")
    print(f"  Train : {X_train.shape[0]}  Test : {X_test.shape[0]}")

    # ======================================================================
    # STEP C: Fit Random Weight Survey BART
    # ======================================================================
    print(f"\n--- Fitting Random Weight Survey BART ---")
    prior = BARTPrior(m=args.m)
    model = SurveyBART(
        prior=prior,
        n_mcmc=args.n_mcmc,
        burn_in=args.burn_in,
        seed=args.seed + 2,
        verbose=args.verbose
    )

    t_start = time.time()
    model.register_test_set(X_test)   # collect posterior draws during MCMC
    model.fit(
        X_train, y_train, df_train,
        strata_col=args.strata_col,
        psu_col=args.psu_col,
        weight_col=args.weight_col
    )
    elapsed = time.time() - t_start

    # ======================================================================
    # STEP D: Posterior predictive inference on test set
    # ======================================================================
    print(f"\n--- Posterior predictive inference ---")
    preds = model.predict()          # uses draws collected during MCMC
    cov95 = model.coverage_95(y_test)
    rmse  = np.sqrt(np.mean((preds["mean"] - y_test) ** 2))
    mab   = np.mean(np.abs(preds["mean"] - y_test))

    print(f"\n{'='*60}")
    print(f"  Survey BART Results")
    print(f"{'='*60}")
    print(f"  RMSE (test)          : {rmse:.4f}")
    print(f"  MAE  (test)          : {mab:.4f}")
    print(f"  95% CI Coverage      : {100*cov95:.1f}%  (target: 95%)")
    print(f"  Mean 95% CI width    : {(preds['upper_95'] - preds['lower_95']).mean():.4f}")
    print(f"  Elapsed              : {elapsed:.1f}s")
    print(f"  Kim-24 OOB rate      : {100*model.oob_mask_draws_.mean():.1f}%")
    print(f"{'='*60}\n")

    # ======================================================================
    # STEP E: Save results
    # ======================================================================
    out_npz = out_stem + "_bart.npz"
    np.savez_compressed(
        out_npz,
        # Predictions
        y_pred_mean        = preds["mean"],
        y_pred_std         = preds["std"],
        y_pred_lower95     = preds["lower_95"],
        y_pred_upper95     = preds["upper_95"],
        y_test             = y_test,
        # In-sample posterior draws
        posterior_draws    = model.posterior_draws_,
        sigma2_draws       = model.sigma2_draws_,
        # Kim-24 weights stored per post-burn-in step
        weight_draws       = model.weight_draws_,
        oob_mask_draws     = model.oob_mask_draws_,
        # Pre-generated B=50 weight matrix
        W_bootstrap        = W,
        oob_mask_bootstrap = oob_mask,
        # RMSE trace
        rmse_train         = np.array(model.rmse_train_),
        # Metrics
        metrics = np.array([rmse, mab, cov95,
                             elapsed,
                             model.oob_mask_draws_.mean()]),
        # Metadata
        feature_names      = np.array(list(df_enc.columns)),
        n_train            = np.array(X_train.shape[0]),
        n_test             = np.array(X_test.shape[0]),
        B                  = np.array(args.B),
        n_mcmc             = np.array(args.n_mcmc),
        burn_in            = np.array(args.burn_in),
        m_trees            = np.array(args.m),
        seed               = np.array(args.seed),
    )

    # --- Summary CSV for easy aggregation across samples ---
    summary = pd.DataFrame([{
        "sample":       stem,
        "n_train":      X_train.shape[0],
        "n_test":       X_test.shape[0],
        "B":            args.B,
        "n_mcmc":       args.n_mcmc,
        "m_trees":      args.m,
        "rmse_test":    rmse,
        "mae_test":     mab,
        "coverage_95":  cov95,
        "ci_width_95":  (preds["upper_95"] - preds["lower_95"]).mean(),
        "oob_rate":     model.oob_mask_draws_.mean(),
        "elapsed_s":    elapsed,
    }])
    out_csv = out_stem + "_bart_summary.csv"
    summary.to_csv(out_csv, index=False)

    print(f"  Saved BART results   : {out_npz}")
    print(f"  Saved summary CSV    : {out_csv}")
    print(f"  Saved weights NPZ    : {npz_wpath}")
    print(f"  Saved weights CSV    : {csv_wpath}")


if __name__ == "__main__":
    main()
