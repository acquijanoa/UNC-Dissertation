#!/usr/bin/env python3
"""
Visualize a Single Survey-Weighted Decision Tree from the Survey Random Forest.

This script loads the first sample, draws a set of Kim-24 bootstrap weights,
fits a shallow decision tree (max_depth=3) using the bootstrap weights,
and saves the tree visualization as a PNG.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from numpy.random import default_rng
import matplotlib.pyplot as plt
from sklearn.tree import DecisionTreeRegressor, plot_tree

# Add the directory containing this script to sys.path
script_dir = Path(__file__).parent.resolve()
sys.path.append(str(script_dir))

try:
    from survey_bart import generate_kim24_weights
except ImportError:
    print("ERROR: Could not import generate_kim24_weights. Ensure survey_bart.py is in the same directory.", file=sys.stderr)
    sys.exit(1)


def main():
    # ---- Load sample data ----
    sample_path = Path("data/samples/sample_0001.csv")
    if not sample_path.exists():
        # Fallback if run from scripts folder
        sample_path = Path("../data/samples/sample_0001.csv")
        
    if not sample_path.exists():
        print(f"ERROR: Sample file not found: sample_0001.csv", file=sys.stderr)
        sys.exit(1)

    print(f"Loading sample data for visualization: {sample_path}")
    df = pd.read_csv(sample_path)

    # Induce MAR missingness at 10% just to split into respondents
    rng = default_rng(20260623)
    x_c2 = df["x_cont_2"].values
    x_b = df["x_bin"].values
    
    # Simple logit threshold for 10% missingness
    p_resp = 1.0 / (1.0 + np.exp(-(0.5 + 0.5 * x_c2 - 0.5 * x_b)))
    R = rng.binomial(1, p_resp).astype(bool)

    # Encode features
    features = ["x_inv", "x_cont_2", "x_bin", "x_cat"]
    available_feats = [f for f in features if f in df.columns]
    df_enc = pd.get_dummies(df[available_feats], drop_first=True)
    feature_names = list(df_enc.columns)
    
    X_resp = df_enc.values[R]
    y_resp = df["y"].values[R]

    # Generate one set of Kim-24 bootstrap weights
    w_star = generate_kim24_weights(
        df, strata_col="strata", psu_col="psuid", weight_col="weight", rng=rng
    )
    w_star_resp = w_star[R]

    # Filter to in-bag units (weight > 0)
    in_bag = w_star_resp > 0
    X_train = X_resp[in_bag]
    y_train = y_resp[in_bag]
    w_train = w_star_resp[in_bag]

    # Fit a shallow tree for visualization
    print("Fitting a shallow survey-weighted decision tree (max_depth=3)...")
    tree = DecisionTreeRegressor(max_depth=3, max_features="sqrt", random_state=42)
    tree.fit(X_train, y_train, sample_weight=w_train)

    # Plot the tree
    fig, ax = plt.subplots(figsize=(16, 10), dpi=300)
    
    # Premium-looking color scheme
    plot_tree(
        tree,
        feature_names=feature_names,
        filled=True,
        rounded=True,
        ax=ax,
        fontsize=10,
        precision=3
    )
    
    ax.set_title(
        "Single Survey-Weighted Decision Tree (Kim-24 Bootstrap Replicate)\n"
        "Node splits are chosen using design-weighted criterion; colors indicate predicted outcome value (y)",
        fontsize=14,
        fontweight="bold",
        pad=20
    )

    # Ensure output directory exists
    out_dir = Path("data/rf_results")
    if not out_dir.exists():
        out_dir = Path("../data/rf_results")
    out_dir.mkdir(parents=True, exist_ok=True)
    
    out_png = out_dir / "survey_tree_visualization.png"
    plt.tight_layout()
    plt.savefig(out_png, bbox_inches="tight")
    plt.close()
    
    print(f"Tree visualization successfully saved to: {out_png.resolve()}")


if __name__ == "__main__":
    main()
