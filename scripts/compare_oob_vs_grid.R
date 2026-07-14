# compare_oob_vs_grid.R
# Compare OOB-selected K results to grid search results

suppressPackageStartupMessages({
  library(data.table)
})

cat("=== Comparing OOB K-Selection vs Grid Search ===\n\n")

# Load OOB results
oob_file <- "data/friedman_results_oob/report_table.csv"
if (!file.exists(oob_file)) {
  stop("OOB results not found: ", oob_file)
}
dt_oob <- fread(oob_file)
dt_oob[, source := "OOB"]

# Load grid results (A2 K-grid)
grid_file <- "data/friedman_K_grid_A2/report_table.csv"
if (!file.exists(grid_file)) {
  cat("Warning: Grid results not found at ", grid_file, "\n")
  cat("Showing OOB results only.\n\n")
  print(dt_oob)
  quit(save = "no", status = 0)
}
dt_grid <- fread(grid_file)
dt_grid[, source := "Grid"]

# Standardize method names for comparison
dt_oob[method == "A_OOB", method_base := "A"]
dt_oob[method == "B_OOB", method_base := "B"]
dt_grid[, method_base := method]

# Key metrics to compare
compare_cols <- c("method_base", "source", "mean_K", "bias", "rel_bias_pct", 
                  "emp_se", "est_se", "se_ratio", "cover_theta", 
                  "rmse_mu", "pp_coverage")

dt_compare <- rbind(
  dt_oob[, ..compare_cols],
  dt_grid[, ..compare_cols],
  fill = TRUE
)

setorder(dt_compare, method_base, source)

cat("Side-by-Side Comparison:\n")
print(dt_compare)

# Highlight key findings
cat("\n=== Key Findings ===\n\n")

for (meth in unique(dt_compare$method_base)) {
  cat(sprintf("** Method %s **\n", meth))
  
  oob_row <- dt_compare[method_base == meth & source == "OOB"]
  grid_rows <- dt_compare[method_base == meth & source == "Grid"]
  
  if (nrow(oob_row) == 0) {
    cat("  No OOB results\n\n")
    next
  }
  
  K_oob <- oob_row$mean_K
  cat(sprintf("  OOB selected K* = %.3f\n", K_oob))
  cat(sprintf("    Coverage(θ) = %.1f%%, RMSE = %.3f, PP(y) = %.1f%%\n",
              oob_row$cover_theta, oob_row$rmse_mu, oob_row$pp_coverage))
  
  if (nrow(grid_rows) > 0) {
    # Find closest grid K to OOB K
    grid_rows[, K_diff := abs(mean_K - K_oob)]
    closest <- grid_rows[which.min(K_diff)]
    
    cat(sprintf("  Closest grid K = %.3f\n", closest$mean_K))
    cat(sprintf("    Coverage(θ) = %.1f%%, RMSE = %.3f, PP(y) = %.1f%%\n",
                closest$cover_theta, closest$rmse_mu, closest$pp_coverage))
    
    # Find optimal coverage and optimal RMSE from grid
    best_cov_idx <- which.min(abs(grid_rows$cover_theta - 95))
    best_rmse_idx <- which.min(grid_rows$rmse_mu)
    
    best_cov <- grid_rows[best_cov_idx]
    best_rmse <- grid_rows[best_rmse_idx]
    
    cat(sprintf("  Grid K for best coverage = %.3f (%.1f%% coverage)\n",
                best_cov$mean_K, best_cov$cover_theta))
    cat(sprintf("  Grid K for best RMSE = %.3f (RMSE = %.3f)\n",
                best_rmse$mean_K, best_rmse$rmse_mu))
    
    if (abs(K_oob - best_rmse$mean_K) < 0.2) {
      cat("  → OOB selected RMSE-optimal K (as expected)\n")
    } else if (abs(K_oob - best_cov$mean_K) < 0.2) {
      cat("  → OOB selected coverage-optimal K (surprising!)\n")
    } else {
      cat("  → OOB selected intermediate K\n")
    }
  }
  cat("\n")
}

cat("=== Interpretation ===\n\n")
cat("OOB loss minimizes predictive MSE, not coverage of θ.\n")
cat("Section 7.3 showed RMSE-optimal K ≠ coverage-optimal K for Method A:\n")
cat("  - K=0.5: minimum RMSE (2.894), coverage 87.8%\n")
cat("  - K=0.1: nominal coverage (94.9%), RMSE 2.948\n")
cat("\nIf OOB selects K closer to 0.5, it confirms the RMSE vs coverage trade-off.\n")
cat("For imputation (predict y_i), OOB-selected K is appropriate.\n")
cat("For mean estimation (E[μ]), coverage-optimal K from grid is preferred.\n")
