# report_friedman_oob_K.R
# Generate Monte Carlo report for OOB K-calibrated Friedman recovery
# Compares OOB-selected K* to grid search results

suppressPackageStartupMessages({
  library(data.table)
})

args_raw <- commandArgs(trailingOnly = TRUE)
args <- list(
  results_dir = "data/friedman_results_oob",
  out_csv = NULL,
  out_md = NULL
)

i <- 1L
while (i <= length(args_raw)) {
  key <- args_raw[i]
  val <- if (i < length(args_raw)) args_raw[i + 1L] else NA_character_
  if (key == "--results_dir") { args$results_dir <- val; i <- i + 2L }
  else if (key == "--out_csv") { args$out_csv <- val; i <- i + 2L }
  else if (key == "--out_md") { args$out_md <- val; i <- i + 2L }
  else { i <- i + 1L }
}

if (is.null(args$out_csv)) args$out_csv <- file.path(args$results_dir, "report_table.csv")
if (is.null(args$out_md)) args$out_md <- file.path(args$results_dir, "report.md")

cat("=== OOB K-Calibration Report ===\n")
cat(sprintf("Results directory: %s\n", args$results_dir))

csv_files <- list.files(args$results_dir, pattern = "\\.csv$", full.names = TRUE)
csv_files <- csv_files[!grepl("report_table\\.csv$", csv_files)]

if (length(csv_files) == 0) {
  stop("No CSV result files found in ", args$results_dir)
}

cat(sprintf("Found %d result files\n", length(csv_files)))

dt <- rbindlist(lapply(csv_files, fread), fill = TRUE)
cat(sprintf("Loaded %d rows\n", nrow(dt)))

if (!"method" %in% names(dt)) dt[, method := "unknown"]
if (!"K" %in% names(dt)) stop("No K column found")

dt[, `:=`(
  bias = as.numeric(bias),
  est_se = as.numeric(est_se),
  cover_mean = as.integer(cover_mean),
  rmse_mu = as.numeric(rmse_mu),
  pp_coverage = as.numeric(pp_coverage),
  f_coverage = as.numeric(f_coverage)
)]

# Compute Monte Carlo summaries by method
summary_dt <- dt[, .(
  n_reps = .N,
  mean_K = mean(K, na.rm = TRUE),
  sd_K = sd(K, na.rm = TRUE),
  bias = mean(bias, na.rm = TRUE),
  rel_bias_pct = 100 * mean(bias / theta, na.rm = TRUE),
  emp_se = sd(theta_hat, na.rm = TRUE),
  est_se = mean(est_se, na.rm = TRUE),
  se_ratio = mean(est_se, na.rm = TRUE) / sd(theta_hat, na.rm = TRUE),
  cover_theta = 100 * mean(cover_mean, na.rm = TRUE),
  rmse_mu = mean(rmse_mu, na.rm = TRUE),
  pp_coverage = 100 * mean(pp_coverage, na.rm = TRUE),
  pointwise_coverage = 100 * mean(f_coverage, na.rm = TRUE),
  mean_oob_mse = mean(K_oob_mse, na.rm = TRUE)
), by = method]

setorder(summary_dt, method)

cat("\n=== Monte Carlo Summary (OOB K-Calibrated) ===\n")
print(summary_dt)

fwrite(summary_dt, args$out_csv)
cat(sprintf("\nCSV report saved to: %s\n", args$out_csv))

# Generate markdown report
md_lines <- c(
  "# Friedman Mean-Surface Recovery: OOB K-Calibration Report",
  "",
  sprintf("**Generated:** %s", Sys.time()),
  sprintf("**Results directory:** `%s`", args$results_dir),
  sprintf("**Total replications:** %d", nrow(dt)),
  "",
  "## Summary",
  "",
  "This report compares OOB-calibrated K selection to grid search.",
  "",
  "### OOB K-Calibration Results",
  ""
)

# Add summary table
header <- sprintf("| %s |", paste(names(summary_dt), collapse = " | "))
sep <- sprintf("|%s|", paste(rep("---:", ncol(summary_dt)), collapse = "|"))
md_lines <- c(md_lines, header, sep)

for (i in seq_len(nrow(summary_dt))) {
  row_vals <- sapply(summary_dt[i, ], function(x) {
    if (is.numeric(x)) sprintf("%.3f", x) else as.character(x)
  })
  md_lines <- c(md_lines, sprintf("| %s |", paste(row_vals, collapse = " | ")))
}

md_lines <- c(md_lines, "", "## Interpretation", "")

# Add method descriptions
md_lines <- c(md_lines,
  "**Methods:**",
  "- **A_OOB**: Kim-Rao BART with `weight_norm`, OOB-selected K*",
  "- **B_OOB**: Kim-Rao BART with raw weights, OOB-selected K*",
  "",
  "**Key Metrics:**",
  "- `mean_K`: Average K selected across replications via OOB loss",
  "- `sd_K`: Standard deviation of selected K (stability measure)",
  "- `cover_theta`: Coverage of mean(mu_f) — target 95%",
  "- `pp_coverage`: Posterior predictive coverage of y_f — target 95%",
  "- `se_ratio`: Est. SE / Emp. SE — target ≈1 for calibrated intervals",
  "- `mean_oob_mse`: Average OOB MSE at selected K*",
  "",
  "**Comparison to Grid Search:**",
  "- Method A grid: K=0.1 gave 94.9% coverage, RMSE=2.948",
  "- Method B grid: K=1.0 gave 86.1% coverage, RMSE=2.912",
  "",
  "**Expected OOB behavior:**",
  "- OOB loss optimizes predictive accuracy (RMSE)",
  "- May not align with coverage-optimal K",
  "- Section 7.3 showed RMSE-optimal K ≠ coverage-optimal K for Method A",
  ""
)

md_text <- paste(md_lines, collapse = "\n")
writeLines(md_text, args$out_md)
cat(sprintf("Markdown report saved to: %s\n", args$out_md))

# Print K distribution statistics
cat("\n=== K Selection Distribution ===\n")
k_stats <- dt[, .(
  min_K = min(K, na.rm = TRUE),
  q25_K = quantile(K, 0.25, na.rm = TRUE),
  median_K = median(K, na.rm = TRUE),
  q75_K = quantile(K, 0.75, na.rm = TRUE),
  max_K = max(K, na.rm = TRUE),
  mean_K = mean(K, na.rm = TRUE),
  sd_K = sd(K, na.rm = TRUE)
), by = method]
print(k_stats)

cat("\nReport generation complete.\n")
