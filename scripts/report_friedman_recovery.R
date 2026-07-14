#!/usr/bin/env Rscript
# report_friedman_recovery.R
# Monte Carlo report: bias, relative bias, empirical SE, estimated SE, coverage
#
# Estimand (per replication, held-out test set):
#   theta     = mean_i mu_f(x_i)
#   theta_hat = mean_i E[f(x_i) | data]
#   bias_r    = theta_hat - theta
#   est_se_r  = posterior SD of mean_i f(x_i)
#   cover_r   = 1{ theta in 95% posterior CI of mean_i f }

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
in_dir <- if (length(args) >= 1) args[[1]] else "data/friedman_results"
out_csv <- if (length(args) >= 2) args[[2]] else file.path(in_dir, "report_table.csv")
out_md  <- if (length(args) >= 3) args[[3]] else file.path(in_dir, "report.md")

files <- list.files(in_dir, pattern = "^result_.*\\.csv$", full.names = TRUE)
if (!length(files)) stop("No result_*.csv files in ", in_dir)

dt <- rbindlist(lapply(files, fread), fill = TRUE)
need <- c("method", "theta", "theta_hat", "bias", "est_se", "cover_mean")
miss <- setdiff(need, names(dt))
if (length(miss)) {
  stop("Results missing columns ", paste(miss, collapse = ", "),
       " — re-run simulation_friedman_recovery.R with the updated metrics.")
}

# Group by method and K if K column exists
by_cols <- if ("K" %in% names(dt)) c("method", "K") else "method"

report <- dt[, {
  bias <- mean(bias)
  theta_bar <- mean(theta)
  .(
    n_rep = .N,
    truth_mean = theta_bar,
    est_mean = mean(theta_hat),
    bias = bias,
    relative_bias = bias / theta_bar,
    emp_se = sd(theta_hat),
    emp_se_bias = sd(bias),
    est_se = mean(est_se),
    se_ratio = mean(est_se) / sd(theta_hat),
    coverage = mean(cover_mean),
    rmse_mu = mean(rmse_mu),
    pointwise_f_coverage = mean(f_coverage),
    pp_coverage = mean(pp_coverage),
    mean_elapsed_sec = mean(elapsed_sec)
  )
}, by = by_cols]

if ("K" %in% names(dt)) {
  setorder(report, method, K)
} else {
  setorder(report, method)
}
fwrite(report, out_csv)

fmt_pct <- function(x) sprintf("%.1f%%", 100 * x)
fmt4 <- function(x) sprintf("%.4f", x)

lines <- c(
  "# Friedman recovery simulation report",
  "",
  "## Estimand",
  "",
  "For each replication, a held-out test set of size `n_test` is drawn from the",
  "population (units not in the survey sample). Define",
  "",
  "- $\\theta = n_{\\mathrm{test}}^{-1}\\sum_i \\mu_f(x_i)$ (true mean regression on the test set)",
  "- $\\hat\\theta = n_{\\mathrm{test}}^{-1}\\sum_i \\widehat{E}[f(x_i)\\mid \\mathrm{data}]$ (posterior mean of $f$)",
  "- Estimated SE: posterior SD of $n_{\\mathrm{test}}^{-1}\\sum_i f^{(s)}(x_i)$",
  "- Coverage: fraction of replications whose 95% posterior interval for $\\theta$ contains $\\theta$",
  "",
  "Methods:",
  "",
  "- **A**: Kim–Rao BART with `weight_norm`, no per-draw $\\sum w^*=n$ (so $K$ scales Term 2)",
  "- **B**: Kim–Rao BART with raw `weight` (original scale), no $\\sum w^*=n`",
  "- **sum1**: Kim–Rao BART with $w_i/\\sum_j w_j$ ($\\sum w=1$), no per-draw renorm",
  "- **chipman**: `BART::wbart` unweighted, covariates = Friedman `x_f` only (no design)",
  "",
  "## Monte Carlo summary",
  ""
)

# Add table header based on whether K column exists
if ("K" %in% names(report)) {
  lines <- c(lines,
    "| Method | K | Reps | Bias | Rel. bias | Emp. SE | Est. SE | SE ratio | Coverage | RMSE($\\mu_f$) | Pointwise cov. |",
    "|--------|--:|-----:|-----:|----------:|--------:|--------:|---------:|---------:|---------------:|---------------:|"
  )
} else {
  lines <- c(lines,
    "| Method | Reps | Bias | Rel. bias | Emp. SE | Est. SE | SE ratio | Coverage | RMSE($\\mu_f$) | Pointwise cov. |",
    "|--------|-----:|-----:|----------:|--------:|--------:|---------:|---------:|---------------:|---------------:|"
  )
}

for (i in seq_len(nrow(report))) {
  r <- report[i]
  if ("K" %in% names(report)) {
    lines <- c(lines, sprintf(
      "| %s | %.1f | %d | %s | %s | %s | %s | %s | %s | %s | %s |",
      r$method, r$K, r$n_rep,
      fmt4(r$bias), fmt_pct(r$relative_bias),
      fmt4(r$emp_se), fmt4(r$est_se), fmt4(r$se_ratio),
      fmt_pct(r$coverage),
      fmt4(r$rmse_mu), fmt_pct(r$pointwise_f_coverage)
    ))
  } else {
    lines <- c(lines, sprintf(
      "| %s | %d | %s | %s | %s | %s | %s | %s | %s | %s |",
      r$method, r$n_rep,
      fmt4(r$bias), fmt_pct(r$relative_bias),
      fmt4(r$emp_se), fmt4(r$est_se), fmt4(r$se_ratio),
      fmt_pct(r$coverage),
      fmt4(r$rmse_mu), fmt_pct(r$pointwise_f_coverage)
    ))
  }
}

lines <- c(
  lines,
  "",
  "## Column notes",
  "",
  "- **Bias**: $\\mathrm{mean}_r(\\hat\\theta_r - \\theta_r)$",
  "- **Rel. bias**: bias / mean($\\theta_r$)",
  "- **Emp. SE**: $\\mathrm{sd}_r(\\hat\\theta_r)$",
  "- **Est. SE**: mean over reps of the posterior SD of $\\hat\\theta$",
  "- **SE ratio**: Est. SE / Emp. SE (near 1 = well calibrated uncertainty)",
  "- **Coverage**: Monte Carlo coverage of the 95% interval for $\\theta$",
  "- **RMSE($\\mu_f$)**: mean over reps of $\\sqrt{\\mathrm{mean}_i(\\bar f_i-\\mu_{f,i})^2}$ (pointwise prediction error)",
  "- **Pointwise cov.**: mean over reps of pointwise 95% coverage of $\\mu_f$",
  "",
  sprintf("Generated from %d result files in `%s`.", length(files), in_dir)
)

writeLines(lines, out_md)
cat("Wrote", out_csv, "\n")
cat("Wrote", out_md, "\n")
print(report)
