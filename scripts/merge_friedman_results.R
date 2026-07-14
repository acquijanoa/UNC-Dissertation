#!/usr/bin/env Rscript
# merge_friedman_results.R — bind per-task CSVs into one summary table
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
in_dir <- if (length(args) >= 1) args[[1]] else "data/friedman_results"
out_file <- if (length(args) >= 2) args[[2]] else file.path(in_dir, "summary.csv")

files <- list.files(in_dir, pattern = "^result_.*\\.csv$", full.names = TRUE)
if (!length(files)) stop("No result_*.csv files found in ", in_dir)

dt <- rbindlist(lapply(files, fread), fill = TRUE)
fwrite(dt, out_file)
cat(sprintf("Wrote %d rows from %d files -> %s\n", nrow(dt), length(files), out_file))

# Group by method and K if K column exists
by_cols <- if ("K" %in% names(dt)) c("method", "K") else "method"

print(dt[, .(
  n = .N,
  bias = if ("bias" %in% names(dt)) mean(bias) else NA_real_,
  est_se = if ("est_se" %in% names(dt)) mean(est_se) else NA_real_,
  coverage = if ("cover_mean" %in% names(dt)) mean(cover_mean) else NA_real_,
  rmse_mu = mean(rmse_mu),
  f_coverage = mean(f_coverage),
  pp_coverage = mean(pp_coverage),
  elapsed_sec = mean(elapsed_sec)
), by = by_cols])
