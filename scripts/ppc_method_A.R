# ppc_method_A.R
# Posterior Predictive Check for Method A (weight_norm + Kim-Rao, no per-draw renorm)
# Uses Kim-Rao OOB units only — no population hold-out required.
#
# Modes:
#   --mode check     : fit at fixed K, report OOB 95% PP coverage of y
#   --mode calibrate : choose K minimizing |OOB-PP-cov - 0.95|, then full-fit check
#
# Example:
#   Rscript scripts/ppc_method_A.R \
#     --single_sample data/samples/sample_0001.csv \
#     --out_file data/friedman_ppc_A/ppc_0001.csv \
#     --mode calibrate --K_grid 0.01,0.1,0.5,1,2

suppressPackageStartupMessages({
  library(data.table)
})

source("scripts/bart_mcmc.R")

args_raw <- commandArgs(trailingOnly = TRUE)
args <- list(
  single_sample = NULL,
  out_file = NULL,
  mode = "check",
  K = 1.0,
  K_grid = "0.01,0.05,0.1,0.2,0.5,1,2",
  target = 0.95,
  level = 0.95,
  response = "y_f",
  n_mcmc = 2000L,
  burn_in = 500L,
  n_trees = 200L,
  calib_trees = 50L,
  calib_iter = 500L,
  calib_burn = 150L,
  seed = 42L
)

i <- 1L
while (i <= length(args_raw)) {
  key <- args_raw[i]
  val <- if (i < length(args_raw)) args_raw[i + 1L] else NA_character_
  if (key == "--single_sample") { args$single_sample <- val; i <- i + 2L }
  else if (key == "--out_file") { args$out_file <- val; i <- i + 2L }
  else if (key == "--mode") { args$mode <- val; i <- i + 2L }
  else if (key == "--K") { args$K <- as.numeric(val); i <- i + 2L }
  else if (key == "--K_grid") { args$K_grid <- val; i <- i + 2L }
  else if (key == "--target") { args$target <- as.numeric(val); i <- i + 2L }
  else if (key == "--level") { args$level <- as.numeric(val); i <- i + 2L }
  else if (key == "--response") { args$response <- val; i <- i + 2L }
  else if (key == "--n_mcmc") { args$n_mcmc <- as.integer(val); i <- i + 2L }
  else if (key == "--burn_in") { args$burn_in <- as.integer(val); i <- i + 2L }
  else if (key == "--n_trees") { args$n_trees <- as.integer(val); i <- i + 2L }
  else if (key == "--calib_trees") { args$calib_trees <- as.integer(val); i <- i + 2L }
  else if (key == "--calib_iter") { args$calib_iter <- as.integer(val); i <- i + 2L }
  else if (key == "--calib_burn") { args$calib_burn <- as.integer(val); i <- i + 2L }
  else if (key == "--seed") { args$seed <- as.integer(val); i <- i + 2L }
  else { i <- i + 1L }
}

MODE <- match.arg(args$mode, c("check", "calibrate"))
if (is.null(args$single_sample)) stop("Must provide --single_sample")
if (is.null(args$out_file)) stop("Must provide --out_file")

FEATURE_COLS <- c("x_f1", "x_f2", "x_f3", "x_f4", "x_f5")
STRATA_COL <- "strata"
PSU_COL <- "psuid"

set.seed(args$seed)
df <- fread(args$single_sample)
need <- c(FEATURE_COLS, args$response, STRATA_COL, PSU_COL, "weight_norm")
miss <- setdiff(need, names(df))
if (length(miss)) stop("sample missing columns: ", paste(miss, collapse = ", "))

X <- as.matrix(df[, ..FEATURE_COLS])
y <- as.numeric(df[[args$response]])
w <- as.numeric(df[["weight_norm"]])
strata <- as.integer(df[[STRATA_COL]])
psu <- as.integer(df[[PSU_COL]])

rep_num <- as.integer(gsub(".*sample_([0-9]+)\\.csv$", "\\1", args$single_sample))
if (is.na(rep_num)) rep_num <- args$seed

K_used <- args$K
calib_gap <- NA_real_
calib_cov <- NA_real_

cat(sprintf(
  "=== Method A OOB PPC | mode=%s | response=%s | rep=%d | n=%d | target=%.0f%% ===\n",
  MODE, args$response, rep_num, length(y), 100 * args$target
))

t0 <- proc.time()[[3]]

if (MODE == "calibrate") {
  K_grid <- as.numeric(strsplit(args$K_grid, ",")[[1]])
  cal <- calibrate_K_ppc(
    X = X, y = y, weights = w, strata = strata, psu = psu,
    K_candidates = K_grid,
    target = args$target,
    level = args$level,
    num_trees = args$calib_trees,
    num_iter = args$calib_iter,
    burn_in = args$calib_burn,
    seed = args$seed
  )
  K_used <- cal$K
  calib_gap <- cal$gap
  calib_cov <- cal$pooled_coverage
  # save calibration grid beside main result if possible
  grid_file <- sub("\\.csv$", "_calib_grid.csv", args$out_file)
  dir.create(dirname(grid_file), showWarnings = FALSE, recursive = TRUE)
  fwrite(cal$grid, grid_file)
}

cat(sprintf("Full Method A fit at K=%.4g (return_weights=TRUE for OOB PPC)...\n", K_used))
fit <- fit_bart_mcmc(
  X = X, y = y, weights = w, strata = strata, psu = psu,
  X_test = NULL,
  num_trees = args$n_trees,
  num_iter = args$n_mcmc,
  burn_in = args$burn_in,
  K = K_used,
  normalize_weights = FALSE,
  return_weights = TRUE,
  seed = args$seed + 1000L
)

ppc <- compute_oob_pp_coverage(fit, y, level = args$level)
elapsed <- proc.time()[[3]] - t0

out <- data.frame(
  method = "A",
  mode = MODE,
  response = args$response,
  rep = rep_num,
  K = K_used,
  level = args$level,
  target_coverage = args$target,
  oob_pp_coverage = ppc$pooled_coverage,
  oob_unit_coverage = ppc$unit_coverage,
  n_oob_pairs = ppc$n_oob_pairs,
  n_units_scored = ppc$n_units_scored,
  oob_mse = ppc$oob_mse,
  mean_pp_width = ppc$mean_interval_width,
  calib_coverage = calib_cov,
  calib_gap = calib_gap,
  n = length(y),
  n_trees = args$n_trees,
  n_mcmc = args$n_mcmc,
  burn_in = args$burn_in,
  seed = args$seed,
  elapsed_sec = elapsed
)

dir.create(dirname(args$out_file), showWarnings = FALSE, recursive = TRUE)
fwrite(out, args$out_file)

cat(sprintf(
  "  K=%.4g | OOB 95%% PP-cov(y)=%.1f%% | unit-cov=%.1f%% | OOB-MSE=%.3f | width=%.3f | %.1fs\n",
  K_used,
  100 * ppc$pooled_coverage,
  100 * (if (is.na(ppc$unit_coverage)) NA_real_ else ppc$unit_coverage),
  ppc$oob_mse,
  ppc$mean_interval_width,
  elapsed
))
if (!is.na(ppc$pooled_coverage)) {
  gap <- abs(ppc$pooled_coverage - args$target)
  status <- if (gap <= 0.02) "NEAR TARGET (+/-2pp)" else if (gap <= 0.05) "ACCEPTABLE (+/-5pp)" else "OFF TARGET"
  cat(sprintf("  PPC status vs %.0f%%: %s (|gap|=%.1fpp)\n",
              100 * args$target, status, 100 * gap))
}
cat("Results saved to:", args$out_file, "\n")
