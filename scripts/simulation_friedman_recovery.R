# simulation_friedman_recovery.R
# Held-out Friedman recovery under INFORMATIONAL PPS: y_f ~ x_f1..x_f5, score vs mu_f
# where mu_f includes design effects (correlated with selection).
#
# Methods:
#   A       — Kim-Rao BART, weight_norm (mean=1, Σw=n), no per-draw Σw*=n
#             so K scales Term 2; E[Σw*]≈n already from weight_norm
#   B       — Kim-Rao BART, raw weight (original scale), no Σw*=n
#   sum1    — Kim-Rao BART, w/Σw (weights sum to 1), no per-draw renorm
#   chipman — BART::wbart unweighted, covariates = x_f only (no design IDs)
#
# Scalar estimand per replication (held-out test set):
#   theta   = mean_i mu_f(x_i)
#   theta_hat = mean_i E[f(x_i) | data]
# Saved so a Monte Carlo report can form bias, rbias, empSE, estSE, coverage.

suppressPackageStartupMessages({
  library(data.table)
})

source("scripts/bart_mcmc.R")

args_raw <- commandArgs(trailingOnly = TRUE)
args <- list(
  single_sample = NULL,
  out_file = NULL,
  method = "A",
  K = 1.0,
  n_test = 2000L,
  n_mcmc = 2000L,
  burn_in = 500L,
  n_trees = 200L,
  seed = 42L
)

i <- 1L
while (i <= length(args_raw)) {
  key <- args_raw[i]
  val <- if (i < length(args_raw)) args_raw[i + 1L] else NA_character_
  if (key == "--single_sample") { args$single_sample <- val; i <- i + 2L }
  else if (key == "--out_file") { args$out_file <- val; i <- i + 2L }
  else if (key == "--method") { args$method <- val; i <- i + 2L }
  else if (key == "--K") { args$K <- as.numeric(val); i <- i + 2L }
  else if (key == "--n_test") { args$n_test <- as.integer(val); i <- i + 2L }
  else if (key == "--n_mcmc") { args$n_mcmc <- as.integer(val); i <- i + 2L }
  else if (key == "--burn_in") { args$burn_in <- as.integer(val); i <- i + 2L }
  else if (key == "--n_trees") { args$n_trees <- as.integer(val); i <- i + 2L }
  else if (key == "--seed") { args$seed <- as.integer(val); i <- i + 2L }
  else { i <- i + 1L }
}

METHOD <- match.arg(args$method, c("A", "B", "sum1", "chipman"))
set.seed(args$seed)

FEATURE_COLS <- c("x_f1", "x_f2", "x_f3", "x_f4", "x_f5")
TARGET_COL <- "y_f"
TRUTH_COL <- "mu_f"
STRATA_COL <- "strata"
PSU_COL <- "psuid"
KEY_COLS <- c("strata", "psuid", "subid")

if (is.null(args$single_sample)) stop("Must provide --single_sample path")
if (is.null(args$out_file)) stop("Must provide --out_file path")

pop <- fread("data/population.csv")
need_pop <- c(KEY_COLS, FEATURE_COLS, TARGET_COL, TRUTH_COL, STRATA_COL, PSU_COL)
miss <- setdiff(need_pop, names(pop))
if (length(miss)) stop("population.csv missing columns: ", paste(miss, collapse = ", "))

pop_keys <- paste(pop[[STRATA_COL]], pop[[PSU_COL]], pop[["subid"]], sep = "_")

build_test_set <- function(pop, df_train, n_test, pop_keys) {
  train_keys <- paste(df_train[[STRATA_COL]], df_train[[PSU_COL]], df_train[["subid"]], sep = "_")
  pool <- pop[!(pop_keys %in% train_keys)]
  n_draw <- min(as.integer(n_test), nrow(pool))
  pool[sample.int(nrow(pool), n_draw)]
}

summarize_draws <- function(f_draws, mu_true, y_obs, sigma2_draws) {
  # f_draws: n_keep x n_test
  post_mean <- colMeans(f_draws)
  rmse_mu <- sqrt(mean((post_mean - mu_true)^2))
  mcmc_sd_f <- mean(apply(f_draws, 2L, sd))

  f_lower <- apply(f_draws, 2L, quantile, probs = 0.025)
  f_upper <- apply(f_draws, 2L, quantile, probs = 0.975)
  f_coverage <- mean(mu_true >= f_lower & mu_true <= f_upper)
  f_ci_width <- mean(f_upper - f_lower)

  n_keep <- nrow(f_draws)
  n_test <- ncol(f_draws)
  pp <- matrix(0, nrow = n_keep, ncol = n_test)
  for (k in seq_len(n_keep)) {
    pp[k, ] <- f_draws[k, ] + rnorm(n_test, 0, sqrt(sigma2_draws[k]))
  }
  pp_lower <- apply(pp, 2L, quantile, probs = 0.025)
  pp_upper <- apply(pp, 2L, quantile, probs = 0.975)
  pp_coverage <- mean(y_obs >= pp_lower & y_obs <= pp_upper)
  pp_width <- mean(pp_upper - pp_lower)

  # Scalar estimand: mean of the regression function on the test set
  theta <- mean(mu_true)
  mean_draws <- rowMeans(f_draws) # posterior of mean_i f(x_i)
  theta_hat <- mean(mean_draws)
  bias <- theta_hat - theta
  est_se <- stats::sd(mean_draws)
  mean_lower <- as.numeric(stats::quantile(mean_draws, 0.025))
  mean_upper <- as.numeric(stats::quantile(mean_draws, 0.975))
  cover_mean <- as.integer(theta >= mean_lower && theta <= mean_upper)

  list(
    rmse_mu = rmse_mu,
    mcmc_sd_f = mcmc_sd_f,
    f_coverage = f_coverage,
    f_ci_width = f_ci_width,
    pp_coverage = pp_coverage,
    pp_width = pp_width,
    theta = theta,
    theta_hat = theta_hat,
    bias = bias,
    est_se = est_se,
    mean_lower = mean_lower,
    mean_upper = mean_upper,
    cover_mean = cover_mean
  )
}

df_train <- fread(args$single_sample)
need_samp <- c(KEY_COLS, FEATURE_COLS, TARGET_COL, STRATA_COL, PSU_COL, "weight", "weight_norm")
miss_s <- setdiff(need_samp, names(df_train))
if (length(miss_s)) stop("sample missing columns: ", paste(miss_s, collapse = ", "))

df_test <- build_test_set(pop, df_train, args$n_test, pop_keys)

X_train <- as.matrix(df_train[, ..FEATURE_COLS])
X_test <- as.matrix(df_test[, ..FEATURE_COLS])
y_train <- df_train[[TARGET_COL]]
y_test <- df_test[[TARGET_COL]]
mu_test <- df_test[[TRUTH_COL]]

rep_num <- as.integer(gsub(".*sample_([0-9]+)\\.csv$", "\\1", args$single_sample))
if (is.na(rep_num)) rep_num <- args$seed

cat(sprintf(
  "=== Friedman recovery | method=%s | K=%.1f | rep=%d | n=%d | n_test=%d ===\n",
  METHOD, args$K, rep_num, nrow(df_train), nrow(df_test)
))
cat(sprintf("MCMC: trees=%d iter=%d burn=%d seed=%d\n",
            args$n_trees, args$n_mcmc, args$burn_in, args$seed))

t0 <- proc.time()[[3]]

if (METHOD == "A") {
  # weight_norm (Σw=n, mean=1) + Kim-Rao bootstrap, NO per-draw re-normalize.
  # Leaving normalize_weights=FALSE lets K inflate/deflate Term 2; with
  # normalize_weights=TRUE, K cancels in the sum-to-n rescaling.
  fit <- fit_bart_mcmc(
    X = X_train, y = y_train,
    weights = df_train[["weight_norm"]],
    strata = df_train[[STRATA_COL]],
    psu = df_train[[PSU_COL]],
    X_test = X_test,
    num_trees = args$n_trees,
    num_iter = args$n_mcmc,
    burn_in = args$burn_in,
    K = args$K,
    normalize_weights = FALSE,
    seed = args$seed
  )
} else if (METHOD == "B") {
  # Original-scale sampling weights + Kim-Rao bootstrap, no Σw*=n
  fit <- fit_bart_mcmc(
    X = X_train, y = y_train,
    weights = df_train[["weight"]],
    strata = df_train[[STRATA_COL]],
    psu = df_train[[PSU_COL]],
    X_test = X_test,
    num_trees = args$n_trees,
    num_iter = args$n_mcmc,
    burn_in = args$burn_in,
    K = args$K,
    normalize_weights = FALSE,
    seed = args$seed
  )
} else if (METHOD == "sum1") {
  # Unit-sum calibrated weights: w_i / Σ_j w_j  (Σw = 1), not mean-1.
  # No per-draw renorm so the Σw*=1 scale (and K) is preserved.
  w_sum1 <- as.numeric(df_train[["weight"]])
  w_sum1 <- w_sum1 / sum(w_sum1)
  fit <- fit_bart_mcmc(
    X = X_train, y = y_train,
    weights = w_sum1,
    strata = df_train[[STRATA_COL]],
    psu = df_train[[PSU_COL]],
    X_test = X_test,
    num_trees = args$n_trees,
    num_iter = args$n_mcmc,
    burn_in = args$burn_in,
    K = args$K,
    normalize_weights = FALSE,
    seed = args$seed
  )
} else {
  fit <- fit_chipman_wbart(
    X_f = X_train, y = y_train,
    X_f_test = X_test,
    include_design = FALSE,
    num_trees = args$n_trees,
    num_iter = args$n_mcmc,
    burn_in = args$burn_in,
    seed = args$seed
  )
}

elapsed <- proc.time()[[3]] - t0
f_draws <- fit$pred_test_samples
if (is.null(f_draws)) stop("No test predictions returned for method ", METHOD)

met <- summarize_draws(f_draws, mu_test, y_test, fit$sigma2_samples)

out <- data.frame(
  method = METHOD,
  rep = rep_num,
  K = args$K,
  theta = met$theta,
  theta_hat = met$theta_hat,
  bias = met$bias,
  est_se = met$est_se,
  mean_lower = met$mean_lower,
  mean_upper = met$mean_upper,
  cover_mean = met$cover_mean,
  rmse_mu = met$rmse_mu,
  mcmc_sd_f = met$mcmc_sd_f,
  f_coverage = met$f_coverage,
  f_ci_width = met$f_ci_width,
  pp_coverage = met$pp_coverage,
  pp_width = met$pp_width,
  n = nrow(df_train),
  n_test = nrow(df_test),
  n_trees = args$n_trees,
  n_mcmc = args$n_mcmc,
  burn_in = args$burn_in,
  seed = args$seed,
  elapsed_sec = elapsed
)

dir.create(dirname(args$out_file), showWarnings = FALSE, recursive = TRUE)
fwrite(out, args$out_file)

cat(sprintf(
  "  method=%s | bias=%.4f | estSE=%.4f | cover=%d | RMSE=%.4f | PP-cov=%.1f%% | %.1fs\n",
  METHOD, met$bias, met$est_se, met$cover_mean,
  met$rmse_mu, 100 * met$pp_coverage, elapsed
))
cat("Results saved to:", args$out_file, "\n")
