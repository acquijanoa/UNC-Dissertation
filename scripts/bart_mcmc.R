# bart_mcmc.R — R wrappers for Kim-Rao bart_mcmc.cpp and Chipman BART::wbart

suppressPackageStartupMessages({
  library(Rcpp)
})

.cpp_loaded <- FALSE
ensure_bart_mcmc_cpp <- function() {
  if (.cpp_loaded) return(invisible(TRUE))
  src <- "scripts/src/bart_mcmc.cpp"
  if (!file.exists(src)) stop("Cannot find ", src, " (run from project root)")
  sourceCpp(src)
  .cpp_loaded <<- TRUE
  invisible(TRUE)
}

#' Scale y to [-0.5, 0.5] like Chipman BART (center + range scaling).
scale_y_chipman <- function(y) {
  y <- as.numeric(y)
  y_min <- min(y)
  y_max <- max(y)
  y_center <- (y_min + y_max) / 2
  y_range <- y_max - y_min
  if (y_range < 1e-12) y_range <- 1.0
  list(
    y_scaled = (y - y_center) / y_range,
    y_center = y_center,
    y_range = y_range
  )
}

#' Fit Kim-Rao weighted BART (one-step random weights each iteration).
#'
#' @param X numeric matrix (n x p)
#' @param y numeric response
#' @param weights design weights (length n)
#' @param strata integer strata ids
#' @param psu integer PSU ids
#' @param X_test optional test matrix (n_test x p)
#' @param K Kim-Rao scale (default 1)
#' @param normalize_weights if TRUE, normalize each bootstrap draw to sum n
fit_bart_mcmc <- function(
    X, y, weights, strata, psu,
    X_test = NULL,
    num_trees = 200,
    num_iter = 2000,
    burn_in = 500,
    min_leaf = 5,
    num_cutpoints = 100,
    alpha = 0.95,
    beta = 2.0,
    K = 1.0,
    normalize_weights = FALSE,
    return_weights = FALSE,
    verbose_every = 0,
    seed = NULL
) {
  ensure_bart_mcmc_cpp()
  if (!is.null(seed)) set.seed(seed)

  X <- as.matrix(X)
  storage.mode(X) <- "double"
  y <- as.numeric(y)
  weights <- as.numeric(weights)
  strata <- as.integer(strata)
  psu <- as.integer(psu)

  sc <- scale_y_chipman(y)
  # Chipman sigma_mu prior: 0.5/(k*sqrt(m)) with k=2
  sigma_mu <- 0.5 / (2 * sqrt(num_trees))

  X_test_arg <- NULL
  if (!is.null(X_test)) {
    X_test_arg <- as.matrix(X_test)
    storage.mode(X_test_arg) <- "double"
  }

  res <- run_bart_mcmc_cpp(
    X = X,
    y_scaled = sc$y_scaled,
    weights = weights,
    strata = strata,
    psu = psu,
    num_trees = as.integer(num_trees),
    num_iter = as.integer(num_iter),
    burn_in = as.integer(burn_in),
    min_leaf = as.integer(min_leaf),
    num_cutpoints = as.integer(num_cutpoints),
    alpha = alpha,
    beta = beta,
    sigma_mu = sigma_mu,
    sigma2_init = 1.0,
    nu = 3.0,
    lambda = 0.1,
    verbose_every = as.integer(verbose_every),
    y_center = sc$y_center,
    y_range = sc$y_range,
    K = as.numeric(K),
    normalize_weights = isTRUE(normalize_weights),
    X_test = X_test_arg,
    return_weights = isTRUE(return_weights)
  )

  list(
    pred_samples = res$pred_samples,
    pred_test_samples = res$pred_test_samples,
    weight_samples = res$weight_samples,
    sigma2_samples = res$sigma2_samples,
    accept_rates = res$accept_rates,
    y_center = res$y_center,
    y_range = res$y_range,
    K = res$K,
    normalize_weights = res$normalize_weights,
    n_mcmc = num_iter,
    burn_in = burn_in
  )
}

#' OOB-calibrate Kim-Rao scale K on training data (unnormalized weights).
#'
#' For each candidate K, run a short MCMC, then compute mean squared error on
#' units with bootstrap weight 0 (Kim-Rao leave-PSU-out draws).
calibrate_K_oob <- function(
    X, y, weights, strata, psu,
    K_candidates = c(1, 5, 15, 30, 50, 80),
    num_trees = 50,
    num_iter = 500,
    burn_in = 150,
    seed = NULL
) {
  cat("\n--- OOB calibration for K (unnormalized weights) ---\n")
  best_loss <- Inf
  best_K <- K_candidates[[1]]
  y <- as.numeric(y)

  for (K_cand in K_candidates) {
    fit <- fit_bart_mcmc(
      X = X, y = y, weights = weights, strata = strata, psu = psu,
      X_test = NULL,
      num_trees = num_trees,
      num_iter = num_iter,
      burn_in = burn_in,
      K = K_cand,
      normalize_weights = FALSE,
      return_weights = TRUE,
      seed = if (is.null(seed)) NULL else as.integer(seed + K_cand)
    )
    W <- fit$weight_samples
    Fmat <- fit$pred_samples
    oob <- (W == 0)
    if (!any(oob)) {
      cat(sprintf("  K = %6.1f | no OOB units — skip\n", K_cand))
      next
    }
    y_mat <- matrix(y, nrow = nrow(Fmat), ncol = length(y), byrow = TRUE)
    oob_loss <- mean((y_mat[oob] - Fmat[oob])^2)
    cat(sprintf("  K = %6.1f | OOB-MSE = %.4f\n", K_cand, oob_loss))
    if (!is.na(oob_loss) && oob_loss < best_loss) {
      best_loss <- oob_loss
      best_K <- K_cand
    }
  }
  cat(sprintf("--- Selected K* = %s (OOB-MSE = %.4f) ---\n\n", best_K, best_loss))
  list(K = best_K, oob_mse = best_loss)
}

#' Posterior predictive check on Kim-Rao OOB units (no population hold-out).
#'
#' For each post-burn-in draw t and unit i with w*_i(t)=0, form the pointwise
#' 95% predictive interval f_i(t) +/- z_{0.975} * sigma(t) and record whether
#' y_i falls inside. Pool over all (t,i) OOB pairs for `pooled_coverage`.
#' Also form a unit-level interval from all OOB predictive draws for that unit.
#'
#' @param fit result of fit_bart_mcmc(..., return_weights = TRUE)
#' @param y observed response (length n)
#' @param level predictive interval level (default 0.95)
compute_oob_pp_coverage <- function(fit, y, level = 0.95) {
  y <- as.numeric(y)
  W <- fit$weight_samples
  Fmat <- fit$pred_samples
  s2 <- as.numeric(fit$sigma2_samples)
  if (is.null(W) || is.null(Fmat) || is.null(s2)) {
    stop("fit must include weight_samples, pred_samples, and sigma2_samples ",
         "(call fit_bart_mcmc with return_weights = TRUE)")
  }
  n_keep <- nrow(Fmat)
  n <- length(y)
  if (nrow(W) != n_keep || ncol(W) != n || ncol(Fmat) != n || length(s2) != n_keep) {
    stop("Dimension mismatch among weight_samples, pred_samples, sigma2_samples, y")
  }

  z <- stats::qnorm(0.5 + level / 2)
  oob <- (W == 0)
  n_oob <- sum(oob)
  if (n_oob == 0L) {
    return(list(
      level = level,
      n_oob_pairs = 0L,
      pooled_coverage = NA_real_,
      unit_coverage = NA_real_,
      n_units_scored = 0L,
      oob_mse = NA_real_,
      mean_interval_width = NA_real_
    ))
  }

  sigma <- sqrt(pmax(s2, 0))
  # pooled: y in [f +/- z*sigma] for every OOB (t,i)
  covered <- logical(n_oob)
  widths <- numeric(n_oob)
  oob_idx <- which(oob, arr.ind = TRUE) # columns: row=t, col=i
  for (k in seq_len(n_oob)) {
    t <- oob_idx[k, 1L]
    i <- oob_idx[k, 2L]
    half <- z * sigma[t]
    lo <- Fmat[t, i] - half
    hi <- Fmat[t, i] + half
    covered[k] <- (y[i] >= lo && y[i] <= hi)
    widths[k] <- hi - lo
  }
  pooled_coverage <- mean(covered)

  y_mat <- matrix(y, nrow = n_keep, ncol = n, byrow = TRUE)
  oob_mse <- mean((y_mat[oob] - Fmat[oob])^2)

  # unit-level: among OOB draws for unit i, predictive draws f + N(0,s2);
  # empirical (level) quantile interval for that unit's predictive distribution
  alpha <- (1 - level) / 2
  unit_cover <- logical(n)
  unit_ok <- logical(n)
  for (i in seq_len(n)) {
    oob_t <- which(W[, i] == 0)
    if (length(oob_t) < 20L) next
    # one predictive replicate per OOB iteration (same RNG seed pattern not needed)
    y_rep <- Fmat[oob_t, i] + stats::rnorm(length(oob_t), 0, sigma[oob_t])
    q <- stats::quantile(y_rep, probs = c(alpha, 1 - alpha), names = FALSE)
    unit_cover[i] <- (y[i] >= q[1L] && y[i] <= q[2L])
    unit_ok[i] <- TRUE
  }
  n_units <- sum(unit_ok)
  unit_coverage <- if (n_units > 0L) mean(unit_cover[unit_ok]) else NA_real_

  list(
    level = level,
    n_oob_pairs = as.integer(n_oob),
    pooled_coverage = pooled_coverage,
    unit_coverage = unit_coverage,
    n_units_scored = as.integer(n_units),
    oob_mse = oob_mse,
    mean_interval_width = mean(widths)
  )
}

#' Calibrate K so OOB posterior predictive coverage of y is near `target`.
#'
#' Uses Method A-style inputs (caller supplies weights, typically weight_norm).
#' Selects K minimizing |pooled_coverage - target|.
calibrate_K_ppc <- function(
    X, y, weights, strata, psu,
    K_candidates = c(0.01, 0.05, 0.1, 0.2, 0.5, 1, 2),
    target = 0.95,
    level = 0.95,
    num_trees = 50,
    num_iter = 500,
    burn_in = 150,
    seed = NULL
) {
  cat(sprintf(
    "\n--- OOB PPC calibration for K (target coverage = %.1f%%) ---\n",
    100 * target
  ))
  y <- as.numeric(y)
  rows <- list()
  best_K <- K_candidates[[1]]
  best_gap <- Inf
  best_cov <- NA_real_

  for (K_cand in K_candidates) {
    fit <- fit_bart_mcmc(
      X = X, y = y, weights = weights, strata = strata, psu = psu,
      X_test = NULL,
      num_trees = num_trees,
      num_iter = num_iter,
      burn_in = burn_in,
      K = K_cand,
      normalize_weights = FALSE,
      return_weights = TRUE,
      seed = if (is.null(seed)) NULL else as.integer(seed + as.integer(1000 * K_cand))
    )
    ppc <- compute_oob_pp_coverage(fit, y, level = level)
    gap <- abs(ppc$pooled_coverage - target)
    cat(sprintf(
      "  K = %6.2f | OOB-PP-cov = %5.1f%% | unit-cov = %5.1f%% | OOB-MSE = %.4f | width = %.3f\n",
      K_cand,
      100 * ppc$pooled_coverage,
      100 * (if (is.na(ppc$unit_coverage)) NA_real_ else ppc$unit_coverage),
      ppc$oob_mse,
      ppc$mean_interval_width
    ))
    rows[[length(rows) + 1L]] <- data.frame(
      K = K_cand,
      pooled_coverage = ppc$pooled_coverage,
      unit_coverage = ppc$unit_coverage,
      oob_mse = ppc$oob_mse,
      mean_interval_width = ppc$mean_interval_width,
      n_oob_pairs = ppc$n_oob_pairs,
      gap = gap
    )
    if (!is.na(gap) && gap < best_gap) {
      best_gap <- gap
      best_K <- K_cand
      best_cov <- ppc$pooled_coverage
    }
  }

  grid <- do.call(rbind, rows)
  cat(sprintf(
    "--- Selected K* = %s (OOB-PP-cov = %.1f%%, |cov-target|=%.3f) ---\n\n",
    best_K, 100 * best_cov, best_gap
  ))
  list(K = best_K, pooled_coverage = best_cov, gap = best_gap, grid = grid)
}

#' Chipman unweighted BART via BART::wbart.
#'
#' Default is a pure non-design baseline (only X_f). Set include_design=TRUE
#' to append strata/psuid as covariates (not the preferred contrast for
#' informative-sampling studies).
fit_chipman_wbart <- function(
    X_f, y,
    strata = NULL, psu = NULL,
    X_f_test = NULL,
    strata_test = NULL,
    psu_test = NULL,
    include_design = FALSE,
    num_trees = 200,
    num_iter = 2000,
    burn_in = 500,
    seed = NULL
) {
  if (!requireNamespace("BART", quietly = TRUE)) {
    stop("Package BART is required for the Chipman baseline. install.packages('BART')")
  }
  if (!is.null(seed)) set.seed(seed)

  X_f <- as.matrix(X_f)
  storage.mode(X_f) <- "double"
  if (isTRUE(include_design)) {
    if (is.null(strata) || is.null(psu))
      stop("strata and psu required when include_design=TRUE")
    x_train <- cbind(X_f, strata = as.numeric(strata), psuid = as.numeric(psu))
  } else {
    x_train <- X_f
  }

  x_test <- NULL
  if (!is.null(X_f_test)) {
    X_f_test <- as.matrix(X_f_test)
    storage.mode(X_f_test) <- "double"
    if (isTRUE(include_design)) {
      if (is.null(strata_test) || is.null(psu_test))
        stop("strata_test and psu_test required when include_design=TRUE")
      x_test <- cbind(
        X_f_test,
        strata = as.numeric(strata_test),
        psuid = as.numeric(psu_test)
      )
    } else {
      x_test <- X_f_test
    }
  }

  ndpost <- as.integer(num_iter - burn_in)
  fit <- BART::wbart(
    x.train = x_train,
    y.train = as.numeric(y),
    x.test = if (is.null(x_test)) matrix(0, 0, ncol(x_train)) else x_test,
    ntree = as.integer(num_trees),
    nskip = as.integer(burn_in),
    ndpost = ndpost,
    printevery = as.integer(max(num_iter + 1, 1))
  )

  pred_test <- NULL
  if (!is.null(x_test) && !is.null(fit$yhat.test)) {
    pred_test <- fit$yhat.test
  }

  list(
    pred_samples = fit$yhat.train,
    pred_test_samples = pred_test,
    sigma2_samples = (fit$sigma[(burn_in + 1):(burn_in + ndpost)])^2,
    n_mcmc = num_iter,
    burn_in = burn_in
  )
}
