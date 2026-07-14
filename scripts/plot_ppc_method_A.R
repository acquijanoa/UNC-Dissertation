# plot_ppc_method_A.R
# Posterior predictive check plots for Method A (OOB, no population hold-out).
# Target: 95% predictive coverage of y.
#
# Produces a multi-panel PNG/PDF:
#   1) Density of y vs y_rep on OOB occasions
#   2) PIT histogram (uniform if well calibrated)
#   3) Unit intervals: observed y vs 95% OOB predictive intervals
#   4) K-grid OOB PP coverage with horizontal 95% target (if --mode calibrate)
#
# Example:
#   Rscript scripts/plot_ppc_method_A.R \
#     --single_sample data/samples/sample_0001.csv \
#     --out_plot figures/ppc_method_A_sample0001.png \
#     --mode calibrate --K 1.0

suppressPackageStartupMessages({
  library(data.table)
})

source("scripts/bart_mcmc.R")

args_raw <- commandArgs(trailingOnly = TRUE)
args <- list(
  single_sample = "data/samples/sample_0001.csv",
  out_plot = "figures/ppc_method_A.png",
  mode = "calibrate",
  K = 1.0,
  K_grid = "0.01,0.05,0.1,0.2,0.5,1,2",
  response = "y_f",
  level = 0.95,
  n_mcmc = 2000L,
  burn_in = 500L,
  n_trees = 200L,
  calib_trees = 50L,
  calib_iter = 500L,
  calib_burn = 150L,
  n_rep_dens = 40L,
  n_units_show = 80L,
  seed = 42L
)

i <- 1L
while (i <= length(args_raw)) {
  key <- args_raw[i]
  val <- if (i < length(args_raw)) args_raw[i + 1L] else NA_character_
  if (key == "--single_sample") { args$single_sample <- val; i <- i + 2L }
  else if (key == "--out_plot") { args$out_plot <- val; i <- i + 2L }
  else if (key == "--mode") { args$mode <- val; i <- i + 2L }
  else if (key == "--K") { args$K <- as.numeric(val); i <- i + 2L }
  else if (key == "--K_grid") { args$K_grid <- val; i <- i + 2L }
  else if (key == "--response") { args$response <- val; i <- i + 2L }
  else if (key == "--level") { args$level <- as.numeric(val); i <- i + 2L }
  else if (key == "--n_mcmc") { args$n_mcmc <- as.integer(val); i <- i + 2L }
  else if (key == "--burn_in") { args$burn_in <- as.integer(val); i <- i + 2L }
  else if (key == "--n_trees") { args$n_trees <- as.integer(val); i <- i + 2L }
  else if (key == "--seed") { args$seed <- as.integer(val); i <- i + 2L }
  else { i <- i + 1L }
}

MODE <- match.arg(args$mode, c("check", "calibrate"))
set.seed(args$seed)

FEATURE_COLS <- c("x_f1", "x_f2", "x_f3", "x_f4", "x_f5")
df <- fread(args$single_sample)
X <- as.matrix(df[, ..FEATURE_COLS])
y <- as.numeric(df[[args$response]])
w <- as.numeric(df[["weight_norm"]])
strata <- as.integer(df[["strata"]])
psu <- as.integer(df[["psuid"]])
n <- length(y)

K_used <- args$K
grid_dt <- NULL

cat("=== Method A PPC plots ===\n")
if (MODE == "calibrate") {
  K_grid <- as.numeric(strsplit(args$K_grid, ",")[[1]])
  cal <- calibrate_K_ppc(
    X = X, y = y, weights = w, strata = strata, psu = psu,
    K_candidates = K_grid,
    target = args$level,
    level = args$level,
    num_trees = args$calib_trees,
    num_iter = args$calib_iter,
    burn_in = args$calib_burn,
    seed = args$seed
  )
  K_used <- cal$K
  grid_dt <- as.data.table(cal$grid)
}

cat(sprintf("Fitting Method A at K=%.4g with return_weights...\n", K_used))
fit <- fit_bart_mcmc(
  X = X, y = y, weights = w, strata = strata, psu = psu,
  X_test = NULL,
  num_trees = args$n_trees,
  num_iter = args$n_mcmc,
  burn_in = args$burn_in,
  K = K_used,
  normalize_weights = FALSE,
  return_weights = TRUE,
  seed = args$seed + 7L
)

ppc <- compute_oob_pp_coverage(fit, y, level = args$level)
cat(sprintf(
  "OOB PP coverage = %.1f%% (target %.0f%%) | unit = %.1f%% | pairs = %d\n",
  100 * ppc$pooled_coverage, 100 * args$level,
  100 * ppc$unit_coverage, ppc$n_oob_pairs
))

W <- fit$weight_samples
Fmat <- fit$pred_samples
s2 <- as.numeric(fit$sigma2_samples)
sigma <- sqrt(pmax(s2, 0))
oob <- (W == 0)
n_keep <- nrow(Fmat)
z <- qnorm(0.5 + args$level / 2)

# ---- build plotting quantities ----
# PIT on OOB pairs: Phi( (y - f) / sigma )
oob_idx <- which(oob, arr.ind = TRUE)
pit <- numeric(nrow(oob_idx))
y_obs_oob <- numeric(nrow(oob_idx))
y_rep_list <- vector("list", min(args$n_rep_dens, n_keep))
# pick evenly spaced iterations for density overlay
rep_iters <- unique(as.integer(round(seq(1, n_keep, length.out = args$n_rep_dens))))

for (k in seq_len(nrow(oob_idx))) {
  t <- oob_idx[k, 1L]
  i <- oob_idx[k, 2L]
  pit[k] <- pnorm(y[i], mean = Fmat[t, i], sd = sigma[t])
  y_obs_oob[k] <- y[i]
}

# dens replicates: for each selected iter, all OOB y_rep at that iter
y_rep_mat <- matrix(NA_real_, nrow = length(rep_iters), ncol = n)
for (j in seq_along(rep_iters)) {
  t <- rep_iters[j]
  oob_i <- which(W[t, ] == 0)
  if (!length(oob_i)) next
  y_rep_mat[j, oob_i] <- Fmat[t, oob_i] + rnorm(length(oob_i), 0, sigma[t])
}

# unit-level intervals for display
unit_lo <- rep(NA_real_, n)
unit_hi <- rep(NA_real_, n)
unit_mid <- rep(NA_real_, n)
unit_cover <- rep(NA, n)
alpha <- (1 - args$level) / 2
for (i in seq_len(n)) {
  oob_t <- which(W[, i] == 0)
  if (length(oob_t) < 20L) next
  y_rep_i <- Fmat[oob_t, i] + rnorm(length(oob_t), 0, sigma[oob_t])
  q <- quantile(y_rep_i, probs = c(alpha, 0.5, 1 - alpha), names = FALSE)
  unit_lo[i] <- q[1]
  unit_mid[i] <- q[2]
  unit_hi[i] <- q[3]
  unit_cover[i] <- (y[i] >= q[1] && y[i] <= q[3])
}
ok <- which(!is.na(unit_lo))
# show a spread of units ordered by observed y
ok <- ok[order(y[ok])]
if (length(ok) > args$n_units_show) {
  keep <- as.integer(round(seq(1, length(ok), length.out = args$n_units_show)))
  ok <- ok[keep]
}

# ---- plot ----
dir.create(dirname(args$out_plot), showWarnings = FALSE, recursive = TRUE)

draw_ppc_panels <- function() {
  n_panels <- if (!is.null(grid_dt)) 4L else 3L
  layout_mat <- if (n_panels == 4L) {
    matrix(c(1, 2, 3, 3, 4, 4), nrow = 2, byrow = TRUE)
  } else {
    matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE)
  }
  layout(layout_mat)
  par(mar = c(4.2, 4.2, 3.2, 1.2), cex = 0.95)

  col_obs <- "#1B4F72"
  col_rep <- adjustcolor("#5DADE2", 0.35)
  col_target <- "#C0392B"
  col_good <- "#1E8449"
  col_bad <- "#922B21"

  # Panel 1: density y vs y_rep (OOB)
  yr <- range(c(y, as.vector(y_rep_mat)), na.rm = TRUE)
  d_obs <- density(y)
  ymax <- max(d_obs$y)
  rep_dens <- list()
  for (j in seq_len(nrow(y_rep_mat))) {
    v <- y_rep_mat[j, ]
    v <- v[!is.na(v)]
    if (length(v) < 30L) next
    d <- density(v, from = yr[1], to = yr[2])
    ymax <- max(ymax, d$y)
    rep_dens[[length(rep_dens) + 1L]] <- d
  }
  plot(NA, xlim = yr, ylim = c(0, 1.05 * ymax), type = "n",
       xlab = expression(y), ylab = "Density",
       main = sprintf("PPC densities (OOB) | Method A, K=%.3g", K_used))
  for (d in rep_dens) lines(d, col = col_rep, lwd = 1)
  lines(d_obs, col = col_obs, lwd = 2.5)
  legend("topright",
         legend = c(expression(y[obs]), expression(y[rep]^{"(OOB)"})),
         col = c(col_obs, "#5DADE2"), lwd = c(2.5, 1.5), bty = "n")

  # Panel 2: PIT histogram
  hist(pit, breaks = seq(0, 1, by = 0.05), col = adjustcolor(col_obs, 0.55),
       border = "white", main = "OOB PIT (uniform if calibrated)",
       xlab = expression(Phi((y - f)/sigma)), ylab = "Frequency",
       xlim = c(0, 1))
  abline(h = length(pit) / 20, col = col_target, lwd = 2, lty = 2)
  legend("topright", legend = "Uniform reference", col = col_target,
         lwd = 2, lty = 2, bty = "n", cex = 0.85)

  # Panel 3: unit 95% predictive intervals
  plot(NA, xlim = c(0.5, length(ok) + 0.5),
       ylim = range(c(unit_lo[ok], unit_hi[ok], y[ok]), na.rm = TRUE),
       xlab = "Units (ordered by observed y)", ylab = expression(y),
       main = sprintf(
         "OOB 95%% predictive intervals | coverage = %.1f%% (target 95%%)",
         100 * mean(unit_cover[ok], na.rm = TRUE)
       ))
  for (j in seq_along(ok)) {
    i <- ok[j]
    col_i <- if (isTRUE(unit_cover[i])) col_good else col_bad
    segments(j, unit_lo[i], j, unit_hi[i], col = adjustcolor(col_i, 0.7), lwd = 1.2)
    points(j, unit_mid[i], pch = 15, cex = 0.35, col = adjustcolor(col_i, 0.9))
    points(j, y[i], pch = 16, cex = 0.55, col = col_obs)
  }
  abline(h = mean(y), col = "grey50", lty = 3)
  legend("topleft",
         legend = c("Observed y", "Covered", "Missed", "95% PI"),
         pch = c(16, 15, 15, NA),
         col = c(col_obs, col_good, col_bad, "grey40"),
         lwd = c(NA, NA, NA, 1.5),
         bty = "n", cex = 0.85)

  # Panel 4: K calibration (optional)
  if (!is.null(grid_dt) && nrow(grid_dt) > 0L) {
    plot(grid_dt$K, 100 * grid_dt$pooled_coverage, type = "b",
         pch = 19, col = col_obs, lwd = 2, log = "x",
         xlab = "K (log scale)", ylab = "OOB PP coverage of y (%)",
         main = "Method A: calibrate K to 95% OOB coverage of y",
         ylim = c(min(80, 100 * min(grid_dt$pooled_coverage, na.rm = TRUE) - 2),
                  max(100, 100 * max(grid_dt$pooled_coverage, na.rm = TRUE) + 2)))
    abline(h = 100 * args$level, col = col_target, lwd = 2, lty = 2)
    points(K_used, 100 * ppc$pooled_coverage, pch = 21, bg = "#F4D03F",
           col = col_target, cex = 1.8, lwd = 2)
    text(K_used, 100 * ppc$pooled_coverage,
         labels = sprintf(" K*=%.3g\n %.1f%%", K_used, 100 * ppc$pooled_coverage),
         pos = 4, cex = 0.9)
    legend("bottomleft",
           legend = c("OOB PP coverage", "95% target", "Selected K* (full fit)"),
           col = c(col_obs, col_target, col_target),
           pch = c(19, NA, 21), pt.bg = c(NA, NA, "#F4D03F"),
           lwd = c(2, 2, NA), lty = c(1, 2, NA), bty = "n", cex = 0.85)
  }

  mtext(sprintf(
    "Method A (weight_norm, Kim-Rao) | n=%d | trees=%d | MCMC=%d | pooled OOB-PP-cov=%.1f%%",
    n, args$n_trees, args$n_mcmc, 100 * ppc$pooled_coverage
  ), side = 1, outer = TRUE, line = -1, cex = 0.75, col = "grey30")
}

# write primary path
ext <- tolower(tools::file_ext(args$out_plot))
if (ext == "pdf") {
  pdf(args$out_plot, width = 11, height = 8.5)
} else {
  png(args$out_plot, width = 1400, height = 1100, res = 140)
}
draw_ppc_panels()
dev.off()
cat("Plot saved to:", args$out_plot, "\n")

# also write the other format
alt <- if (ext == "pdf") {
  sub("\\.pdf$", ".png", args$out_plot, ignore.case = TRUE)
} else {
  sub("\\.png$", ".pdf", args$out_plot, ignore.case = TRUE)
}
if (identical(alt, args$out_plot) || !nzchar(tools::file_ext(alt))) {
  alt <- paste0(tools::file_path_sans_ext(args$out_plot), ".pdf")
}
if (tolower(tools::file_ext(alt)) == "pdf") {
  pdf(alt, width = 11, height = 8.5)
} else {
  png(alt, width = 1400, height = 1100, res = 140)
}
draw_ppc_panels()
dev.off()
cat("Also saved:", alt, "\n")

metrics_file <- sub("\\.(png|pdf)$", "_metrics.csv", args$out_plot, ignore.case = TRUE)
fwrite(data.frame(
  method = "A",
  K = K_used,
  oob_pp_coverage = ppc$pooled_coverage,
  oob_unit_coverage = ppc$unit_coverage,
  n_oob_pairs = ppc$n_oob_pairs,
  oob_mse = ppc$oob_mse,
  mean_pp_width = ppc$mean_interval_width,
  sample = args$single_sample,
  out_plot = args$out_plot
), metrics_file)
cat("Metrics saved to:", metrics_file, "\n")
