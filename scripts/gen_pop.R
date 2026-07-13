##
## Population generation — Cross-sectional
##
##  Based on Kim & Rao (2024) Supplemental S9 design (stratified two-stage PPS)
##
## Author: Álvaro Quijano
##
## Creation date: 17 jun 26

# Set library
library(data.table)

# Set parameters
set.seed(170626)

design <- list(
  H  = 30, # number of strata
  N1 = 50, # minimum clusters per stratum (floor)
  N2 = 20 # minimum individuals per cluster (floor)
)

# Latent class parameters
# Class 1: low        Class 2: moderate        Class 3: high
K <- 3L
cp <- list(
  pi      = c(0.50, 0.30, 0.20), # mixing proportions
  alpha   = c(5.0, 10.0, 15.0), # class intercepts
  beta1   = c(0.8, 1.0, 1.2), # slope on x_inv
  sigma_u = c(0.3, 0.5, 0.7) # within-class individual SD
)

# Outcome model coefficients
coef_pop <- list(
  delta_cont2 = 0.50, # x_cont_2
  delta_bin = 0.80, # x_bin
  delta_cat = c(0, 0.30, 0.60, 0.90), # categorical (cat 1 = reference)
  delta_int_cont2_bin = -0.40, # interaction: x_cont_2 * x_bin
  delta_int_inv_bin = -0.20 # interaction: x_inv * x_bin
)

# Initialize the dataset with strata information
dt <- data.table(
  strata   = 1:design$H,
  lambda_h = rexp(design$H, rate = 1) # stratum effect
)

# Number of clusters per stratum
dt[, n_clust := 8L * rpois(.N, lambda = lambda_h) + design$N1]

# Expand the dataset to include clusters
dt <- dt[rep(1:.N, n_clust)][, n_clust := NULL]

# Initialize PSUs and cluster-level effects
dt[, psuid := seq_len(.N), by = strata]
dt[, eta_i := rexp(rate = 1, n = 1), by = .(strata, psuid)]

# Generate cluster sizes (informative PPS)
dt[, M_hi := 8L * rpois(.N, lambda = lambda_h + eta_i) + design$N2]

# Stratum total: denominator for PPS inclusion probability
dt[, N_h := sum(M_hi), by = strata]

# Expand to individuals
dt <- dt[rep(1:.N, M_hi)]
dt[, subid := seq_len(.N), by = .(strata, psuid)]

# Generate latent class and class-specific random intercept
dt[, class := sample.int(K, size = .N, replace = TRUE, prob = cp$pi)]
dt[, u_ij := rnorm(.N, mean = 0, sd = cp$sigma_u[class])]

# Generate covariates

# Continuous 1: x_inv ~ Uniform(0, 20)
dt[, x_inv := runif(.N, min = 0, max = 20)]

# Continuous 2: x_cont_2 ~ N(0, 1)
dt[, x_cont_2 := rnorm(.N, mean = 0, sd = 1)]

# Binary: x_bin ~ Bernoulli(p(x_cont_2))
dt[, x_bin := rbinom(.N, size = 1L, prob = plogis(-0.3 + 0.4 * x_cont_2))]

# Categorical: x_cat ∈ {1, 2, 3, 4}  [proportions: 30 / 30 / 20 / 20 %]
dt[, x_cat := sample.int(4L,
  size = .N, replace = TRUE,
  prob = c(0.30, 0.30, 0.20, 0.20)
)]

# Friedman #1 inputs (independent of the survey DGP covariates)
# Classic MARS benchmark: x_fk ~ Uniform(0, 1), k = 1,...,5
dt[, `:=`(
  x_f1 = runif(.N),
  x_f2 = runif(.N),
  x_f3 = runif(.N),
  x_f4 = runif(.N),
  x_f5 = runif(.N)
)]

# Generate residual
dt[, eps := rnorm(.N)]

# Generate outcome (informative PPS / latent-class DGP)
dt[, y := cp$alpha[class] +
  cp$beta1[class] * x_inv +
  coef_pop$delta_cont2 * x_cont_2 +
  coef_pop$delta_bin * x_bin +
  coef_pop$delta_cat[x_cat] +
  coef_pop$delta_int_cont2_bin * (x_cont_2 * x_bin) +
  coef_pop$delta_int_inv_bin * (x_inv * x_bin) +
  2.0 * eta_i +
  2.0 * lambda_h +
  u_ij + eps]

# Friedman #1 at the population level (Friedman, 1991), plus the same
# stratum/cluster effects as y (informative PPS linkage):
#   mu_f = 10 sin(π x1 x2) + 20 (x3 - 0.5)^2 + 10 x4 + 5 x5
#          + 2 λ_h + 2 η_i
dt[, mu_f := 10 * sin(pi * x_f1 * x_f2) +
  20 * (x_f3 - 0.5)^2 +
  10 * x_f4 +
  5 * x_f5 +
  2.0 * lambda_h +
  2.0 * eta_i]
dt[, eps_f := rnorm(.N)]
dt[, y_f := mu_f + eps_f]

# Final key and column order
setkey(dt, strata, psuid, subid)

setcolorder(dt, c(
  # --- identifiers ---
  "strata", "psuid", "subid",
  # --- design effects (informative PPS) ---
  "lambda_h", "eta_i", "M_hi", "N_h",
  # --- covariates (survey DGP) ---
  "x_inv", "x_cont_2", "x_bin", "x_cat",
  # --- covariates (Friedman #1) ---
  "x_f1", "x_f2", "x_f3", "x_f4", "x_f5",
  # --- latent class and random effects (truth) ---
  "class", "u_ij", "eps", "eps_f",
  # --- outcomes ---
  "y", "mu_f", "y_f"
))


# Save population
data_dir <- here::here("data")
fwrite(dt, file = file.path(data_dir, "population.csv"))

message(
  "Population saved to ", data_dir,
  "\n  rows : ", nrow(dt),
  "\n  cols : ", ncol(dt),
  "\n  files: population.csv"
)

# Column guide
#
# IDENTIFIERS
#   strata   : stratum index h = 1, ..., H
#   psuid    : cluster index i within stratum
#   subid    : individual index j within cluster
#
# DESIGN (informative PPS — do not pass to analysis models)
#   lambda_h : stratum Exp(1) effect → 0.3*lambda_h in y; drives M_h
#   eta_i    : cluster Exp(1) effect → 0.3*eta_i in y; drives M_hi
#   M_hi     : cluster size (PPS selection prob: p_hi = M_hi / N_h)
#   N_h      : stratum total individuals (PPS weight denominator)
#
# COVARIATES (survey DGP)
#   x_inv    : Uniform(0, 20)        — continuous
#   x_cont_2 : N(0, 1)              — continuous
#   x_bin    : Bernoulli(plogis(-0.3 + 0.4*x_cont_2))  — binary
#   x_cat    : Categorical {1,2,3,4} probs (0.30, 0.30, 0.20, 0.20)
#
# COVARIATES (Friedman #1 — independent U(0,1) inputs)
#   x_f1..x_f5 : Uniform(0, 1)
#
# LATENT / RANDOM EFFECTS (truth — do not pass to analysis models)
#   class    : ground-truth latent class z_i ∈ {1, 2, 3}
#   u_ij     : N(0, sigma_u[k]^2) class-specific individual random intercept
#   eps      : N(0, 1) residual for y
#   eps_f    : N(0, 1) residual for y_f
#
# OUTCOMES
#   y        : informative PPS / latent-class DGP outcome
#   mu_f     : Friedman #1 mean + 2*lambda_h + 2*eta_i (no noise; for scoring)
#   y_f      : mu_f + eps_f  (observed Friedman outcome)
#
