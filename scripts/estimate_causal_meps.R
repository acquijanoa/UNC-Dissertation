# estimate_causal_meps.R
library(data.table)

source("scripts/survey_bart_rcpp.R")

args_raw <- commandArgs(trailingOnly=TRUE)
args <- list(
    n_mcmc = 1000,
    burn_in = 500,
    n_trees = 200,
    outcome = "totexp21",
    chain_id = 1
)

i <- 1
while(i <= length(args_raw)) {
    if (args_raw[i] == "--n_mcmc") { args$n_mcmc <- as.integer(args_raw[i+1]); i <- i+2 }
    else if (args_raw[i] == "--burn_in") { args$burn_in <- as.integer(args_raw[i+1]); i <- i+2 }
    else if (args_raw[i] == "--n_trees") { args$n_trees <- as.integer(args_raw[i+1]); i <- i+2 }
    else if (args_raw[i] == "--outcome") { args$outcome <- args_raw[i+1]; i <- i+2 }
    else if (args_raw[i] == "--chain_id") { args$chain_id <- as.integer(args_raw[i+1]); i <- i+2 }
    else { i <- i+1 }
}

# Seed based on chain to ensure independent exploration
set.seed(20260700 + args$chain_id)

# Load MEPS data
df <- fread("data/meps_2021_hc233.csv")

# Preprocessing
cat_cols <- c('sex', 'racethx', 'povcat21', 'rthlth42', 'mnhlth42', 'insurc21')
cat_cols <- intersect(cat_cols, colnames(df))
for (col in cat_cols) {
    df <- df[get(col) >= 0]
    df[[col]] <- as.factor(df[[col]])
}

dx_cols <- c('hibpdx', 'chddx', 'angidx', 'midx', 'ohrtdx', 'strkdx', 
             'emphdx', 'chbrondx', 'choldx', 'cancerdx', 'diabdx', 'diabdx_m18',
             'jntwkdx', 'arthdx', 'asthdx', 'adhdaddx', 'pregnt31')
dx_cols <- intersect(dx_cols, colnames(df))
for (col in dx_cols) {
    df[[col]] <- as.integer(df[[col]] == 1)
}

df <- df[haveus42 >= 0]
df$Z <- as.integer(df$haveus42 == 1)

if (args$outcome == "totslf21") {
    df$Y <- df$totslf21
} else {
    df$Y <- df$totexp21
}

# Formula
formula_str <- paste("Y ~ Z + age21x +", paste(dx_cols, collapse = " + "), "+", paste(cat_cols, collapse = " + "))
f <- as.formula(formula_str)

# Extract model matrix to create counterfactuals
mf <- model.frame(f, as.data.frame(df), na.action = na.pass)
X_train <- model.matrix(f, mf)
if("(Intercept)" %in% colnames(X_train)) X_train <- X_train[, colnames(X_train) != "(Intercept)"]

# Counterfactuals
X_treated <- X_train
X_treated[, "Z"] <- 1
X_control <- X_train
X_control[, "Z"] <- 0
X_test <- rbind(X_treated, X_control)

# Calibration of K (One-step weight adjustment)
w_orig <- df$perwt21f
n <- nrow(df)
K_used <- max(1.0, n / (mean(w_orig) * 5))
w_scaled <- w_orig * K_used

cat(sprintf("Fitting Two-Part Causal Survey BART (n=%d, K=%.2f, outcome=%s, chain=%d)...\n", n, K_used, args$outcome, args$chain_id))

fit <- survey_bart(
    formula = f,
    data = as.data.frame(df),
    family = "hurdle",
    weights = w_scaled,
    strata = df$varstr,
    psu = df$varpsu,
    n_trees = args$n_trees,
    n_mcmc = args$n_mcmc,
    burn_in = args$burn_in,
    X_test = X_test
)

# Extract predictions
n_post <- args$n_mcmc - args$burn_in
y_hat_treated <- fit$y_hat_test[, 1:n]
y_hat_control <- fit$y_hat_test[, (n+1):(2*n)]

causal_diff <- y_hat_treated - y_hat_control

# SATE
sate_draws <- rowMeans(causal_diff)

# CATEs by Insurance
insurc_idx <- as.integer(as.character(df$insurc21))
cate_priv_draws <- rowMeans(causal_diff[, insurc_idx == 1])
cate_pub_draws  <- rowMeans(causal_diff[, insurc_idx == 2])
cate_unins_draws <- rowMeans(causal_diff[, insurc_idx == 3])

# CATEs by Poverty Status (1=Poor, 2=Near Poor, 3=Low, 4=Middle, 5=High)
povcat_idx <- as.integer(as.character(df$povcat21))
cate_poor_draws <- rowMeans(causal_diff[, povcat_idx == 1])
cate_nearpoor_draws <- rowMeans(causal_diff[, povcat_idx == 2])
cate_lowinc_draws <- rowMeans(causal_diff[, povcat_idx == 3])
cate_midinc_draws <- rowMeans(causal_diff[, povcat_idx == 4])
cate_highinc_draws <- rowMeans(causal_diff[, povcat_idx == 5])

res <- list(
    chain_id = args$chain_id,
    outcome = args$outcome,
    sate_draws = sate_draws,
    cate_priv_draws = cate_priv_draws,
    cate_pub_draws = cate_pub_draws,
    cate_unins_draws = cate_unins_draws,
    cate_poor_draws = cate_poor_draws,
    cate_nearpoor_draws = cate_nearpoor_draws,
    cate_lowinc_draws = cate_lowinc_draws,
    cate_midinc_draws = cate_midinc_draws,
    cate_highinc_draws = cate_highinc_draws
)

out_file <- sprintf("data/meps_chain_%s_%d.rds", args$outcome, args$chain_id)
saveRDS(res, out_file)
cat(sprintf("Saved %d draws to %s\n", n_post, out_file))
