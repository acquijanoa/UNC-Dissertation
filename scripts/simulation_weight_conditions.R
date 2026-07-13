# simulation_weight_conditions.R
library(data.table)

source("scripts/survey_bart_rcpp.R")

args_raw <- commandArgs(trailingOnly=TRUE)
args <- list(
    n_samples = 50, single_sample = NULL, out_file = NULL, condition = NULL,
    K_grid = "1,5,15,30,50,80", n_mcmc = 800, burn_in = 250, n_trees = 50, seed = 42
)
i <- 1
while(i <= length(args_raw)) {
    if (args_raw[i] == "--n_samples") { args$n_samples <- as.integer(args_raw[i+1]); i <- i+2 }
    else if (args_raw[i] == "--single_sample") { args$single_sample <- args_raw[i+1]; i <- i+2 }
    else if (args_raw[i] == "--out_file") { args$out_file <- args_raw[i+1]; i <- i+2 }
    else if (args_raw[i] == "--condition") { args$condition <- args_raw[i+1]; i <- i+2 }
    else if (args_raw[i] == "--K_grid") { args$K_grid <- args_raw[i+1]; i <- i+2 }
    else if (args_raw[i] == "--n_mcmc") { args$n_mcmc <- as.integer(args_raw[i+1]); i <- i+2 }
    else if (args_raw[i] == "--burn_in") { args$burn_in <- as.integer(args_raw[i+1]); i <- i+2 }
    else if (args_raw[i] == "--n_trees") { args$n_trees <- as.integer(args_raw[i+1]); i <- i+2 }
    else if (args_raw[i] == "--seed") { args$seed <- as.integer(args_raw[i+1]); i <- i+2 }
    else { i <- i+1 }
}
set.seed(args$seed)

FEATURE_COLS <- c("x_inv", "x_cont_2", "x_bin", "x_cat")
TARGET_COL <- "y"
STRATA_COL <- "strata"
PSU_COL <- "psuid"
WEIGHT_COL <- "weight"
N_TEST <- 2000

design_effect <- function(weights) {
    n <- length(weights)
    return((n * sum(weights^2)) / (sum(weights)^2))
}

# Load population
pop <- fread("data/population.csv")
pop_keys <- paste(pop[[STRATA_COL]], pop[[PSU_COL]], pop[["subid"]], sep="_")

build_test_set <- function(pop, df_train, n_test, pop_keys) {
    train_keys <- paste(df_train[[STRATA_COL]], df_train[[PSU_COL]], df_train[["subid"]], sep="_")
    mask_out <- !(pop_keys %in% train_keys)
    pool <- pop[mask_out, ]
    n_draw <- min(n_test, nrow(pool))
    return(pool[sample(nrow(pool), n_draw), ])
}

# OOB Calibration
calibrate_K_oob <- function(df_train, X_test_mat, y_test, K_candidates) {
    cat("\n--- Starting Fast OOB Calibration for K ---\n")
    best_loss <- Inf
    best_K <- K_candidates[1]
    
    for (K_cand in K_candidates) {
        # Quick MCMC
        fit <- survey_bart(
            formula = y ~ x_inv + x_cont_2 + x_bin + x_cat,
            data = as.data.frame(df_train),
            family = "continuous",
            weights = df_train[[WEIGHT_COL]] * K_cand,
            strata = df_train[[STRATA_COL]],
            psu = df_train[[PSU_COL]],
            n_trees = 50, n_mcmc = 500, burn_in = 150
        )
        
        # Calculate OOB MSE
        burn_idx <- 151:500
        n_post <- 350
        w_post <- fit$W_mat[burn_idx, ]
        oob_mask <- (w_post == 0)
        
        y_train <- df_train[[TARGET_COL]]
        y_mat <- matrix(y_train, nrow = n_post, ncol = length(y_train), byrow = TRUE)
        err_sq <- (y_mat[oob_mask] - fit$y_hat_train[oob_mask])^2
        oob_loss <- mean(err_sq)
        
        cat(sprintf("  Candidate K = %4.1f | OOB-MSE = %.4f\n", K_cand, oob_loss))
        
        if (!is.na(oob_loss) && oob_loss < best_loss) {
            best_loss <- oob_loss
            best_K <- K_cand
        }
    }
    cat(sprintf("--- Calibration Complete: Optimal K* = %s (OOB-MSE = %.4f) ---\n\n", best_K, best_loss))
    return(best_K)
}


# Main Loop setup
sample_files <- list.files("data/samples", pattern="sample_.*\\.csv", full.names=TRUE)
if (!is.null(args$single_sample)) {
    sample_files <- args$single_sample
} else {
    sample_files <- head(sort(sample_files), args$n_samples)
}

conditions <- if (!is.null(args$condition)) args$condition else c("A", "B", "C")
K_candidates <- as.numeric(strsplit(args$K_grid, ",")[[1]])

results <- list()
idx_count <- 1

for (sample_path in sample_files) {
    df_train <- fread(sample_path)
    df_test <- build_test_set(pop, df_train, N_TEST, pop_keys)
    X_test_mat <- as.matrix(df_test[, ..FEATURE_COLS])
    y_test <- df_test[[TARGET_COL]]
    
    rep_num <- as.numeric(gsub(".*sample_([0-9]+)\\.csv", "\\1", sample_path))
    if (is.na(rep_num)) rep_num <- 1
    
    for (cond in conditions) {
        n <- nrow(df_train)
        w_orig <- df_train[[WEIGHT_COL]]
        
        if (cond == "A") {
            w_use <- w_orig * n / sum(w_orig)
            K_used <- NA
        } else if (cond == "B") {
            w_use <- w_orig
            K_used <- 1
        } else if (cond == "C") {
            if (length(K_candidates) > 0) {
                K_used <- calibrate_K_oob(df_train, X_test_mat, y_test, K_candidates)
            } else {
                K_used <- max(1.0, n / (mean(w_orig) * 5))
            }
            w_use <- w_orig * K_used
        }
        
        fit <- survey_bart(
            formula = y ~ x_inv + x_cont_2 + x_bin + x_cat,
            data = as.data.frame(df_train),
            family = "continuous",
            weights = w_use,
            strata = df_train[[STRATA_COL]],
            psu = df_train[[PSU_COL]],
            n_trees = args$n_trees,
            n_mcmc = args$n_mcmc,
            burn_in = args$burn_in,
            X_test = X_test_mat
        )
        
        # Calculate OOB loss for reporting
        n_post <- args$n_mcmc - args$burn_in
        burn_idx <- (args$burn_in + 1):args$n_mcmc
        w_post <- fit$W_mat[burn_idx, ]
        oob_mask <- (w_post == 0)
        y_train <- df_train[[TARGET_COL]]
        y_mat <- matrix(y_train, nrow = n_post, ncol = length(y_train), byrow = TRUE)
        oob_loss <- mean((y_mat[oob_mask] - fit$y_hat_train[oob_mask])^2)
        
        # Metrics on X_test
        f_draws <- fit$y_hat_test
        pp_draws <- matrix(0, nrow = n_post, ncol = length(y_test))
        for (i in 1:n_post) {
            pp_draws[i, ] <- f_draws[i, ] + rnorm(length(y_test), 0, sqrt(fit$sigma2_draws[i]))
        }
        
        post_mean <- colMeans(f_draws)
        rmse <- sqrt(mean((post_mean - y_test)^2))
        
        mcmc_sd_f <- mean(apply(f_draws, 2, sd))
        
        f_lower <- apply(f_draws, 2, quantile, probs=0.025)
        f_upper <- apply(f_draws, 2, quantile, probs=0.975)
        f_ci_width <- mean(f_upper - f_lower)
        
        pp_lower <- apply(pp_draws, 2, quantile, probs=0.025)
        pp_upper <- apply(pp_draws, 2, quantile, probs=0.975)
        pp_coverage <- mean(y_test >= pp_lower & y_test <= pp_upper)
        pp_width <- mean(pp_upper - pp_lower)
        
        deff <- design_effect(w_orig)
        
        results[[idx_count]] <- data.frame(
            condition = cond,
            rmse = rmse,
            mcmc_sd_f = mcmc_sd_f,
            f_ci_width = f_ci_width,
            pp_coverage = pp_coverage,
            pp_width = pp_width,
            deff = deff,
            K = K_used,
            N_hat = sum(w_orig),
            n = n,
            oob_loss = oob_loss,
            rep = rep_num
        )
        
        cat(sprintf("  Rep %3d | Cond %s | PP-Cov=%.1f%%  f-SD=%.3f  f-CI=%.3f  RMSE=%.4f  Deff=%.3f  K=%.1f\n", 
                    rep_num, cond, 100*pp_coverage, mcmc_sd_f, f_ci_width, rmse, deff, K_used))
        
        idx_count <- idx_count + 1
    }
}

df_res <- rbindlist(results)
if (!is.null(args$out_file)) {
    out_path <- args$out_file
} else {
    out_path <- "data/simulation_weight_conditions.csv"
}

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
fwrite(df_res, out_path)
cat(sprintf("\nResults saved to: %s\n", out_path))
