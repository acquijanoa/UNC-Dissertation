# survey_bart_rcpp.R
# R wrapper for the C++ Rcpp implementation of Survey BART

library(Rcpp)
library(RcppArmadillo)

# Source the C++ code
sourceCpp("scripts/src/survey_bart.cpp")

#' Generate Kim-Rao Stratified Multinomial Bootstrap Weights
#' @keywords internal
generate_kim_rao_weights <- function(data, strata_col, psu_col, weight_col, B) {
    n <- nrow(data)
    W_mat <- matrix(0, nrow = B, ncol = n)
    
    # Split by strata
    strata_list <- split(seq_len(n), data[[strata_col]])
    
    for (b in seq_len(B)) {
        w_star <- numeric(n)
        for (idx in strata_list) {
            stratum_data <- data[idx, , drop = FALSE]
            psus <- unique(stratum_data[[psu_col]])
            n_h <- length(psus)
            
            if (n_h == 1) {
                w_star[idx] <- stratum_data[[weight_col]]
                next
            }
            
            # Multinomial draw
            r_h <- as.numeric(rmultinom(1, n_h - 1, rep(1/n_h, n_h)))
            k_h <- n_h / (n_h - 1)
            
            # Broadcast to units
            psu_map <- setNames(r_h, psus)
            r_expanded <- psu_map[as.character(stratum_data[[psu_col]])]
            
            w_star[idx] <- k_h * r_expanded * stratum_data[[weight_col]]
        }
        W_mat[b, ] <- w_star
    }
    return(W_mat)
}


#' Survey BART Model
#'
#' @param formula Model formula (e.g., y ~ x1 + x2)
#' @param data Data frame
#' @param family "continuous", "probit", or "hurdle"
#' @param weights Vector of sampling weights (default 1s)
#' @param strata Vector of strata IDs (optional)
#' @param psu Vector of PSU IDs (optional)
#' @param n_trees Number of trees (default 200)
#' @param n_mcmc Total MCMC iterations (default 2000)
#' @param burn_in Number of burn-in iterations (default 500)
#' @param X_test Optional matrix for out-of-sample predictions (e.g., G-computation)
#' @return A list with MCMC results
survey_bart <- function(formula, data, family = "continuous", 
                        weights = NULL, strata = NULL, psu = NULL,
                        n_trees = 200, n_mcmc = 2000, burn_in = 500,
                        X_test = NULL) {
    
    # 1. Parse formula to get X and y
    mf <- model.frame(formula, data, na.action = na.pass)
    y <- model.response(mf)
    X <- model.matrix(formula, mf)
    
    if("(Intercept)" %in% colnames(X)) {
        X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
    }
    
    # 2. Survey Weights Generation
    n <- length(y)
    if (is.null(weights)) {
        # i.i.d. case
        W_mat <- matrix(1, nrow = n_mcmc, ncol = n)
    } else if (!is.null(strata) && !is.null(psu)) {
        cat("Generating Kim-Rao bootstrap weight matrix...\n")
        # Temporarily append for easy generation
        tmp_df <- data.frame(weight = weights, strata = strata, psu = psu)
        W_mat <- generate_kim_rao_weights(tmp_df, "strata", "psu", "weight", n_mcmc)
    } else {
        # Random Exponential weights (Fallback if no design provided)
        cat("Warning: Design missing. Using Exp(1) random weights.\n")
        W_mat <- matrix(rexp(n_mcmc * n), nrow = n_mcmc, ncol = n)
        W_mat <- sweep(W_mat, 2, weights, `*`)
    }
    
    if (!family %in% c("continuous", "probit", "hurdle")) {
        stop("family must be one of: continuous, probit, hurdle")
    }
    
    # 3. Call C++ Backend
    cat(sprintf("Starting Survey BART C++ backend... [%s]\n", family))
    
    time_start <- Sys.time()
    result <- fit_bart_cpp(
        X = as.matrix(X), 
        y = as.numeric(y), 
        W_mat = as.matrix(W_mat), 
        n_trees = as.integer(n_trees), 
        n_mcmc = as.integer(n_mcmc), 
        burn_in = as.integer(burn_in), 
        family = as.character(family),
        X_test = if (is.null(X_test)) NULL else as.matrix(X_test)
    )
    time_end <- Sys.time()
    
    cat("Finished in", round(difftime(time_end, time_start, units = "secs"), 2), "seconds.\n")
    
    out <- list(
        y_hat_train = result$y_hat_train,
        y_hat_test = result$y_hat_test,
        sigma2_draws = result$sigma2_draws,
        W_mat = W_mat,
        formula = formula,
        family = family,
        n_mcmc = n_mcmc,
        burn_in = burn_in
    )
    class(out) <- "survey_bart"
    return(out)
}
