# merge_meps_chains.R

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
    outcome <- args[1]
} else {
    outcome <- "totslf21"
}

cat(sprintf("\nMerging 4 MCMC Chains for %s...\n", outcome))

sate_all <- numeric(0)
cate_priv_all <- numeric(0)
cate_pub_all <- numeric(0)
cate_unins_all <- numeric(0)
cate_poor_all <- numeric(0)
cate_nearpoor_all <- numeric(0)
cate_lowinc_all <- numeric(0)
cate_midinc_all <- numeric(0)
cate_highinc_all <- numeric(0)

for (i in 1:4) {
    f <- sprintf("data/meps_chain_%s_%d.rds", outcome, i)
    if (!file.exists(f)) {
        stop(sprintf("File %s not found!", f))
    }
    res <- readRDS(f)
    sate_all <- c(sate_all, res$sate_draws)
    cate_priv_all <- c(cate_priv_all, res$cate_priv_draws)
    cate_pub_all <- c(cate_pub_all, res$cate_pub_draws)
    cate_unins_all <- c(cate_unins_all, res$cate_unins_draws)
    
    if (!is.null(res$cate_poor_draws)) {
        cate_poor_all <- c(cate_poor_all, res$cate_poor_draws)
        cate_nearpoor_all <- c(cate_nearpoor_all, res$cate_nearpoor_draws)
        cate_lowinc_all <- c(cate_lowinc_all, res$cate_lowinc_draws)
        cate_midinc_all <- c(cate_midinc_all, res$cate_midinc_draws)
        cate_highinc_all <- c(cate_highinc_all, res$cate_highinc_draws)
    }
}

cat("\n========================================\n")
cat(sprintf("CAUSAL INFERENCE RESULTS (MEPS 2021)\n"))
cat("========================================\n")
cat(sprintf("Outcome: %s\n", outcome))
cat(sprintf("Total Posterior Draws: %d\n", length(sate_all)))
cat("\n--- OVERALL SATE ---\n")
cat(sprintf("Mean:      $%.2f\n", mean(sate_all)))
cat(sprintf("Std Err:   $%.2f\n", sd(sate_all)))
ci <- quantile(sate_all, probs=c(0.025, 0.975))
cat(sprintf("95%% Cred:  [$%.2f, $%.2f]\n", ci[1], ci[2]))

cat("\n--- CATE: Private Insurance ---\n")
cat(sprintf("Mean:      $%.2f\n", mean(cate_priv_all)))
ci <- quantile(cate_priv_all, probs=c(0.025, 0.975))
cat(sprintf("95%% Cred:  [$%.2f, $%.2f]\n", ci[1], ci[2]))

cat("\n--- CATE: Public Only ---\n")
cat(sprintf("Mean:      $%.2f\n", mean(cate_pub_all)))
ci <- quantile(cate_pub_all, probs=c(0.025, 0.975))
cat(sprintf("95%% Cred:  [$%.2f, $%.2f]\n", ci[1], ci[2]))

cat("\n--- CATE: Uninsured ---\n")
cat(sprintf("Mean:      $%.2f\n", mean(cate_unins_all)))
ci <- quantile(cate_unins_all, probs=c(0.025, 0.975))
cat(sprintf("95%% Cred:  [$%.2f, $%.2f]\n", ci[1], ci[2]))

if (length(cate_poor_all) > 0) {
    cat("\n--- CATE: Poor ---\n")
    cat(sprintf("Mean:      $%.2f\n", mean(cate_poor_all)))
    ci <- quantile(cate_poor_all, probs=c(0.025, 0.975))
    cat(sprintf("95%% Cred:  [$%.2f, $%.2f]\n", ci[1], ci[2]))

    cat("\n--- CATE: Near Poor ---\n")
    cat(sprintf("Mean:      $%.2f\n", mean(cate_nearpoor_all)))
    ci <- quantile(cate_nearpoor_all, probs=c(0.025, 0.975))
    cat(sprintf("95%% Cred:  [$%.2f, $%.2f]\n", ci[1], ci[2]))

    cat("\n--- CATE: Low Income ---\n")
    cat(sprintf("Mean:      $%.2f\n", mean(cate_lowinc_all)))
    ci <- quantile(cate_lowinc_all, probs=c(0.025, 0.975))
    cat(sprintf("95%% Cred:  [$%.2f, $%.2f]\n", ci[1], ci[2]))

    cat("\n--- CATE: Middle Income ---\n")
    cat(sprintf("Mean:      $%.2f\n", mean(cate_midinc_all)))
    ci <- quantile(cate_midinc_all, probs=c(0.025, 0.975))
    cat(sprintf("95%% Cred:  [$%.2f, $%.2f]\n", ci[1], ci[2]))

    cat("\n--- CATE: High Income ---\n")
    cat(sprintf("Mean:      $%.2f\n", mean(cate_highinc_all)))
    ci <- quantile(cate_highinc_all, probs=c(0.025, 0.975))
    cat(sprintf("95%% Cred:  [$%.2f, $%.2f]\n", ci[1], ci[2]))
}

cat("========================================\n")
