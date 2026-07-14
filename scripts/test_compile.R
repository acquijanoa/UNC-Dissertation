# test_compile.R
# Run from project root: Rscript scripts/test_compile.R
source("scripts/survey_bart_rcpp.R")

cat("Testing dummy fit...\n")
set.seed(42)
n <- 100
df <- data.frame(
    y = rnorm(n),
    y_bin = rbinom(n, 1, 0.5),
    y_hurdle = ifelse(runif(n) > 0.5, 0, rlnorm(n)),
    x1 = rnorm(n),
    x2 = runif(n)
)

cat("--- CONTINUOUS ---\n")
fit1 <- survey_bart(y ~ x1 + x2, data = df, family = "continuous", n_mcmc = 10, burn_in = 5)
print(dim(fit1$y_hat_train))

cat("--- PROBIT ---\n")
fit2 <- survey_bart(y_bin ~ x1 + x2, data = df, family = "probit", n_mcmc = 10, burn_in = 5)
print(dim(fit2$y_hat_train))

cat("--- HURDLE ---\n")
fit3 <- survey_bart(y_hurdle ~ x1 + x2, data = df, family = "hurdle", n_mcmc = 10, burn_in = 5)
print(dim(fit3$y_hat_train))

cat("Success!\n")
