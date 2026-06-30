# evaluate_coverage.R
# Estimates the mean of y and evaluates coverage across 1,000 samples using the survey package

suppressPackageStartupMessages({
  library(survey)
  library(data.table)
  library(parallel)
})

# Adjust for lonely PSUs if any stratum has only one PSU sampled
options(survey.lonely.psu = "adjust")

# 1. Load population and calculate true mean
pop_file <- here::here("data", "population.csv")
if (!file.exists(pop_file)) {
  stop("Population file not found at: ", pop_file, ". Please run scripts/gen_pop.R first.")
}

message("Loading population data...")
pop <- fread(pop_file)
true_mean_y <- mean(pop$y)
message(sprintf("True Population Mean of y: %.6f", true_mean_y))

# 2. Get list of all generated sample files
samples_dir <- here::here("data", "samples")
sample_files <- list.files(samples_dir, pattern = "^sample_\\d{4}\\.csv$", full.names = TRUE)
n_samples <- length(sample_files)

if (n_samples == 0) {
  stop("No sample files found in: ", samples_dir, ". Please run the sample generation first.")
}

message(sprintf("Found %d sample files. Processing in parallel...", n_samples))

# Detect available cores (leave one core free)
num_cores <- max(1, detectCores() - 1)

# 3. Process each sample
process_sample <- function(file_path) {
  # Extract sample index from filename
  file_name <- basename(file_path)
  idx <- as.integer(gsub("[^0-9]", "", file_name))
  
  # Read sample data
  sample_dt <- fread(file_path)
  
  # --- 1. Weighted (Survey Design) Estimation ---
  # Define the survey design taking into account strata, psuid (clusters), and weights
  design <- svydesign(
    ids     = ~psuid,
    strata  = ~strata,
    weights = ~weight_norm,
    data    = sample_dt,
    nest    = TRUE
  )
  
  est_w <- svymean(~y, design)
  mean_w <- coef(est_w)[1]
  se_w <- SE(est_w)[1]
  ci_w <- confint(est_w, level = 0.95)
  covered_w <- (true_mean_y >= ci_w[1, 1]) && (true_mean_y <= ci_w[1, 2])
  
  # --- 2. Unweighted (Naive) Estimation ---
  mean_u <- mean(sample_dt$y)
  se_u <- sd(sample_dt$y) / sqrt(nrow(sample_dt))
  ci_u_low <- mean_u - 1.96 * se_u
  ci_u_high <- mean_u + 1.96 * se_u
  covered_u <- (true_mean_y >= ci_u_low) && (true_mean_y <= ci_u_high)
  
  return(data.table(
    sample_index = idx,
    mean_weighted = mean_w,
    se_weighted = se_w,
    covered_weighted = as.integer(covered_w),
    mean_unweighted = mean_u,
    se_unweighted = se_u,
    covered_unweighted = as.integer(covered_u)
  ))
}

# Run loop in parallel
results_list <- mclapply(sample_files, process_sample, mc.cores = num_cores)

# Bind results into a single data.table
results <- rbindlist(results_list)
setorder(results, sample_index)

# 4. Save detailed results
results_file <- here::here("data", "coverage_results.csv")
fwrite(results, results_file)
message("Detailed simulation results saved to: ", results_file)

# 5. Summarize findings
summary_stats <- data.table(
  Metric = c(
    "Number of Samples",
    "True Population Mean",
    "Average Estimate",
    "Empirical Bias",
    "Relative Bias (%)",
    "Average Estimated SE",
    "Empirical SD of Estimates",
    "95% CI Coverage Rate (%)"
  ),
  Weighted = c(
    n_samples,
    true_mean_y,
    mean(results$mean_weighted),
    mean(results$mean_weighted) - true_mean_y,
    ((mean(results$mean_weighted) - true_mean_y) / true_mean_y) * 100,
    mean(results$se_weighted),
    sd(results$mean_weighted),
    mean(results$covered_weighted) * 100
  ),
  Unweighted = c(
    n_samples,
    true_mean_y,
    mean(results$mean_unweighted),
    mean(results$mean_unweighted) - true_mean_y,
    ((mean(results$mean_unweighted) - true_mean_y) / true_mean_y) * 100,
    mean(results$se_unweighted),
    sd(results$mean_unweighted),
    mean(results$covered_unweighted) * 100
  )
)

# Print Summary Table
cat("\n==========================================================\n")
cat("          SIMULATION STUDY SUMMARY STATISTICS             \n")
cat("==========================================================\n")
print(summary_stats, digits = 4)
cat("==========================================================\n")
