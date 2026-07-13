# Set library
suppressPackageStartupMessages({
  library(survey)
  library(data.table)
})

# Parse command-line argument for sample index (defaults to 1 if not provided)
args <- commandArgs(trailingOnly = TRUE)
sample_index <- if (length(args) > 0) as.integer(args[1]) else 1

# Set seed using the sample index to guarantee uniqueness and reproducibility
set.seed(06222026 + sample_index)

# Load population
data_dir <- here::here("data")
dt <- fread(file.path(data_dir, "population.csv"))

# Sample size
n <- 6100
avg_clust_size <- round(36.72)
n_clusters <- n/avg_clust_size

# Number of clusters per strata and their sizes/populations
data <- unique(dt[, .(strata, psuid, M_hi, N_h)])
data[, C_h := .N, by = strata]

# Calculate total number of unique clusters in the population
total_clusters <- nrow(data)

# Proportion of clusters to sample per strata
data[, prop := C_h / total_clusters]
data[, sample_size_h := ceiling(prop * (n / avg_clust_size))]

# Sample clusters per stratum using SRS (equal probability without replacement)
data[, psu_sample := psuid %in% sample(psuid, size = sample_size_h[1], replace = FALSE), by = strata]

# First-stage weight under SRS: w1 = C_h / sample_size_h
data[, w1 := C_h / sample_size_h]

# Keep only sampled clusters using a composite key join (since psuid is nested within strata)
sampled_clusters <- data[psu_sample == TRUE, .(strata, psuid, w1)]
sample_dt <- dt[sampled_clusters, on = .(strata, psuid)]

# compute cluster_size to sample subjects (use avg_clust_size from line 11)
n2 <- avg_clust_size

# Sample subjects within each cluster (SRS without replacement)
sample_dt <- sample_dt[, .SD[sample(.N, size = min(.N, n2), replace = FALSE)], by = .(strata, psuid)]

# Compute second-stage weight and final weight
sample_dt[, n_sampled := .N, by = .(strata, psuid)]
sample_dt[, w2 := M_hi / n_sampled]
sample_dt[, weight := w1 * w2]

# Compute normalized weights that sum to the final sample size (number of individuals)
sample_dt[, weight_norm := weight * (nrow(sample_dt) / sum(weight))]

# Keep only columns needed for fitting / Kim-Rao (drop DGP truth & intermediates)
sample_dt <- sample_dt[, .(
  strata, psuid, subid,
  weight, weight_norm,
  x_inv, x_cont_2, x_bin, x_cat,
  x_f1, x_f2, x_f3, x_f4, x_f5,
  y, y_f
)]

# Create output folder and save sample
output_dir <- "data/samples"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
file_name <- file.path(output_dir, sprintf("sample_%04d.csv", sample_index))
fwrite(sample_dt, file = file_name)
writeLines(sprintf("Sample %d saved to: %s", sample_index, file_name))


