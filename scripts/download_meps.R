# scripts/download_meps.R
# ==============================================================================
# Script to install MEPS and download 2021 Full Year Consolidated data
# ==============================================================================

# Install MEPS from GitHub (R defaults to the user's personal library if needed)
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("MEPS", quietly = TRUE)) {
  remotes::install_github("e-mitchell/meps_r_pkg/MEPS", upgrade = "never")
}

library(MEPS)
library(dplyr)

# Download 2021 Household Component FYC data (HC-233)
meps_data <- read_MEPS(year = 2021, type = "FYC")
names(meps_data) <- tolower(names(meps_data))

# Save output to data/
dir.create("data", showWarnings = FALSE, recursive = TRUE)
write.csv(meps_data, file = "data/meps_2021_hc233.csv", row.names = FALSE)
cat("MEPS data downloaded successfully to data/meps_2021_hc233.csv!\n")
