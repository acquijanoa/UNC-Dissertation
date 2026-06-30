# Path to Rscript (fallback to /usr/local/bin/Rscript if not in PATH)
RSCRIPT ?= $(shell which Rscript 2>/dev/null || echo /usr/local/bin/Rscript)

.PHONY: all generate_samples evaluate_coverage clean help

# Default target
all: evaluate_coverage

# Generate 1,000 samples in parallel using xargs
generate_samples:
	@echo "Generating 1,000 samples..."
	seq 1 1000 | xargs -I {} -P 4 $(RSCRIPT) scripts/sample.R {}

# Evaluate survey design coverage across the 1,000 samples
evaluate_coverage:
	@echo "Evaluating coverage across samples..."
	$(RSCRIPT) scripts/evaluate_coverage.R

# Clean up generated samples and results
clean:
	@echo "Cleaning up generated files..."
	rm -rf data/samples/*.csv
	rm -f data/coverage_results.csv

# Display help message
help:
	@echo "Available Makefile targets:"
	@echo "  make generate_samples  - Generate 1,000 samples in parallel"
	@echo "  make evaluate_coverage - Estimate mean of y and calculate coverage rates"
	@echo "  make clean             - Remove all generated sample CSV files and results"
	@echo "  make all               - Run evaluate_coverage (default)"
