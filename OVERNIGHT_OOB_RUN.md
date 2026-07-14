# OOB K-Calibration: Overnight Run Summary

**Submitted:** $(date)  
**User:** aquijano

## Jobs Submitted

| Job Type | Job ID | Array Size | Method | Status |
|:---------|:-------|:-----------|:-------|:-------|
| Method A (OOB) | 58911044 | 1000 samples | A (weight_norm) | Running |
| Method B (OOB) | 58911045 | 1000 samples | B (raw weights) | Running |
| Report | 58911046 | Single | Comparison | Pending (depends on A+B) |

## What's Running

### Approach
For each of 100 samples, the script:
1. **Runs OOB calibration** on a K grid using short MCMC (50 trees, 500 iterations)
   - Method A: K ∈ {0.01, 0.05, 0.1, 0.2, 0.5, 1.0}
   - Method B: K ∈ {0.5, 1.0, 2.0, 5.0, 10.0}
2. **Selects K*** that minimizes average OOB MSE across post-burn-in iterations
3. **Runs full MCMC** with selected K* (200 trees, 2000 iterations)
4. **Saves results** including: K*, OOB MSE, bias, coverage, RMSE, PP coverage

### Expected Completion
- Each sample: ~30-60 minutes
- Total runtime: ~1-2 hours (array jobs run in parallel)
- Report generation: ~5 minutes after jobs complete

## Output Locations

```
data/friedman_results_oob/
├── result_0001_A_OOB.csv
├── result_0001_B_OOB.csv
├── ...
├── result_0100_A_OOB.csv
├── result_0100_B_OOB.csv
├── report_table.csv          # Auto-generated after completion
└── report.md                 # Auto-generated after completion
```

## Monitoring

```bash
# Check job status
squeue -u aquijano

# Check specific job details
scontrol show job 58910821

# View logs (while running or after completion)
tail -f logs/friedman_oob_K_1.out
tail -f logs/report_oob_K.out

# Check how many results are complete
ls data/friedman_results_oob/*.csv | wc -l
```

## What to Expect Tomorrow Morning

The report will show:

1. **Mean K selected** across 100 replications for each method
2. **Standard deviation of K** (stability of OOB selection)
3. **Coverage of θ** (mean surface) — compare to grid results
4. **RMSE** — should be near-optimal since OOB minimizes MSE
5. **PP coverage** — should be near-nominal

### Key Questions to Answer

1. **Does OOB select the RMSE-optimal K from the grid?**
   - For Method A: Grid showed K=0.5 minimizes RMSE (2.894)
   - For Method A: Grid showed K=0.1 gives nominal coverage (94.9%)
   - **Hypothesis:** OOB will select K ≈ 0.5 (RMSE-optimal, not coverage-optimal)

2. **Is OOB K-selection stable across samples?**
   - Low sd(K) → stable, reliable selection
   - High sd(K) → noisy, sensitive to sample variation

3. **Trade-off confirmation:**
   - If OOB K gives better RMSE but worse coverage than grid K=0.1, confirms the RMSE vs coverage trade-off
   - Validates the interpretation: **OOB optimizes prediction (RMSE), not inference (coverage)**

## Comparison Script

After jobs complete, run:

```bash
Rscript scripts/compare_oob_vs_grid.R
```

This will produce a side-by-side comparison:
- OOB-selected K* vs grid search results
- Highlight whether OOB selected RMSE-optimal or coverage-optimal K
- Confirm theoretical predictions about the trade-off

## Next Steps

Tomorrow morning:
1. Check that all 200 result files exist (100 per method)
2. Review `data/friedman_results_oob/report.md`
3. Run comparison script to validate against grid
4. Update Section 7 of `weight_calibration_proposal.md` with OOB findings
5. Interpret whether OOB is appropriate for imputation vs mean estimation

## Theoretical Implications

**Expected finding:** OOB will select K that optimizes predictive accuracy (RMSE), not coverage of population means (θ).

**Why this matters:**
- **For imputation:** OOB-selected K is the right choice (predicting individual y_i)
- **For mean estimation:** Coverage-optimal K from grid is the right choice (estimating E[μ])
- **For your dissertation:** Demonstrates that the research question determines the calibration criterion

This overnight run will provide empirical validation of the RMSE vs coverage distinction emphasized in the recent theoretical interpretation.
