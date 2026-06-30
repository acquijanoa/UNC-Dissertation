# Population Generation and Sample Size Estimation

**Author:** Álvaro Quijano  
**Date:** June 2026  
**Scripts:** `scripts/gen_pop.R`, `scripts/sample_size.py`

This report documents the finite population data-generating process (DGP) and the empirical sample size calculation used to calibrate the simulation sampling design.

---

## 1. Overview

The simulation study follows Kim & Rao (2024) Supplemental S9: a **stratified two-stage cluster design** with **informative probability-proportional-to-size (PPS)** cluster selection. The cross-sectional population script (`gen_pop.R`) generates a single finite population saved to `data/population.csv`. The sample size script (`sample_size.py`) reads that population, estimates the intraclass correlation and design effect from the realized outcome `y`, and derives the target individual-level sample size used in `scripts/sample.R` ($n = 6{,}100$).

```mermaid
flowchart LR
  A[gen_pop.R] --> B[population.csv]
  B --> C[sample_size.py]
  C --> D["n ≈ 6,100"]
  B --> E[sample.R]
  D --> E
```

---

## 2. Population Generation (`gen_pop.R`)

### 2.1 Design rationale

The DGP embeds a **four-level hierarchy** — strata → clusters (PSUs) → individuals — with random effects that enter both cluster sizes and the outcome. Because cluster size $M_{hi}$ depends on $\lambda_h$ and $\eta_i$, which also appear in the outcome equation, PPS sampling proportional to $M_{hi}$ is **informative** with respect to $y$. This mirrors the inflated Type I error scenario documented in Kim & Rao (2024).

The population also includes **$K = 3$ latent classes** at the individual level. Ground-truth class labels (`class`) are retained for post-hoc evaluation but are not passed to samplers or estimators.

### 2.2 Structural parameters

| Symbol | Role | Value |
|:---|:---|:---|
| $H$ | Number of strata | 30 |
| $N_1$ | Minimum clusters per stratum (floor) | 50 |
| $N_2$ | Minimum individuals per cluster (floor) | 20 |
| $K$ | Number of latent classes | 3 |
| Seed | Reproducibility | `170626` |

### 2.3 Hierarchical data-generating process

**Level 1 — Strata.** For each stratum $h = 1, \ldots, H$, draw a stratum effect:

$$
\lambda_h \sim \mathrm{Exp}(1)
$$

The number of clusters in stratum $h$ is:

$$
M_h = 8 \cdot \mathrm{Poisson}(\lambda_h) + N_1
$$

**Level 2 — Clusters.** Within stratum $h$, for each cluster $i$, draw a cluster effect:

$$
\eta_i \sim \mathrm{Exp}(1)
$$

The number of individuals in cluster $(h,i)$ is:

$$
M_{hi} = 8 \cdot \mathrm{Poisson}(\lambda_h + \eta_i) + N_2
$$

The stratum total $N_h = \sum_i M_{hi}$ serves as the denominator for PPS inclusion probabilities ($p_{hi} = M_{hi}/N_h$).

**Level 3 — Individuals.** Within each cluster, for each individual $j$:

1. Draw latent class $z_{ij} \sim \mathrm{Categorical}(\pi_1, \pi_2, \pi_3)$.
2. Draw class-specific random intercept $u_{ij} \sim N(0, \sigma_{u,z_{ij}}^2)$.
3. Draw covariates and residual $\varepsilon_{ij} \sim N(0,1)$.

### 2.4 Latent class parameters

| Class $k$ | $\pi_k$ | $\alpha_k$ (intercept) | $\beta_{1k}$ (slope on $x_{\mathrm{inv}}$) | $\sigma_{u,k}$ |
|:---:|:---:|:---:|:---:|:---:|
| 1 — low | 0.50 | 5.0 | 0.8 | 0.3 |
| 2 — moderate | 0.30 | 10.0 | 1.0 | 0.5 |
| 3 — high | 0.20 | 15.0 | 1.2 | 0.7 |

### 2.5 Covariates

| Variable | Distribution | Notes |
|:---|:---|:---|
| $x_{\mathrm{inv}}$ | $\mathrm{Uniform}(0, 20)$ | Continuous, time-invariant analogue |
| $x_{\mathrm{cont\_2}}$ | $N(0, 1)$ | Continuous |
| $x_{\mathrm{bin}}$ | $\mathrm{Bernoulli}(p)$, $p = \mathrm{logit}^{-1}(-0.3 + 0.4\, x_{\mathrm{cont\_2}})$ | Binary, depends on $x_{\mathrm{cont\_2}}$ |
| $x_{\mathrm{cat}}$ | $\{1,2,3,4\}$ with probs $(0.30, 0.30, 0.20, 0.20)$ | Categorical; level 1 is reference |

### 2.6 Outcome model

The continuous outcome for individual $j$ in cluster $(h,i)$ belonging to class $k$ is:

$$
\begin{aligned}
y_{hij} =\;& \alpha_k + \beta_{1k}\, x_{\mathrm{inv}} \\
&+ \delta_{\mathrm{cont2}}\, x_{\mathrm{cont\_2}} + \delta_{\mathrm{bin}}\, x_{\mathrm{bin}} + \delta_{\mathrm{cat}}[x_{\mathrm{cat}}] \\
&+ \delta_{\mathrm{int,cont2\_bin}}\, (x_{\mathrm{cont\_2}} \cdot x_{\mathrm{bin}}) + \delta_{\mathrm{int,inv\_bin}}\, (x_{\mathrm{inv}} \cdot x_{\mathrm{bin}}) \\
&+ 2.0\,\eta_i + 2.0\,\lambda_h + u_{ij} + \varepsilon_{ij}
\end{aligned}
$$

**Fixed population coefficients:**

| Coefficient | Value |
|:---|:---:|
| $\delta_{\mathrm{cont2}}$ | 0.50 |
| $\delta_{\mathrm{bin}}$ | 0.80 |
| $\delta_{\mathrm{cat}}$ | $(0,\; 0.30,\; 0.60,\; 0.90)$ |
| $\delta_{\mathrm{int,cont2\_bin}}$ | $-0.40$ |
| $\delta_{\mathrm{int,inv\_bin}}$ | $-0.20$ |

The coefficients on $\eta_i$ and $\lambda_h$ are set to 2.0 (rather than 0.3 in the longitudinal extension documented in `simulation_setup.md`), strengthening the informativeness of the size–outcome linkage in this cross-sectional arm.

### 2.7 Output structure

The generated dataset is saved to `data/population.csv` with the following column groups:

| Group | Columns | Purpose |
|:---|:---|:---|
| Identifiers | `strata`, `psuid`, `subid` | Nested indexing |
| Design (informative PPS) | `lambda_h`, `eta_i`, `M_hi`, `N_h` | Size mechanism; **not** for analysis models |
| Covariates | `x_inv`, `x_cont_2`, `x_bin`, `x_cat` | Observed predictors |
| Latent truth | `class`, `u_ij`, `eps` | Ground truth; **not** for analysis models |
| Outcome | `y` | Target variable |

### 2.8 Realized population (current run)

| Quantity | Value |
|:---|---:|
| Total individuals ($N$) | 65,360 |
| Strata ($H$) | 30 |
| Clusters (PSUs) | 1,780 |
| Average cluster size ($\bar{m}$) | 36.72 |

---

## 3. Sample Size Estimation (`sample_size.py`)

### 3.1 Objective

Given the realized finite population, estimate the **individual-level sample size** required to achieve a margin of error $e = 0.05$ on the population mean of $y$, accounting for **cluster-induced correlation** via Lohr's design effect framework.

### 3.2 Population summary

The script constructs a unique cluster identifier `cluster_id = strata + '_' + psuid` and computes:

$$
N = 65{,}360, \quad H = 30, \quad n_{\mathrm{clust}} = 1{,}780, \quad \bar{m} = N / n_{\mathrm{clust}} \approx 36.72
$$

### 3.3 Lohr's adjusted $R^2_\alpha$ (ICC proxy)

Following Lohr (*Sampling: Design and Analysis*), the intraclass correlation is estimated from the analysis of variance decomposition of $y$:

1. **Total variance:** $S^2 = \mathrm{Var}(y)$
2. **Within-cluster mean square:**

$$
\mathrm{MSW} = \frac{\sum_{j}(y_j - \bar{y}_{c(j)})^2}{N - n_{\mathrm{clust}}}
$$

where $\bar{y}_{c(j)}$ is the cluster mean for the cluster containing individual $j$.

3. **Lohr's adjusted $R^2_\alpha$:**

$$
R^2_\alpha = 1 - \frac{\mathrm{MSW}}{S^2}
$$

This quantity acts as an ICC-like measure of the proportion of total variance attributable to between-cluster differences.

**Empirical estimates from `population.csv`:**

| Statistic | Value |
|:---|---:|
| $S^2$ | 65.15 |
| $\mathrm{MSW}$ | 57.37 |
| $R^2_\alpha$ | 0.119 |

### 3.4 Design effect

The design effect for cluster sampling is:

$$
\mathrm{DEFF} = 1 + R^2_\alpha\,(\bar{m} - 1)
$$

**Estimated DEFF:** $5.27$

This implies that a cluster sample needs roughly 5.3 times as many individuals as a simple random sample (SRS) to achieve the same precision for the population mean of $y$.

### 3.5 SRS baseline sample size

Assuming a 95% confidence interval for the population mean with margin of error $e = 0.05$ on the scale of $y$:

$$
n_{\mathrm{SRS}} = \left(\frac{z_{0.975}\, s}{e}\right)^2 = \left(\frac{1.96 \times 8.07}{0.5}\right)^2 \approx 1{,}001
$$

where $s = \sqrt{S^2} \approx 8.07$. The finite population correction (FPC) is intentionally omitted because the target design is a complex cluster sample, not an SRS from the finite population.

### 3.6 Final sample size

The adjusted sample size incorporates the design effect and a **15% conservative buffer** on the within-cluster design effect ($\mathrm{DEFF}_w$):

$$
n = n_{\mathrm{SRS}} \times \mathrm{DEFF} \times 1.15 \approx 6{,}060 \;\Rightarrow\; \boxed{n = 6{,}100}
$$

(rounding to the nearest hundred).

This value is used directly in `scripts/sample.R` as the target number of sampled individuals.

### 3.7 Implied sampling fractions

| Level | Implied quantity | Approximate value |
|:---|:---|---:|
| Individual | $n / N$ | 9.3% |
| Cluster | $n / \bar{m}$ | $\approx 166$ clusters |
| Cluster | $166 / 1{,}780$ | 9.3% |

---

## 4. Connection between the two scripts

| Step | Script | Output |
|:---|:---|:---|
| 1. Generate superpopulation draw | `gen_pop.R` | `data/population.csv` (65,360 rows) |
| 2. Estimate DEFF and $n$ | `sample_size.py` | $n = 6{,}100$ |
| 3. Draw replicate samples | `sample.R` | `data/samples/sample_XXXX.csv` |

The sample size calculation is **empirical**: it conditions on the specific finite population realization produced by `gen_pop.R`. Re-running `gen_pop.R` with a different seed would change $R^2_\alpha$, DEFF, and the recommended $n$.

---

## 5. Key assumptions and limitations

1. **Cross-sectional DGP.** `gen_pop.R` implements a single-wave outcome. The longitudinal extension (panel waves, AR(1) errors, time trends $\gamma_k$) is specified separately in `simulation_setup.md`.

2. **Informative design by construction.** Cluster and stratum effects enter both $M_{hi}$ and $y$ with coefficient 2.0, ensuring PPS selection is informative. Naïve estimators that ignore the design are expected to be biased.

3. **Sample size targets the population mean of $y$.** The DEFF calculation does not account for additional precision requirements for regression coefficients, latent class recovery, or multivariate DPMM estimation. The 15% buffer provides partial insurance against underestimation.

4. **Sampling implementation note.** `sample.R` currently implements **SRS cluster selection** (equal probability) rather than PPS with replacement as in Kim & Rao S9. The sample size was calibrated under the realized population correlation structure; the sampling mechanism in `sample.R` may warrant alignment with the intended S9 design in future iterations.

---

## References

- Kim, J. K., Rao, J. N. K., & Wang, Z. (2024). Hypotheses testing from complex survey data using bootstrap weights: A unified approach. *Journal of the American Statistical Association*, 119(546), 1229–1239.
- Lohr, S. L. *Sampling: Design and Analysis* (2nd ed.). Chapman & Hall/CRC.
