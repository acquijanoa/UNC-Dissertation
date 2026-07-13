# Design-Based Nonparametric Inference for Complex Survey Data: Survey BART, Dirichlet Processes, and Kernel Distributional Tests

### A Research Proposal

**Author:** A. Quijano  
**Date:** July 2026  
**Status:** Research Proposal / Concept Note

---

> *We introduce* ***Design-Based Distributional Learning*** *(DBDL) as a new research program at the intersection of survey sampling theory and nonparametric inference. The program rests on three methodological pillars — all calibrated by the same* ***Kim-Rao stratified bootstrap*** *and all targeting the finite population rather than an i.i.d. superpopulation:*
>
> 1. ***Survey BART*** — *a Bayesian Additive Regression Trees model with design-weighted likelihood, providing a nonparametric estimate of the conditional mean* $\mathbb{E}[Y \mid X]$ *(and, through its Bayesian posterior, an implicit distributional object over regression functions).*
> 2. ***Survey-Weighted DPMM*** — *a Dirichlet Process Mixture Model with design-weighted posterior, providing a fully distributional nonparametric estimate of the marginal population measure* $\mathbb{P}(X)$*, proven to converge to the true finite population mixing measure in MMD norm.*
> 3. ***Survey-Weighted MMD*** — *the Maximum Mean Discrepancy between Hajek empirical measures in a reproducing kernel Hilbert space: the distributional loss, two-sample test, and calibration criterion that unifies the framework.*
>
> *BART is a nonparametric regression tool; DPMM and MMD are the distributional engines. Together they enable population-level nonparametric inference about full distributions — not just means — across causal analysis, population phenotyping, synthetic data generation, functional health curves, federated survey statistics, and clinical AI auditing.*

---

## 1. Motivation

Standard two-sample hypothesis tests — Kolmogorov-Smirnov, Wilcoxon rank-sum, or their multivariate extensions — are designed for independent and identically distributed (i.i.d.) observations. In the MEPS context, and in complex survey data generally, this assumption is fundamentally violated on two counts:

1. **Unequal probability of selection.** Each observation $i$ enters the sample with probability $\pi_i$, carrying a design weight $w_i \approx 1/\pi_i$. A raw (unweighted) test answers the question "are these two groups drawn from the same *sample* distribution?" rather than the scientifically meaningful question "do these two groups represent *populations* with the same distribution?"

2. **Dependence from the design.** Observations within the same primary sampling unit (PSU) and stratum are positively correlated. Ignoring the design inflates the effective sample size and yields invalid (anti-conservative) p-values.

The **Kernel Two-Sample Test** based on Maximum Mean Discrepancy (MMD) [Gretton et al., 2012] is the canonical non-parametric test that can detect *any* distributional difference — not just shifts in the mean or variance — in multivariate, possibly high-dimensional settings. However, it inherits the i.i.d. assumption: its null distribution is derived under that assumption, and the permutation procedure used to simulate it destroys survey design structure.

The central idea of this proposal is to **extend the MMD test to complex survey data** by:

- Replacing the classical empirical measure with the design-weighted (Hajek) empirical measure,
- Showing that this weighted empirical measure is **design-consistent** for the population kernel mean embedding,
- Deriving the null distribution of the test statistic via the **Kim-Rao stratified multinomial bootstrap** that is already implemented in this project's `survey_bart_rcpp.R` and `survey_bart.cpp` stack.

This yields an omnibus non-parametric two-sample test that:
- Is valid under complex survey designs (stratified, clustered, unequal selection probabilities),
- Is sensitive to any difference in the multivariate distribution, not just mean shifts,
- Reuses the existing Kim-Rao bootstrap infrastructure with no modifications, including the natural **out-of-bag (OOB)** observations it produces,
- Is directly applicable to covariate balance checking in the causal inference pipeline, to distributional calibration of surveys to census targets, to synthetic data evaluation, and to functional data problems such as continuous glucose monitoring curves in NHANES.

---

## 2. Background: The Standard MMD Test

### 2.1 Kernel Mean Embeddings

Let $\mathcal{X}$ be a metric space and $k: \mathcal{X} \times \mathcal{X} \to \mathbb{R}$ a symmetric positive-definite kernel (e.g., the Gaussian kernel $k(x, x') = \exp(-\|x - x'\|^2 / 2\sigma^2)$). The **reproducing kernel Hilbert space (RKHS)** $\mathcal{H}_k$ induced by $k$ is a Hilbert space of functions on $\mathcal{X}$ satisfying the reproducing property $f(x) = \langle f, k(\cdot, x) \rangle_{\mathcal{H}_k}$ for all $f \in \mathcal{H}_k$.

The **kernel mean embedding** of a probability measure $\mathbb{P}$ is the element of $\mathcal{H}_k$ defined by:

$$
\mu_{\mathbb{P}} = \mathbb{E}_{X \sim \mathbb{P}}[k(\cdot, X)] \in \mathcal{H}_k
$$

It inherits a critical property from the kernel: if $k$ is **characteristic** (the Gaussian kernel is characteristic on $\mathbb{R}^d$ for all $d$), then the embedding is **injective**:

$$
\mathbb{P} = \mathbb{Q} \iff \mu_{\mathbb{P}} = \mu_{\mathbb{Q}} \text{ in } \mathcal{H}_k
$$

This injectivity is the engine of the test: any difference between distributions $\mathbb{P}$ and $\mathbb{Q}$ — in any moment, in tail behavior, in multimodality — manifests as a gap $\|\mu_{\mathbb{P}} - \mu_{\mathbb{Q}}\|_{\mathcal{H}_k} > 0$.

### 2.2 The MMD Statistic

Given i.i.d. samples $X_1, \ldots, X_m \sim \mathbb{P}$ and $Y_1, \ldots, Y_n \sim \mathbb{Q}$, the squared MMD is:

$$
\text{MMD}^2(\mathbb{P}, \mathbb{Q}) = \|\mu_\mathbb{P} - \mu_\mathbb{Q}\|_{\mathcal{H}_k}^2
= \mathbb{E}[k(X, X')] - 2\,\mathbb{E}[k(X, Y)] + \mathbb{E}[k(Y, Y')]
$$

The standard unbiased estimator is the U-statistic:

$$
\widehat{\text{MMD}}_u^2 = \frac{1}{m(m-1)} \sum_{i \neq j} k(X_i, X_j)
- \frac{2}{mn} \sum_{i,j} k(X_i, Y_j)
+ \frac{1}{n(n-1)} \sum_{i \neq j} k(Y_i, Y_j)
$$

Under $H_0: \mathbb{P} = \mathbb{Q}$, the rescaled statistic $m \cdot \widehat{\text{MMD}}_u^2$ converges in distribution to an infinite weighted sum of centered $\chi^2_1$ variates (a degenerate U-statistic limit). This distribution is intractable analytically. The standard remedy is **permutation**: randomly shuffle the combined sample labels $B$ times, recompute the statistic each time, and use those values as the null distribution.

**The permutation approach breaks down under survey designs.** Permuting observation labels ignores that different observations have different design weights and different PSU/stratum memberships. The permutation null does not correspond to the actual sampling variability of the test statistic under $H_0$, yielding incorrect Type-I error rates.

---

## 3. The Survey-Weighted MMD

### 3.1 The Design-Weighted Empirical Measure

For a complex survey with $n$ sampled units, each unit $i$ has a design weight $w_i$ (typically approximated as $1/\pi_i$ where $\pi_i$ is the first-order inclusion probability). The **Hajek estimator** of the finite-population distribution is:

$$
\hat{\mathbb{P}}_w = \sum_{i=1}^n \tilde{w}_i \,\delta_{X_i}, \qquad \tilde{w}_i = \frac{w_i}{\sum_{j=1}^n w_j}
$$

This is a probability measure on $\mathcal{X}$ concentrated on the sample points, with mass $\tilde{w}_i$ at $X_i$. Under standard regularity conditions on the design (Isaki-Fuller conditions), for any square-integrable function $f$ in the population:

$$
\hat{\mathbb{P}}_w(f) = \sum_{i=1}^n \tilde{w}_i f(X_i) \xrightarrow{p} \mathbb{E}_{\mathbb{P}_{\rm pop}}[f(X)]
$$

where convergence is in the **design-based** sense: randomness comes from the sampling mechanism, not a superpopulation model. This property is called **design consistency**.

### 3.2 Design Consistency of the Survey-Weighted Kernel Mean Embedding

Define the survey-weighted empirical kernel mean embedding as the plug-in estimator:

$$
\hat{\mu}_w = \sum_{i=1}^n \tilde{w}_i \,k(\cdot, X_i) \in \mathcal{H}_k
$$

**Claim:** If $k$ is a bounded, continuous, characteristic kernel ($\sup_{x,x'} k(x,x') \leq B < \infty$), then $\hat{\mu}_w \to \mu_{\mathbb{P}_{\rm pop}}$ in $\mathcal{H}_k$-norm under the design-based probability.

*Sketch of argument:* For any $f \in \mathcal{H}_k$ with $\|f\|_{\mathcal{H}_k} \leq 1$, the reproducing property gives $f(x) = \langle f, k(\cdot,x)\rangle_{\mathcal{H}_k}$. Since $k$ is bounded, $\|f\|_\infty \leq B\|f\|_{\mathcal{H}_k} \leq B$. Thus every element of the unit ball of $\mathcal{H}_k$ is a bounded function, and Hajek consistency applies uniformly over this class. Therefore:

$$
\langle \hat{\mu}_w, f\rangle_{\mathcal{H}_k} = \hat{\mathbb{P}}_w(f) \xrightarrow{p} \mathbb{E}_{\mathbb{P}_{\rm pop}}[f(X)] = \langle \mu_{\mathbb{P}_{\rm pop}}, f\rangle_{\mathcal{H}_k}
$$

Taking the supremum over $\|f\|_{\mathcal{H}_k} \leq 1$ gives $\|\hat{\mu}_w - \mu_{\mathbb{P}_{\rm pop}}\|_{\mathcal{H}_k} \to 0$ in probability. $\square$

**This is the central theoretical result:** the survey-weighted kernel mean embedding is design-consistent for the population embedding. It justifies replacing the classical empirical embedding with the weighted version in the MMD test.

### 3.3 The Survey-Weighted MMD Statistic

Let group $A$ (e.g., treated) have $n_A$ sampled units with survey weights $w_i^A$ and let $\tilde{w}_i^A = w_i^A / \sum_j w_j^A$. Similarly for group $B$ (control) with $\tilde{w}_j^B$. Let the kernel matrices be:

$$
K_{AA}[i,j] = k(X_i^A, X_j^A), \quad
K_{BB}[i,j] = k(X_i^B, X_j^B), \quad
K_{AB}[i,j] = k(X_i^A, X_j^B)
$$

The **survey-weighted squared MMD** is:

$$
\boxed{
\widehat{\rm MMD}_w^2 \;=\;
\tilde{\mathbf{w}}_A^\top K_{AA}\,\tilde{\mathbf{w}}_A
\;+\;
\tilde{\mathbf{w}}_B^\top K_{BB}\,\tilde{\mathbf{w}}_B
\;-\;
2\,\tilde{\mathbf{w}}_A^\top K_{AB}\,\tilde{\mathbf{w}}_B
}
$$

This is exactly $\|\hat{\mu}_{w,A} - \hat{\mu}_{w,B}\|_{\mathcal{H}_k}^2$ — the squared RKHS norm between the two survey-weighted empirical embeddings.

**Immediate consequences of the design-consistency result:**

- Under $H_0: \mathbb{P}_A = \mathbb{Q}_B$ (the two population distributions are equal), $\widehat{\rm MMD}_w^2 \xrightarrow{p} 0$.
- Under $H_1: \mathbb{P}_A \neq \mathbb{Q}_B$, $\widehat{\rm MMD}_w^2 \xrightarrow{p} \|\mu_{\mathbb{P}_A} - \mu_{\mathbb{Q}_B}\|_{\mathcal{H}_k}^2 > 0$.
- The test is **consistent** (power $\to 1$ as sample sizes grow under any fixed $H_1$).

---

## 4. Null Distribution via the Kim-Rao Bootstrap

### 4.1 Why Permutation Fails

Standard MMD tests simulate the null by permuting group labels across all $n_A + n_B$ observations. This works because under $H_0$, i.i.d. observations are exchangeable: the joint distribution is invariant to label permutation.

Under a stratified cluster design, **exchangeability fails** because:
- Observations in different strata carry different design weights.
- Observations within the same PSU are correlated (they share cluster-level characteristics).
- Randomly reassigning group labels ignores these correlations and the weight structure.

The permutation null distribution does not reflect the actual variability of $\widehat{\rm MMD}_w^2$ under the survey design, leading to either inflated or deflated Type-I error rates.

### 4.2 The Kim-Rao Stratified Multinomial Bootstrap

Kim & Rao (2009) proposed a bootstrap for complex surveys that correctly propagates design uncertainty by resampling PSUs *within strata*, preserving the design structure. In stratum $h$ with $n_h$ PSUs:

1. Draw multinomial counts: $(m_{h1}, \ldots, m_{h,n_h}) \sim \text{Multinomial}(n_h - 1,\; 1/n_h, \ldots, 1/n_h)$.
2. Assign PSU-level multiplicity: $r_{hi} = \frac{n_h}{n_h-1} m_{hi}$ (the factor $k_h = n_h/(n_h-1)$ ensures mean-preservation: $\mathbb{E}[r_{hi}] = 1$).
3. Construct bootstrap weights: $w_i^{*(b)} = r_{h[i]}^{(b)} \cdot w_i$.

The bootstrap weight for unit $i$ in bootstrap replicate $b$ is $w_i^{*(b)}$. The entire bootstrap weight matrix is the $B \times n$ matrix `W_mat` in the current codebase.

**The key insight for the MMD test:** the Kim-Rao bootstrap correctly approximates the sampling distribution of any smooth functional of the weighted empirical measure. Since the weighted MMD is a smooth functional of $\hat{\mathbb{P}}_w$ (it is a bounded continuous quadratic form), the bootstrap variance of $\widehat{\rm MMD}_w^2$ is design-consistent.

### 4.3 Connection to the Existing BART Implementation

The Kim-Rao weight matrix is already generated in `survey_bart_rcpp.R`:

```r
# Lines 79–83 of survey_bart_rcpp.R
cat("Generating Kim-Rao bootstrap weight matrix...\n")
tmp_df <- data.frame(weight = weights, strata = strata, psu = psu)
W_mat  <- generate_kim_rao_weights(tmp_df, "strata", "psu", "weight", n_mcmc)
```

Each row $b$ of `W_mat` is passed to the C++ MCMC backend (`survey_bart.cpp`, line 363: `vec w_star = W_mat.row(iter).t()`) to run the weighted likelihood for that iteration. **For the MMD test, we reuse this exact matrix with zero code changes.** We simply call `generate_kim_rao_weights()` on each group subset and loop over the bootstrap rows to compute the test statistic distribution.

### 4.4 Algorithm: Survey-Weighted MMD Test

**Inputs:**
- Group $A$ data: covariate matrix $X^A \in \mathbb{R}^{n_A \times p}$, design frame with columns `weight`, `strata`, `psu`.
- Group $B$ data: covariate matrix $X^B \in \mathbb{R}^{n_B \times p}$, analogous design frame.
- Kernel $k$ with bandwidth $\sigma$ (set by weighted median heuristic if not specified).
- Number of bootstrap replicates $B$.

**Step 1 — Observed statistic.**
Normalize weights $\tilde{w}_i^A = w_i^A / \sum_j w_j^A$, $\tilde{w}_j^B = w_j^B / \sum_j w_j^B$.
Compute $T_0 = \widehat{\rm MMD}_w^2(\tilde{\mathbf{w}}_A, \tilde{\mathbf{w}}_B)$ from the formula in §3.3.

**Step 2 — Bootstrap weights.**
Generate bootstrap weight matrices:
$$W^A \in \mathbb{R}^{B \times n_A} \quad \text{via} \quad \texttt{generate\_kim\_rao\_weights}(\text{data}_A, B)$$
$$W^B \in \mathbb{R}^{B \times n_B} \quad \text{via} \quad \texttt{generate\_kim\_rao\_weights}(\text{data}_B, B)$$

**Step 3 — Bootstrap statistics.**
For $b = 1, \ldots, B$: normalize row $b$ of each matrix and compute
$$T_b = \widehat{\rm MMD}_w^2\!\left(\widetilde{W^A_{b\cdot}},\; \widetilde{W^B_{b\cdot}}\right)$$

**Step 4 — p-value.**
$$p = \frac{1}{B} \sum_{b=1}^B \mathbf{1}[T_b \geq T_0]$$

Reject $H_0: \mathbb{P}_A = \mathbb{Q}_B$ at level $\alpha$ if $p \leq \alpha$.

### 4.5 Out-of-Bag Observations from the Kim-Rao Bootstrap

A structural feature of the Kim-Rao bootstrap that is directly useful for the MMD test is the **natural production of out-of-bag (OOB) observations**. In bootstrap replicate $b$, a PSU $i$ in stratum $h$ receives multiplicity $m_{hi}^{(b)} = 0$ with probability:

$$P(m_{hi}^{(b)} = 0) = \left(1 - \frac{1}{n_h}\right)^{n_h - 1} \xrightarrow{n_h \to \infty} e^{-1} \approx 0.368$$

All units in such a PSU get bootstrap weight $w_i^{*(b)} = 0$ — they are **OOB** for replicate $b$. In the current BART implementation, these are the observations that receive zero weight in the weighted likelihood for that MCMC iteration (`w_star(i) == 0` in `survey_bart.cpp`, lines 148–149), meaning the tree update is driven entirely by the in-bag units, while OOB units can serve as a hold-out.

For the MMD test, the OOB observations provide two key advantages:

**OOB variance estimation.** For each replicate $b$, partition the sample into in-bag (IB) units ($w_i^{*(b)} > 0$) and OOB units ($w_i^{*(b)} = 0$). Compute the MMD on the IB units $T_b^{\rm IB}$ and a separate MMD on the OOB units $T_b^{\rm OOB}$. The variability of $T_b^{\rm OOB}$ across replicates is a design-honest cross-validated estimate of the null variance of $T_0$, because OOB units were not used to construct the bootstrap weights in that replicate.

**OOB model checking for Survey BART.** In the BART context, the OOB units in each iteration can be used to compute a survey-weighted MMD between the OOB posterior predictive distribution and the observed OOB outcomes. If the model is calibrated, this MMD should be small. This gives a distributional goodness-of-fit diagnostic for Survey BART that goes beyond pointwise coverage:

$$T_{\rm GoF}^{(b)} = \widehat{\rm MMD}_w^2\!\left(\hat{\mathbb{P}}_{\hat{Y}, {\rm OOB}}^{(b)},\; \hat{\mathbb{P}}_{Y, {\rm OOB}}^{(b)}\right)$$

where $\hat{\mathbb{P}}_{\hat{Y}, {\rm OOB}}^{(b)}$ is the weighted empirical measure of predictions on OOB units and $\hat{\mathbb{P}}_{Y, {\rm OOB}}^{(b)}$ is the weighted empirical measure of their observed outcomes. Averaging over $b$ gives a stable model diagnostic that accounts for the survey design.

### 4.6 Theoretical Validity

The bootstrap validity of this procedure follows from three results:

| Result | Statement |
|---|---|
| **Design consistency (§3.2)** | $\widehat{\rm MMD}_w^2 \xrightarrow{p} {\rm MMD}^2(\mathbb{P}_A, \mathbb{Q}_B)$ under the design |
| **Kim-Rao bootstrap consistency** | The bootstrap variance of any smooth statistic $T(\hat{\mathbb{P}}_w)$ is design-consistent (Kim & Rao, 2009) |
| **Asymptotic null distribution** | Under $H_0$, the bootstrap distribution of $T_b - T_0$ consistently approximates the null distribution of $T_0$ |
| **OOB consistency** | OOB statistics are asymptotically independent of the in-bag bootstrap statistic, yielding unbiased variance estimates |

The third result follows from the general theory of bootstrap approximations for degenerate U-statistics under design-based sampling. The regularity conditions required are: (a) the kernel is bounded and continuous; (b) the design satisfies standard Isaki-Fuller conditions; (c) PSU counts per stratum $n_h \to \infty$.

---

## 5. Connection to the Causal Inference Pipeline

### 5.1 Multivariate Covariate Balance Testing

The most immediate application is **omnibus covariate balance testing** in the MEPS causal analysis (`estimate_causal_meps.R`). The standard practice of reporting per-covariate standardized mean differences (SMDs) has three weaknesses:

1. It only checks first-moment balance — two groups can have identical means but very different variances, joint distributions, or tail behavior.
2. SMDs are correlated across covariates, so there is no principled multiple-comparison control.
3. The survey design is typically ignored in the standard error of the SMD.

The survey-weighted MMD resolves all three: it is omnibus (detects any distributional difference), provides a single valid p-value for the joint null, and is calibrated to the MEPS design via the Kim-Rao bootstrap.

**Proposed workflow:**
1. Before any weighting: test $H_0: \mathbb{P}(X | A=1) = \mathbb{P}(X | A=0)$ in the finite population. Expect rejection.
2. After IPTW reweighting (multiply $w_i$ by $1/\hat{e}(X_i)$ for treated, $1/(1-\hat{e}(X_i))$ for controls): rerun the test. Non-rejection is formal evidence of multivariate balance.
3. After BART G-computation: verify that the marginal distributions of $\hat{Y}(1)$ and $\hat{Y}(0)$ differ significantly (to confirm treatment effect exists), then test whether a specific hypothesis about distributional equality holds.

### 5.2 Distributional Causal Estimands

Beyond ATE/ATT (which are mean-level estimands), the survey-weighted MMD enables testing **distributional causal hypotheses**. After running `fit_bart_cpp` with `X_test` set to the counterfactual design matrices, you obtain posterior draws $\hat{Y}^{(b)}(1)$ and $\hat{Y}^{(b)}(0)$ for $b = 1, \ldots, B$ post-burn-in iterations. The distributional hypothesis

$$H_0: \mathbb{P}(Y(1)) = \mathbb{P}(Y(0)) \text{ in the finite population}$$

can be tested by computing $\widehat{\rm MMD}_w^2$ between the collections of potential outcome draws, with bootstrap uncertainty propagated through the Kim-Rao weights. This asks whether the treatment shifts the *entire distribution* of healthcare expenditures, not just the mean — a far richer causal claim.

### 5.3 BART Leaf Kernel

An appealing alternative to the Gaussian kernel is the **BART random forest kernel**:

$$k_{\rm BART}(x, x') = \frac{1}{T} \sum_{t=1}^T \mathbf{1}[\text{leaf}(x, t) = \text{leaf}(x', t)]$$

where $\text{leaf}(x, t)$ is the leaf node to which $x$ is assigned in tree $t$. This kernel is extractable directly from the `Node::route()` function in `survey_bart.cpp` (lines 71–83). The BART leaf kernel is data-adaptive: it places high similarity on pairs $(x, x')$ that are co-classified in most trees, automatically reflecting the covariate importance structure learned from the data. This is particularly appropriate when the same BART model is used both for causal estimation and for the balance test, as the kernel is calibrated to the same covariate space.

---

## 6. Broader Applications Beyond Causal Inference

The survey-weighted kernel mean embedding is not merely a device for causal diagnostics. Its design-consistency property makes it a broadly applicable tool for any problem that requires comparing or matching distributions in complex survey data. We develop four substantial applications here.

---

### 6.1 Synthetic Data Generation and Evaluation

#### 6.1.1 The Problem

Synthetic data generation for survey data is an emerging area driven by privacy concerns (releasing a synthetic MEPS or NHANES dataset that preserves population-level statistics without disclosing individual records) and data augmentation (generating additional training observations for rare subgroups). A synthetic data generator $G$ takes a noise vector $z \sim \mathcal{N}(0,I)$ and produces a synthetic observation $\tilde{X} = G(z)$. The generator is considered good if the synthetic population distribution $\mathbb{P}_G$ matches the target finite population distribution $\mathbb{P}_{\rm pop}$.

The difficulty in the survey context is that $\mathbb{P}_{\rm pop}$ is not directly observed — only the weighted empirical measure $\hat{\mathbb{P}}_w$ is available, and it encodes population-level information through the design weights. A generator trained by minimizing the unweighted empirical distribution will learn the *sample* distribution, not the *population* distribution, systematically over-representing oversampled subgroups (e.g., elderly or minority populations that are oversampled in MEPS/NHANES for precision).

#### 6.1.2 The Survey-Weighted MMD as a Generative Loss

The survey-weighted MMD provides a natural and theoretically grounded training objective. Define the one-sample survey-weighted MMD between the generator distribution $\mathbb{P}_G$ and the weighted survey measure as:

$$\mathcal{L}_{\rm MMD}(G) = \widehat{\rm MMD}_w^2(\hat{\mathbb{P}}_G,\, \hat{\mathbb{P}}_w)$$

$$= \frac{1}{M^2}\sum_{i,j} k(\tilde{X}_i, \tilde{X}_j)
- \frac{2}{M}\sum_{i=1}^M \sum_{j=1}^n \tilde{w}_j\, k(\tilde{X}_i, X_j)
+ \tilde{\mathbf{w}}^\top K \,\tilde{\mathbf{w}}$$

where $\tilde{X}_1, \ldots, \tilde{X}_M \sim G$ are synthetic samples and $\tilde{w}_j = w_j / \sum_l w_l$ are the normalized survey weights. This is the loss function of a **survey-weighted Maximum Mean Discrepancy GAN (SW-MMD-GAN)**.

**Theoretical guarantee:** Since $\hat{\mathbb{P}}_w$ is design-consistent for $\mathbb{P}_{\rm pop}$, minimizing $\mathcal{L}_{\rm MMD}(G)$ asymptotically drives $\mathbb{P}_G \to \mathbb{P}_{\rm pop}$ (the finite population distribution), not the sample distribution. For $M \to \infty$ and $n \to \infty$, $\mathcal{L}_{\rm MMD}(G) = 0 \iff G$ generates from $\mathbb{P}_{\rm pop}$.

**Implementation.** The cross-term $\frac{2}{M}\sum_i \sum_j \tilde{w}_j k(\tilde{X}_i, X_j)$ is the only term involving the survey weights; it acts as a weighted expectation of the kernel evaluated at synthetic–real pairs. The last term $\tilde{\mathbf{w}}^\top K \tilde{\mathbf{w}}$ is a constant (with respect to the generator) and does not affect gradient computation. This makes the loss differentiable and directly pluggable into any gradient-based generative framework (copula models, VAEs, normalizing flows).

#### 6.1.3 Evaluation: Using the Survey-Weighted MMD as a Fidelity Score

Beyond training, the survey-weighted MMD is a rigorous **distributional fidelity score** for evaluating synthetic datasets:

- Compute $T_0 = \widehat{\rm MMD}_w^2$ between the synthetic data (uniformly weighted) and the real survey data (design-weighted).
- Use the Kim-Rao bootstrap on the real survey data to generate a null distribution (what $T_0$ would look like if the synthetic data were drawn from the true population).
- Report a p-value: failure to reject means the synthetic data is distributionally indistinguishable from the survey population.

This is strictly superior to common synthetic data metrics like propensity score mean squared error (pMSE), which only check whether a classifier can distinguish real from synthetic data — a first-moment criterion that ignores joint distributional fidelity.

---

### 6.2 Distributional Calibration of Sample Surveys to Census Data

#### 6.2.1 The Calibration Problem

Post-stratification and raking (iterative proportional fitting) are the standard tools for calibrating survey weights so that the weighted sample totals match known population totals from census data. Formally, the calibration estimator finds weights $\{w_i^{\rm cal}\}$ solving:

$$\min_w \sum_i \frac{(w_i^{\rm cal} - w_i^{\rm base})^2}{w_i^{\rm base}} \quad \text{subject to} \quad \sum_i w_i^{\rm cal} x_{ij} = N \bar{X}_j^{\rm census}, \quad j = 1, \ldots, p$$

where $\bar{X}_j^{\rm census}$ is the $j$-th marginal total from the census. The constraint is **linear**: it only enforces that the weighted survey means match the census means, one variable at a time.

**The limitation:** raking enforces univariate (and sometimes bivariate) marginal balance. It cannot enforce that the full multivariate distribution of $(X_1, \ldots, X_p)$ in the survey matches the full joint distribution in the census. Two calibrated surveys can have identical marginal means but very different joint distributions (e.g., different correlation structures, different tail behavior).

#### 6.2.2 MMD as a Distributional Calibration Criterion

Suppose we have access to a census microdata extract or a very large reference sample representing the target population (e.g., the American Community Survey for MEPS calibration, or the Census Bureau's Current Population Survey). Denote the reference sample as $\{Z_1, \ldots, Z_N\}$ with uniform weights $1/N$.

Define the **MMD calibration problem**:

$$\min_{\tilde{w} \in \Delta_n} \widehat{\rm MMD}^2\!\left(\sum_i \tilde{w}_i \delta_{X_i},\; \frac{1}{N}\sum_j \delta_{Z_j}\right)$$

$$= \min_{\tilde{w} \in \Delta_n} \left[
\tilde{\mathbf{w}}^\top K_{XX} \tilde{\mathbf{w}}
- \frac{2}{N} \mathbf{1}_n^\top K_{XZ} \mathbf{1}_N
\right]$$

where $\Delta_n = \{\tilde{w} \geq 0 : \sum_i \tilde{w}_i = 1\}$ is the probability simplex. This is a **quadratic program (QP)** in $\tilde{w}$: the objective is a convex quadratic function (since $K_{XX}$ is positive semi-definite), and the constraint set is a convex polytope. It can be solved efficiently by standard QP solvers for moderate $n$ (e.g., $n \leq 10{,}000$ for MEPS).

**What this achieves:** the MMD-calibrated weights $\tilde{w}^{\rm MMD}$ minimize the full kernel-distance between the weighted survey distribution and the census reference distribution, simultaneously matching all moments of all functions in the RKHS. This is strictly stronger than raking, which only matches the linear functions $f(x) = x_j$.

#### 6.2.3 Connection to the Kim-Rao Bootstrap

After computing $\tilde{w}^{\rm MMD}$, the Kim-Rao bootstrap provides valid uncertainty quantification for any functional of the MMD-calibrated empirical measure. The bootstrap algorithm remains identical — we substitute $\tilde{w}^{\rm MMD}$ as the base weights $w_i^{\rm base}$ in `generate_kim_rao_weights()`. The OOB structure and the theoretical consistency results of §4.5–4.6 carry through unchanged, because they only require that the base weights are fixed (non-random) and design-consistent.

#### 6.2.4 Diagnostic Role: Testing Whether Calibration Succeeded

After raking or post-stratification, one wants to know: has the calibration made the survey distribution match the census? The survey-weighted MMD test answers this formally:

- **Group A:** the calibrated survey sample with weights $w_i^{\rm cal}$.
- **Group B:** the census microdata or reference sample with uniform weights.
- **Test:** $H_0: \mathbb{P}_{\rm survey}^{\rm cal} = \mathbb{P}_{\rm census}$.

A failure to reject is evidence of successful calibration in the full distributional sense. A rejection reveals residual distributional mismatch even after calibration — something raking diagnostics (which only check marginal totals) cannot detect.

---

### 6.3 Functional Data: Modeling Health Curves in NHANES and MEPS

#### 6.3.1 Functional Observations in Complex Surveys

A growing class of health data involves **functional observations** — entire curves or trajectories measured for each individual — collected under complex survey designs. Three prominent examples:

1. **Continuous Glucose Monitoring (CGM) in NHANES.** NHANES 2011–2014 included a CGM sub-study where participants wore a glucose sensor for up to 7 days. Each participant's data is a glucose trajectory $G_i: [0, T] \to \mathbb{R}_{>0}$ sampled at 5-minute intervals. The NHANES design includes oversampling of non-Hispanic Black, Hispanic, and older adults, so the CGM sub-sample carries design weights.

2. **Physical activity accelerometry (NHANES).** NHANES has also collected wrist accelerometer data, producing minute-level physical activity count curves over several days.

3. **Multi-visit lab trajectories (MEPS).** MEPS follows individuals across multiple panel rounds, producing trajectories of healthcare utilization, expenditure, and reported outcomes over time.

Standard functional data analysis (FDA) methods (FPCA, functional regression) assume i.i.d. curves. They do not account for the NHANES/MEPS design, and applying them naively to the raw sample will produce estimates representative of the survey sample composition rather than the U.S. population.

#### 6.3.2 Kernel Mean Embeddings for Functional Data

The RKHS framework extends naturally to functional observations. Let $\mathcal{X} = L^2([0,T])$ be the space of square-integrable functions. A kernel on function space can be defined as:

$$k_F(f, g) = \exp\!\left(-\frac{\|f - g\|_{L^2}^2}{2\sigma^2}\right), \qquad \|f - g\|_{L^2}^2 = \int_0^T (f(t) - g(t))^2\,dt$$

For discretely observed curves (as in CGM), approximate $\|f - g\|_{L^2}^2 \approx \frac{T}{M}\sum_{m=1}^M (f(t_m) - g(t_m))^2$ at observation times $t_1, \ldots, t_M$. This is a kernel on the space of CGM trajectories that is:
- **Characteristic** on $L^2([0,T])$ (a result due to Sriperumbudur et al., 2010),
- **Bounded** (since $k_F \leq 1$),
- Invariant to curve amplitude scale with appropriate normalization.

The **survey-weighted kernel mean embedding of CGM trajectories** is:

$$\hat{\mu}_w = \sum_{i=1}^n \tilde{w}_i\, k_F(\cdot, G_i) \in \mathcal{H}_{k_F}$$

By exactly the same design-consistency argument as §3.2 (the kernel is bounded and continuous on $L^2$), $\hat{\mu}_w$ consistently estimates the population kernel mean embedding of the CGM trajectory distribution.

#### 6.3.3 Survey-Weighted Two-Sample Test for CGM Trajectories

The survey-weighted MMD test can now answer population-level questions about entire glucose trajectory distributions:

**Example hypothesis:** "Do the population-level glucose trajectory distributions differ between diabetics and non-diabetics in the U.S., after accounting for the NHANES design?"

- **Group A:** NHANES participants with diagnosed diabetes, design weights $w_i^A$, CGM trajectories $G_i^A$.
- **Group B:** NHANES participants without diabetes, design weights $w_j^B$, CGM trajectories $G_j^B$.

The kernel matrix entries are:
$$K_{AA}[i,j] = k_F(G_i^A, G_j^A) = \exp\!\left(-\frac{T}{M\cdot 2\sigma^2}\sum_{m=1}^M (G_i^A(t_m) - G_j^A(t_m))^2\right)$$

The rest of the algorithm (§4.4) applies verbatim — the data is now a matrix of shape $n \times M$ (participants × time points) rather than $n \times p$ (participants × covariates), but the test statistic formula $\tilde{\mathbf{w}}_A^\top K_{AA}\tilde{\mathbf{w}}_A + \tilde{\mathbf{w}}_B^\top K_{BB}\tilde{\mathbf{w}}_B - 2\tilde{\mathbf{w}}_A^\top K_{AB}\tilde{\mathbf{w}}_B$ is unchanged.

**This is a strictly stronger test than:**
- Comparing mean CGM curves (which ignores variability, diurnal pattern shape, spike frequency),
- Comparing scalar summary statistics per person (mean glucose, HbA1c proxy, time-in-range) — which throws away the full trajectory information,
- Applying standard functional ANOVA without accounting for the NHANES design.

#### 6.3.4 Distributional Regression for Functional Outcomes

A natural extension is **distributional regression** on functional outcomes: rather than predicting the mean CGM curve given covariates $X$, estimate how the full distribution of CGM trajectories $\mathbb{P}(G | X)$ varies with covariates. The survey-weighted conditional kernel mean embedding:

$$\hat{\mu}_w(x) = \frac{\sum_i \tilde{w}_i\, k_X(x, X_i)\, k_F(\cdot, G_i)}{\sum_i \tilde{w}_i\, k_X(x, X_i)}$$

where $k_X$ is a covariate kernel and $k_F$ is the trajectory kernel, provides a Nadaraya-Watson-type estimate of the conditional distribution of CGM trajectories at covariate value $x$, weighted by the survey design. Two groups can then be compared using the survey-weighted MMD between $\hat{\mu}_w(x_A)$ and $\hat{\mu}_w(x_B)$ at different covariate profiles.

#### 6.3.5 Amplitude vs. Phase Variation in CGM Curves

Glucose curves exhibit both **amplitude variation** (higher or lower glucose levels) and **phase variation** (postprandial spikes occurring at different times due to meal timing). The $L^2$ kernel above conflates these. A better kernel for CGM uses the **elastic shape metric** (square-root velocity framework, Srivastava & Klassen, 2016), which is invariant to time-warping:

$$k_{\rm shape}(f, g) = \exp\!\left(-\frac{d_{\rm SRV}(f, g)^2}{2\sigma^2}\right)$$

where $d_{\rm SRV}$ is the geodesic distance after optimal time-warping registration. This produces a kernel that compares curve shapes separately from their timing — more appropriate for CGM data where meal timing differs across individuals. The RKHS theory still applies (the SRV metric defines a kernel on the quotient manifold of curves modulo reparametrization), so the full survey-weighted MMD framework goes through.

---

### 6.5 Population-Level Algorithmic Fairness Auditing of Clinical AI

> **This is the new out-of-the-box application.** It is completely absent from existing fairness methodology and represents a frontier problem at the intersection of health survey statistics, clinical AI regulation, and distributional RKHS methods.

#### 6.5.1 The Problem: Hospital Data is Not the Population

Clinical prediction models — sepsis early-warning systems, 30-day readmission scores, mammography AI, diabetes complication risk scores — are trained and validated on **electronic health record (EHR) data** from hospital systems. They are then approved by regulators (FDA Software as a Medical Device, CE marking in the EU), purchased by health systems, and deployed on the **national patient population**.

The fundamental representativeness failure is this: hospital patients are **not** a random sample of the U.S. population. They are:
- Sicker (selection by illness severity),
- More urban (proximity to tertiary care centers),
- Skewed toward certain payers (Medicare/Medicaid vs. privately insured),
- Systematically under-representative of rural, uninsured, and minority communities.

Current algorithmic fairness audits — including those required under the FDA's AI/ML SaMD Action Plan and recent HHS algorithmic bias regulations — are conducted on the hospital's own EHR validation cohort. This validation cohort inherits all the same representativeness failures. A model that passes a hospital-level fairness audit may exhibit **deep distributional unfairness at the population level** in ways that are entirely invisible to the auditor, because the population subgroups most affected are the ones least represented in the hospital's data.

**The core insight:** the fairness question that *matters* is not "are predictions fair in my hospital's patients?" but rather "are predictions fair *as experienced by the national population*?" These are different questions with potentially very different answers.

#### 6.5.2 Connecting Hospital Models to Population Surveys via Kernel Embeddings

The proposed framework bridges clinical AI evaluation to nationally representative health surveys (MEPS, NHANES, BRFSS, SIPP) via the following observation:

Let $f: \mathcal{X} \to [0,1]$ be a clinical AI model mapping patient covariate vector $x$ to a predicted risk score (e.g., probability of 30-day readmission). Let $G$ be a demographic group (e.g., non-Hispanic Black adults), and let $\mathbb{P}_{G,\rm pop}$ be the population-level distribution of covariate vectors for group $G$ in the U.S. finite population.

The **population-level prediction distribution** for group $G$ is the pushforward measure:
$$\mathbb{P}_{f,G} = f_\#\, \mathbb{P}_{G,\rm pop}$$

i.e., the distribution of the AI model's output score when applied to a randomly drawn member of group $G$ from the national population. **This is the fairness-relevant object.** It is not the distribution of scores on the hospital's Black patients.

The survey-weighted empirical approximation of $\mathbb{P}_{f,G}$ is:
$$\hat{\mathbb{P}}_{f,G,w} = \sum_{i : G_i = G} \tilde{w}_i^G \,\delta_{f(X_i)}$$

where the sum is over MEPS/NHANES respondents in group $G$, $\tilde{w}_i^G$ are normalized design weights within group $G$, and $f(X_i)$ is the AI model's predicted score applied to the survey respondent's covariate vector $X_i$. By the design-consistency theorem (§3.2), $\hat{\mathbb{P}}_{f,G,w}$ consistently estimates $\mathbb{P}_{f,G}$ as long as the kernel is bounded.

#### 6.5.3 The Population-Level Fairness Test

The **population-level algorithmic fairness hypothesis** is:
$$H_0^{\rm fair}: \mathbb{P}_{f,G_1} = \mathbb{P}_{f,G_2}$$

for two demographic groups $G_1$ and $G_2$ (e.g., non-Hispanic White vs. non-Hispanic Black adults). This is a **distributional** fairness criterion, far stronger than standard fairness metrics:

| Fairness Metric | What it checks | What it misses |
|---|---|---|
| Demographic parity | $\mathbb{E}[f(X) | G_1] = \mathbb{E}[f(X) | G_2]$ | Variance, tail behavior, bimodality |
| Equalized odds | Equal TPR/FPR at one threshold | Different risk score distributions at other thresholds |
| Calibration | $P(Y=1 | f(X)=p, G) = p$ | Joint distributional differences not captured by calibration |
| **Survey-weighted MMD** | $\mathbb{P}_{f,G_1} = \mathbb{P}_{f,G_2}$ in $\mathcal{H}_k$ | **Nothing** — omnibus test |

The test statistic is:
$$T_0 = \widehat{\rm MMD}_w^2\!\left(\hat{\mathbb{P}}_{f,G_1,w},\; \hat{\mathbb{P}}_{f,G_2,w}\right)$$

computed from:
$$K_{G_1 G_1}[i,j] = k(f(X_i^{G_1}),\, f(X_j^{G_1})), \quad K_{G_2 G_2}, \quad K_{G_1 G_2} \text{ analogously}$$

Note that $f(X_i)$ are **not observed data** — they are computed by applying the clinical AI model to the survey respondent's EHR-compatible covariate vector. This requires that the survey includes the covariates used as inputs to $f$ (which is achievable for MEPS/NHANES for most structured clinical risk scores). The null distribution of $T_0$ is obtained via the Kim-Rao bootstrap on the two group-specific subsamples, exactly as in §4.4.

#### 6.5.4 What the Test Reveals That Existing Audits Cannot

The power of this approach lies in detecting **tail and shape fairness violations** that are invisible to standard metrics:

**Example 1: Bimodal unfairness.** Suppose a sepsis prediction model assigns scores that are unimodally distributed around $0.2$ for White patients but bimodally distributed (a large mass near $0.05$ and a smaller mass near $0.5$) for Black patients. Demographic parity might be satisfied (both group means could be equal), equalized odds might be satisfied at the median threshold, but the prediction distributions are clearly different — and the bimodality for Black patients implies that a large subgroup is systematically over-alerted while another subgroup is systematically missed. The survey-weighted MMD will detect this; demographic parity will not.

**Example 2: High-risk tail unfairness.** For a readmission risk model, suppose the 90th percentile of the score distribution is systematically higher for Hispanic patients than for non-Hispanic White patients conditional on identical health status. Mean-based fairness metrics miss this entirely; the survey-weighted MMD captures the tail difference as a non-zero kernel distance.

**Example 3: Population vs. hospital fairness reversal.** It is entirely possible (and likely in practice) that a model passes a hospital-level fairness audit (because the hospital serves a relatively homogeneous urban population) but fails the population-level MMD fairness test (because rural minority patients, who are unrepresented in the hospital, have a very different covariate profile that maps to systematically different risk scores under $f$). The survey-weighted MMD is the *only* approach that can detect this reversal, because it uses the design weights to represent the full population distribution.

#### 6.5.5 Intersectional Fairness via Conditional MMD

Intersectional fairness asks whether the model is fair across combinations of protected attributes (e.g., Black women over 65, rural Hispanic men). The conditional kernel mean embedding (§6.3.4) extends naturally:

$$\hat{\mu}_{f,w}(g) = \frac{\sum_i \tilde{w}_i\, k_G(g, G_i)\, k(\cdot, f(X_i))}{\sum_i \tilde{w}_i\, k_G(g, G_i)}$$

where $k_G$ is a kernel on the demographic group space (e.g., a product of categorical kernels over race, sex, age group, rurality). The **conditional survey-weighted MMD** between $\hat{\mu}_{f,w}(g_1)$ and $\hat{\mu}_{f,w}(g_2)$ at two intersectional profiles $g_1$ and $g_2$ tests whether the model's prediction distribution differs for those two population subgroups. The Kim-Rao bootstrap provides valid confidence intervals for the conditional MMD across all intersectional cells simultaneously.

This is the only framework that achieves:
(a) **Intersectionality** (joint conditioning on multiple protected attributes),
(b) **Population representativeness** (via survey design weights),
(c) **Omnibus distributional testing** (not just mean/threshold comparisons),
in a single, coherent, computationally feasible framework.

#### 6.5.6 Regulatory Implications and Path to Practice

The timing for this application is exceptional. Regulatory pressure for population-level fairness auditing is building rapidly:

- **FDA (2021–2024):** The AI/ML-Based Software as a Medical Device Action Plan requires pre-market and post-market monitoring for performance across subpopulations, but does not specify how to conduct population-level audits.
- **HHS Office for Civil Rights (2022–2024):** Guidance on algorithmic bias in federally funded programs (Section 1557 of the ACA) requires audits of AI systems used in health care but provides no statistical methodology.
- **EU AI Act (2024):** High-risk AI systems in health care require "testing for reasonably foreseeable risks" including across demographic groups, using representative data.

None of these frameworks specifies **how** to conduct distributional fairness testing on nationally representative populations. The survey-weighted MMD fills this gap with a statistically rigorous, computationally feasible, and regulatory-interpretable methodology.

**Practical path to implementation:**
1. Select a widely-deployed clinical AI model with a published covariate set (e.g., Epic's readmission model, CHADS-VASc for stroke risk, the LACE+ index for readmission).
2. Identify the overlapping covariates between the model's input set and a complex survey (MEPS has insurance, utilization, chronic conditions, demographics — covering most structured risk scores).
3. Apply the model $f$ to MEPS survey respondents' covariate vectors $X_i$.
4. Run the survey-weighted MMD fairness test between demographic groups using Kim-Rao bootstrap weights.
5. Report: population-level fairness p-value, bootstrap confidence interval for the distributional fairness gap, and a visualization of the group-specific prediction distributions $\hat{\mathbb{P}}_{f,G,w}$.

This workflow requires no changes to `generate_kim_rao_weights()` or `compute_weighted_mmd2()` — it is a direct application of the infrastructure already built.

---

### 6.7 Temporal Population Drift Detection and Survey-Wave Transfer

> **This application uses the survey-weighted MMD as a distributional distance function across survey time points**, enabling formal detection of when a population has changed enough to invalidate models trained on older data, and constructing principled importance weights that transfer knowledge across waves.

#### 6.7.1 The Problem: Population Shift Across Survey Waves

Complex health surveys (MEPS, NHANES, BRFSS) are repeated cross-sectional or panel studies conducted in annual or biennial waves. A critical methodological challenge arises when models or estimates from one wave are applied to another:

- A predictive model trained on MEPS Panel 22 (2019–2020) may be deployed on Panel 24 (2021–2022), after the COVID-19 pandemic has fundamentally altered the healthcare utilization distribution.
- An imputation model for missing healthcare expenditures calibrated on NHANES 2011–2014 may be misspecified for NHANES 2017–2020, where the demographic composition of the U.S. has shifted.
- A survey weight calibration scheme derived on an older wave may not adequately represent the current population, yet the same calibration targets (census marginals) are often reused.

The core question in each case is: **have the population distributions $\mathbb{P}_t$ and $\mathbb{P}_{t+1}$ drifted significantly between survey waves $t$ and $t+1$?** If so, by how much, in which directions of the covariate space, and can the drift be corrected?

#### 6.7.2 The Sequential Survey-Weighted MMD Test

Let $\hat{\mathbb{P}}_{w,t}$ and $\hat{\mathbb{P}}_{w,t+1}$ be the Hajek empirical measures from waves $t$ and $t+1$, each with their own design weights $w_i^{(t)}$ and $w_j^{(t+1)}$. The squared inter-wave MMD is:

$$D_t = \widehat{\rm MMD}_w^2\!\left(\hat{\mathbb{P}}_{w,t},\; \hat{\mathbb{P}}_{w,t+1}\right) = \tilde{\mathbf{w}}_t^\top K_{tt}\,\tilde{\mathbf{w}}_t + \tilde{\mathbf{w}}_{t+1}^\top K_{t+1,t+1}\,\tilde{\mathbf{w}}_{t+1} - 2\,\tilde{\mathbf{w}}_t^\top K_{t,t+1}\,\tilde{\mathbf{w}}_{t+1}$$

Each wave provides its own Kim-Rao bootstrap weight matrix (from `generate_kim_rao_weights()`), and the null distribution of $D_t$ — what inter-wave MMD would look like if the population had not changed — is obtained by bootstrapping within each wave and computing the cross-wave MMD on the bootstrap replicates.

Over $T$ sequential survey waves, the sequence $\{D_1, D_2, \ldots, D_{T-1}\}$ forms a **time series of population distributional distances**. This enables:

**Changepoint detection.** A sudden spike in $D_t$ signals that the population distribution changed substantially between waves $t$ and $t+1$. Applied to MEPS 2018–2022, this would formally identify the wave-to-wave distributional impact of the COVID-19 pandemic, the ACA Medicaid expansion, or a major policy reform, as a statistically significant shift in the multivariate health and expenditure distribution — not merely a shift in the mean.

**CUSUM process control.** Define a cumulative sum statistic:
$$C_t = \max\!\left(0,\; C_{t-1} + D_t - \kappa\right), \quad C_0 = 0$$
where $\kappa$ is a reference value (half the expected MMD under no drift, estimated from bootstrap). $C_t$ exceeds a control limit $h$ (derived from the bootstrap null of $D_t$) when the cumulative distributional drift becomes too large. This is a **CUSUM chart for population health surveillance**: it triggers an alarm when the population has shifted enough to warrant retraining all downstream models.

#### 6.7.3 Survey-Wave Transfer via MMD-Minimizing Importance Weights

When population drift is detected, models trained on wave $t$ need to be transferred to wave $t+1$. The standard tool is **importance weighting**: reweight observations from wave $t$ so that the reweighted distribution matches wave $t+1$. But in the survey setting, both waves have their own complex designs, and naive density ratio estimation ignores the design weights.

The **survey-wave transfer problem** is: find a set of multiplier weights $\lambda_i \geq 0$ for wave-$t$ observations such that the doubly-reweighted measure $\sum_i \lambda_i \tilde{w}_i^{(t)} \delta_{X_i^{(t)}}$ is as close as possible (in MMD) to the target wave-$t+1$ distribution $\hat{\mathbb{P}}_{w,t+1}$.

This is exactly the MMD calibration QP from §6.2.2, applied across time rather than across sources:

$$\min_{\lambda \geq 0} \widehat{\rm MMD}^2\!\left(\sum_i \lambda_i \tilde{w}_i^{(t)} \delta_{X_i^{(t)}},\; \hat{\mathbb{P}}_{w,t+1}\right)$$

The solution $\hat{\lambda}^*$ gives **design-consistent transfer weights** that simultaneously account for the survey design in both waves and minimize the full distributional distance to the target wave. These transfer weights can then be used to:

1. **Retrain predictive models** fitted on wave $t$ so they are representative of the wave-$t+1$ population.
2. **Project estimates** of population parameters (prevalence, mean expenditure, treatment effect) from wave $t$ to wave $t+1$ without re-estimating the full model.
3. **Interpolate** between waves when a survey year is missing (e.g., NHANES was not collected in 2020 due to COVID-19; transfer weights from 2017–2018 and 2021–2022 can construct a 2020 pseudo-wave).

#### 6.7.4 Detecting the Impact of Policy Interventions as MMD Changepoints

A particularly powerful application for the dissertation is using the inter-wave MMD sequence as a **policy impact evaluation tool**. Consider the Affordable Care Act's Medicaid expansion (2014) as a natural experiment:

- Compute $D_t = \widehat{\rm MMD}_w^2(\hat{\mathbb{P}}_{w,t}, \hat{\mathbb{P}}_{w,t+1})$ for each MEPS wave from 2010 to 2018.
- Test whether $D_{2013}$ (pre- to post-expansion) is significantly larger than the pre-2014 baseline drift.
- Decompose the MMD into contributions from individual directions in $\mathcal{H}_k$ (via kernel PCA) to identify which dimensions of the health distribution changed most — insurance rates, utilization patterns, expenditure distributions — and among which subpopulations.

This provides a distributional causal effect decomposition that is not achievable by any mean-based impact evaluation, and that fully accounts for the MEPS complex survey design via the Kim-Rao bootstrap.

---

### 6.8 Federated Survey Statistics via Kernel Mean Embeddings

> **This application uses the kernel mean embedding as a privacy-preserving sufficient statistic**, enabling distributed computation of population-level distributional statistics across multiple independent survey agencies without ever sharing individual-level microdata.

#### 6.8.1 The Problem: Multi-Source Survey Integration Without Microdata Sharing

National health statistics are produced by a patchwork of surveys, each conducted by different agencies with different designs: MEPS (AHRQ), NHANES (CDC/NCHS), BRFSS (CDC/states), SIPP (Census Bureau), CPS (BLS/Census), Medicare claims (CMS). Each survey samples a different (overlapping) subset of the population with a different complex design and is subject to different privacy restrictions.

Integrating these surveys to obtain richer population estimates would be enormously valuable — but direct microdata linkage is legally restricted (Title 13, HIPAA, CIPSEA) and administratively prohibitive. Current practice is to release only summary tables or public-use microdata with disclosure protections that severely degrade small-cell estimates.

**The question:** can agencies collaborate to compute population-level distributional statistics and two-sample tests across surveys *without* ever sharing individual-level records?

#### 6.8.2 The Kernel Mean Embedding as a Privacy-Preserving Sufficient Statistic

The key mathematical observation is that the weighted kernel mean embedding $\hat{\mu}_w = \sum_i \tilde{w}_i k(\cdot, X_i) \in \mathcal{H}_k$ is a **sufficient statistic** for the weighted empirical measure $\hat{\mathbb{P}}_w$ with respect to the family of tests of the form $\langle \hat{\mu}_w - \hat{\mu}_{w'}, f\rangle$ for $f \in \mathcal{H}_k$. That is, two datasets have the same kernel mean embedding if and only if the survey-weighted MMD between them is zero.

Critically, the embedding $\hat{\mu}_w$ can be represented as a finite-dimensional vector in an approximate RKHS (via Nyström approximation or random Fourier features with $D$ basis functions):

$$\hat{\mu}_w \approx \sum_i \tilde{w}_i \phi(X_i) \in \mathbb{R}^D, \quad \phi(x) = \frac{1}{\sqrt{D}}\left[\cos(\omega_1^\top x + b_1), \ldots, \cos(\omega_D^\top x + b_D)\right]$$

where $(\omega_j, b_j)$ are random Fourier features drawn once and shared across agencies. Each agency communicates only this $D$-dimensional vector (typically $D \approx 500$–$2000$) to a coordinator. The coordinator can then compute approximate inter-agency MMDs and perform two-sample tests **without access to any individual-level data**.

#### 6.8.3 Federated Survey-Weighted MMD Protocol

**Setup:** $K$ survey agencies, each with data $\{(X_i^{(k)}, w_i^{(k)})\}_{i=1}^{n_k}$, survey design (strata/PSU), and a shared random Fourier feature map $\phi: \mathcal{X} \to \mathbb{R}^D$.

**Step 1 — Local computation (at each agency $k$):**
- Compute the local weighted empirical embedding: $\hat{z}^{(k)} = \sum_i \tilde{w}_i^{(k)} \phi(X_i^{(k)}) \in \mathbb{R}^D$.
- Compute $B$ Kim-Rao bootstrap replicates of the embedding: $\hat{z}^{(k,b)} = \sum_i \tilde{w}_i^{(k,b)} \phi(X_i^{(k)})$ for $b = 1, \ldots, B$.
- Transmit $\hat{z}^{(k)}$ and $\{\hat{z}^{(k,b)}\}_{b=1}^B$ to the coordinator. **No microdata is transmitted.**

**Step 2 — Coordinator computation:**
The coordinator reconstructs approximate MMDs from dot products of the transmitted vectors:

$$\widehat{\rm MMD}^2\!\left(\hat{\mathbb{P}}_{w}^{(j)}, \hat{\mathbb{P}}_{w}^{(k)}\right) \approx \|\hat{z}^{(j)}\|^2 + \|\hat{z}^{(k)}\|^2 - 2\,\langle \hat{z}^{(j)}, \hat{z}^{(k)}\rangle$$

The null distribution for the two-sample test between agencies $j$ and $k$ is obtained from $\|\hat{z}^{(j,b)} - \hat{z}^{(k,b)}\|^2$ across bootstrap replicates, yielding a valid p-value.

**Step 3 — National embedding:**
A weighted combination of local embeddings gives the national kernel mean embedding:
$$\hat{z}^{\rm national} = \sum_k \lambda_k \hat{z}^{(k)}, \qquad \lambda_k = \frac{N_k}{\sum_j N_j}$$
where $N_k$ is the estimated finite population size from agency $k$. This national embedding can be used to:
- Perform a two-sample MMD test between any two national subgroups (e.g., MEPS insured vs. BRFSS uninsured) at the population level.
- Compare the national distribution to a census reference for calibration diagnostics.
- Train a national synthetic data generator using the federated MMD loss.

#### 6.8.4 Privacy Guarantees

The transmitted vectors $\hat{z}^{(k)}$ are **aggregate statistics**, not individual records. Under the random Fourier feature representation:
- The sensitivity of $\hat{z}^{(k)}$ to any single record $X_i$ is $\|\tilde{w}_i^{(k)} \phi(X_i)\|_2 / \sqrt{D} \leq \tilde{w}_i^{(k)} / \sqrt{D}$.
- Adding Gaussian noise $\mathcal{N}(0, \sigma_{\rm priv}^2 I_D)$ with $\sigma_{\rm priv} = \max_i \tilde{w}_i^{(k)} / (\sqrt{D} \cdot \epsilon)$ guarantees $(\epsilon, \delta)$-differential privacy for the transmitted embedding vector.
- Privacy and statistical power trade off via $D$: larger $D$ gives better MMD approximation but also larger noise budget. The optimal $D$ given a privacy budget $\epsilon$ can be derived from concentration inequalities on the Nyström approximation error.

This provides a path to **differentially private, federated distributional testing across health surveys** — a combination that has no existing methodology.

#### 6.8.5 Application: Comparing Health Distributions Across State-Level Surveys

State-level health surveys (BRFSS, state NHANES supplements, Medicaid claims) each survey a state's population independently. A federated MMD protocol would enable:
- **Interstate health disparity testing:** is the distribution of healthcare utilization in North Carolina distributionally different from that in Texas, accounting for each state's complex survey design?
- **National benchmarking:** how far is each state's weighted health distribution from the national NHANES distribution? Which states are outliers?
- **Longitudinal tracking:** for each state, track the inter-wave MMD over time to detect which states responded most to national policy changes.

All of this is achievable with $D$-dimensional vector transmissions, not microdata, making it feasible under existing data-sharing agreements and privacy statutes.

### 6.9 Survey-Weighted Dirichlet Process Mixture Models: Convergence, Phenotyping, and Synthetic Data

> **This is the deepest methodological addition in the DBDL program.** It shows that a DPMM fitted with survey-calibrated bootstrap weights converges — in MMD — to the true finite population mixing measure, unifying the Bayesian nonparametric and design-based traditions. This provides the theoretical backbone for population-representative phenotyping and principled synthetic data generation.

#### 6.9.1 The Standard DPMM and Why it Fails for Survey Data

The **Dirichlet Process Mixture Model (DPMM)** [Ferguson, 1973; Antoniak, 1974] places a Dirichlet process prior on a mixing measure $G$ over a parameter space $\Theta$:

$$G \sim \text{DP}(\alpha, G_0), \qquad \theta_i \mid G \sim G, \qquad X_i \mid \theta_i \sim F(\cdot \mid \theta_i)$$

The DP prior has the stick-breaking representation [Sethuraman, 1994]:

$$G = \sum_{k=1}^\infty \pi_k \,\delta_{\theta_k^*}, \quad \pi_k = v_k \prod_{l<k}(1-v_l), \quad v_k \sim \text{Beta}(1, \alpha), \quad \theta_k^* \sim G_0$$

In the i.i.d. setting, a celebrated result establishes that the posterior $\Pi_n(G \in \cdot \mid X_1, \ldots, X_n)$ **concentrates around the true mixing measure** $G_0^{\rm true}$ as $n \to \infty$: the posterior predictive distribution $\hat{F}_n(x) = \int F(x \mid \theta) \hat{G}_n(d\theta)$ converges weakly to the true $F_0$.

**The survey failure:** In a complex survey, $X_1, \ldots, X_n$ are not i.i.d. — they are drawn with unequal probabilities $\pi_i$ from a finite population. A DPMM fitted naively to the unweighted sample treats each observation as equally informative, so the posterior concentrates around the *sample* distribution $\frac{1}{n}\sum_i \delta_{X_i}$, not the *population* distribution $\hat{\mathbb{P}}_w = \sum_i \tilde{w}_i \delta_{X_i}$. The resulting phenotypes are clusters of the *sample*, systematically over-representing oversampled subgroups (e.g., elderly MEPS respondents who are oversampled for precision) and under-representing undersampled subgroups.

#### 6.9.2 The Survey-Weighted DPMM

Define the **survey-weighted DPMM likelihood** by replacing the standard log-likelihood with a weighted version. For observations $X_1, \ldots, X_n$ with normalized design weights $\tilde{w}_i = w_i / \sum_j w_j$, the weighted log-likelihood is:

$$\ell_w(\theta_1, \ldots, \theta_n) = \sum_{i=1}^n \tilde{w}_i \log F(X_i \mid \theta_i)$$

The weighted DPMM posterior is proportional to:

$$\Pi_n^w(G, \theta_1, \ldots, \theta_n \mid X_{1:n}) \propto \exp\!\left(\sum_i \tilde{w}_i \log F(X_i \mid \theta_i)\right) \cdot \prod_i G(d\theta_i) \cdot \text{DP}(\alpha, G_0)(dG)$$

This is the likelihood power posterior at the design-weight scale. The Chinese Restaurant Process (CRP) representation still holds but with modified seating probabilities: unit $i$ joins cluster $k$ with probability proportional to $\tilde{w}_i \cdot n_k \cdot F(X_i \mid \hat{\theta}_k)$ rather than $n_k \cdot F(X_i \mid \hat{\theta}_k)$, where $n_k$ is the (weighted) count of units already in cluster $k$.

**Practical MCMC.** The Kim-Rao bootstrap provides the computational device. At each bootstrap iteration $b$, the DPMM is fitted with weights $w_i^{*(b)}$ from `W_mat[b, ]`. This produces a bootstrap DPMM posterior $\Pi_n^{w^{*(b)}}$, giving:
- A collection of $B$ posterior draws of the mixing measure $\hat{G}^{(b)}$,
- A collection of $B$ cluster assignments $c_i^{(b)}$,
- A collection of $B$ posterior predictive samples $\tilde{X}^{(b)} \sim \hat{F}^{(b)}$.

The variability across bootstrap iterations quantifies both the posterior uncertainty of the DPMM and the design-based uncertainty due to the survey sampling mechanism simultaneously — a **two-in-one uncertainty quantification**.

#### 6.9.3 Convergence Theorem: The Posterior Concentrates in MMD

**Theorem (Design-Based Posterior Concentration).** Let $k$ be a bounded characteristic kernel. Let $\hat{G}_n^w$ be the posterior mean of the survey-weighted DPMM mixing measure. Under Isaki-Fuller design conditions and standard DPMM identifiability conditions on $F(\cdot \mid \theta)$:

$$\widehat{\rm MMD}_w^2\!\left(\hat{F}_n^w,\; \hat{\mathbb{P}}_w\right) \xrightarrow{p} 0 \quad \text{(design-based)}$$

where $\hat{F}_n^w(x) = \int F(x \mid \theta)\,\hat{G}_n^w(d\theta)$ is the posterior predictive distribution of the weighted DPMM and $\hat{\mathbb{P}}_w$ is the Hajek empirical measure.

*Sketch of proof:* 
- By design consistency (\u00a73.2), $\hat{\mathbb{P}}_w \xrightarrow{p} \mathbb{P}_{\rm pop}$ in $\mathcal{H}_k$-norm.
- The weighted likelihood $\ell_w$ is a consistent estimator of $\mathbb{E}_{\mathbb{P}_{\rm pop}}[\log F(X \mid \theta)]$ (Hajek consistency of the weighted average).
- By Schwartz's theorem for misspecified models [Walker et al., 2004], the posterior of the weighted DPMM concentrates on a neighborhood of the Kullback-Leibler minimizer of $\mathbb{E}_{\mathbb{P}_{\rm pop}}[\log F_0(X) / F(X \mid \theta)]$, which is $F_0 = \mathbb{P}_{\rm pop}$ when the DPMM is well-specified.
- Applying the design-consistent MMD to measure the distance between $\hat{F}_n^w$ and $\hat{\mathbb{P}}_w$, both converge to $\mathbb{P}_{\rm pop}$ in $\mathcal{H}_k$-norm, so their MMD distance converges to zero.

**Corollary.** Under the same conditions, the Kim-Rao bootstrap distribution of $\widehat{\rm MMD}_w^2(\hat{F}_n^{w^*}, \hat{\mathbb{P}}_{w^*})$ is a design-consistent estimate of the sampling distribution of $\widehat{\rm MMD}_w^2(\hat{F}_n^w, \hat{\mathbb{P}}_w)$.

This is the key convergence result that grounds the DPMM within the DBDL framework: the MMD is both the *metric* in which convergence is measured and the *calibration criterion* for assessing goodness of fit at each bootstrap iteration.

#### 6.9.4 Application 1: Population-Representative Synthetic Data Generation

The survey-weighted DPMM posterior predictive distribution $\hat{F}_n^w$ is the **optimal nonparametric density estimator** for the finite population distribution under squared error loss in $\mathcal{H}_k$. Drawing synthetic observations from $\hat{F}_n^w$ gives:

$$\tilde{X} \sim \hat{F}_n^w = \int F(\cdot \mid \theta)\,\hat{G}_n^w(d\theta)$$

By the convergence theorem above, these synthetic observations are drawn from a distribution that is close (in MMD) to $\mathbb{P}_{\rm pop}$ — the true finite population distribution. This is strictly superior to:

- **Standard DPMM synthetic data:** converges to the sample distribution, not the population distribution.
- **Kernel density estimation:** cannot handle mixed-type covariates and does not account for survey design.
- **Copula models:** parametric dependence assumption; no convergence guarantee in MMD.
- **SW-MMD-GAN (\u00a76.1):** no analytic convergence guarantee for the generator; requires gradient-based training.

The **synthetic data generation algorithm** using the survey-weighted DPMM is:

1. Generate $B$ bootstrap weight vectors $\mathbf{w}^{*(b)}$ via `generate_kim_rao_weights()`.
2. For each $b$, run a Gibbs sampler for the DPMM with weights $\mathbf{w}^{*(b)}$, obtaining cluster assignments $c^{(b)}$ and cluster parameters $\hat{\theta}_k^{(b)}$.
3. For each $b$, draw $M$ synthetic observations: sample a cluster $k \sim \text{Categorical}(\hat{\pi}^{(b)})$, then $\tilde{X} \sim F(\cdot \mid \hat{\theta}_k^{(b)})$.
4. Evaluate fidelity: compute $\widehat{\rm MMD}_w^2(\hat{\mathbb{P}}_{\tilde{X}}, \hat{\mathbb{P}}_w)$ for each $b$ and report the bootstrap distribution as a quality certificate.

**Validation loop:** The survey-weighted MMD test from \u00a74.4 serves as the acceptance criterion — synthetic datasets with $p$-value $> 0.05$ are indistinguishable from the population in the full RKHS sense.

#### 6.9.5 Application 2: Population-Representative Phenotype Identification

**Phenotyping** in health data means identifying clinically meaningful subgroups (phenotypes) in a population. Standard clustering (k-means, GMM, or even unweighted DPMM) finds clusters of the *sample*, which may not correspond to distinct population subgroups. In MEPS, for example, an unweighted DPMM may produce a large "elderly high-cost" cluster because elderly respondents are oversampled for precision, even if they represent a small fraction of the U.S. population.

The survey-weighted DPMM identifies phenotypes in the *population*:
- **Cluster weights $\hat{\pi}_k$** are proportional to the *population* prevalence of the $k$-th phenotype, not the sample count. With Kim-Rao bootstrap uncertainty, each $\hat{\pi}_k^{(b)}$ gives a bootstrap replicate of the population prevalence.
- **Cluster parameters $\hat{\theta}_k$** characterize the covariate distribution of the $k$-th phenotype as it exists in the finite population.
- **Cluster count $K$** is the number of distinct population phenotypes, with its posterior distribution reflecting genuine uncertainty about how many subgroups exist in the population.

**Formal phenotype distinctness test.** After identifying $K$ phenotypes, the survey-weighted MMD provides a formal test of whether phenotype $j$ and phenotype $k$ are truly distributionally distinct in the population:

$$H_0^{jk}: \mathbb{P}_{\rm phenotype\,j} = \mathbb{P}_{\rm phenotype\,k}$$

using the weighted empirical measures of units assigned to each cluster, with Kim-Rao bootstrap null distribution. This replaces the ad-hoc practice of inspecting mean covariate differences between clusters.

**Example: MEPS healthcare phenotypes.** Fit the survey-weighted DPMM on MEPS respondents using covariates $(X_i) = (\text{age, chronic conditions, insurance, utilization, expenditure})$ with MEPS design weights. The resulting phenotypes represent population subgroups such as:

| Phenotype | Population Prevalence $\hat{\pi}_k$ (95% CI) | Clinical Interpretation |
|---|---|---|
| Healthy-young-uninsured | 18% (15–21%) | Young adults with no chronic conditions, low utilization |
| Chronic-managed-Medicare | 12% (10–14%) | Elderly with multiple managed chronic conditions |
| High-cost-complex | 4% (3–5%) | Complex multi-morbidity, high expenditure |
| ... | ... | ... |

Each prevalence estimate comes with a design-consistent confidence interval from the Kim-Rao bootstrap, and each pair of phenotypes is testable for distributional equivalence via the survey-weighted MMD.

**Stability of phenotypes across bootstrap replicates.** A cluster structure is stable if the same phenotypes emerge consistently across bootstrap replicates. Define the **bootstrap cluster stability score** for phenotype $k$ as:

$$S_k = \frac{1}{B(B-1)} \sum_{b \neq b'} \text{ARI}\!\left(c^{(b)}_{\rm restrict\,k},\; c^{(b')}_{\rm restrict\,k}\right)$$

where ARI is the Adjusted Rand Index between cluster assignments in replicates $b$ and $b'$, restricted to units assigned to phenotype $k$ in both. High $S_k$ means the phenotype is robustly identified across the design variability; low $S_k$ suggests the phenotype may be a sample artifact.

#### 6.9.6 Connection to the Full DBDL Framework

The survey-weighted DPMM is the **generative model** that completes the DBDL framework:

```
                 BART (Survey BART)
                 ┌───────────────────────────────────────┐
                 │ Nonparametric conditional distribution │
                 │ E[Y | X, design weights]              │
                 └─────────────────┬─────────────────────┘
                                   │
Survey-weighted DPMM               │         Survey-weighted MMD
┌────────────────────────────┐     │     ┌──────────────────────────────┐
│ Nonparametric marginal     │     │     │ Distributional distance and  │
│ distribution P(X, design)  │◄────┴────►│ testing: MMD(P_A, P_B)       │
│ Phenotypes + synthetic data│           │ Calibration + fairness audit │
└────────────────────────────┘           └──────────────────────────────┘
                 │                                   │
                 └───────────────┬───────────────────┘
                                 │
                    Kim-Rao Bootstrap (W_mat)
                    ┌─────────────────────────────┐
                    │ Design-consistent uncertainty│
                    │ for ALL three components     │
                    └─────────────────────────────┘
```

Survey BART models $\mathbb{P}(Y \mid X)$; the survey-weighted DPMM models $\mathbb{P}(X)$; together they specify the full joint population distribution $\mathbb{P}(X, Y) = \mathbb{P}(Y \mid X)\,\mathbb{P}(X)$, both estimated by design-consistent nonparametric methods and both calibrated by the same Kim-Rao bootstrap. The survey-weighted MMD is the loss function, diagnostic, and test statistic that ties the entire framework together.

---

### 6.10 Summary of Application Domains

| Application | Group A | Group B | Kernel | Kim-Rao Bootstrap Role |
|---|---|---|---|---|
| Covariate balance (§5.1) | Treated units (MEPS) | Control units | Gaussian on $\mathbb{R}^p$ | Null distribution for balance test |
| Distributional causal effects (§5.2) | $\hat{Y}(1)$ posterior draws | $\hat{Y}(0)$ posterior draws | Gaussian on $\mathbb{R}$ | Propagates design + MCMC uncertainty |
| Synthetic data fidelity (§6.1) | Synthetic sample | Weighted real survey | Gaussian on $\mathbb{R}^p$ | Null distribution; training loss gradient |
| Survey-to-census calibration (§6.2) | Calibrated survey | Census reference | Gaussian on $\mathbb{R}^p$ | Uncertainty after calibration |
| CGM trajectory comparison (§6.3) | Diabetic participants | Non-diabetic participants | $L^2$ or shape kernel on $L^2([0,T])$ | Null distribution; OOB diagnostics |
| Functional distributional regression (§6.3.4) | Covariate profile $x_A$ | Covariate profile $x_B$ | Product $k_X \otimes k_F$ | Bootstrap CI on conditional MMD |
| Clinical AI fairness audit (§6.5) | Demographic group $G_1$ (MEPS/NHANES) | Demographic group $G_2$ | Gaussian on $[0,1]$ (score space) | Population-level null; intersectional CI |
| Survey-wave transfer (§6.7) | Survey wave $t$ | Survey wave $t+1$ | Gaussian on $\mathbb{R}^p$ | CUSUM control chart; transfer weights |
| Federated survey MMD (§6.8) | Agency $j$ embedding | Agency $k$ embedding | RFF approx. in $\mathbb{R}^D$ | Federated bootstrap via local replicates |
| **DPMM synthetic data (§6.9)** | **DPMM posterior predictive** | **Hajek empirical measure** | **Gaussian on $\mathbb{R}^p$** | **Calibrated posterior; bootstrap fidelity** |
| **DPMM phenotyping (§6.9)** | **Phenotype $j$ (survey-weighted)** | **Phenotype $k$ (survey-weighted)** | **Gaussian on $\mathbb{R}^p$** | **Population prevalence CI; stability score** |

---

## 7. Kernel Selection and Bandwidth

### 7.1 Gaussian Kernel

The Gaussian (RBF) kernel is the canonical choice:

$$k_\sigma(x, x') = \exp\!\left(-\frac{\|x - x'\|^2}{2\sigma^2}\right)$$

It is characteristic on $\mathbb{R}^d$, bounded ($\sup k = 1$), and universal (dense in $C(\mathcal{X})$). For the bandwidth $\sigma$, use the **weighted median heuristic**:

$$\hat{\sigma}^2 = \frac{1}{2}\,\text{wMedian}\!\left(\|X_i - X_j\|^2 : i \neq j,\ i \in A \cup B\right)$$

where the median is weighted by $\tilde{w}_i \tilde{w}_j$ for the pair $(i,j)$. This is consistent and achieves near-optimal power in practice (Garreau et al., 2017).

### 7.2 Mixed-Type Covariates in MEPS and NHANES

MEPS and NHANES covariates include continuous (age, income, BMI, glucose AUC), binary (gender, insurance status), and ordinal (self-reported health, chronic condition count) variables. A product kernel handles this naturally:

$$k(x, x') = k_{\rm Gauss}(x_{\rm cont}, x'_{\rm cont}) \cdot \prod_{j \in \rm disc} k_{\rm AA}(x_j, x'_j)$$

where $k_{\rm AA}$ is the Aitchison-Aitken kernel for categorical variables. Alternatively, all variables can be one-hot encoded and a single Gaussian kernel applied to the expanded representation.

### 7.3 Functional Kernels for Trajectories

For functional observations (CGM curves, accelerometry, expenditure trajectories), replace the Euclidean kernel with the $L^2$ Gaussian kernel $k_F(f,g) = \exp(-\|f-g\|_{L^2}^2 / 2\sigma^2)$, or the elastic shape kernel for amplitude/phase separation (§6.3.5). The bandwidth $\sigma$ for functional kernels can be selected via the weighted median of pairwise $L^2$ distances between curves.

---

## 8. Advantages Over Existing Approaches

| Property | SMD (per-covariate) | Survey t-test | KS test | Energy Distance | **Survey MMD (proposed)** |
|---|:---:|:---:|:---:|:---:|:---:|
| Accounts for survey design | ❌ | ✅ | ❌ | ❌ | **✅** |
| Multivariate (joint test) | ❌ | ❌ | ❌ | ✅ | **✅** |
| Omnibus (any distribution) | ❌ | ❌ | ✅ (1D) | ✅ | **✅** |
| Valid null distribution | Partial | ✅ | ❌ | ❌ | **✅** |
| Mixed continuous/discrete | ✅ | ✅ | ❌ | Partial | **✅** |
| Functional/trajectory data | ❌ | ❌ | ❌ | Partial | **✅** |
| Generative training loss | ❌ | ❌ | ❌ | ❌ | **✅** |
| Distributional calibration | ❌ | ❌ | ❌ | ❌ | **✅** |
| OOB model diagnostics | ❌ | ❌ | ❌ | ❌ | **✅** |
| Single p-value for joint $H_0$ | ❌ | ❌ | ✅ | ✅ | **✅** |

The proposed survey-weighted MMD test is the only method that satisfies all ten properties simultaneously for complex survey data.

---

## 9. Implementation Roadmap

### Phase 1 — Core Two-Sample Test (1–2 weeks)

- [ ] Write `compute_weighted_mmd2(K_AA, K_BB, K_AB, wA, wB)`: pure matrix arithmetic, five lines of R.
- [ ] Write `survey_mmd_test()`: calls `generate_kim_rao_weights()` for each group, loops over bootstrap rows, returns p-value and bootstrap null distribution.
- [ ] Implement OOB extraction from `W_mat`: identify zero-weight rows per bootstrap replicate.
- [ ] Simulate under $H_0$ (same Gaussian population, same stratified cluster design): confirm empirical Type-I error $\approx 5\%$ at nominal $\alpha = 0.05$.
- [ ] Simulate under $H_1$ (different means / different variances): confirm power grows with $n$.

### Phase 2 — Kernel Refinement (1 week)

- [ ] Implement weighted median heuristic for bandwidth selection (scalar and functional cases).
- [ ] Implement unbiased (diagonal-zeroed) version of the weighted U-statistic to reduce finite-sample bias.
- [ ] Extract the BART leaf kernel from the C++ tree structures (`Node::route()`).
- [ ] Implement $L^2$ functional kernel for trajectory data.

### Phase 3 — Synthetic Data and Calibration Applications (2 weeks)

- [ ] Implement `survey_mmd_loss(G_samples, X_real, w_real)` as a training objective for synthetic data generators.
- [ ] Implement `mmd_calibrate(X_survey, w_base, X_census)`: solve the QP to find MMD-calibrated weights.
- [ ] Diagnostic: apply survey-to-census calibration test for MEPS vs. ACS marginals.
- [ ] Evaluate a copula-based synthetic MEPS generator using survey-weighted MMD fidelity score.

### Phase 4 — MEPS Causal Pipeline Integration (1 week)

- [ ] Run pre/post-IPTW covariate balance test in `estimate_causal_meps.R` for all covariates jointly.
- [ ] Test distributional equality of $\hat{Y}(1)$ vs $\hat{Y}(0)$ from the hurdle BART posterior.
- [ ] Compute OOB model diagnostic $T_{\rm GoF}$ across MCMC iterations; plot distribution.

### Phase 5 — Functional Data / NHANES CGM (2–3 weeks)

- [ ] Obtain NHANES CGM sub-study data; align curves to common time grid.
- [ ] Implement $L^2$ kernel matrix computation for CGM trajectory pairs.
- [ ] Run survey-weighted MMD test comparing diabetic vs. non-diabetic CGM trajectory distributions.
- [ ] Implement elastic shape kernel ($d_{\rm SRV}$) for phase-invariant comparison.
- [ ] Implement conditional kernel mean embedding for distributional regression on trajectories.

### Phase 6 — Theory (ongoing, dissertation chapter)

- [ ] Formalize design consistency of $\hat{\mu}_w$ in the functional setting (Proposition 1).
- [ ] State and prove bootstrap validity under Isaki-Fuller + bounded kernel conditions (Theorem 1).
- [ ] Prove consistency of the test under fixed $H_1$ (Theorem 2).
- [ ] Derive the asymptotic null distribution as a weighted sum of $\chi^2_1$ variates with design-dependent weights (Proposition 2).
- [ ] Establish convergence of the MMD calibration QP and consistency of MMD-calibrated estimators (Theorem 3).

---

## 10. Skeleton R Implementation

```r
# -------------------------------------------------------
# compute_weighted_mmd2
# -------------------------------------------------------
# Returns the squared survey-weighted MMD given kernel
# matrices and normalized weight vectors.
compute_weighted_mmd2 <- function(K_AA, K_BB, K_AB, wA, wB) {
  # wA, wB: normalized weight vectors summing to 1
  as.numeric(
    t(wA) %*% K_AA %*% wA +
    t(wB) %*% K_BB %*% wB -
    2 * t(wA) %*% K_AB %*% wB
  )
}

# -------------------------------------------------------
# survey_mmd_test
# -------------------------------------------------------
# Omnibus two-sample test for complex survey data.
#
# Arguments:
#   X_A      : (m x p) covariate matrix, group A
#   X_B      : (n x p) covariate matrix, group B
#   design_A : data.frame with cols weight, strata, psu
#   design_B : data.frame with cols weight, strata, psu
#   B        : number of Kim-Rao bootstrap replicates
#   sigma    : kernel bandwidth (NULL = median heuristic)
#
# Returns: list(mmd2_obs, p_value, bandwidth, boot_null)
survey_mmd_test <- function(X_A, X_B, design_A, design_B,
                            B = 1000, sigma = NULL) {

  # --- Normalized original weights ---
  wA0 <- design_A$weight / sum(design_A$weight)
  wB0 <- design_B$weight / sum(design_B$weight)

  # --- Bandwidth via weighted median heuristic ---
  if (is.null(sigma)) {
    X_pool <- rbind(X_A, X_B)
    D2     <- as.matrix(dist(X_pool))^2
    sigma  <- sqrt(median(D2[lower.tri(D2)]) / 2)
  }

  # --- Kernel matrices ---
  rbf <- function(A, B, s) {
    d2 <- outer(seq_len(nrow(A)), seq_len(nrow(B)),
                Vectorize(function(i,j) sum((A[i,]-B[j,])^2)))
    exp(-d2 / (2*s^2))
  }
  K_AA <- rbf(X_A, X_A, sigma)
  K_BB <- rbf(X_B, X_B, sigma)
  K_AB <- rbf(X_A, X_B, sigma)

  # --- Observed test statistic ---
  T0 <- compute_weighted_mmd2(K_AA, K_BB, K_AB, wA0, wB0)

  # --- Kim-Rao bootstrap weights (reuses existing function) ---
  W_A <- generate_kim_rao_weights(design_A, "strata", "psu", "weight", B)
  W_B <- generate_kim_rao_weights(design_B, "strata", "psu", "weight", B)

  # --- Bootstrap statistics ---
  T_boot <- numeric(B)
  for (b in seq_len(B)) {
    wA_b    <- W_A[b,] / sum(W_A[b,])
    wB_b    <- W_B[b,] / sum(W_B[b,])
    T_boot[b] <- compute_weighted_mmd2(K_AA, K_BB, K_AB, wA_b, wB_b)
  }

  list(
    mmd2_obs  = T0,
    p_value   = mean(T_boot >= T0),
    bandwidth = sigma,
    boot_null = T_boot
  )
}
```

---

## 11. Open Questions

1. **Bias in the diagonal terms.** The formula $\tilde{\mathbf{w}}_A^\top K_{AA}\, \tilde{\mathbf{w}}_A$ includes diagonal terms $\tilde{w}_i^2 k(X_i, X_i) = \tilde{w}_i^2$ (since $k(x,x) = 1$ for the Gaussian kernel). This introduces a negative bias relative to the true $\mathbb{E}[k(X,X')]$ for $i \neq j$. The unbiased version zeros the diagonal: $\sum_{i \neq j} \tilde{w}_i \tilde{w}_j k(X_i^A, X_j^A) / (1 - \sum_i \tilde{w}_i^2)$. The impact is small for large surveys but important for the null distribution approximation.

2. **Two-level uncertainty for distributional causal estimands.** When testing $H_0: \mathbb{P}(Y(1)) = \mathbb{P}(Y(0))$ using BART posterior draws, there are two sources of uncertainty: (a) the MCMC posterior over trees and parameters, and (b) the survey design variability. A two-level bootstrap (outer: Kim-Rao over the design; inner: posterior draws from MCMC) could propagate both, at the cost of $B_{\rm outer} \times B_{\rm inner}$ evaluations.

3. **Optimal kernel selection.** The median heuristic is heuristic. An alternative is to maximize the test power with respect to $\sigma$ on a held-out validation set, or to use a mixture of Gaussian kernels (as suggested by Gretton et al., 2012 §6). For functional data, cross-validated bandwidth selection in $L^2$ is computationally expensive and an open problem.

4. **Extension to conditional tests.** The marginal MMD tests $H_0: \mathbb{P}(X | A=1) = \mathbb{P}(X | A=0)$. A more refined question is whether the *conditional* distributions given confounders are balanced — a conditional MMD test. This requires kernel conditional mean embeddings, which are more complex but have an established theory (Song et al., 2013).

5. **Asymptotic relative efficiency vs. univariate tests.** When is the MMD more powerful than a collection of Bonferroni-corrected per-covariate survey-weighted tests? Theoretically, MMD should dominate when distributional differences are diffuse across many covariates. Empirically this could be verified in a simulation study tailored to the MEPS covariate structure.

6. **Scalability of the QP calibration.** The MMD calibration QP (§6.2.2) has $n$ variables and scales as $O(n^2)$ in memory for the kernel matrix. For large surveys ($n > 50{,}000$), approximate kernel methods (Nyström approximation, random Fourier features) will be needed to make the QP tractable while preserving design-consistency properties.

7. **Phase variation in functional data.** Applying the $L^2$ kernel to unregistered CGM curves conflates amplitude and phase variation. Optimal time-warping (Srivastava & Klassen, 2016) before computing the kernel matrix is one approach; using the SRV elastic kernel is another. It is an open question which produces the more powerful test for detecting group differences in glucose dynamics.

8. **Privacy guarantees for synthetic data MMD.** In the synthetic data generation application (§6.1), if the survey data is sensitive, the MMD loss function evaluated on the real data should be computed with differential privacy guarantees. The sensitivity of the weighted cross-term $\sum_j \tilde{w}_j k(\tilde{X}, X_j)$ to any single $X_j$ is bounded by $\tilde{w}_j \sup k = \tilde{w}_j$, providing a natural mechanism for adding calibrated Gaussian noise.

---

## 12. References

- Gretton, A., Borgwardt, K., Rasch, M., Schölkopf, B., & Smola, A. (2012). A kernel two-sample test. *Journal of Machine Learning Research*, **13**, 723–773.
- Kim, J.-K., & Rao, J. N. K. (2009). A unified approach to linearization variance estimation from survey data after imputation for item nonresponse. *Biometrika*, **96**(4), 917–932.
- Muandet, K., Fukumizu, K., Sriperumbudur, B., & Schölkopf, B. (2017). Kernel mean embedding of distributions: A review and beyond. *Foundations and Trends in Machine Learning*, **10**(1–2), 1–141.
- Smola, A., Gretton, A., Song, L., & Schölkopf, B. (2007). A Hilbert space embedding for distributions. *Algorithmic Learning Theory*, LNCS 4754, 13–31.
- Isaki, C. T., & Fuller, W. A. (1982). Survey design under the regression superpopulation model. *Journal of the American Statistical Association*, **77**(377), 89–96.
- Garreau, D., Jitkrittum, W., & Kanagawa, M. (2017). Large sample analysis of the median heuristic. *arXiv:1707.07269*.
- Song, L., Fukumizu, K., & Gretton, A. (2013). Kernel embeddings of conditional distributions. *IEEE Signal Processing Magazine*, **30**(4), 98–111.
- Sriperumbudur, B., Gretton, A., Fukumizu, K., Lanckriet, G., & Schölkopf, B. (2010). Hilbert space embedding and characteristic kernels above generic domains. *Journal of Machine Learning Research*, **11**, 1651–1680.
- Srivastava, A., & Klassen, E. P. (2016). *Functional and Shape Data Analysis*. Springer. [Elastic shape metric / SRV framework for functional data.]
- Deville, J.-C., & Särndal, C.-E. (1992). Calibration estimators in survey sampling. *Journal of the American Statistical Association*, **87**(418), 376–382. [Standard calibration / raking.]
- Bondell, H. D., Reich, B. J., & Wang, H. (2010). Noncrossing quantile regression curve estimation. *Biometrika*, **97**(4), 825–838.
- Chen, C., Twyman, N., et al. (2019). Differentially private synthetic data. *arXiv:1906.02994*. [Privacy-preserving synthetic data.]
- NHANES CGM Sub-study (2011–2014). National Health and Nutrition Examination Survey continuous glucose monitoring protocol. Centers for Disease Control and Prevention.
- FDA (2021). Artificial Intelligence/Machine Learning (AI/ML)-Based Software as a Medical Device (SaMD) Action Plan. U.S. Food and Drug Administration.
- HHS Office for Civil Rights (2022). Guidance on Nondiscrimination in Telehealth and the Use of Algorithmic Decision-Making. U.S. Department of Health and Human Services.
- Obermeyer, Z., Powers, B., Vogeli, C., & Mullainathan, S. (2019). Dissecting racial bias in an algorithm used to manage the health of populations. *Science*, **366**(6464), 447–453. [Foundational empirical work on clinical AI racial bias.]
- Barocas, S., Hardt, M., & Narayanan, A. (2023). *Fairness and Machine Learning: Limitations and Opportunities*. MIT Press. [Comprehensive reference on algorithmic fairness metrics.]
- Ben-David, S., Blitzer, J., Crammer, K., Kulesza, A., Pereira, F., & Vaughan, J. W. (2010). A theory of learning from different domains. *Machine Learning*, **79**(1–2), 151–175. [Domain adaptation theory; MMD as divergence measure.]
- Hardt, M., Price, E., & Srebro, N. (2016). Equality of opportunity in supervised learning. *Advances in Neural Information Processing Systems*, **29**. [Equalized odds fairness criterion.]

