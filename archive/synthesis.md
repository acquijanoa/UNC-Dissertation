# Literature Synthesis: Frequentist and Bayesian Paradigms for Complex Survey Data and Infinite-Dimensional Models

## 1. Introduction and Overarching Themes

A foundational assumption in classical statistical theory is that observations are independent and identically distributed (i.i.d.) draws from some underlying superpopulation. In practice, however, data from large-scale social, health, and economic surveys are collected using **complex survey designs**. These designs employ features such as stratification, clustering, multi-stage sampling, and unequal selection probabilities to maximize efficiency under budget constraints. Applying standard i.i.d.-based statistical methods to such data leads to:
1. **Bias in Parameter Estimation:** If the sampling design is *informative* (i.e., inclusion probabilities are correlated with the outcome variable of interest even conditional on covariates), ignoring design weights leads to severe estimation bias.
2. **Underestimated Variance (Overconfidence):** Clustering induces positive correlation among units within the same primary sampling unit (PSU). Ignoring this correlation understates the standard errors, leading to artificially narrow confidence or credible intervals and inflated Type I error rates in hypothesis testing.

This literature synthesis examines **ten core references** bridging the gap between complex survey sampling, semiparametric/nonparametric estimation, and asymptotic theory. These works fall into two broad statistical traditions:
* **The Frequentist Paradigm:** Relying on design-based randomization, Horvitz–Thompson (HT) estimators, inverse probability weighting (IPW), sandwich variance estimators, empirical processes, and rescaled bootstrap calibration.
* **The Bayesian Paradigm:** Employing pseudo-posterior likelihoods, Bayesian finite population predictive inference (treating non-sampled units as missing data), synthetic population generation via the Finite Population Bayesian Bootstrap (FPBB), and Weighted Dirichlet Process Mixtures (WDPM).

A core narrative emerging from this synthesis is the convergence of the two paradigms: Frequentist methods are adopting computationally intensive resampling (bootstraps) to avoid intractable variance derivations, while Bayesian methods are increasingly relying on design-based corrections (pseudo-likelihoods and covariance projections) to guarantee frequentist coverage.

---

## 2. Taxonomy and Intellectual Lineage

The papers can be organized not just by topic, but by their intellectual lineage. The mathematical foundations establish the limits of infinite-dimensional models, which are then applied to the Frequentist and Bayesian branches.

| Citation | Core Methodology | Key Contribution & Lineage |
| :--- | :--- | :--- |
| **Ghosal, Ghosh & van der Vaart (2000)** | Nonparametric Bayes | Establishes general rates of posterior convergence in infinite-dimensional models based on entropy and prior mass. *(Foundational theory)* |
| **Ghosal & van der Vaart (2007)** | Nonparametric Bayes (non-i.i.d.) | Extends the above theory to non-i.i.d. observations. *(Crucial for adapting nonparametric Bayes to clustered survey data)* |
| **Han & Wellner (2021)** | Empirical Processes | Proves uniform limit theorems for the HT empirical process. *(Foundational for proving properties of complex survey estimators)* |
| **Breslow & Wellner (2007)** | Weighted Likelihood / Semiparametric | Derives asymptotics for Cox regression in two-phase sampling. *(Applies empirical process theory to specific semiparametric models)* |
| **Kim, Rao & Wang (2024)** | Bootstrap Calibration | Unified bootstrap calibration for hypothesis tests. *(Provides a computational solution to the intractable variances found in models like Breslow & Wellner)* |
| **Savitsky & Toth (2016)** | Pseudo-Posterior | Proves consistency of pseudo-posteriors under informative sampling. *(The Bayesian analogue to Breslow & Wellner's weighted likelihood)* |
| **Williams & Savitsky (2021)** | MCMC post-processing | "Scale and shape projection" of MCMC draws. *(Corrects the covariance mismatch in Savitsky & Toth to match Frequentist properties)* |
| **Elliott & Xia (2021)** | Weighted DPM / Predictive Inference | Integrates Weighted DPMs for regression. *(Applies Ghosal-style infinite-dimensional models to the pseudo-posterior framework)* |
| **Dong, Elliott & Raghunathan (2014)** | Finite Population Bayesian Bootstrap | FPBB to generate synthetic SRS from complex samples. *(Alternative to pseudo-posteriors: invert the design first)* |
| **Zhou, Elliott & Raghunathan (2016)** | Double-Weighted FPBB | Synthetic MI for multi-stage designs. *(Scales Dong et al. to hierarchical clustered data)* |

---

## 3. Pillar I: Mathematical Foundations & Asymptotics

To build rigorous models for complex data, one must establish asymptotic consistency and rate of convergence results. The fundamental challenge is that complex surveys break the i.i.d. assumption that underlies classical asymptotics.

### Posterior Convergence in Infinite-Dimensional Models
**Ghosal, Ghosh, and van der Vaart (2000)** establish the modern framework for the rate of convergence of posterior distributions. In infinite-dimensional settings, the posterior distribution $\Pi_n(B \mid X_1, \dots, X_n)$ contracts around the true parameter $\theta_0$ at a rate $\epsilon_n$ (with $\epsilon_n \to 0$) if the posterior mass outside an $\epsilon_n$-neighborhood of $\theta_0$ vanishes in probability.
The two key conditions for a rate $\epsilon_n$ are:
1. **Entropy (Covering Numbers):** The size of the parameter space must not be too large, preventing the model from overfitting. Measured by $\log N(\epsilon_n, \Theta_n, d) \le n \epsilon_n^2$.
2. **Prior Mass:** The prior must assign sufficient mass to a Kullback-Leibler (KL) neighborhood around the true density.

**Ghosal and van der Vaart (2007)** extend this theory to **non-i.i.d. observations**, which is the critical bridge for survey data where clustering and stratification create dependence. They show that if testing conditions and prior mass conditions are satisfied uniformly over the non-i.i.d. structure, the posterior rate of convergence is preserved. This theory is what eventually allows models like Elliott and Xia's (2021) WDPM to possess theoretical guarantees when applied to clustered samples.

### Horvitz–Thompson Empirical Processes
In the frequentist design-based setting, M- and Z-estimators are analyzed using empirical process theory. **Han and Wellner (2021)** provide a comprehensive treatment of the **Horvitz–Thompson (HT) empirical process** under complex sampling designs:
$$
\mathbb{G}_N^{\pi}(f) = \frac{1}{\sqrt{N}} \sum_{i=1}^N \left( \frac{\xi_i}{\pi_i} - 1 \right) f(X_i)
$$
By establishing that $\mathbb{G}_N^{\pi}$ converges weakly to a Gaussian process across diverse designs (stratified, rejective, Poisson), Han and Wellner provide the unified theoretical foundation needed to prove the consistency and asymptotic normality of downstream complex survey estimators, freeing analysts from proving basic consistency anew for every proposed model.

---

## 4. Pillar II: Frequentist Semiparametrics & Bootstrap Calibration

Frequentist methods adjust for complex survey designs by incorporating sampling weights directly into the estimating equations (Inverse Probability Weighting) and calibrating the variances to match design-based properties (Sandwich Estimators).

### Semiparametric Weighted Likelihood
**Breslow and Wellner (2007)** apply empirical process principles to semiparametric models under two-phase stratified sampling (e.g., Cox regression). They solve the inverse probability weighted (IPW) likelihood equations:
$$
\Psi_{N}^{\pi}(\theta, \eta) = \mathbb{P}_N^{\pi} \dot{\ell}_{\theta, \eta} = 0
$$
They mathematically prove that the parameter estimator $\hat{\theta}$ is $\sqrt{N}$-consistent and asymptotically normal. A profound insight in this work is **Pierce’s Correction**: estimating the sampling weights $\pi_i$ via a parametric model based on phase 1 data actually *improves* the efficiency of the estimator compared to using the true, known design weights. The reduction in variance occurs because estimating the weights projects the influence function onto the score space of the auxiliary variables.

### Unified Bootstrap Hypothesis Testing
While Breslow and Wellner show how to derive the asymptotic variance, the derivations are mathematically complex and specific to the model. For general hypothesis testing, the classical likelihood ratio test statistic $W(\theta_0)$ under complex designs converges to a weighted sum of Chi-squares: $W(\theta_0) \xrightarrow{d} \sum_{i=1}^p c_i Z_i^2$, where the weights $c_i$ depend on intricate design effects.

**Kim, Rao, and Wang (2024)** resolve this computational bottleneck by proposing a **unified bootstrap weight calibration approach**. Instead of analytically deriving the covariance matrices, they generate rescaled bootstrap weights $w_i^*$ satisfying:
$$
E_*\{S_w^*(\theta)\} = S_w(\theta), \quad \text{and} \quad V_*\{S_w^*(\theta)\} = \hat{V}\{S_w(\theta) \mid \mathcal{F}_N\}
$$
By simply bootstrapping the sample and running the unadjusted likelihood ratio test on the bootstrapped weights, the empirical distribution of the test statistic provides asymptotically valid p-values. This democratizes frequentist complex survey inference, removing the need for specialized design-effect software.

---

## 5. Pillar III: Bayesian & Pseudo-Bayesian Paradigms

Classical Bayesian inference struggles with complex surveys because the survey weights do not map cleanly to a generative probability model. Modern Bayesian approaches either distort the likelihood with weights or build highly flexible imputation models.

### The Pseudo-Posterior Approach
**Savitsky and Toth (2016)** introduce the **pseudo-posterior** by raising individual likelihood contributions to the power of their normalized survey weights $w_i$:
$$
\pi(\theta \mid Y_{obs}, w) \propto \left[ \prod_{i \in A} f(y_i \mid \theta)^{w_i} \right] \pi(\theta)
$$
They establish that this pseudo-posterior asymptotically contracts around the superpopulation parameter. This approach allows Bayesians to utilize complex survey data while retaining MCMC computational engines.

### MCMC Uncertainty Calibration (The Covariance Mismatch)
A fundamental issue with the pseudo-posterior is that its covariance matrix does not match the frequentist asymptotic variance. The pseudo-posterior variance captures model uncertainty ($H^{-1}$) but ignores the design-induced variance ($J$) captured by the frequentist sandwich estimator $\Sigma = H^{-1} J H^{-1}$. In clustered samples ($J > H$), pseudo-posteriors yield artificially narrow credible intervals.

**Williams and Savitsky (2021)** resolve this by proposing a **scale and shape projection** of the MCMC draws. By centering the raw pseudo-posterior draws around their mean and stretching them using $R = \hat{\Sigma}^{1/2} \hat{H}^{1/2}$, they force the MCMC empirical covariance to match the design-based sandwich covariance, restoring nominal frequentist coverage without altering the original sampling algorithm.

### Weighted Dirichlet Process Mixture (WDPM) Models
An alternative Bayesian approach avoids altering the likelihood and instead uses **Bayesian finite population inference**, treating the unsampled population as missing data to be imputed.

**Elliott and Xia (2021)** utilize a **Weighted Dirichlet Process Mixture (WDPM)** model for this task. The WDPM assigns weights to locations in the covariate space, allowing mixture components to adapt locally to informative sampling probabilities. The model induces a data-driven bias-variance tradeoff:
* If the true model is linear and weights are non-informative, the model shrinks to a single mixture component, behaving like an efficient unweighted MLE.
* If there is non-linearity or informative sampling, it expands to multiple components to correct for bias, behaving like a fully-weighted estimator.
This is grounded in the infinite-dimensional asymptotics established by Ghosal et al., ensuring that the flexible prior space contracts properly around the true data generating mechanism.

---

## 6. Pillar IV: Synthetic Populations & Multiple Imputation

When item-level missing data are present in complex surveys, analysts face a dual challenge: modeling the missingness mechanism while simultaneously accounting for the complex design features.

### Nonparametric Synthetic Populations
**Dong, Elliott, and Raghunathan (2014)** propose a nonparametric method using the **Finite Population Bayesian Bootstrap (FPBB)**. Given a sample of size $n$, the FPBB generates a complete synthetic population of size $N$ by drawing units sequentially with probabilities adjusted by their survey weights. 
By generating synthetic populations, the complex survey design features are "inverted" into the synthetic data. The synthetic population acts as a **simple random sample (SRS)** from the superpopulation. Secondary analysts can thus apply standard i.i.d. statistical methods (e.g., standard regression) to the synthetic datasets without knowing the complex design or weights.

### Multiple Imputation in Multi-Stage Designs
**Zhou, Elliott, and Raghunathan (2016)** extend the FPBB to address item-level missing data in stratified multi-stage designs. Under standard Multiple Imputation (MI), the analyst must include dummy variables representing all strata and PSUs to avoid misspecification, leading to over-parameterized models.
They develop a **double-weighted FPBB procedure (SYN1)**:
1. Re-sample PSUs using FPBB adjusted by PSU-level weights.
2. Re-sample elements within selected PSUs using conditional element-level weights.
3. Impute missing values on the generated synthetic population using an SRS assumption.
This method elegantly preserves the complex design structure while drastically reducing the operational and computational burden for downstream analysts.

---

## 7. Comparative Analysis: Frequentist vs. Bayesian Approaches

| Dimension | Frequentist Paradigm (IPW / Bootstrap) | Bayesian Paradigm (Pseudo-Likelihood / Predictive) |
| :--- | :--- | :--- |
| **Philosophical Focus** | Randomization-based; values in population are fixed, sampling is random. | Model-based; population values are random variables from a superpopulation. |
| **Handling of Weights** | Weights ($w_i = \pi_i^{-1}$) adjust the estimating equations (IPW) directly. | Weights exponentiate the likelihood (pseudo-posterior) or adjust resampling probabilities (FPBB). |
| **Model Misspecification** | Robust; consistently estimates the census parameter even if the model is wrong. | Mitigated via highly flexible infinite-dimensional nonparametric priors (e.g. DPM). |
| **Variance Estimation** | Empirical sandwich estimators ($\Sigma = H^{-1} J H^{-1}$) or rescaled bootstraps. | Adjusted via scale/shape projection of MCMC draws or predictive population imputation. |
| **Hypothesis Testing** | Bootstrap weight calibration to match design score moments (Kim et al., 2024). | Posterior credible intervals and Bayes factors derived from calibrated posteriors. |

---

## 8. Core Dissertation Proposal: Privacy-Preserving Design Inversion via Dirichlet Process Bootstraps

The culmination of this literature synthesis points directly to a novel, high-impact dissertation proposal. By synthesizing frequentist empirical processes, Bayesian nonparametrics, and bootstrap calibration, we can solve a critical modern problem: generating synthetic populations that correct for informative sampling without violating data privacy.

### 8.1 The Methodological Gap
Under complex survey designs, secondary analysts struggle to fit models because they must correctly specify the sampling weights and the design structure (strata, clusters) to avoid bias and underestimated variances. 
**Design Inversion** via the Finite Population Bayesian Bootstrap (FPBB; Dong et al., 2014; Zhou et al., 2016) brilliantly solves this: it generates a synthetic population of size $N$ by resampling the observed units with probabilities proportional to their survey weights. The resulting synthetic population acts as a Simple Random Sample (SRS), allowing analysts to use standard i.i.d. software without knowing the complex design.

**The Privacy Flaw:** The FPBB relies on the Polya posterior, which is fundamentally an empirical step-function distribution. It generates synthetic populations by making exact copies of the observed data $y_i$. In an era of strict privacy regulations, releasing datasets containing exact replicas of human subjects is often prohibited due to re-identification risks.

### 8.2 A Critical Regularity Constraint: Kim's Conditions and the DPMM Incompatibility

A fundamental theoretical tension must be confronted at the outset. Kim et al. (2024) derive the validity of their calibrated bootstrap under a set of regularity conditions on the quasi-log-likelihood $l_w(\theta)$, which include:
* **Concavity:** $l_w(\theta)$ must be globally concave in $\theta$, ensuring a unique quasi-MLE $\hat{\theta}$.
* **Unimodality of the Score:** The weighted score function $S_w(\theta) = \nabla_\theta l_w(\theta)$ must have a unique root, allowing the bootstrap score moments to be well-defined and matched at a single point.

These conditions are satisfied for standard models in the exponential family (logistic regression, Cox regression, Poisson GLMs), where the log-likelihood is a sum of individual concave contributions. **They are violated under a DPMM.** The DPMM log-likelihood takes the form:
$$
l_{DPMM}(\boldsymbol{\theta}, G \mid \mathbf{y}) = \sum_{i \in A} w_i \log \int K(y_i; \phi) \, dG(\phi)
$$
where $G \sim DP(\alpha, G_0)$ is the Dirichlet Process measure and $K(\cdot; \phi)$ is a continuous kernel (e.g., Gaussian). Because $\log \int K(y_i; \phi) \, dG(\phi)$ is a log-sum-exp of kernel evaluations, it is **not globally concave** and induces a likelihood surface with multiple local modes — one for each plausible clustering configuration of the data. Directly applying Kim's bootstrap calibration to the DPMM pseudo-likelihood is therefore theoretically invalid, since the bootstrap score conditions cannot be matched at a unique parameter value.

### 8.3 The Proposed Resolution: Two-Stage Decoupling

To preserve both the formal validity of Kim's bootstrap *and* the flexibility of the DPMM, the dissertation proposes a **two-stage decoupling** architecture that cleanly separates the role of bootstrap calibration from the role of nonparametric density estimation.

**Stage 1: Bootstrap Weight Calibration on a Working Parametric Model**

In the first stage, Kim's calibrated bootstrap is applied to a **working parametric model** $l_w^{par}(\theta)$ — for instance, a weighted linear model or a GLM — chosen to be globally concave in $\theta$. Kim's regularity conditions are fully satisfied here. The outputs of this stage are calibrated bootstrap weight vectors $\{w_i^{*(b)}\}_{b=1}^B$ for each bootstrap replicate $b$, satisfying:
$$
E_*\{S_w^{par,*}(\theta)\} = S_w^{par}(\hat{\theta}), \quad V_*\{S_w^{par,*}(\hat{\theta})\} = \hat{V}\{S_w^{par}(\hat{\theta})\}
$$
Critically, these calibrated weights encode the **full design variance structure** — the clustering, stratification, and unequal probabilities — without requiring the working model to be correctly specified for the outcome distribution. The weights are design-consistent objects, not model-specific ones.

**Stage 2: Weighted DPMM Synthetic Population Generation**

In the second stage, the calibrated weight vector $w_i^{*(b)}$ from each bootstrap replicate $b$ is treated as a fixed, known input — a design object — and is fed into the pseudo-likelihood of the DPMM:
$$
l_{DPMM}^{*(b)}(G \mid \mathbf{y}) = \sum_{i \in A} w_i^{*(b)} \log \int K(y_i; \phi) \, dG(\phi)
$$
For each bootstrap replicate, the DPMM is fitted via MCMC, and the unobserved population units $Y_{nob}^{*(b)}$ are drawn from the posterior predictive distribution:
$$
p^*(y_{new}) = \int K(y_{new}; \phi) \, d\hat{G}^{*(b)}(\phi)
$$
Because this distribution is an integral against a continuous kernel $K$, the synthetic values are drawn from a **strictly continuous distribution** — they are not exact copies of any observed $y_i$, satisfying privacy requirements. Because the weights $w_i^{*(b)}$ are calibrated to the complex design, the fitted DPMM correctly accounts for informative selection.

**The key decoupling insight:** Kim's concavity conditions apply only to the Stage 1 working model, which is a concave parametric object. The DPMM in Stage 2 receives the calibrated weights as external, fixed inputs and is never itself subject to Kim's score-matching conditions. The validity of the two stages is thus maintained separately and cleanly.

### 8.4 The Asymptotic Proof Strategy: Weak Convergence of Empirical Processes

The core theoretical contribution of the dissertation will be to prove formally that, under this two-stage decoupling, the empirical process of the DPMM-generated synthetic populations converges weakly to the HT empirical process of the true finite population.

Let $\hat{F}_N^{*(b)}$ be the empirical distribution of the synthetic population $\{Y_{nob}^{*(b)}, \mathbf{y}\}$ generated under bootstrap replicate $b$. We define the DPMM empirical process:
$$
\mathbb{G}_N^{DPMM,*}(f) = \frac{1}{\sqrt{N}} \sum_{j=1}^N f\!\left(y_j^{syn}\right), \quad y_j^{syn} \sim p^*(\cdot)
$$
and compare it against the HT empirical process:
$$
\mathbb{G}_N^{\pi}(f) = \frac{1}{\sqrt{N}} \sum_{i=1}^N \left( \frac{\xi_i}{\pi_i} - 1 \right) f(X_i)
$$
The proof proceeds in three steps:

1. **Stage 1 — Design Variance Transfer:** By construction of the Kim (2024) bootstrap weights, for any fixed $\theta$:
$$
\sqrt{N}\left(\frac{1}{B}\sum_{b=1}^B S_w^{par,*}(\theta)^{(b)} - S_w^{par}(\theta)\right) \xrightarrow{d} \mathbb{G}^{\pi}
$$
This shows that the distribution of the calibrated weight vectors across bootstrap replicates converges weakly to the HT Gaussian process.

2. **Stage 2 — Posterior Contraction of the DPMM:** Using **Ghosal, Ghosh, and van der Vaart (2000, 2007)**, for any fixed calibrated weight vector $w^{*(b)}$, the DPMM posterior $\Pi_n(\cdot \mid \mathbf{y}, w^{*(b)})$ contracts around the true weighted density $f_0^{\pi}$ at rate $\epsilon_n$, satisfying:
$$
\Pi_n\!\left(\left\|\hat{f}^{DPMM} - f_0^{\pi}\right\|_{L_2} > M_n \epsilon_n \,\Big|\, \mathbf{y}, w^{*(b)}\right) \to 0
$$
The non-i.i.d. extension (Ghosal & van der Vaart, 2007) is essential here since the observations are not identically distributed under informative sampling.

3. **Functional Equivalence & Weak Convergence:** The composition of Steps 1 and 2 yields — via a continuous mapping argument in $l^\infty(\mathcal{F})$ — that $\mathbb{G}_N^{DPMM,*}(f)$ and $\mathbb{G}_N^{\pi}(f)$ share the same asymptotic distribution uniformly over $f \in \mathcal{F}$ (a Donsker class). Formally:
$$
\sup_{f \in \mathcal{F}} \left| \mathbb{G}_N^{DPMM,*}(f) - \mathbb{G}_N^{\pi}(f) \right| \xrightarrow{p} 0
$$
This establishes that the DPMM synthetic populations asymptotically target the true finite population in the sense of empirical process weak convergence, validating all downstream inference performed on the synthetic data as asymptotically equivalent to design-based inference on the true population.

---

## 9. Simulation Design: Longitudinal Multistage Complex Survey

To empirically validate the two-stage DPMM proposal, we extend Kim et al.'s (2024) Simulation 3 (stratified two-stage cluster sampling) to a **panel longitudinal setting**, where the same sampled units are followed over $T$ time waves. This design introduces an additional source of within-unit temporal correlation that must be correctly propagated through both the bootstrap and the DPMM.

### 9.1 Design Taxonomy

Three longitudinal survey designs are distinguished by the relationship between the sampling mechanism and the time dimension:

| Design | Sampling Unit | Followed Over Time? | Bootstrap Unit |
| :--- | :--- | :--- | :--- |
| **Panel (cohort)** | Clusters sampled once at baseline | Yes — same units at every $t$ | Cluster |
| **Rotating panel** | Units enter for $k$ waves, then replaced | Partially | Cluster-wave block |
| **Repeated cross-section** | Independent samples at each $t$ | No | Individual per wave |

The **panel design** is the most methodologically demanding and is the primary target of this simulation, as it concentrates three distinct correlation sources: stratum-level, cluster-level, and individual temporal (AR(1)).

### 9.2 Population Generation

The finite population has a **four-level nested structure**: strata $h$ → clusters $i$ → individuals $j$ → time waves $t$.

**Level 1 — Strata.** For each stratum $h = 1, \ldots, H$, draw a stratum effect $\lambda_s \sim \text{Exp}(1)$.

**Level 2 — Clusters.** Within stratum $h$, the number of clusters is $M_h \sim \text{Poisson}(\lambda_h) \times 5 + N_1$. For each cluster $i$, draw a cluster effect $\eta_i \sim \text{Exp}(1)$.

**Level 3 — Individuals.** Within cluster $i$ in stratum $h$, the number of individuals is $M_{hi} \sim \text{Poisson}(\lambda_h + \eta_i) \times 5 + N_2$. Each individual $j$ has a time-invariant covariate $x_{ij}^{inv} \sim \text{Uniform}(0, 20)$ and a random intercept $u_{ij} \sim N(0, \sigma_u^2)$.

**Level 4 — Time.** For each individual $j$, observations over $T$ waves are generated under an AR(1) error structure with autocorrelation $\rho$. The $T \times T$ covariance matrix is $\Sigma_{AR(1)} = [\rho^{|t - t'|}]_{t,t'}$. The outcome model is:
$$
y_{ijt} = x_{ij}^{inv} + x_t^{var} + 0.3\,\eta_i + 0.3\,\lambda_h + u_{ij} + \varepsilon_{ijt}
$$
where $x_t^{var} \sim \text{Uniform}(0, 5)$ is a common time-varying covariate, and $(\varepsilon_{ij1}, \ldots, \varepsilon_{ijT})^\top \sim N(\mathbf{0}, \Sigma_{AR(1)})$.

```r
library(MASS)

gen.popu.longitudinal <- function(N.strata = 50, N.1 = 30, N.2 = 10, T = 4, rho = 0.6) {
  stratum.effect <- rexp(N.strata, rate = 1)
  N.cluster      <- rpois(N.strata, lambda = stratum.effect) * 5 + N.1
  Sigma.ar1      <- rho^abs(outer(1:T, 1:T, "-"))  # AR(1) covariance matrix

  popu <- NULL
  for (s in 1:N.strata) {
    cluster.effect <- rexp(N.cluster[s])
    N.ele <- rpois(N.cluster[s], lambda = stratum.effect[s] + cluster.effect) * 5 + N.2

    for (i in 1:N.cluster[s]) {
      n_i    <- N.ele[i]
      x_inv  <- runif(n_i, 0, 20)         # time-invariant covariate
      u_ij   <- rnorm(n_i, 0, sd = 0.5)  # individual random intercepts
      x_time <- runif(T, 0, 5)           # common time-varying covariate

      for (j in 1:n_i) {
        eps_j <- mvrnorm(1, mu = rep(0, T), Sigma = Sigma.ar1)
        y_jt  <- x_inv[j] + x_time +
                 0.3 * cluster.effect[i] + 0.3 * stratum.effect[s] +
                 u_ij[j] + eps_j

        popu <- rbind(popu, data.frame(
          stratum = s, cluster = i, unit = j, time = 1:T,
          x_inv = x_inv[j], x_time = x_time, y = y_jt,
          M.h = N.cluster[s], M.hi = n_i, N.h = sum(N.ele)
        ))
      }
    }
  }
  popu$N.units <- nrow(popu) / T
  return(popu)
}
```

### 9.3 Panel Sampling

Within each stratum $h$:
1. **Stage 1:** Sample $n_1$ clusters **with replacement** using PPS (probability $\propto M_{hi}$).
2. **Stage 2:** Sample $n_2$ individuals **without replacement** within each selected cluster.
3. **Observation:** All $T$ time waves are observed for every sampled individual.

The sampling weights factor as:
$$
w_{hij} = w_{hi}^{(1)} \times w_{ij}^{(2)} = \left(\frac{n_1 M_{hi}}{N_h}\right)^{-1} \times \frac{M_{hi}}{n_2}
$$

```r
sample.f.longitudinal <- function(popu, n.1 = 5, n.2 = 5) {
  sample.result <- NULL
  for (s in unique(popu$stratum)) {
    popu.s    <- subset(popu, stratum == s & time == 1)
    summary.s <- aggregate(M.hi ~ cluster, data = popu.s, FUN = unique)

    sampled.clus <- sample(summary.s$cluster, size = n.1,
                           prob = summary.s$M.hi, replace = TRUE)

    for (ii in 1:n.1) {
      popu.ii    <- subset(popu.s, cluster == sampled.clus[ii])
      sel.units  <- sample(popu.ii$unit, size = n.2, replace = FALSE)

      # Pull ALL T waves for the selected individuals
      obs.ii     <- subset(popu, stratum == s & cluster == sampled.clus[ii] &
                           unit %in% sel.units)
      M.hi.val   <- unique(obs.ii$M.hi)
      N.h.val    <- unique(obs.ii$N.h)
      obs.ii$w1  <- (n.1 * M.hi.val / N.h.val)^(-1)
      obs.ii$w2  <- M.hi.val / n.2
      obs.ii$weight     <- obs.ii$w1 * obs.ii$w2
      obs.ii$clus.index <- ii    # bootstrap cluster index within stratum
      sample.result <- rbind(sample.result, obs.ii)
    }
  }
  return(sample.result)
}
```

### 9.4 Bootstrap Weight Calibration for Longitudinal Data

The Kim (2024) bootstrap must be adapted to the longitudinal structure. The **fundamental rule** is: resample at the cluster level and propagate the same multiplier to all $n_2$ individuals **and** all $T$ time waves within each selected cluster. This preserves three correlation structures simultaneously: between-stratum, between-cluster, and within-individual temporal.

For each bootstrap replicate $b$ and stratum $h$, draw:
$$
(r_{h1}^{*(b)}, \ldots, r_{hn_1}^{*(b)}) \sim \text{Multinomial}\!\left(n_1 - 1,\, \tfrac{1}{n_1}\mathbf{1}_{n_1}\right)
$$
and set the bootstrap weight for all rows belonging to cluster $i$ in stratum $h$ as:
$$
w_{hijt}^{*(b)} = k_h \cdot r_{hi}^{*(b)} \cdot w_{hij}, \quad k_h = \frac{n_1}{n_1 - 1}
$$
The expansion factor is `each = n.2 * T` (not `each = n.2` as in the cross-sectional case), because the bootstrap multiplier must be broadcast across all individuals within the cluster **and** all time points within each individual.

```r
bootstrap.iter.longitudinal <- function(sample.result, n.1, n.2, T, theta.full) {
  S.j <- sample.result
  k.h <- n.1 / (n.1 - 1)

  for (s in unique(S.j$stratum)) {
    # One Multinomial draw per stratum for the n.1 sampled clusters
    r.h <- as.vector(rmultinom(1, n.1 - 1, rep(1/n.1, n.1)))

    # Expand to n.1 clusters × n.2 individuals × T time waves
    r.expanded       <- rep(r.h, each = n.2 * T)
    idx              <- S.j$stratum == s
    S.j$weight[idx]  <- S.j$weight[idx] * k.h * r.expanded
  }

  # Fit working model (e.g., weighted GEE) on the bootstrap-weighted sample
  estimation.f.long(sample.result = S.j, beta1 = theta.full[2])
}
```

### 9.5 Key Structural Differences from Cross-Sectional Simulation 3

| Dimension | Kim Simulation 3 (Cross-sectional) | Longitudinal Extension |
| :--- | :--- | :--- |
| **Population levels** | Strata → Clusters → Individuals | Strata → Clusters → Individuals → Time waves |
| **Correlation sources** | Stratum effect + cluster effect | Stratum + cluster + individual intercept + AR(1) |
| **Bootstrap expansion** | `rep(r_hi, each = n.2)` | `rep(r_hi, each = n.2 * T)` |
| **Estimation model** | Weighted OLS | Weighted GEE (Liang-Zeger, 1986) |
| **Score function** | Cross-sectional weighted score $S_w(\beta)$ | Weighted GEE score: $\sum_i w_i D_i^\top V_i^{-1}(y_i - \mu_i)$ |
| **Kim concavity** | Satisfied (WLS) | Satisfied (GEE with canonical working correlation) |
| **DPMM Stage 2** | Univariate DPMM per observation | **Multivariate DPMM** over $(y_{ij1}, \ldots, y_{ijT})^\top$ |

### 9.6 Simulation Grid

Mimicking Kim's grid structure, the following parameter configurations are proposed:

| Parameter | Values | Role |
| :--- | :--- | :--- |
| `N.strata` | 20, 50 | Few vs. many strata |
| `N.1` | 20 | Clusters per stratum (population) |
| `N.2` | 30 | Individuals per cluster (population) |
| `n.1` | 5 | Sampled clusters per stratum |
| `n.2` | 10 | Sampled individuals per cluster |
| `T` | 4 | Number of time waves |
| `rho` | 0.3, 0.6 | Weak vs. strong temporal autocorrelation |
| `beta1` | 1.00, 1.01 | Null (size) vs. alternative (power) |

Each cell: 1,000 Monte Carlo replicates × 1,000 bootstrap iterations. Output columns track **NQS** (naïve quasi-score), **BQS** (bootstrap quasi-score), and, for the DPMM extension, **DSYN** (DPMM synthetic population rejection rate).

---

## 10. Conclusion

The literature reveals a convergent evolution between frequentist and Bayesian methods for complex survey data. Frequentists have moved away from rigid parametric sandwich corrections toward flexible, software-friendly rescaled bootstraps (`Kim-24`), democratizing access to complex inference. Meanwhile, Bayesians have transitioned from simple parametric models to pseudo-likelihoods (`Savitsky-16`) and non-parametric predictive frameworks (`Elliot-21`, `Zhou-16`) that naturally incorporate weights and correct for informative selection. 

By unifying these threads—using empirical process theory to bridge bootstrap calibration and non-parametric Bayesian DPMMs—we open a profound methodological frontier. The proposed dissertation framework not only resolves the covariance mismatch between paradigms but provides a vital, privacy-safe mechanism for publishing complex survey data in the modern era.

