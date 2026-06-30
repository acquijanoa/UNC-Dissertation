# Simulation Setup

## Overview

This document specifies the simulation study used to empirically validate the dissertation's
proposed two-stage methodology:

- **Stage 1.** Bootstrap weight calibration via the Kim, Rao & Wang (2024) multinomial
  rescaling bootstrap applied to a concave working parametric model.
- **Stage 2.** Weighted Dirichlet Process Mixture Model (DPMM) fitted under the calibrated
  bootstrap weights to generate privacy-preserving synthetic populations.

The design is directly grounded in Kim & Rao's (2024) Simulation Study (Sections 6.1 and
Supplemental S9 — stratified two-stage cluster sampling), which we extend to a **longitudinal
panel setting** augmented with a DPMM synthetic population generation step.

Critically, the population DGP embeds **$K$ latent trajectory classes** at the individual
level. Each individual $j$ belongs to one of $K$ unobserved classes, each defined by its
own intercept, time trend, covariate slope, and within-class variance. This serves two
purposes simultaneously:

1. **Bootstrap/DPMM validation** — ground-truth class labels allow direct assessment of
   whether the DPMM posterior correctly recovers the latent structure under complex sampling.
2. **Latent trajectory model benchmark** — future Group-Based Trajectory Models (GBTM),
   Latent Class Growth Analysis (LCGA), or Growth Mixture Models (GMM) fitted to
   the simulated data have a known truth to evaluate against, controlling for the bias
   introduced by informative sampling.

Throughout, we deliberately mirror Kim & Rao's notation to make the extension transparent
and the results directly comparable.

---

## 1. Superpopulation Model

### 1.1 Rationale

Kim & Rao (2024) generate their finite population as an i.i.d. draw from a
superpopulation — the standard assumption under which their CLT holds (cf. Corollary
1.3.2.1 of Fuller, 2009). We retain this assumption. The population has a **four-level
nested structure** — strata → clusters → individuals → time waves — which is the natural
extension of Kim & Rao's stratified two-stage design (Supplemental S9) to a panel setting.

### 1.2 Structural Parameters

| Symbol | Role | Value(s) |
|:---|:---|:---|
| $H$ | Number of strata | 20, 50 |
| $M_h$ | Number of clusters in stratum $h$ | $\sim \text{Poisson}(\lambda_h) \times 5 + N_1$ |
| $M_{hi}$ | Number of individuals in cluster $(h,i)$ | $\sim \text{Poisson}(\lambda_h + \eta_i) \times 5 + N_2$ |
| $N_1$ | Minimum clusters per stratum | 30 |
| $N_2$ | Minimum individuals per cluster | 10 |
| $T$ | Number of time waves (panel) | 4 |
| $\rho$ | AR(1) temporal autocorrelation | 0.3, 0.6 |
| $K$ | Number of latent trajectory classes | 3 (default), also 2, 4 |
| $\boldsymbol{\pi}$ | Class mixing proportions $\pi_k = P(z_{ij} = k)$ | See §1.4 |
| $\alpha_k$ | Class-specific intercept | See §1.4 |
| $\gamma_k$ | Class-specific linear time trend | See §1.4 |
| $\sigma_{u,k}^2$ | Class-specific within-individual variance | See §1.4 |

### 1.3 Data-Generating Process

**Level 1 — Strata.** For each stratum $h$, draw a stratum random effect:
$$
\lambda_h \sim \text{Exp}(1)
$$
$\lambda_h$ enters the outcome model (coefficient 0.3) and also serves as the Poisson
lambda for cluster count generation. The number of clusters in stratum $h$ is:
$$
M_h \sim \text{Poisson}(\lambda_h) \times 5 + N_1
$$
The $\times 5$ scaling ensures meaningful size variation around the floor $N_1$;
$E[M_h] \approx 5 + N_1$. Because $M_h$ depends on $\lambda_h$ (which enters the
outcome), stratum size is correlated with the outcome — mirroring Kim's S9 design.

**Level 2 — Clusters.** Within stratum $h$, for each cluster $i$, draw a cluster
random effect:
$$
\eta_i \sim \text{Exp}(1)
$$
$\eta_i$ enters the outcome model (coefficient 0.3). The number of individuals in
cluster $(h, i)$ is:
$$
M_{hi} \sim \text{Poisson}(\lambda_h + \eta_i) \times 5 + N_2
$$
Using $\lambda_h + \eta_i$ as the Poisson lambda (mean $\approx 2$ under Exp(1))
makes $M_{hi}$ depend on both $\lambda_h$ and $\eta_i$, both of which enter the
outcome — so **PPS sampling proportional to $M_{hi}$ is informative** with respect
to the outcome distribution, exactly as in Kim & Rao's S9.

**Level 3 — Individuals.** Within cluster $(h, i)$, each individual $j$ has:
- A **time-invariant covariate** $x_{ij}^{inv} \sim \text{Uniform}(0, 20)$
  (mirroring Kim & Rao's $x_i \sim U(-5, 5)$ for the cross-sectional case).
- A **latent class assignment** drawn once per individual:
$$
z_{ij} \sim \text{Categorical}(\pi_1, \ldots, \pi_K), \quad \sum_{k=1}^K \pi_k = 1
$$
- A **class-specific individual random intercept**
  $u_{ij} \sim N(0,\, \sigma_{u, z_{ij}}^2)$, capturing within-class heterogeneity
  (set $\sigma_{u,k}^2 = 0$ for the LCGA/GBTM variant; nonzero for GMM).

**Level 4 — Time.** For each individual $j$ in class $k = z_{ij}$, observed across $T$
waves, draw a **common time-varying covariate** $x_t^{var} \sim \text{Uniform}(0, 5)$
and class-specific AR(1) errors:
$$
(\varepsilon_{ij1}, \ldots, \varepsilon_{ijT})^\top \sim N\!\left(\mathbf{0},\; \Sigma_{AR(1)}\right),
\quad [\Sigma_{AR(1)}]_{t,t'} = \rho^{|t - t'|}
$$

The **class-conditional outcome model** is:
$$
y_{ijt} \mid z_{ij} = k \;=\;
  \alpha_k
  + \beta_{1k}\, x_{ij}^{inv}
  + \gamma_k \cdot t
  + 0.3\,\eta_i + 0.3\,\lambda_h
  + u_{ij} + \varepsilon_{ijt}
$$

where $\alpha_k$ is the class intercept, $\beta_{1k}$ is the class-specific slope on the
time-invariant covariate, and $\gamma_k$ is the class-specific **linear time trend** —
the primary parameter distinguishing trajectory shapes across classes.

For hypothesis testing purposes, the **marginal** slope on $x_{ij}^{inv}$ averaged across
classes is:
$$
\bar{\beta}_1 = \sum_{k=1}^K \pi_k \beta_{1k}
$$
and the null of interest remains $H_0: \bar{\beta}_1 = \bar{\beta}_1^{(0)}$, which is
estimated from the weighted GEE in Stage 1 (marginalizing over latent classes).

The cluster and stratum effects enter additively with coefficient 0.3, making the design
**informative**: PPS probabilities based on $M_{hi}$ correlate with $\eta_i$, which enters
the outcome — replicating the inflated-Type-I-error scenario of Kim & Rao (2024, Table 1).
Importantly, since $\eta_i$ is also correlated with which class individuals in cluster $i$
tend to belong to (through the cluster-size mechanism), the sampling is informative
**with respect to both the marginal mean and the latent class distribution**.

### 1.4 Default Latent Class Parameters ($K = 3$)

The following class structure is designed to produce clearly separated, interpretable
trajectory shapes that mirror patterns found in child BMI trajectory research (e.g.,
low-stable, moderate-rising, and high-declining subgroups):

| Class $k$ | $\pi_k$ | $\alpha_k$ | $\beta_{1k}$ | $\gamma_k$ (time trend/wave) | $\sigma_{u,k}$ | Shape |
|:---:|:---:|:---:|:---:|:---:|:---:|:---|
| 1 | 0.50 | 5.0 | 0.8 | +0.2 | 0.3 | Low-stable, slight rise |
| 2 | 0.30 | 10.0 | 1.0 | +0.8 | 0.5 | Moderate-rising |
| 3 | 0.20 | 15.0 | 1.2 | −0.4 | 0.7 | High-declining |

> **Separation heuristic.** The class intercepts are spaced at least $2\sigma$ apart
> (approximately), ensuring the DPMM has sufficient signal to recover the $K$ components
> under moderate sample sizes. Class separation can be deliberately reduced in a sensitivity
> arm (e.g., $\alpha_k \in \{8, 10, 12\}$) to study DPMM behavior under weak separation.

For the **$K = 2$ variant**: merge classes 1 and 2 into a single "non-high" group
($\pi_1 = 0.80, \alpha_1 = 7, \gamma_1 = +0.4$) and retain class 3 as the "high-declining"
group ($\pi_2 = 0.20$). For the **$K = 4$ variant**: split the moderate class into
"moderate-stable" and "moderate-accelerating" subgroups.

---

## 2. Sampling Design

The sampling design follows Kim & Rao's Supplemental S9 (stratified two-stage cluster
sampling with clusters selected with replacement by PPS).

### 2.1 Stage 1 — Cluster Sampling (PPS, with replacement)

Within stratum $h$, select $n_1 = 5$ clusters **with replacement** using probability
proportional to size (PPS), with selection probability:
$$
p_{hi} = \frac{M_{hi}}{\sum_{i'=1}^{M_h} M_{hi'}}
$$

This makes the design informative because $M_{hi}$ is a deterministic function of $\eta_i$
and $\lambda_h$, both of which appear in the outcome model.

### 2.2 Stage 2 — Individual Sampling (SRSWOR within cluster)

Within each selected cluster $(h, i)$, sample $n_2 = 10$ individuals **without
replacement** using simple random sampling.

### 2.3 Observation Rule

All $T = 4$ time waves are observed for every sampled individual (pure panel design; no
attrition). This concentrates three correlation structures in the observed sample:
stratum-level, cluster-level, and individual-level temporal (AR(1)).

### 2.4 Sampling Weights

The design weight for individual $j$ in cluster $(h, i)$ across all $T$ time points
factors as:
$$
w_{hij} = w_{hi}^{(1)} \times w_{ij}^{(2)}
        = \left(\frac{n_1 \cdot M_{hi}}{M_h^{tot}}\right)^{-1} \times \frac{M_{hi}}{n_2}
$$
where $M_h^{tot} = \sum_{i'} M_{hi'}$ is the total population size in stratum $h$. The
same scalar weight $w_{hij}$ is broadcast to all $T$ observations for that individual.

---

## 3. Bootstrap Weight Calibration (Stage 1 of Proposed Method)

### 3.1 Procedure

Following Kim & Rao (2024, Example 1 and Supplemental S4), for each bootstrap replicate
$b = 1, \ldots, B$ and each stratum $h$:

1. Draw a rescaling vector from a Multinomial distribution with $n_1 - 1$ trials and
   equal probability $n_1^{-1}$:
$$
(r_{h1}^{*(b)}, \ldots, r_{hn_1}^{*(b)}) \sim \text{Multinomial}\!\left(n_1 - 1,\; \tfrac{1}{n_1}\mathbf{1}_{n_1}\right)
$$

2. Set the bootstrap weight for all rows belonging to cluster $i$ in stratum $h$ as:
$$
w_{hijt}^{*(b)} = k_h \cdot r_{hi}^{*(b)} \cdot w_{hij}, \quad k_h = \frac{n_1}{n_1 - 1}
$$

The rescaling factor $k_h$ ensures that:
$$
E_*\!\left\{\widehat{S}_w^{*(b)}(\theta)\right\} = \widehat{S}_w(\theta)
\quad \text{and} \quad
V_*\!\left\{\widehat{S}_w^{*(b)}(\theta)\right\} = \widehat{V}\!\left\{\widehat{S}_w(\theta) \mid \mathcal{F}_N\right\}
$$
which are Kim & Rao's Conditions C6–C7, guaranteeing bootstrap CLT validity.

### 3.2 Longitudinal Adaptation

The key structural difference from the cross-sectional case is how the bootstrap multiplier
is broadcast. Because we sample clusters and then observe all $n_2 \times T$ rows for each
cluster, the expansion must cover **all individuals within the cluster across all time
waves**:

```r
r.expanded <- rep(r.h, each = n.2 * T)   # longitudinal case
# vs.
r.expanded <- rep(r.h, each = n.2)        # cross-sectional Kim & Rao S9
```

This preserves the within-individual temporal correlation in the bootstrap distribution.

### 3.3 Working Model for Stage 1

Kim's regularity conditions (C1: concavity; C5: uniform convergence of the information)
require the Stage 1 working model to be globally concave. We use a **weighted GEE
(Generalized Estimating Equations)** model with a working AR(1) correlation structure and
identity link:

$$
\widehat{S}_w^{GEE}(\beta) = \sum_{j \in A} w_j \, D_j^\top V_j^{-1}(y_j - X_j \beta) = \mathbf{0}
$$

where $D_j = \partial \mu_j / \partial \beta^\top$, $V_j$ is the working covariance matrix
(AR(1) with $T \times T$ structure), and $y_j = (y_{ij1}, \ldots, y_{ijT})^\top$. Under a
canonical working correlation and identity link, the GEE score is linear in $\beta$ and
hence globally concave — satisfying Kim's Condition C1.

---

## 4. DPMM Synthetic Population Generation (Stage 2 of Proposed Method)

### 4.1 Decoupling Rationale

As established in the dissertation theory (Synthesis I, §8.2–8.3), the DPMM
log-likelihood is not globally concave and therefore Kim's bootstrap conditions cannot be
applied to it directly. The two-stage decoupling resolves this: the calibrated bootstrap
weight vectors $\{w_{hijt}^{*(b)}\}$ from Stage 1 are treated as **fixed external inputs**
to the DPMM. The DPMM itself is never subjected to score-matching conditions.

### 4.2 Model

For each bootstrap replicate $b$, fit a **multivariate DPMM** to the $T$-dimensional
outcome vectors $\mathbf{y}_j = (y_{ij1}, \ldots, y_{ijT})^\top$ using the calibrated
pseudo-likelihood:

$$
l_{DPMM}^{*(b)}(G \mid \mathbf{y}) = \sum_{j \in A} w_{j}^{*(b)} \log \int K(\mathbf{y}_j;\, \phi)\, dG(\phi)
$$

where $G \sim DP(\alpha, G_0)$, $G_0 = \text{NIW}(\mu_0, \kappa_0, \Psi_0, \nu_0)$
(Normal-Inverse-Wishart base measure), and $K(\cdot; \phi)$ is a multivariate Gaussian
kernel with mean $\mu$ and covariance $\Sigma$ (so $\phi = (\mu, \Sigma)$).

### 4.3 Synthetic Population Draw

For each bootstrap replicate $b$, after fitting the DPMM via MCMC (truncated stick-breaking
with $L = 20$ components, 2,000 burn-in, 1,000 post-burn-in draws), draw the unobserved
population units from the posterior predictive distribution:

$$
p^*(y_{new}) = \int K(y_{new};\, \phi)\, d\hat{G}^{*(b)}(\phi)
$$

Because $K$ is a continuous Gaussian kernel, the synthetic values are **not exact copies**
of any observed unit, satisfying privacy requirements.

### 4.4 Evaluation Criterion

The quality of each synthetic population is evaluated via the **HT-MMD** (Synthesis II,
§5.2):

$$
\widehat{\text{MMD}}^2_w = \left\|\hat{\mu}_w^{DPMM} - \hat{\mu}_w^{true}\right\|_{\mathcal{H}}^2
$$

where $\hat{\mu}_w^{DPMM}$ is the HT-weighted kernel mean embedding of the synthetic
population and $\hat{\mu}_w^{true}$ is the HT-weighted kernel mean embedding of the full
finite population. Rejection of the null $H_0: \text{MMD}^2 = 0$ indicates the synthetic
population does not match the true population distribution.

### 4.5 DPMM as a Latent Trajectory Recovery Engine

The $K$-class DGP creates a direct theoretical link between the DPMM Stage 2 and
latent trajectory/class models:

- **DPMM ↔ GMM.** The DPMM with Gaussian kernels is the Bayesian nonparametric
  generalization of the $K$-component Gaussian Mixture Model. Under the $K$-class DGP,
  the posterior of the DPMM should concentrate on $K$ active mixture components (one per
  class), with component weights converging to $(\pi_1, \ldots, \pi_K)$ and component
  means converging to the class-specific trajectory means.

- **DPMM ↔ GBTM/LCGA.** When $\sigma_{u,k}^2 = 0$ for all $k$ (no within-class
  individual heterogeneity), the DGP reduces to a GBTM/LCGA structure. The DPMM
  recovers this by collapsing each mixture component to a degenerate Gaussian with
  very small variance — the posterior will show $K$ tight clusters.

- **Design correction matters.** Without bootstrap weight calibration, the DPMM
  pseudo-likelihood is distorted by informative sampling: units from large clusters
  (high $\eta_i$, high $M_{hi}$, high $w_i$) are over-represented in the effective
  likelihood, systematically biasing the class proportion estimates toward whichever
  class is more prevalent in large clusters. The Stage 1 calibrated weights correct
  this by down-weighting over-sampled units — a bias that standard GBTM/GMM software
  (e.g., `lcmm`, `flexmix`, `mplus`) does not address.

**Class recovery metrics** added to §7:
- **ARI** (Adjusted Rand Index): agreement between DPMM-assigned classes and ground-truth
  $z_{ij}$, averaged over posterior draws.
- **$\hat{K}$**: effective number of DPMM components with posterior weight $> 0.01$;
  should concentrate at $K$ under the default parameters.
- **$\hat{\pi}_k$ bias**: $|\hat{\pi}_k - \pi_k|$ for each class; should shrink to zero
  under correct weight calibration and increase under naïve (unweighted) DPMM.

---

## 5. R Implementation

### 5.0 Correspondence with Kim's Reference Code

Kim's single-stage PPS code (`gen.popu` / `one.iter`) is the canonical reference. The
table below maps each piece to our longitudinal extension:

| Kim (Section 6.1, single-stage) | Our Extension (Stratified two-stage, longitudinal) |
|:---|:---|
| `popu$y = 1 + x + ε` | `y_ijt = α_k + β_{1k} x_inv + γ_k t + 0.3η_i + 0.3λ_h + u_ij + ε_ijt` |
| `cov_xy = 6*I + 1` → `sel.prob` | PPS ∝ $M_{hi}$, where $M_{hi}$ depends on `cluster.effect` (outcome-entering) |
| `weight = 1/(n * sel.prob)` | `weight = w1 * w2` (two-stage factor weights) |
| `wei.boo = k.n * rmultinom(1, n-1, rep(1/n, n))` | Same draw; broadcast `each = n.2 * T` within each cluster-stratum block |
| `sample.result.boo$weight = weight * wei.boo` | `S.j$weight[idx] = weight * k.h * r.expanded` |
| `estimation.f`: WLS on `(y ~ x)` | `estimation.f.long`: weighted GEE on `(y_t ~ x_inv + t)` with AR(1) working correlation |
| `W.n.ori > quantile(result.i[,1], 0.95)` | Same; bootstrap distribution replaces $\chi^2(1)$ |

**Critical structural difference:** Kim's informative design creates selection–residual
dependence through `(y - mean_y) * x > 0`. In the two-stage case the mechanism is
replaced by **cluster-size informativeness**: because $M_{hi}$ depends on `cluster.effect`
(which enters $y$), PPS selection ∝ $M_{hi}$ over-samples clusters with large random
effects — the direct analogue of Kim's `cov_xy` construction.

### 5.1 Population Generation

```r
library(MASS)

# Default K = 3 latent class parameters (see §1.4)
make.class.params <- function(K = 3) {
  if (K == 3) {
    list(
      pi    = c(0.50, 0.30, 0.20),          # mixing proportions
      alpha = c(5.0,  10.0, 15.0),          # class intercepts
      beta1 = c(0.8,   1.0,  1.2),          # slopes on x_inv
      gamma = c(0.2,   0.8, -0.4),          # linear time trends
      sigma.u = c(0.3,  0.5,  0.7)          # within-class SD (GMM)
      # set sigma.u = rep(0, K) for LCGA/GBTM variant
    )
  } else if (K == 2) {
    list(
      pi    = c(0.80, 0.20),
      alpha = c(7.0,  15.0),
      beta1 = c(0.9,   1.2),
      gamma = c(0.4,  -0.4),
      sigma.u = c(0.4,  0.7)
    )
  } else {
    stop("Define class.params for K = ", K)
  }
}

gen.popu.longitudinal <- function(N.strata = 50, N.1 = 30, N.2 = 10,
                                   T = 4, rho = 0.6, K = 3,
                                   cp = make.class.params(K)) {
  # Mirrors Kim & Rao S9 exactly: Exp(1) effects, rpois() * 5 + floor
  stratum.effect <- rexp(N.strata, rate = 1)
  N.cluster      <- rpois(N.strata, lambda = stratum.effect) * 5 + N.1
  Sigma.ar1      <- rho^abs(outer(1:T, 1:T, "-"))  # T×T AR(1) covariance
  t.vec          <- 1:T

  popu <- NULL
  for (s in 1:N.strata) {
    cluster.effect <- rexp(N.cluster[s], rate = 1)
    N.ele          <- rpois(N.cluster[s],
                            lambda = stratum.effect[s] + cluster.effect) * 5 + N.2

    for (i in 1:N.cluster[s]) {
      n_i    <- N.ele[i]
      x_inv  <- runif(n_i, 0, 20)           # time-invariant covariate
      x_time <- runif(T, 0, 5)              # common time-varying covariate

      # ---- Draw latent class for each individual (once per person) --------
      z_ij   <- sample(1:K, size = n_i, replace = TRUE, prob = cp$pi)

      for (j in 1:n_i) {
        k     <- z_ij[j]                    # individual j's latent class
        u_ij  <- rnorm(1, 0, sd = cp$sigma.u[k])  # class-specific intercept
        eps_j <- mvrnorm(1, mu = rep(0, T), Sigma = Sigma.ar1)

        # Class-conditional trajectory
        y_jt  <- cp$alpha[k] +
                 cp$beta1[k] * x_inv[j] +
                 cp$gamma[k] * t.vec +      # linear time trend per class
                 0.3 * cluster.effect[i] +
                 0.3 * stratum.effect[s] +
                 u_ij + eps_j

        popu <- rbind(popu, data.frame(
          stratum = s, cluster = i, unit = j, time = t.vec,
          x_inv   = x_inv[j], x_time = x_time, y = y_jt,
          class   = k,                      # ground-truth class label
          M.h     = N.cluster[s], M.hi = n_i,
          N.h     = sum(N.ele)
        ))
      }
    }
  }
  popu$N.units <- nrow(popu) / T
  return(popu)
}
```

> **Note on `class` column.** The `class` column stores the ground-truth latent class
> label $z_{ij}$ for every row. This column is **not passed to the sampler or the
> estimator** — it is retained only in the full population object (`popu`) for
> post-hoc evaluation of class recovery. The sample object (`samp`) returned by
> `sample.f.longitudinal()` will contain `class` as an observed column **only if**
> you choose to retain it; in the primary analysis it should be dropped to simulate
> the latent (unobserved) scenario.

### 5.2 Panel Sampling

```r
sample.f.longitudinal <- function(popu, n.1 = 5, n.2 = 10) {
  sample.result <- NULL
  for (s in unique(popu$stratum)) {
    popu.s    <- subset(popu, stratum == s & time == 1)
    summary.s <- aggregate(M.hi ~ cluster, data = popu.s, FUN = unique)

    # Stage 1: PPS cluster sampling with replacement
    sampled.clus <- sample(summary.s$cluster, size = n.1,
                           prob = summary.s$M.hi, replace = TRUE)

    for (ii in 1:n.1) {
      popu.ii   <- subset(popu.s, cluster == sampled.clus[ii])
      sel.units <- sample(popu.ii$unit, size = n.2, replace = FALSE)

      # Pull all T waves for selected individuals
      obs.ii    <- subset(popu, stratum == s &
                            cluster == sampled.clus[ii] &
                            unit %in% sel.units)

      M.hi.val  <- unique(obs.ii$M.hi)
      N.h.val   <- unique(obs.ii$N.h)
      n.h.tot   <- n.1 * M.hi.val / N.h.val    # first-stage sampling fraction

      obs.ii$w1         <- 1 / n.h.tot           # first-stage weight
      obs.ii$w2         <- M.hi.val / n.2         # second-stage weight
      obs.ii$weight     <- obs.ii$w1 * obs.ii$w2
      obs.ii$clus.index <- ii

      sample.result <- rbind(sample.result, obs.ii)
    }
  }
  return(sample.result)
}
```

### 5.3 Bootstrap Weight Calibration

```r
bootstrap.iter.longitudinal <- function(sample.result, n.1, n.2, T,
                                         theta.full) {
  S.j <- sample.result
  k.h <- n.1 / (n.1 - 1)          # rescaling factor per Kim & Rao

  for (s in unique(S.j$stratum)) {
    # Multinomial draw: n.1 - 1 trials, equal probability 1/n.1
    r.h <- as.vector(rmultinom(1, n.1 - 1, rep(1 / n.1, n.1)))

    # Broadcast to all n.2 individuals × T time waves per cluster
    r.expanded      <- rep(r.h, each = n.2 * T)
    idx             <- S.j$stratum == s
    S.j$weight[idx] <- S.j$weight[idx] * k.h * r.expanded
  }

  # Fit working GEE model on bootstrap-reweighted sample
  estimation.f.long(sample.result = S.j, beta1.ref = theta.full[2])
}
```

### 5.4 Monte Carlo Loop

```r
run.simulation <- function(N.strata, rho, beta1.null, n.MC = 1000, B = 1000,
                            n.1 = 5, n.2 = 10, T = 4) {
  results <- matrix(NA, nrow = n.MC, ncol = 3,
                    dimnames = list(NULL, c("NQS", "BQS", "DSYN")))

  for (m in 1:n.MC) {
    popu          <- gen.popu.longitudinal(N.strata = N.strata, T = T, rho = rho)
    samp          <- sample.f.longitudinal(popu, n.1 = n.1, n.2 = n.2)
    theta.full    <- fit.gee.weighted(samp)        # full-sample estimate

    # Naïve quasi-score test (no design correction)
    results[m, "NQS"] <- naive.qs.test(samp, beta1.null)

    # Bootstrap quasi-score test (Kim & Rao calibration)
    boot.stats        <- replicate(B, bootstrap.iter.longitudinal(
                           samp, n.1, n.2, T, theta.full))
    results[m, "BQS"] <- mean(boot.stats > qs.test(samp, beta1.null))

    # DPMM synthetic population test (HT-MMD)
    results[m, "DSYN"] <- dpmm.htmmd.test(samp, popu, B.boot = B)
  }

  return(results)
}
```

---

## 6. Simulation Grid

The grid mirrors Kim & Rao (2024, Table 1) and extends it to the longitudinal setting.
Each cell is run with $n_{MC} = 1{,}000$ Monte Carlo replicates and $B = 1{,}000$ bootstrap
iterations. The nominal significance level is $\alpha = 0.05$.

### 6.1 Primary Grid (Hypothesis Testing Focus)

| Parameter | Values | Role |
|:---|:---|:---|
| `N.strata` ($H$) | 20, 50 | Few vs. many strata |
| `N.1` | 30 | Minimum clusters per stratum |
| `N.2` | 10 | Minimum individuals per cluster |
| `n.1` | 5 | Sampled clusters per stratum |
| `n.2` | 10 | Sampled individuals per cluster |
| `T` | 4 | Time waves (panel) |
| `rho` ($\rho$) | 0.3, 0.6 | Weak vs. strong temporal autocorrelation |
| `beta1.null` ($\bar{\beta}_1^{(0)}$) | 1.00, 1.01 | Size ($H_0$ true) vs. power ($H_0$ false) |
| `K` | **3** (fixed) | Latent classes (default structure from §1.4) |

This yields $2 \times 2 \times 2 = 8$ cells. Total: $8{,}000$ outer replicates $\times$
$1{,}000$ bootstrap draws each.

### 6.2 Latent Class Recovery Grid (Class Model Focus)

Run separately with `beta1.null` fixed at the true $\bar{\beta}_1$ (size scenario only),
varying the parameters that affect DPMM class recovery:

| Parameter | Values | Role |
|:---|:---|:---|
| `K` | 2, 3, 4 | Number of true latent classes |
| `rho` ($\rho$) | 0.3, 0.6 | Temporal correlation (affects within-person information) |
| `N.strata` ($H$) | 20, 50 | Sample size proxy |
| **Separation** | Full (§1.4 defaults), Reduced ($\Delta\alpha = 2$) | Tests DPMM under weak vs. strong class separation |

This yields $3 \times 2 \times 2 \times 2 = 24$ additional cells.
Metrics reported: ARI, $\hat{K}$, $\hat{\pi}_k$ bias (see §7).

---

## 7. Output Metrics

### 7.1 Hypothesis Testing Metrics (Primary Grid)

| Metric | Symbol | Description |
|:---|:---|:---|
| Naïve quasi-score rejection rate | **NQS** | Type I error / power of naïve GEE test, $\chi^2(1)$ reference |
| Bootstrap quasi-score rejection rate | **BQS** | Type I error / power of Kim & Rao calibrated bootstrap |
| DPMM synthetic rejection rate | **DSYN** | HT-MMD rejection rate: synthetic population vs. true population |
| Design effect | **DEFF** | Estimated $\widehat{V}_{design} / \widehat{V}_{SRS}$ for $\hat{\bar{\beta}}_1$ |
| Synthetic coverage | **COV** | Coverage of 95% credible intervals from DPMM posteriors |

The primary comparison is **NQS vs. BQS**: NQS should show severely inflated Type I error
under $\bar{\beta}_1^{(0)} = 1.00$ (replicating Kim & Rao, Table 1), while BQS should
control it at $\approx 0.05$. **DSYN** quantifies whether the DPMM synthetic populations
are distributionally equivalent to the true finite population under the HT-MMD.

### 7.2 Latent Class Recovery Metrics (Class Recovery Grid)

| Metric | Symbol | Description |
|:---|:---|:---|
| Adjusted Rand Index | **ARI** | Agreement between DPMM posterior class assignments and ground-truth $z_{ij}$; 0 = chance, 1 = perfect |
| Effective component count | **$\hat{K}$** | Number of DPMM components with posterior weight $> 0.01$; should concentrate at $K$ |
| Class proportion bias | **$\Delta\hat{\pi}_k$** | $\lvert \hat{\pi}_k - \pi_k \rvert$ per class; calibrated DPMM should reduce vs. naïve DPMM |
| Class mean RMSE | **RMSE$_\mu$** | $\sqrt{\frac{1}{K}\sum_k (\hat{\alpha}_k - \alpha_k)^2}$; trajectory intercept recovery |
| Class trend RMSE | **RMSE$_\gamma$** | $\sqrt{\frac{1}{K}\sum_k (\hat{\gamma}_k - \gamma_k)^2}$; trajectory slope recovery |

**Key comparison:** $\Delta\hat{\pi}_k$ and RMSE$_\mu$ under the **calibrated** DPMM
(using Stage 1 bootstrap weights) vs. the **naïve** DPMM (using raw survey weights or
no weights). The informative sampling hypothesis predicts that naïve DPMM will
systematically over-estimate $\pi_k$ for classes over-represented in large clusters,
while the calibrated DPMM corrects this.

---

## 8. Connection to FLOR Study

The simulation design was calibrated to resemble key features of the FLOR/HCHS-SOL
study described in the methods document:

| FLOR Feature | Simulation Analogue |
|:---|:---|
| Four U.S. field centers (Bronx, Chicago, Miami, San Diego) | $H = 4$ strata (small $H$ scenario) |
| Mother–child dyads nested within field centers | Individuals nested within clusters within strata |
| Repeated HCHS/SOL baseline + FLOR follow-up visits | $T = 2$ wave variant of the longitudinal design |
| Unequal sampling across centers (oversampling of specific Latina subgroups) | PPS design inducing informative sampling |
| Multiple imputation of covariates (10 imputed datasets) | DPMM posterior draws naturally produce multiply-imputed synthetic completions |
| BMI z-score as continuous primary outcome | Continuous $y_{ijt}$ generated from class-conditional linear model |
| Latent subgroups of child BMI trajectories (clinical intuition: lean-stable, overweight-rising, obese-declining) | $K = 3$ latent trajectory classes (§1.4) with distinct $\alpha_k$, $\gamma_k$ |
| Future LCGA / GMM of child BMI trajectories | Class recovery grid (§6.2) provides power / ARI curves as a function of $K$, $\rho$, and separation |

This parallel motivates interpreting the simulation's **BQS** results as a benchmark for
what hypothesis testing performance one should expect when applying the same bootstrap
calibration approach to FLOR model coefficients, and **DSYN** as a benchmark for the
quality of DPMM-generated synthetic FLOR data released under privacy constraints.

---

## References

- Kim, J. K., Rao, J. N. K., & Wang, Z. (2024). Hypotheses testing from complex survey
  data using bootstrap weights: A unified approach. *Journal of the American Statistical
  Association*, 119(546), 1229–1239.
- Han, Q., & Wellner, J. A. (2021). Complex sampling designs: Uniform limit theorems and
  applications. *Annals of Statistics*, 49(1), 459–485.
- Ghosal, S., Ghosh, J. K., & van der Vaart, A. W. (2000). Convergence rates of posterior
  distributions. *Annals of Statistics*, 28(2), 500–531.
- Ghosal, S., & van der Vaart, A. W. (2007). Convergence rates of posterior distributions
  for non-i.i.d. observations. *Annals of Statistics*, 35(1), 192–223.
- Savitsky, T. D., & Toth, D. (2016). Bayesian estimation under informative sampling.
  *Electronic Journal of Statistics*, 10(1), 1677–1708.
- Williams, M. R., & Savitsky, T. D. (2021). Uncertainty estimation for pseudo-Bayesian
  inference under complex sampling. *International Statistical Review*, 89(1), 72–107.
- Elliott, M. R., & Xia, F. (2021). A flexible two-phase sampling design for estimation
  of cumulative incidence functions. *Biometrics*, 77(1), 28–41.
- Dong, Q., Elliott, M. R., & Raghunathan, T. E. (2014). A nonparametric method to
  generate synthetic populations to adjust for complex sampling design features.
  *Survey Methodology*, 40(1), 29–46.
- Zhou, H., Elliott, M. R., & Raghunathan, T. E. (2016). Multiple imputation in
  two-stage cluster samples using the weighted finite population Bayesian bootstrap.
  *Journal of Survey Statistics and Methodology*, 4(2), 139–170.
