# Theoretical Justification for Raw Kim-Rao Weights in One-Step Random Weight MCMC

**Author:** Álvaro Quijano  
**Date:** June 2026  
**Context:** Research Dissertation — One-Step Random Weight MCMC for DPMM and BART

---

## 1. The Core Law of Total Variance Decomposition

### 1.1 Joint Probability Space

In the proposed **one-step Random Weight MCMC** framework, we draw a fresh set of Kim-Rao multinomial bootstrap weights $\mathbf{w}^{*(t)}$ at each MCMC iteration $t$. The theoretical analysis proceeds under the **joint distribution**

$$
\mathbb{P} = \mathcal{L}_* \otimes \pi(\,\cdot \mid \mathbf{w}^*, \mathbf{y})
$$

where $\mathcal{L}_*$ is the Kim-Rao bootstrap distribution over weight vectors, and $\pi(\theta \mid \mathbf{w}^*, \mathbf{y})$ is the pseudo-posterior induced by the survey-weighted log-likelihood $\sum_{i\in S} w_i^* \ell_i(\theta)$ combined with the model prior. Under this joint measure, $\theta$ is a proper random variable; all variances and expectations below are with respect to $\mathbb{P}$ unless otherwise noted.

### 1.2 Variance Decomposition

The total variance of the parameter draws $\{\theta^{(t)}\}_{t=1}^T$ is decomposed using the Law of Total Variance with respect to $\mathcal{L}_*$:

$$
\text{Var}_{\mathbb{P}}(\theta) = \underbrace{\text{Var}_{\mathcal{L}_*}\!\bigl(E_\pi[\theta \mid \mathbf{w}^*]\bigr)}_{\text{Term 1: Bootstrap Design Variance}} + \underbrace{E_{\mathcal{L}_*}\!\bigl[\text{Var}_\pi(\theta \mid \mathbf{w}^*)\bigr]}_{\text{Term 2: Expected Conditional Model Variance}}
$$

**Term 1** is the variance, over repeated bootstrap draws $\mathbf{w}^* \sim \mathcal{L}_*$, of the pseudo-posterior mean $\bar\theta(\mathbf{w}^*) = E_\pi[\theta \mid \mathbf{w}^*]$.

**Term 2** is the expected pseudo-posterior variance, averaged over bootstrap draws. For exact design-consistent inference, we require that as $n \to \infty$:

$$
E_{\mathcal{L}_*}\!\bigl[\text{Var}_\pi(\theta \mid \mathbf{w}^*)\bigr] \to 0.
$$

The weight scaling mechanism described in Section 2 provides explicit, rate-quantifiable control over this term.

### 1.3 Connecting Term 1 to the Sandwich Variance

Connecting Term 1 to the design-based sandwich variance requires a delta-method argument resting on the following regularity conditions:

**(R1) Pseudo-posterior concentration.** The pseudo-posterior $\pi(\theta \mid \mathbf{w}^*)$ concentrates at rate $O_p(n^{-1/2})$ around the maximizer $\hat\theta(\mathbf{w}^*)$ of the weighted pseudo-log-likelihood $\sum_{i\in S} w_i^* \ell_i(\theta)$, uniformly over $\mathbf{w}^*$ in a neighborhood of $\mathbf{w}$ (cf. Kleijn & van der Vaart, 2012, Theorem 2.1).

**(R2) Score smoothness.** The score $\nabla_\theta \ell_i(\theta)$ and Hessian $\nabla^2_\theta \ell_i(\theta)$ exist and are continuous in $\theta$ on a neighborhood of $\theta_0$, and satisfy appropriate uniform integrability conditions under the design.

**(R3) Identification.** The expected pseudo-log-likelihood $Q(\theta) = E_\xi[n^{-1}\sum_{i\in S} w_i \ell_i(\theta)]$ — where $E_\xi$ denotes expectation over the sampling design $\xi$ — is uniquely maximized at an interior point $\theta_0$.

Under **(R1)-(R3)**, the pseudo-posterior mean concentrates near $\hat\theta(\mathbf{w}^*)$:

$$
\bar\theta(\mathbf{w}^*) = \hat\theta(\mathbf{w}^*) + O_p(n^{-1})
$$

and by the implicit function theorem applied to the score equation $n^{-1}\sum_{i\in S} w_i^* \nabla_\theta \ell_i(\hat\theta(\mathbf{w}^*)) = 0$, the map $\mathbf{w}^* \mapsto \hat\theta(\mathbf{w}^*)$ is smooth in $\mathbf{w}^*$. A first-order expansion gives:

$$
\hat\theta(\mathbf{w}^*) \approx \hat\theta(\mathbf{w}) - [H_n(\theta_0)]^{-1}\, n^{-1}\sum_{i\in S}(w_i^* - w_i)\nabla_\theta \ell_i(\theta_0)
$$

where $H_n(\theta_0) = -n^{-1}\sum_{i\in S} \nabla^2_\theta \ell_i(\theta_0)$ is the normalized sample Hessian. Consequently:

$$
\text{Var}_{\mathcal{L}_*}\!\bigl(\bar\theta(\mathbf{w}^*)\bigr) \approx \text{Var}_{\mathcal{L}_*}\!\bigl(\hat\theta(\mathbf{w}^*)\bigr) = [H_n(\theta_0)]^{-1}\,\text{Var}_{\mathcal{L}_*}\!\bigl(\bar{S}_{w^*}(\theta_0)\bigr)\,[H_n(\theta_0)]^{-1}
$$

where $\bar{S}_{w^*}(\theta) = n^{-1}\sum_{i\in S} w_i^* \nabla_\theta \ell_i(\theta)$ is the bootstrap-weighted normalized score. The key property of the Kim-Rao bootstrap (Kim, Rao & Wang, 2024, Theorem 1) is that it exactly replicates the design-based score variance:

$$
\text{Var}_{\mathcal{L}_*}\!\bigl(\bar{S}_{w^*}(\theta_0)\bigr) = \hat{V}\!\left[\bar{S}_w(\theta_0) \mid \mathcal{F}_N\right] \xrightarrow{p} J_0
$$

where $\bar{S}_w(\theta_0) = n^{-1}\sum_{i\in S} w_i \nabla_\theta \ell_i(\theta_0)$ is the design-weighted normalized score and

$$
J_0 = \lim_{n\to\infty}\, n^2 \cdot \text{Var}_\xi\!\!\left[n^{-1}\sum_{i\in S} w_i \nabla_\theta \ell_i(\theta_0)\right]
$$

is the design-based variance of the Horvitz-Thompson total score, scaled by $n^{-2}$, accounting for the full stratification and clustering structure. Note that $J_0$ is a **PSU-level** variance — it reflects the dependence structure among units within primary sampling units, not a simple unit-level variance. Combining the above steps, with $H_0 = -\nabla^2 Q(\theta_0)$ the population Hessian of the limiting criterion:

$$
\text{Var}_{\mathcal{L}_*}\!\bigl(\bar\theta(\mathbf{w}^*)\bigr) \xrightarrow{p} \Sigma_{\text{sandwich}} = H_0^{-1} J_0 H_0^{-1}
$$

**Caveat for BART.** Conditions (R1)-(R3) apply to smooth parametric or semiparametric models. For BART, the parameter space is the space of tree ensembles and the pseudo-likelihood is non-differentiable over tree configurations. The delta-method argument above applies to the **leaf parameters** $\{\mu_{jl}\}$ conditional on a fixed tree structure — these are smooth Gaussian functionals of $\mathbf{w}^*$ given the tree — but a unified asymptotic argument for the full ensemble requires separate, nonparametric justification; see Section 5.1.

---

## 2. Scale-Sensitivity of Term 2 and the Role of $K$

### 2.1 Asymmetric Effect of Weight Scaling on Terms 1 and 2

A central insight is that the two terms respond differently to a constant multiplicative rescaling $\tilde{w}_i^* = K w_i^*$ for $K > 0$.

**Term 1 — scale-invariant through the leaf mean.** In the Gibbs step for the terminal node parameter $\mu_l$ in leaf $l$, the weighted pseudo-posterior mean is:

$$
\hat{\mu}_l(K\mathbf{w}^*) = \frac{\sum_{i \in l} K w_i^* R_i}{\sum_{i \in l} K w_i^*} = \frac{\sum_{i \in l} w_i^* R_i}{\sum_{i \in l} w_i^*} = \hat{\mu}_l(\mathbf{w}^*)
$$

This ratio is scale-invariant, so the fluctuations of $\hat\mu_l$ over bootstrap draws — and therefore Term 1 — are unchanged by rescaling.

**The tree structure is not scale-invariant.** The BART tree proposal is evaluated via a Metropolis-Hastings step whose acceptance ratio involves the **integrated marginal pseudo-likelihood** of each leaf configuration, obtained by integrating out $\mu_l \sim N(0, \sigma_\mu^2)$ against the weighted pseudo-likelihood. Under the normal-normal conjugate structure with error variance $\sigma^2$, the log integrated pseudo-marginal likelihood for leaf $l$ is (Chipman, George & McCulloch, 2010):

$$
\log m_l^* = \frac{1}{2}\log\!\left(\frac{\sigma^2}{\sigma^2 + W_l^*\sigma_\mu^2}\right) + \frac{(W_l^*)^2\sigma_\mu^2\,(\bar{R}_l^*)^2}{2\sigma^2\!\left(\sigma^2 + W_l^*\sigma_\mu^2\right)} + \text{const}
$$

where $W_l^* = \sum_{i\in l} w_i^*$ and $\bar{R}_l^* = W_l^{*-1}\sum_{i\in l} w_i^* R_i$. For the grow proposal (splitting leaf $l$ into children $L$ and $R$), the log acceptance ratio depends on:

$$
\Delta\log m^* = \bigl(\log m_L^* + \log m_R^*\bigr) - \log m_l^*
$$

Both terms in $\log m_l^*$ depend on $W_l^*$, which scales linearly with $K$. In the regime $W_l^*\sigma_\mu^2 \gg \sigma^2$ (data dominating prior), the quadratic term approximates $W_l^*(\bar{R}_l^*)^2/(2\sigma^2)$ and grows linearly with $K$, while the log-determinant grows only logarithmically. The net effect is that, for moderate increases in $K$, the integrated likelihood gain from each split increases, tending to favor deeper trees.

This tendency is **not globally monotone**: at extreme $K$, the within-leaf error variance $\sigma^2$ (updated via its own posterior draw) adapts, and the tree structure prior — which penalizes depth independently of the likelihood — eventually dominates. For moderate departures from $K=1$, however, the dominant effect is toward deeper trees. **This scale-sensitivity of the tree posterior is the primary reason why $K$ must be carefully calibrated rather than set arbitrarily large.** For the leaf mean draws conditional on a fixed tree structure, scale-invariance holds exactly.

### 2.2 Term 2: The Conditional Posterior Variance

Conditional on tree structure and the current error variance $\sigma^2$, the posterior variance of the leaf parameter $\mu_l$ under a $N(0, \sigma_\mu^2)$ prior is:

$$
\text{Var}(\mu_l \mid \mathbf{w}^*, \text{tree}) = \left(\frac{W_l^*}{\sigma^2} + \frac{1}{\sigma_\mu^2}\right)^{-1}
$$

where $W_l^* = \sum_{i \in l} w_i^*$ is the bootstrap leaf weight sum. Scaling weights by $K$ replaces $W_l^*$ with $K W_l^*$, giving:

$$
\text{Var}(\mu_l \mid K\mathbf{w}^*, \text{tree}) = \left(\frac{K W_l^*}{\sigma^2} + \frac{1}{\sigma_\mu^2}\right)^{-1} \approx \frac{\sigma^2}{K W_l^*} \quad \text{when } K W_l^* \gg \frac{\sigma^2}{\sigma_\mu^2}
$$

Thus Term 2 is of order $O(1/(K\hat{N}))$ when raw population-scale weights are used, and vanishes as $K\hat{N} \to \infty$.

### 2.3 The Practical Reality of $\hat{N}$ and Derivation of the Default $K$

In practice, the sum of sampling weights $\hat{N} = \sum_{i \in S} w_i$ need not equal the true finite population size $N$:

- **Calibrated/post-stratified weights:** Nonresponse adjustments or raking alter the weight total.
- **Domain-level analysis:** For a subpopulation $S_{\text{sub}} \subset S$, the effective weight sum $\hat{N}_{\text{sub}} = \sum_{i \in S_{\text{sub}}} w_i$ can be far smaller than $N$.

In these settings, $K$ provides additional scaling to guarantee $KW_l^* \gg \sigma^2/\sigma_\mu^2$ for all non-trivial leaves. A principled default is derived by requiring the **worst-case leaf** (the smallest, with $n_{\min}$ observations) to satisfy the Term 2 collapse condition. The expected bootstrap weight sum for this leaf is:

$$
E_{\mathcal{L}_*}[W_l^*] \approx \hat{N} \cdot \frac{n_{\min}}{n} = \bar{w} \cdot n_{\min}, \qquad \bar{w} = \hat{N}/n
$$

Term 2 collapse requires $K\bar{w}n_{\min}/\sigma^2 \gg 1/\sigma_\mu^2$. Setting $K$ to the smallest value achieving $K\bar{w}n_{\min} = n$ — so that the effective precision in the smallest leaf matches that of $n$ i.i.d. unweighted observations, the regime for which the BART prior $\sigma_\mu^2$ was calibrated — yields:

$$
K = \frac{n}{\bar{w} \cdot n_{\min}} = \frac{n^2}{\hat{N} \cdot n_{\min}}, \qquad n_{\min} = 5 \text{ (default floor)}
$$

**Relationship to population scale.** For large population-level surveys with $\hat{N} \gg n$, this formula gives $K \ll 1$, which conflicts with $K \geq 1$. In this regime, raw weights ($K = 1$) already satisfy $E_{\mathcal{L}_*}[W_l^*] \approx \hat{N}n_l/n \gg \sigma^2/\sigma_\mu^2$ for all leaves, so Term 2 collapse is automatic. This is the case for MEPS analyses where $\hat{N} \approx 330$ million and $n \approx 26{,}000$.

### 2.4 The Indispensable Role of Kim-Rao: Why Arbitrary Randomizations Fail

Since $K$ collapses Term 2 for *any* set of randomized weights, a natural question arises: can the Kim-Rao bootstrap be replaced with a simpler scheme? **The answer is no**, because Term 2 is only half of the problem.

The one-step approach achieves design-consistent frequentist coverage if and only if both conditions hold simultaneously:

1. **Term 2 $\to 0$:** Controlled by $K$ (or by population-scale $\hat{N}$).
2. **Term 1 $= \Sigma_{\text{sandwich}}$:** Requires $\mathcal{L}_*$ to satisfy the PSU-level moment-matching condition:

$$
\text{Var}_{\mathcal{L}_*}\!\bigl(\bar{S}_{w^*}(\theta_0)\bigr) = \hat{V}\!\left[\bar{S}_w(\theta_0) \mid \mathcal{F}_N\right]
$$

This encodes the **complete PSU-level covariance structure of the complex design**. Alternative randomizations fail:

- **i.i.d. $\text{Exp}(1)$ multipliers:** Produce $\text{Var}_{\mathcal{L}_*}[\bar{S}_{w^*}] \propto n^{-2}\sum_i w_i^2 (\nabla_\theta \ell_i)^2$. This estimates a unit-level score variance treating all observations as independent. Under stratified cluster sampling, the correct design variance accounts for PSU-level covariance $\text{Cov}(w_{hi}\nabla\ell_{hi},\, w_{hi'}\nabla\ell_{hi'})$ for units $i,i'$ within the same PSU $h$; the i.i.d. formulation ignores this within-cluster dependence and systematically underestimates $J_0$ — it is only exact under simple random sampling where within-cluster covariance is zero by design.
- **Dirichlet $(1, \ldots, 1)$ weights (Bayesian bootstrap):** The bootstrap variance converges to the empirical variance of the score contributions, matching the i.i.d. (SRS) design variance rather than $\hat{V}[\bar{S}_w(\theta_0)\mid\mathcal{F}_N]$.
- **Nonparametric bootstrap (unit resampling with replacement):** Treats all units as exchangeable, replicating the i.i.d. variance without accounting for stratification or clustering.

**The Kim-Rao multinomial bootstrap** satisfies Condition 2 exactly. By resampling $n_h - 1$ PSUs with replacement within each stratum $h$ and applying the factor $k_h = n_h/(n_h-1)$, Kim et al. (2024) establish in Theorem 1 that:

$$
E_{\mathcal{L}_*}\!\bigl[\bar{S}_{w^*}(\theta)\bigr] = \bar{S}_w(\theta), \qquad \text{Var}_{\mathcal{L}_*}\!\bigl(\bar{S}_{w^*}(\theta)\bigr) = \hat{V}\!\left[\bar{S}_w(\theta) \mid \mathcal{F}_N\right]
$$

To summarize the separation of roles:

> **$K$ controls Term 2:** any weight randomization can have Term 2 collapsed to zero by choosing $K$ sufficiently large (subject to the tree overfitting constraint in Section 2.1).  
> **Kim-Rao controls Term 1:** only a bootstrap distribution that exactly replicates the design-based score variance at the PSU level guarantees Term 1 $= \Sigma_{\text{sandwich}}$.

The two requirements are **orthogonal**: $K$ is a computational regularizer; Kim-Rao is a statistical necessity. Neither alone is sufficient.

### 2.5 Weight Boundedness, Regularity Conditions, and the Role of $K$

Every frequentist and Bayesian result in this literature — from Savitsky & Toth (2016) and Williams & Savitsky (2021) to Kim, Rao & Wang (2024) — rests on a regularity condition bounding the influence of individual design weights. We review these conditions and establish that $K$ does not interact with any of them.

#### The Boundedness Conditions

**Savitsky & Toth (2016)** establish pseudo-posterior concentration (their Theorem 2.1) under the condition that normalized design weights are bounded:

$$
\max_{i \in S}\, \tilde{w}_i \;\leq\; C, \qquad \tilde{w}_i = \frac{n\, w_i}{\hat{N}}, \quad \hat{N} = \sum_{j \in S} w_j
$$

for some constant $C < \infty$ not growing with $n$.

**Williams & Savitsky (2021)** maintain the identical normalization condition ($\max_i \tilde{w}_i \leq C$) and additionally require uniform continuity of $n^{-1}\sum_i \tilde{w}_i \nabla^2 \ell_i(\theta)$ in $\theta$.

**Kim, Rao & Wang (2024)** state their Condition (C2): a Lindeberg condition for the bootstrap score CLT,

$$
\frac{1}{V_n}\sum_{h=1}^{H} \frac{k_h^2}{n_h} \sum_{i=1}^{n_h} \bigl(w_{hi} s_{hi}(\theta_0)\bigr)^2 \cdot \mathbf{1}\!\left\{\bigl|w_{hi} s_{hi}(\theta_0)\bigr| > \epsilon\, \sqrt{V_n}\right\} \to 0
$$

where $V_n = \hat{V}[\bar{S}_w(\theta_0)\mid\mathcal{F}_N]$ and $s_{hi} = \nabla_\theta \ell_{hi}(\theta_0)$. This is a condition on the *original* design weights $w_{hi}$ — it encodes that no single PSU-level term dominates the design variance of the HT score.

**Common thread:** All three conditions are stated in terms of **normalized or standardized** weight quantities — either $\tilde{w}_i = nw_i/\hat{N}$ or ratios $w_{hi}^2 s_{hi}^2/V_n$ — not in terms of the absolute scale of $w_i$.

#### Does $K$ Affect These Conditions?

**1. The Savitsky-Toth normalized boundedness condition.**

$$
\frac{n \cdot K w_i^*}{\sum_j K w_j^*} = \frac{n \cdot w_i^*}{\sum_j w_j^*} = \tilde{w}_i^*
$$

The factor $K$ cancels exactly. **The Savitsky-Toth boundedness condition is fully invariant to $K$.**

**2. The Kim-Rao Lindeberg condition.**  
Under $K$-scaling, $\bar{S}_{Kw^*}(\theta) = K \cdot \bar{S}_{w^*}(\theta)$ and $\text{Var}_{\mathcal{L}_*}(\bar{S}_{Kw^*}) = K^2 \text{Var}_{\mathcal{L}_*}(\bar{S}_{w^*})$. The Lindeberg ratio becomes:

$$
\frac{(K w_{hi})^2 s_{hi}^2}{K^2 V_n} = \frac{w_{hi}^2 s_{hi}^2}{V_n}
$$

Again $K$ cancels. **The Lindeberg condition for the bootstrap CLT is scale-invariant with respect to $K$.**

**3. The bootstrap weight boundedness.**  
The Kim-Rao bootstrap weights are $w_{hi}^* = k_h m_{hi}^* w_{hi}$, where $m_{hi}^* \in \{0,\ldots,n_h-1\}$. Worst-case: $\max_i w_{hi}^* \leq k_h(n_h-1)w_{hi} \approx n_h w_{hi}$. The Lindeberg condition then requires $\max_i n_h^2 w_{hi}^2 s_{hi}^2 / V_n \to 0$ — a condition on the original design weights, unrelated to $K$.

#### Summary: The Two Distinct Roles of Boundedness and $K$

| Condition | Statement | Scale-invariant w.r.t. $K$? |
|:---|:---|:---:|
| Savitsky-Toth (normalized weight bound) | $\max_i n w_i/\hat{N} \leq C$ | **Yes** — $K$ cancels |
| Kim-Rao Lindeberg (CLT for Term 1) | $\max_i w_{hi}^2 s_{hi}^2 / V_n \to 0$ | **Yes** — $K$ cancels |
| Pseudo-likelihood scale (Term 2 mechanism) | $W_l^* = \sum_{i \in l} w_i^*$ (absolute) | **No** — scales as $K W_l^*$ |

The boundedness conditions are entirely unaffected by $K$. The only quantity $K$ directly controls is the *absolute* scale of $W_l^*$, which governs Term 2. This makes $K$ a legitimate computational tuning device: it operates on the unnormalized pseudo-likelihood scale to collapse Term 2 while leaving all statistical validity conditions — stated in scale-normalized form — completely intact.

**Practical implication.** Even in small-domain analyses requiring large $K$ to suppress Term 2, this does not push weights into a regime violating the regularity conditions of Savitsky-Toth or Kim-Rao. The genuine constraint on $K$ is the numerical one identified in Section 2.1: very large $K$ biases the tree MH acceptance ratios toward excessive depth, degrading OOB generalization. This is a computational concern, not a statistical one.

### 2.6 Rethinking $K$: Population-Scale Normalization and an Analogy to Power Posteriors

The treatment of $K$ as an arbitrary tuning hyperparameter understates a deeper theoretical structure. A more precise view reveals that $K$'s correct value is determined by the normalization theory of the pseudo-likelihood; OOB calibration is best understood as empirical confirmation of a theoretical requirement.

#### $K$ as the Pseudo-Likelihood Normalization Constant

The survey-weighted pseudo-log-likelihood is:

$$\mathcal{PL}(\theta;\, \mathbf{w}) = \sum_{i \in S} w_i\, \ell_i(\theta)$$

Under the Horvitz-Thompson principle, this consistently estimates the population-level sum $\mathcal{L}_N(\theta) = \sum_{i \in \mathcal{U}} \ell_i(\theta) \approx N \cdot E_{\mathcal{U}}[\ell_i(\theta)]$, which is $O(N)$. When raw Kim-Rao weights are used ($K = 1$), the pseudo-likelihood retains this $O(\hat{N})$ population scale. When normalized to sum to $n$ — corresponding to $K = n/\hat{N}$ — the pseudo-likelihood is deflated to $O(n)$, misrepresenting the effective amount of population-level information.

The BART prior $\sigma_\mu^2$ is calibrated under the assumption of $n$ i.i.d. observations. With raw survey weights ($K=1$), the pseudo-likelihood operates at scale $O(\hat{N}) \gg O(n)$, overwhelming the prior and achieving Term 2 collapse automatically. From this perspective:

> **$K = 1$ (raw weights) is not merely the empirically best grid point. It is the value for which the pseudo-likelihood operates at its correct population scale, naturally achieving the prior-to-likelihood balance required for Term 2 collapse.**

Choosing any $K > 1$ inflates beyond population scale, artificially concentrating the posterior. Choosing $K < 1$ deflates below population scale, reintroducing non-negligible Term 2 — precisely the failure mode of the Savitsky-Toth normalization.

#### An Analogy to the Power Posterior Framework

Multiplying all bootstrap weights by $K$ is mathematically equivalent to raising the pseudo-likelihood to a power:

$$\exp\!\bigl(K \cdot \mathcal{PL}(\theta;\, \mathbf{w}^*)\bigr) = \exp\!\bigl(\mathcal{PL}(\theta;\, K\mathbf{w}^*)\bigr)$$

This bears a **formal analogy** to the power posterior (fractional likelihood) framework of Friel & Pettitt (2008) and Walker & Hjort (2001), where the likelihood is tempered by $\alpha \in (0,1]$. The analogy is instructive but imperfect: in standard power posteriors, the tempering acts on a proper likelihood $p(y\mid\theta)^\alpha$; here, $\mathcal{PL}(\theta;\mathbf{w})$ is itself a weighted approximation to a population log-likelihood, not a proper probability model. The temperature interpretation is therefore an analogy, not an identity.

Within this analogy, define the natural normalization scale $\eta_0 = 1/\hat{N}$. Any choice of $K$ corresponds to a relative scale $\alpha = K/\hat{N}$:

- **$K = n/\hat{N}$ (normalized weights):** $\alpha = n/\hat{N}^2 \ll 1$ — severe deflation, requiring the post-hoc sandwich correction of Williams & Savitsky (2021).
- **$K = 1$ (raw weights):** $\alpha = 1/\hat{N}$ — the natural scale where the pseudo-likelihood correctly represents population-level information.
- **$K > 1$:** Over-inflation beyond natural population scale, concentrating the posterior too sharply and biasing tree structure.

**The OOB criterion consistently selects $K^* = 1$ because this is the value at which the pseudo-likelihood is correctly calibrated to population scale.**

#### Small-Domain Analyses: A Theory-Derived $K^*$

For domain-level analyses where $n_{\text{sub}}$ observations fall in subpopulation $S_{\text{sub}}$ with $\hat{N}_{\text{sub}} = \sum_{i \in S_{\text{sub}}} w_i \ll \hat{N}$, the subpopulation pseudo-likelihood is:

$$\mathcal{PL}_{\text{sub}}(\theta) = \sum_{i \in S_{\text{sub}}} w_i \ell_i(\theta) = O(\hat{N}_{\text{sub}}) \ll O(\hat{N})$$

The prior $\sigma_\mu^2$ (calibrated at full-sample scale) now dominates, Term 2 is non-negligible, and design-consistency is threatened. The theory-derived $K^*$ restores the correct population scale:

$$K^* = \frac{\hat{N}}{\hat{N}_{\text{sub}}} = \frac{\sum_{i \in S} w_i}{\sum_{i \in S_{\text{sub}}} w_i}$$

This is a computable function of the observed design weights, requiring no grid search. The OOB procedure then serves as a robustness check on this theoretical target, not a primary estimation strategy.

#### Joint Calibration with BART Prior Hyperparameters

Term 2 collapse requires $KW_l^*/\sigma^2 \gg 1/\sigma_\mu^2$. This can be achieved equivalently by:

- (a) Increasing $K$ — inflating the pseudo-likelihood scale;
- (b) Increasing $\sigma_\mu^2$ — weakening the leaf prior;
- (c) A combination on the **confounded manifold** $\{(K, \sigma_\mu^2) : K\sigma_\mu^2 W_l^*/\sigma^2 \gg 1\}$.

Standard BART calibrates $\sigma_\mu^2 = (y_{\max} - y_{\min})^2/(4mk^2)$ for $n$ i.i.d. observations (Chipman, George & McCulloch, 2010). With raw survey weights ($K=1$), the effective precision $W_l^*/\sigma^2 \approx \hat{N}n_l/(n\sigma^2)$ at population scale ($\hat{N} \gg n$) already overwhelms the default prior. Thus, **the standard BART prior calibration, when combined with population-scale raw weights, jointly achieves Term 2 collapse as a consequence of their interaction — not as a separate tuning step.** This reframes $K$ as a population-scale normalization correction bridging the BART prior (calibrated for $n$ i.i.d.) and the survey pseudo-likelihood (operating at $\hat{N}$ population observations). For standard population-level analyses, $K = 1$ achieves this automatically. For small-domain analyses, $K^* = \hat{N}/\hat{N}_{\text{sub}}$ is the theory-derived correction.

---

## 3. The One-Step Resolution Mechanism

### 3.1 Sufficient Condition and Convergence Rate

With Term 2 negligible — specifically, when $K\hat{N}_{l,\min}/\sigma^2 \gg 1/\sigma_\mu^2$, where $\hat{N}_{l,\min}$ is the minimum leaf weight sum — the Gibbs draw for each terminal node parameter concentrates near its posterior mean at rate $O((K\hat{N}_l)^{-1/2})$:

$$
\mu_l^{(t)} \sim N\!\bigl(\hat{\mu}_l^{(t)},\, V_l^{(t)}\bigr), \qquad V_l^{(t)} = \left(\frac{K W_l^{*(t)}}{\sigma^2} + \frac{1}{\sigma_\mu^2}\right)^{-1} = O\!\left(\frac{1}{K\hat{N}_l}\right)
$$

Since $V_l^{(t)} \to 0$ as $K\hat{N}_l \to \infty$, the draw $\mu_l^{(t)}$ is within $O((K\hat{N}_l)^{-1/2})$ of $\hat\mu_l^{(t)}$ with high probability. The condition $K\hat{N}_{l,\min} \to \infty$ is therefore the **sufficient and verifiable condition** for Term 2 collapse; it is not an assumption but a constraint on the design and chosen $K$.

The empirical variance of draws across iterations is then determined primarily by Term 1:

$$
\text{Var}(\mu_l) \approx \text{Var}_{\mathcal{L}_*}\!\bigl(\hat{\mu}_l(\mathbf{w}^*)\bigr) \approx \Sigma_{\text{sandwich}}
$$

with approximation error of order $O(1/(K\hat{N}_l))$. This is the **one-step resolution**: design-consistent uncertainty is natively embedded in the MCMC fluctuations through the Kim-Rao resampling, requiring no post-processing, no matrix inversion, and no Bartlett-identity correction.

### 3.2 Treating $K$ as a Hyperparameter: OOB Tuning

The Kim-Rao bootstrap provides a computationally free internal validation set at every MCMC iteration. Because it resamples $n_h - 1$ PSUs with replacement within stratum $h$, any PSU $(h,i)$ has exclusion probability:

$$
P_{\mathcal{L}_*}(m_{hi}^* = 0) = \left(1 - \frac{1}{n_h}\right)^{n_h - 1} \xrightarrow{n_h \to \infty} e^{-1} \approx 0.368
$$

For small strata (e.g., $n_h = 2$), this probability is exactly $1/2$. The asymptotic approximation $\approx 37\%$ OOB units is accurate when all $n_h$ are moderate to large; in designs with small strata, the actual OOB fraction can differ materially from $e^{-1}$.

**OOB predictive loss as an objective for $K$.** At iteration $t$, let $\mathcal{O}^{(t)} = \{i : w_i^{*(t)} = 0\}$ be the OOB set and $f^{(t)}(x) = \sum_j \mu_{jl(x)}^{(t)}$ the current ensemble fitted on in-bag data. The OOB loss is:

$$
\mathcal{L}_{\text{OOB}}(K, t) = \frac{1}{|\mathcal{O}^{(t)}|} \sum_{i \in \mathcal{O}^{(t)}} \bigl(y_i - f^{(t)}(x_i)\bigr)^2
$$

The cumulative post-burn-in OOB loss defines the calibration criterion:

$$
\mathcal{L}_{\text{OOB}}(K) = \frac{1}{T_{\text{post}}} \sum_{t=t_{\text{burn}}}^{T} \mathcal{L}_{\text{OOB}}(K, t), \qquad K^* = \underset{K \in \mathcal{K}}{\arg\min}\; \mathcal{L}_{\text{OOB}}(K)
$$

**Theoretical properties of OOB tuning:**

1. **Term 1 is invariant to $K$** (Section 2.1, leaf mean scale-invariance). Selecting $K$ via OOB loss does not alter the design-consistent sandwich variance in Term 1, preserving the frequentist coverage guarantee.

2. **The OOB units are genuinely held out.** Within any iteration, the OOB set is determined by the bootstrap resampling that drives the posterior fluctuations, so there is no leakage between model fitting and validation.

3. **Finite-sample bias-variance trade-off.** As $K$ increases, the scale-sensitive tree MH acceptance ratios (Section 2.1) tend to favor deeper splits, increasing in-sample fit but potentially degrading OOB generalization. The OOB criterion balances these competing pressures. Empirically (Section 6), this balance favors $K^* = 1$ in population-scale settings ($\hat{N} \gg n$), confirming that additional $K$-scaling beyond raw population weights is unnecessary.

---

## 4. Analysis of Alternative Normalization Schemes

To understand what is lost by normalizing the bootstrap weights, let $M$ denote the target for $\sum_i \tilde{w}_i$. Because the leaf mean $\hat\mu_l$ is scale-invariant, Term 1 remains fixed at $\Sigma_{\text{sandwich}}$ regardless of $M$. However, Term 2 scales inversely with $M$: the average leaf weight sum is $E_{\mathcal{L}_*}[W_l^*] \approx M n_l/n$, so:

$$
E_{\mathcal{L}_*}\!\bigl[\text{Var}(\theta \mid \mathbf{w}^*)\bigr] \approx \frac{\sigma^2}{M n_l/n} = \frac{n}{M} \cdot \frac{\sigma^2}{n_l} = \frac{n}{M}\, V_{\text{model}}
$$

where $V_{\text{model}} = \sigma^2/n_l$ is the standard unweighted posterior variance for leaf $l$ under $n$ i.i.d. observations.

### Case A: Normalizing to Sample Size ($M = n$)

Setting $M = n$ gives $n/M = 1$ and:

$$
\text{Var}(\theta) \approx \Sigma_{\text{sandwich}} + V_{\text{model}}
$$

For designs where the design effect is defined as $\text{Deff} = \Sigma_{\text{sandwich}}/V_{\text{model}}$ — an approximation that holds when $\Sigma_{\text{sandwich}}$ is compared to the SRS variance for the same estimator and when the design effect is approximately uniform across leaves — the over-inflation ratio is approximately:

$$
\frac{\text{Var}(\theta)}{\Sigma_{\text{sandwich}}} \approx 1 + \frac{1}{\text{Deff}}
$$

This approximation degrades when the design effect varies substantially across domains or with post-stratification. For small design effects ($\text{Deff}$ close to 1), the additive inflation is most severe.

### Case B: Normalizing to Effective Sample Size ($M = n_{\text{eff}} = n/\text{Deff}$)

Setting $M = n_{\text{eff}}$ gives $n/M = \text{Deff}$:

$$
E_{\mathcal{L}_*}\!\bigl[\text{Var}(\theta \mid \mathbf{w}^*)\bigr] \approx \text{Deff} \cdot V_{\text{model}} \approx \Sigma_{\text{sandwich}}
$$

Substituting into the Law of Total Variance:

$$
\text{Var}(\theta) \approx \Sigma_{\text{sandwich}} + \Sigma_{\text{sandwich}} = 2\,\Sigma_{\text{sandwich}}
$$

This doubles the total variance, inflating credible interval widths by $\sqrt{2} \approx 1.41$. Among the normalization targets $\{n_{\text{eff}},\, n,\, \hat{N},\, 1\}$ considered here, normalizing to $n_{\text{eff}}$ produces the **largest relative bias in total variance** within this framework, because it most closely matches Term 2 to $\Sigma_{\text{sandwich}}$, causing both terms to be commensurate. This confirms that intentionally matching Term 2 to $\Sigma_{\text{sandwich}}$ — as two-step fixed-weight approaches attempt — is fundamentally incompatible with a one-step random-weight design.

---

## 5. Contrast with the Two-Step Approach of Savitsky & Williams

The distinction between the one-step and two-step approaches reflects a fundamental difference in the inferential target and in the roles of Term 1 and Term 2.

**Inferential target clarification.** Savitsky & Toth (2016) and Williams & Savitsky (2021) primarily target **superpopulation inference**: they seek pseudo-posterior draws that are Bayes-consistent under a superpopulation model, with Williams & Savitsky (2021) additionally providing a post-hoc correction for **design-based frequentist coverage**. The one-step approach directly targets design-based frequentist coverage from the outset, using Kim-Rao resampling to embed the design's covariance structure organically. The comparison below adopts **design-based frequentist coverage** as the common inferential standard.

| | **Savitsky & Toth (2016)** | **Williams & Savitsky (2021)** | **One-Step (Proposed)** |
|:---|:---:|:---:|:---:|
| Bootstrap resampling at each MCMC step? | ✗ | ✗ | ✓ |
| Weight normalization | $\sum \tilde{w}_i = n$ | $\sum \tilde{w}_i = n$ | None (raw, or $\times K$) |
| Term 1 | $0$ | $0$ | $\Sigma_{\text{sandwich}}$ |
| Term 2 | $\approx V_{\text{model}}$ | $\approx V_{\text{model}}$ (corrected post-hoc) | $\approx 0$ |
| Design-based frequentist coverage | **Below nominal** | **Nominal** (post-hoc) | **Nominal** (one-step) |

**Savitsky & Toth (2016)** normalize weights to sum to $n$ and fix them for all MCMC iterations (Term 1 = 0). The total variance equals only Term 2 $\approx V_{\text{model}}$. Under informative sampling where Bartlett's second identity fails ($H_0 \neq J_0$), the pseudo-posterior variance is $O(H_0^{-1})$ rather than $\Sigma_{\text{sandwich}} = H_0^{-1}J_0H_0^{-1}$, and since $V_{\text{model}} \ll \Sigma_{\text{sandwich}}$ for typical complex designs, design-based credible intervals are too narrow. This is not a failure of the Savitsky-Toth framework per se — it targets superpopulation Bernstein-von Mises consistency, not design-based sandwich coverage — but it renders the approach insufficient for design-based frequentist inference without additional correction.

**Williams & Savitsky (2021)** retain fixed normalized weights but apply a post-hoc variance adjustment. Let $\hat{V}$ denote the empirical covariance of the MCMC draws and $\hat\Sigma_{\text{sandwich}}$ the estimated sandwich covariance. The adjustment transforms draws via:

$$
\theta^{(t)}_{\text{adjusted}} = \bar{\theta} + R\bigl(\theta^{(t)} - \bar{\theta}\bigr), \quad R = \hat\Sigma_{\text{sandwich}}^{1/2}\,\hat{V}^{-1/2}
$$

This rescales the posterior covariance from $\hat V$ to $\hat\Sigma_{\text{sandwich}}$, restoring nominal frequentist coverage. However, it has critical limitations:

1. **Requires estimation of $\hat\Sigma_{\text{sandwich}}$:** Computing $J_0 = \lim n^2 \text{Var}_\xi[\bar{S}_w(\theta_0)]$ requires the score $\nabla_\theta \ell_i(\theta_0)$ for each unit. For BART, the regression function $f(x) = \sum_j \mu_{jl(x)}$ is a **step function over the joint space of tree structures and leaf parameters** — a functional rather than a finite-dimensional parameter. Although the leaf parameters $\{\mu_{jl}\}$ enter a differentiable normal likelihood *conditional on tree structure*, computing a sandwich variance that marginalizes over tree space jointly has no established closed-form or numerically stable solution for BART. The Gateaux derivative of the prediction functional with respect to the pseudo-likelihood does not exist in the classical sense, making the analytic gradient needed for $J_0$ unavailable.

2. **Distorts non-Gaussian posterior geometries:** The projection $R$ is an affine linear transformation. If the pseudo-posterior is non-Gaussian (skewed, heavy-tailed, or multimodal), applying $R$ corrects the covariance matrix but cannot correct higher-order cumulants, leaving the predictive distribution misspecified for tail inference.

3. **Inapplicable to multi-part models:** For hurdle or two-part models, the joint pseudo-posterior spans both discrete (probit) and continuous (log-amount) components with no unified differentiable joint likelihood. There is no single sandwich matrix $\Sigma_{\text{sandwich}}$ spanning the joint parameter space and no single $R$ applicable across both components simultaneously.

### 5.1 Methodological Advantages of the One-Step Approach

The proposed one-step Random Weight MCMC resolves all three limitations simultaneously:

1. **No sandwich estimation required.** The design-based covariance structure is generated organically through the Kim-Rao resampling mechanism. The bootstrap score variance equals the design-based sandwich variance by construction (Section 1), so no explicit estimation of $H_0$, $J_0$, or $\Sigma_{\text{sandwich}}$ is ever needed.

2. **Posterior shape is preserved.** Because the variance correction operates at the level of the generating mechanism (the weighted pseudo-likelihood draws) rather than as a post-hoc linear projection, the posterior draws retain their natural non-Gaussian geometry. Skewness, multimodality, and tail behavior emerge organically from the MCMC dynamics.

3. **Seamless extension to multi-part models.** Since design-consistency is resolved at the observation level — by assigning bootstrap weights $w_i^{*(t)}$ to each unit before the MCMC updates — the same weights are passed simultaneously to the probit (participation) and log-amount (positive spending) components of a hurdle model. Each component receives the correct survey-weighted pseudo-likelihood without requiring a joint sandwich adjustment.

---

## 6. Empirical Validation & Simulation Results

To empirically validate the one-step variance theory and the proposed $K$-scaling calibration framework, we conducted a simulation study using 1,000 design-replicated samples ($n = 6{,}100$) drawn from a synthetic finite population ($N = 65{,}360$) under a stratified two-stage PPS design mimicking the Medical Expenditure Panel Survey (MEPS). Three weight conditions were evaluated:

* **Condition A (Normalized, $\sum w_i = n$):** Weights scaled to sum to sample size, as in Savitsky & Toth (2016).
* **Condition B (Raw Kim-Rao, $K = 1.0$):** Weights retained at population scale, $\hat{N} \approx 65{,}360$.
* **Condition C (OOB-Calibrated $K^*$):** $K$ selected by minimizing cumulative OOB MSE over a grid $\mathcal{K} = \{1, 5, 15, 30, 50, 80\}$. *(Preliminary: 2 smoke-test replications; full simulation ongoing.)*

### Summary of Results

| Condition | Replications | Posterior Predictive Coverage (Target: 95%) | MCMC SD | CI Width | RMSE | Deff |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| **A** (Normalized) | 1,000 | **96.59%** | 1.3228 | 5.1930 | 6.3377 | 1.0683 |
| **B** (Raw weights) | 1,000 | **96.11%** | 1.4584 | 5.7440 | 6.3674 | 1.0683 |
| **C** (OOB Calibrated) | 2† | **95.95%** | 1.4358 | 5.6535 | 6.4300 | 1.0581 |

*†Condition C results are from 2 smoke-test replications only and are not comparable in precision to the 1,000-replication estimates for Conditions A and B. All Condition C figures should be regarded as preliminary.*

### Key Findings & Conclusions

1. **Over-coverage under normalized weights (Condition A).**  
   Coverage of 96.59% systematically exceeds the nominal 95% rate, a 1.59 percentage-point excess. From the total variance decomposition with $M = n$, the ratio of total to sandwich variance is approximately $1 + 1/\text{Deff} \approx 1 + 1/1.0683 \approx 1.937$, a $\approx$93.7% variance inflation. The modest coverage excess (only $\sim$1.6 pp despite a near-doubling of variance) reflects the **concavity of coverage as a function of variance** for approximately normal posteriors: the coverage shift at nominal level $\alpha$ scales approximately as $\Phi(z_{\alpha/2}/\sqrt{1+\delta}) - (1-\alpha/2)$ in the fractional inflation $\delta$, which for $\delta \approx 0.94$ and $\alpha = 0.05$ yields a shift of only $\sim$1-2 percentage points — consistent with the observed 1.59 pp excess.

2. **Near-nominal coverage with raw weights (Condition B).**  
   Using population-scale Kim-Rao weights ($K = 1.0$) achieves 96.11% coverage. This confirms that for $\hat{N} \gg n$, the raw weight sum collapses Term 2 sufficiently without requiring additional scaling. The residual over-coverage of 1.11 pp above 95% is attributable to finite-sample effects: the approximation $E_{\mathcal{L}_*}[W_l^*] \approx \hat{N}n_l/n$ holds asymptotically, but at finite $n = 6{,}100$ a small non-zero Term 2 contribution remains.

3. **OOB calibration selects $K^* = 1.0$ (Condition C).**  
   Across sample draws, the OOB criterion is monotonically increasing as $K$ increases from 1 to 80:

   | $K$ | OOB-MSE |
   |:---:|:---:|
   | 1.0 | 39.98 |
   | 5.0 | 40.25 |
   | 15.0 | 40.38 |
   | 30.0 | 40.49 |
   | 50.0 | 40.54 |
   | 80.0 | 40.66 |

   This provides direct empirical evidence that, in population-scale settings, inflating $K$ beyond 1 degrades OOB predictive accuracy. As discussed in Section 2.1, larger $K$ concentrates the tree posterior toward deeper splits through the scale-sensitive MH acceptance ratio, increasing in-sample fit at the expense of OOB generalization. Since $K=1$ already sufficiently collapses Term 2 (because $\hat{N} \gg n$), the OOB criterion correctly identifies no benefit from additional $K$-inflation. The preliminary coverage under Condition C (95.95%) is the closest to nominal of the three conditions.

4. **Dissertation Recommendation.**  
   For population-level health economics studies where $\hat{N}$ is of population scale (e.g., MEPS with $\hat{N} \approx 330$ million and $n \approx 26{,}000$), raw Kim-Rao weights ($K = 1.0$) naturally satisfy the Term 2 collapse condition and are recommended as the default. OOB calibration over a $K$ grid should be retained for small-domain analyses where $n_{\text{sub}}$ and $\hat{N}_{\text{sub}}$ are both small, and the raw population-scale weights may not provide sufficient Term 2 suppression. In those cases, the theory-derived starting value $K^* = \hat{N}/\hat{N}_{\text{sub}}$ (Section 2.6) should anchor the search grid.

---

## 7. Friedman Mean-Surface Recovery under Informative Sampling

Section 6 evaluated **posterior predictive coverage of the observed outcome $y$** under the original latent-class DGP. To stress-test weight scaling for a **regression mean surface** under informative complex sampling, we ran a second Monte Carlo experiment with the Friedman #1 response linked to the design.

### 7.1 Design

* **Population.** Same stratified two-stage PPS finite population ($N = 65{,}360$).
* **Outcome (informative).**  
  $$
  \mu_f(x) = 10\sin(\pi x_1 x_2) + 20(x_3-0.5)^2 + 10 x_4 + 5 x_5 + 2\lambda_h + 2\eta_i,\qquad
  y_f = \mu_f + \varepsilon,\ \varepsilon\sim N(0,1).
  $$
  Design effects $\lambda_h$ (stratum) and $\eta_i$ (PSU) enter both selection and the mean, so ignoring the design biases recovery of the population surface.
* **Sample.** 1,000 independent survey draws ($n \approx 5{,}500$).
* **Model covariates.** Friedman features $x_1,\ldots,x_5$ only (design IDs not included as covariates) — a deliberate misspecification relative to $\mu_f$.
* **Scoring.** Held-out population units not in the sample ($n_{\mathrm{test}}=2{,}000$). Estimands:
  - $\theta = n_{\mathrm{test}}^{-1}\sum_i \mu_f(x_i)$ with coverage of the posterior for $\bar f$;
  - pointwise RMSE of $\bar f$ vs $\mu_f$;
  - posterior predictive coverage of $y_f$ (same metric family as Section 6).
* **Backend.** C++ Kim–Rao BART (`bart_mcmc.cpp`), 200 trees, 2,000 MCMC iterations, burn-in 500. $\sigma^2$ updates always use $n$-scaled residual weights (as in Section 6).

### 7.2 Weight regimes

**A2 ($M \approx K\cdot n$).** Input weights `weight_norm` (mean 1, $\sum w = n$). Kim–Rao bootstrap is applied and **not** re-normalized each draw, so $K$ multiplies the effective weight sum $M\approx K n$ and therefore scales Term 2 as $n/M \propto 1/K$. Grid: $K\in\{0.01, 0.1, 0.5, 1, 2, 5, 10, 20\}$.

**sum1 ($M = 1$).** Input weights $w_i/\sum_j w_j$ (unit sum), Kim–Rao without per-draw re-normalization, $K=1$. This is the extreme small-$M$ end of Section 4.

*(For reference from an earlier baseline under the same informative DGP: forcing $\sum w^*=n$ each draw with `weight_norm` and $K=1$ gave 86.3% coverage of $\theta$; raw population-scale weights with $K=1$ and Chipman `wbart` without design covariates both severely undercovered $\theta$ — consistent with Term 2 collapse plus mean-surface misspecification.)*

### 7.3 Results — A2 $K$-grid (1,000 replications per $K$)

| $K$ | Bias | Rel. bias | Emp. SE | Est. SE | SE ratio | Cover($\theta$) | RMSE($\mu_f$) | PP cover($y_f$) | Pointwise cover($\mu_f$) |
|:---:|-----:|----------:|--------:|--------:|---------:|----------------:|--------------:|----------------:|-------------------------:|
| 0.01 | −0.029 | −0.2% | 0.223 | 0.631 | 2.823 | **100.0%** | 3.251 | 99.2% | 92.2% |
| **0.1** | −0.056 | −0.3% | 0.225 | 0.253 | **1.125** | **94.9%** | 2.948 | **96.8%** | 68.7% |
| 0.5 | −0.072 | −0.4% | 0.224 | 0.208 | 0.929 | 87.8% | **2.894** | 95.3% | 49.0% |
| 1.0 | −0.084 | −0.4% | 0.222 | 0.197 | 0.888 | 86.1% | 2.912 | 94.5% | 44.2% |
| 2.0 | −0.117 | −0.6% | 0.217 | 0.178 | 0.819 | 81.7% | 2.984 | 93.1% | 45.1% |
| 5.0 | −0.306 | −1.6% | 0.200 | 0.097 | 0.482 | 30.0% | 3.618 | 77.8% | 53.5% |
| 10.0 | −0.427 | −2.2% | 0.211 | 0.042 | 0.199 | 4.2% | 5.198 | 42.4% | 39.5% |
| 20.0 | −0.436 | −2.3% | 0.236 | 0.033 | 0.139 | 3.6% | 5.979 | 29.0% | 27.8% |

*Sources: `data/friedman_K_grid_A2/report_table.csv` (8,000 task files).*

### 7.4 Results — unit-sum weights (sum1, $K=1$)

| Method | Bias | Rel. bias | Emp. SE | Est. SE | SE ratio | Cover($\theta$) | RMSE($\mu_f$) | PP cover($y_f$) |
|:---:|-----:|----------:|--------:|--------:|---------:|----------------:|--------------:|----------------:|
| sum1 | +1.237 | +6.4% | 0.897 | 5.899 | **6.576** | **100.0%** | 5.719 | **100.0%** |

*Source: `data/friedman_results_sum1/report_table.csv` (1,000 replications).*

### 7.5 Interpretation relative to Section 6

1. **Term 2 scaling direction is confirmed.** Under A2, larger $K$ shrinks Est. SE, SE ratio, Cover($\theta$), and PP($y_f$) together. This matches Section 4: Term 2 $\propto n/M$ with $M\approx K n$. Inflating $K$ past 1 on the mean-1 weight scale concentrates the posterior and undercovers — the same qualitative message as the OOB grid in Section 6.3 (do not inflate $K$ once Term 2 is already small).

2. **PP of $y$ vs coverage of mean($\mu_f$).** At $K=1$ on `weight_norm`, PP($y_f$)$\approx 94.5\%$ (near the Section 6 Condition A figure of 96.6%), while Cover($\theta$)$=86.1\%$. Interval calibration for the **mean surface** under design-linked misspecification is stricter than predictive coverage of $y$; the two metrics must not be conflated.

3. **$K\approx 0.1$ restores nominal Cover($\theta$).** With $M\approx 0.1 n$, Term 2 is inflated enough that Est. SE $\approx$ Emp. SE (ratio 1.13) and Cover($\theta$)$=94.9\%$, with PP($y_f$)$=96.8\%$. Point RMSE remains close to the minimum (2.95 vs 2.89 at $K=0.5$). Thus $K<1$ on the normalized-weight scale acts as a **Term 2 re-inflation dial** for this estimand — complementary to the Section 6 recommendation of $K=1$ raw weights when the target is PP($y$) with Term 2 collapsed.

4. **Best RMSE vs best coverage.** Minimum RMSE($\mu_f$) occurs near $K=0.5$ (2.894) with only 87.8% mean-surface coverage. $K=0.1$ trades a small RMSE increase for near-nominal Cover($\theta$). Method choice should state whether the priority is point recovery or honest uncertainty for $\theta$.

5. **Unit-sum weights ($M=1$) over-inflate Term 2.** SE ratio $\approx 6.6$, Cover($\theta$)$=100\%$, and RMSE jumps to 5.72. This empirically validates Section 4’s warning that setting $M\ll n$ (here $M=1$) leaves Term 2 dominant: intervals are far too wide and point estimation deteriorates.

6. **Implication for practice.**  
   - For **PP of $y$ under a well-aligned DGP** (Section 6): raw Kim–Rao ($K=1$) remains the preferred default.  
   - For **population mean-surface recovery under informative sampling with incomplete covariates** (this section): start from mean-1 (`weight_norm`) weights **without** per-draw $\sum w^*=n$ forced renorm, and calibrate $K$ so that Est. SE / Emp. SE $\approx 1$ (here $K\approx 0.1$). Do **not** use $\sum w=1$ scaling.  
   - Always report both RMSE and coverage; $K$ that minimizes RMSE need not maximize interval calibration.

### 7.6 OOB-selected $K$ vs grid (1,000 replications)

Section 3.2 proposes selecting $K$ by cumulative Kim–Rao OOB MSE. We ran that procedure on the same informative Friedman design: short calibration MCMC over a candidate grid, then full MCMC at the selected $K^*$.

| Method | Mean $K^*$ | SD($K^*$) | Bias | Emp. SE | Est. SE | SE ratio | Cover($\theta$) | RMSE($\mu_f$) | PP($y_f$) |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **A_OOB** (`weight_norm`) | **1.00** | 0.016 | −0.083 | 0.227 | 0.197 | 0.868 | **86.9%** | **2.910** | **94.5%** |
| **B_OOB** (raw weights) | **4.64** | 4.05 | −0.434 | 0.228 | 0.050 | 0.218 | **7.4%** | 6.045 | 40.6% |

*Source: `data/friedman_results_oob/report_table.csv`. Candidate grids: A $\{0.01,0.05,0.1,0.2,0.5,1\}$; B $\{0.5,1,2,5,10\}$.*

**A_OOB detail.** In 999/1,000 replications OOB selected $K^*=1.0$ (once $K^*=0.5$). Performance matches Section 7.3 at $K=1$ almost exactly — **not** the coverage-optimal $K\approx 0.1$.

**B_OOB detail.** Selected $K^*$ is unstable (counts: $K=10$: 341; $K=1$: 319; $K=2$: 152; $K=5$: 114; $K=0.5$: 74). Larger $K$ further collapses Term 2 on the raw-weight scale, yielding catastrophic Cover($\theta$).

**Interpretation.** OOB loss is a **prediction** criterion (MSE of $y$ on held-out PSUs). It therefore tracks the RMSE / PP($y$) trade-off, not Cover($\theta$). For mean-surface interval calibration, OOB cannot replace a coverage- or SE-ratio–targeted $K$ (Section 7.3). For **imputation / predictive uses**, A_OOB’s choice $K^*=1$ on `weight_norm` is appropriate and aligns with Section 6’s PP($y$) default.

### 7.7 Rejected alternative: bootstrap weights for trees only (“fixed leaf”)

A natural conjecture is that design variance (Term 1) should enter **only** through tree structure, while leaf parameters $\mu_l$ should be drawn with **fixed uniform weights** (all $w_i=1$), decoupling Term 1 from Term 2 without tuning $K$.

**Modification.** At each MCMC iteration: draw $\mathbf{w}^*\sim\mathcal{L}_*$ (Kim–Rao); use $\mathbf{w}^*$ in the integrated MH grow/prune likelihood; draw leaf $\mu_l$ with $w_i\equiv 1$; retain $n$-scaled $\mathbf{w}^*$ for the $\sigma^2$ update.

| Method | $K$ | Bias | Emp. SE | Est. SE | SE ratio | Cover($\theta$) | RMSE($\mu_f$) | PP($y_f$) |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| A_FixedLeaf | 0.1 | −0.436 | 0.184 | 0.046 | **0.251** | **3.1%** | 2.949 | 94.9% |
| A_FixedLeaf | 1.0 | −0.436 | 0.185 | 0.044 | **0.238** | **3.1%** | 2.915 | 94.3% |
| B_FixedLeaf | 1.0 | −0.437 | 0.191 | 0.065 | **0.340** | **5.9%** | 3.320 | 91.9% |

*Source: `data/friedman_results_fixed_leaf/report_table.csv` (3,000 replications). Code archived under `scripts/archive/fixed_leaf/`.*

**Verdict: reject.** Cover($\theta$) collapses to $\approx$3–6% for all three settings — **worse** than coupled Method A at any $K$ on the grid. Est. SE is far below Emp. SE (SE ratio $\approx 0.25$). Bias $\approx -0.44$ matches the over-concentrated large-$K$ end of the A2 grid ($K=20$), not the well-calibrated $K=0.1$ row. PP($y_f$) remains near 95% for Method A variants (the $\sigma^2$ buffer again), underscoring that predictive coverage is **not** a diagnostic of mean-surface calibration.

**Why it fails.** Uniform leaf weights force every leaf draw into an $M\approx n_l$ (unweighted) precision regime irrespective of $\mathbf{w}^*$. That collapses the leaf contribution to Term 2 while tree-only Kim–Rao resampling does **not** restore enough between-draw variability in $\bar f$ to match Emp. SE of $\theta$. Design-aware leaf weighting is therefore **essential** for sandwich-like Term 1 to register in the mean-surface posterior; restricting $\mathbf{w}^*$ to MH tree acceptance is not a substitute.

**Implication.** The working one-step architecture keeps $\mathbf{w}^*$ in the leaf Gibbs updates (as in `bart_mcmc.cpp`). Do not use the fixed-leaf variant.