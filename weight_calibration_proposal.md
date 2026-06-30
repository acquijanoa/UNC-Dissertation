# Theoretical Justification for Raw Kim-Rao Weights in One-Step Random Weight MCMC

**Author:** Álvaro Quijano  
**Date:** June 2026  
**Context:** Research Dissertation — One-Step Random Weight MCMC for DPMM and BART

---

## 1. The Core Law of Total Variance Decomposition

In the proposed **one-step Random Weight MCMC** framework, we draw a fresh set of Kim-Rao multinomial bootstrap weights $\mathbf{w}^{*(t)}$ at each MCMC iteration $t$. The total variance of the resulting parameter draws $\{\theta^{(t)}\}_{t=1}^T$ is decomposed using the Law of Total Variance with respect to the bootstrap distribution $\mathcal{L}_*$:

$$
\text{Var}(\theta) = \underbrace{\text{Var}_{\mathbf{w}^*} \bigl( E[\theta \mid \mathbf{w}^*] \bigr)}_{\text{Term 1: Bootstrap Design Variance}} + \underbrace{E_{\mathbf{w}^*} \bigl[ \text{Var}(\theta \mid \mathbf{w}^*) \bigr]}_{\text{Term 2: Conditional Posterior Variance}}
$$

**Term 1** is the variance of the weighted point estimator $\hat{\theta}(\mathbf{w}^*)$ over repeated bootstrap draws. By construction of the Kim and Rao (2024) multinomial bootstrap — which produces rescaling factors $r_i^*$ with $E_*[r_i^*] = 1$ and $\text{Var}_*(S_w^*(\theta)) = \hat{V}[S_w(\theta) \mid \mathcal{F}_N]$ — this term converges in probability to the design-consistent sandwich covariance:

$$
\text{Var}_{\mathbf{w}^*} \bigl( E[\theta \mid \mathbf{w}^*] \bigr) \xrightarrow{p} \Sigma_{\text{sandwich}} = H_0^{-1} J_0 H_0^{-1}
$$

where $H_0 = -E[\nabla^2 \ell(\theta_0)]$ is the population Hessian and $J_0 = \text{Var}[\nabla \ell(\theta_0)]$ is the design-based score variance.

**Term 2** is the expected model-based posterior variance, averaged over bootstrap draws.

To achieve exact one-step variance resolution — where the empirical variance of the MCMC draws matches $\Sigma_{\text{sandwich}}$ — we require Term 2 to be asymptotically negligible:

$$
E_{\mathbf{w}^*} \bigl[ \text{Var}(\theta \mid \mathbf{w}^*) \bigr] \to 0
$$

---

## 2. Scale-Invariance of Term 1 and the Role of $K$

### 2.1 Why Scaling Weights Does Not Affect Term 1

A critical property is that **Term 1 is invariant to any constant multiplicative rescaling of the weights**. This follows from the fact that the weighted point estimator is scale-invariant. In BART, for example, the weighted mean residual in leaf $l$ is:

$$
\hat{\mu}_l(\mathbf{w}^*) = \frac{\sum_{i \in l} w_i^* R_i}{\sum_{i \in l} w_i^*}
$$

Replacing $w_i^* \mapsto K w_i^*$ for any constant $K > 0$ leaves $\hat{\mu}_l$ unchanged. Consequently, $\text{Var}_{\mathbf{w}^*}(E[\theta \mid K\mathbf{w}^*]) = \text{Var}_{\mathbf{w}^*}(E[\theta \mid \mathbf{w}^*]) \approx \Sigma_{\text{sandwich}}$ for any $K$. This means $K$ is a free parameter that can be used to control Term 2 *without* distorting Term 1.

### 2.2 Term 2 and the Role of the Weight Sum $\hat{N}$

In traditional fixed-weight pseudo-Bayesian methods (Savitsky & Toth, 2016), the weights are fixed and normalized to $\sum_i \tilde{w}_i = n$. In that setting, the posterior variance is the sole carrier of uncertainty and must be calibrated to match $\Sigma_{\text{sandwich}}$. Under the one-step random weight approach, the logic is inverted: Term 1 already captures $\Sigma_{\text{sandwich}}$, so we want Term 2 to vanish.

In the BART Gibbs step, the conditional posterior variance of the leaf parameter $\mu_l$ is:

$$
V_l^{(t)} = \frac{1}{\dfrac{W_l^{*(t)}}{\sigma^2} + \dfrac{1}{\sigma_\mu^2}}
$$

where $W_l^{*(t)} = \sum_{i \in l} w_i^{*(t)}$ is the sum of bootstrap weights in leaf $l$. The expected conditional posterior variance is therefore:

$$
E_{\mathbf{w}^*}\bigl[\text{Var}(\mu_l \mid \mathbf{w}^*)\bigr] = E_{\mathbf{w}^*}[V_l^{(t)}] \approx \frac{\sigma^2}{E_*[W_l^*] + \sigma^2/\sigma_\mu^2}
$$

This vanishes if and only if $E_*[W_l^*] \gg \sigma^2/\sigma_\mu^2$, i.e., the expected leaf weight sum is large relative to the prior variance ratio.

### 2.3 The Practical Reality of $\hat{N}$

In practice, the sum of sampling weights $\hat{N} = \sum_{i \in S} w_i$ is a random variable that is **not guaranteed to equal the true finite population size $N$**. In particular:

- **Calibrated/post-stratified weights:** Weights may be adjusted by nonresponse correction or raking, altering their total.
- **Domain-level analysis:** For a subpopulation $S_{\text{sub}} \subset S$, the effective leaf weight sum $\hat{N}_{\text{sub}} = \sum_{i \in S_{\text{sub}}} w_i$ can be far smaller than $N$, making Term 2 non-negligible.

In these cases, the conditional posterior variance is of order:

$$
E_{\mathbf{w}^*}\bigl[\text{Var}(\theta \mid \mathbf{w}^*)\bigr] \approx \frac{n}{\hat{N}} V_{\text{model}} > 0
$$

which adds spurious variance inflation (over-coverage) on top of Term 1.

### 2.4 Using $K$ as a Variance-Collapsing Guarantee

To guarantee that Term 2 collapses to zero regardless of the realized scale of $\hat{N}$, multiply all bootstrap weights by a constant scaling factor $K \gg 1$:

$$
w_{hij}^{*(t)} = K \cdot r_{hi}^{*(t)} \, w_{hij}
$$

By §2.1, this does not alter Term 1. It scales the leaf weight sum to $K \cdot \hat{N}_l$, so:

$$
V_l^{(t)} \approx \frac{\sigma^2}{K \cdot \hat{N}_l / n \cdot n_l} \xrightarrow{K \to \infty} 0
$$

and therefore:

$$
E_{\mathbf{w}^*}\bigl[\text{Var}(\theta \mid K\mathbf{w}^*)\bigr] \approx O\!\left(\frac{1}{K \cdot \hat{N}}\right) \approx 0
$$

**$K$ is thus not a design-effect tuning parameter** (as in two-step post-hoc adjustments). It is a purely computational device that ensures the conditional posterior variance is dominated by the prior only in the trivial zero-weight (OOB) case.

### 2.5 Practical Rule for Choosing $K$

A sensible default is to choose $K$ so that the expected leaf weight sum exceeds the effective sample size $n_{\text{eff}} = n / \widehat{\text{Deff}}$ by a comfortable margin. Concretely, set:

$$
K = \frac{n}{\bar{w} \cdot n_{\min}}
$$

where $\bar{w} = \hat{N}/n$ is the average design weight and $n_{\min}$ is the smallest expected leaf count (a tuning floor, e.g., $n_{\min} = 5$). Under this choice, the leaf weight sum $K \cdot \hat{N}_l / n \cdot n_l \geq n_l / n_{\min} \geq 1$ for all non-trivial leaves, ensuring Term 2 is negligible relative to $\Sigma_{\text{sandwich}}$.

### 2.6 The Indispensable Role of Kim-Rao: Why Arbitrary Randomizations Fail

The argument above raises a natural question: since $K$ collapses Term 2 for *any* set of randomized weights $\mathbf{w}^{*(t)}$, can we replace the Kim-Rao bootstrap with a simpler or more convenient randomization scheme — for example, multiplying the design weights $w_i$ by i.i.d. $\text{Exp}(1)$ draws, or drawing Dirichlet weights, or using a simple nonparametric bootstrap?

**The answer is no.** $K$ controls Term 2 regardless of how the randomization is structured. But Term 2 is only half of the Law of Total Variance. The one-step resolution works if and only if *both* conditions hold simultaneously:
1. Term 2 $\to 0$ (controlled by $K$).
2. Term 1 $= \Sigma_{\text{sandwich}}$ (controlled exclusively by the choice of weight distribution $\mathcal{L}_*$).

Condition 2 is a precise requirement on the bootstrap distribution: the variance of the bootstrap-weighted score function $S_{w^*}(\theta) = N^{-1} \sum_{i \in S} w_i^{*(t)} \nabla_\theta \ell(\theta; x_i, y_i)$ must satisfy:

$$
\text{Var}_*\bigl(S_{w^*}(\theta_0)\bigr) = \hat{V}\bigl[S_w(\theta_0) \mid \mathcal{F}_N\bigr] \xrightarrow{p} J_0
$$

where $J_0 = \text{Var}_{\mathcal{F}_N}[\nabla_\theta \ell(\theta_0; x_i, y_i)]$ is the true design-based score variance under the complex survey.

**This condition is not free.** It encodes the entire covariance structure of the complex sampling design — stratification, clustering, and PPS selection. An arbitrary randomization fails this requirement:

- **i.i.d.\ $\text{Exp}(1)$ multipliers:** Produce $\text{Var}_*[S_{w^*}] \propto \sum_i w_i^2 S_i(\theta_0)^2$, which ignores within-cluster correlation and equals the design-based score variance only under simple random sampling. Under stratified cluster sampling, this systematically underestimates $J_0$, making Term 1 smaller than $\Sigma_{\text{sandwich}}$ and causing undercoverage.

- **Dirichlet$(\alpha \mathbf{1}_n)$ weights:** The bootstrap variance of the score depends on $\alpha$ and converges to the empirical score variance under i.i.d. sampling, not to $\hat{V}[S_w(\theta_0) \mid \mathcal{F}_N]$.

- **Nonparametric bootstrap (resampling units with replacement):** Treats all units as exchangeable, ignoring the stratified cluster structure. The resulting bootstrap score variance is the i.i.d. score variance, not the design-based sandwich variance.

**The Kim-Rao multinomial bootstrap is designed precisely to satisfy Condition 2.** By resampling $n_h - 1$ clusters with replacement within each stratum $h$ and using the rescaling factor $r_{hi}^* = k_h m_{hi}^*$ (where $k_h = n_h/(n_h - 1)$), Kim et al. (2024) prove in their Lemma S7 that:

$$
E_*\bigl[S_{w^*}(\theta)\bigr] = S_w(\theta), \qquad \text{Var}_*\bigl(S_{w^*}(\theta)\bigr) = \hat{V}\bigl[S_w(\theta) \mid \mathcal{F}_N\bigr]
$$

It is this exact moment-matching property — not any particular distributional form — that makes the Kim-Rao weights uniquely suited for the one-step approach. To summarize the separation of roles:

> **$K$ controls Term 2:** any weight randomization can have Term 2 collapsed to zero by choosing $K$ large enough.  
> **Kim-Rao controls Term 1:** only a bootstrap distribution that exactly replicates the design-based score variance guarantees Term 1 $= \Sigma_{\text{sandwich}}$.

The two requirements are therefore **orthogonal**: $K$ is a computational device; Kim-Rao is a statistical requirement. Neither alone is sufficient.

---


## 3. The One-Step Resolution Mechanism

With Term 2 collapsed — whether by the population-scale of $\hat{N}$ alone or by the additional scaling factor $K$ — the Gibbs draw for each terminal node parameter degenerates to its posterior mean:

$$
\mu_l^{(t)} \sim N\!\bigl(\hat{\mu}_l^{(t)},\, V_l^{(t)}\bigr) \xrightarrow{K\hat{N} \to \infty} \hat{\mu}_l^{(t)} = \frac{\sum_{i \in l} w_i^{*(t)} R_i}{\sum_{i \in l} w_i^{*(t)}}
$$

The MCMC draws at different iterations $t$ differ only because the bootstrap weights $\mathbf{w}^{*(t)}$ differ. The empirical variance of the draws is then entirely determined by Term 1:

$$
\text{Var}(\mu_l) \approx \text{Var}_*\!\bigl(\hat{\mu}_l(\mathbf{w}^*)\bigr) \approx \Sigma_{\text{sandwich}}
$$

This is the **one-step resolution**: no post-processing, no matrix inversion, no Bartlett-identity correction. The design-consistent uncertainty is natively embedded in the MCMC fluctuations.

### 3.1 Treating $K$ as a Bayesian Hyperparameter: OOB Tuning

In the classical BART framework, hyperparameters such as the prior variance $\sigma_\mu^2$ and the error variance $\sigma^2$ are fixed by cross-validation or empirical Bayes prior to fitting. The scaling factor $K$ can be incorporated into the same framework. However, the Kim-Rao bootstrap offers something more powerful: a **computationally free internal validation set** at every MCMC iteration.

**The OOB mechanism.** Because the Kim-Rao multinomial resampling draws $n_h - 1$ clusters with replacement within each stratum $h$, any individual PSU $(h, i)$ has probability $(1 - 1/n_h)^{n_h - 1} \approx e^{-1} \approx 0.37$ of receiving a zero count $m_{hi}^* = 0$. Units in such PSUs receive bootstrap weight $w_{hij}^* = 0$ and are entirely excluded from the weighted pseudo-likelihood at iteration $t$. These $\approx 37\%$ OOB units are **not in-bag observations** — they constitute a genuine held-out validation set that changes at every MCMC step.

**OOB predictive loss as an objective for $K$.** At iteration $t$, let $\mathcal{I}^{(t)} = \{i : w_i^{*(t)} > 0\}$ be the in-bag set and $\mathcal{O}^{(t)} = \{i : w_i^{*(t)} = 0\}$ be the OOB set. Given the current MCMC state $f^{(t)}(x) = \sum_j \mu_{jl(x)}^{(t)}$ fitted on the in-bag data, the OOB predictive loss is:

$$
\mathcal{L}_{\text{OOB}}(K, t) = \frac{1}{|\mathcal{O}^{(t)}|} \sum_{i \in \mathcal{O}^{(t)}} \bigl(y_i - f^{(t)}(x_i)\bigr)^2
$$

Note that this loss depends on $K$ through the tree structure and leaf parameters, which are fitted using the $K$-scaled in-bag weights. Define the cumulative OOB loss across post-burn-in iterations:

$$
\mathcal{L}_{\text{OOB}}(K) = \frac{1}{T_{\text{post}}} \sum_{t=t_{\text{burn}}}^{T} \mathcal{L}_{\text{OOB}}(K, t)
$$

The optimal hyperparameter is then:

$$
K^* = \underset{K \in \mathcal{K}}{\arg\min} \; \mathcal{L}_{\text{OOB}}(K)
$$

**Theoretical justification.** This tuning strategy has three important properties:

1. **Term 1 is invariant to $K$** (by §2.1), so selecting $K$ via OOB loss does not alter the design-consistent sandwich variance captured in Term 1. The frequentist coverage guarantee of the one-step approach is preserved regardless of the chosen $K$.

2. **The OOB units are truly held out.** Because the Kim-Rao bootstrap is the resampling mechanism inside the MCMC loop, the OOB set is determined by the same randomization that drives the posterior fluctuations. There is no leakage between the model-fitting step and the validation step within any given iteration.

3. **Large-sample convergence.** As $n \to \infty$, the optimal $K^*$ from OOB tuning grows without bound, since larger $K$ drives the conditional posterior variance (Term 2) to zero, producing sharper in-bag fits with smaller residuals and consequently smaller OOB loss. In finite samples, the OOB criterion implicitly selects the $K$ that best balances in-sample fit against leave-cluster-out generalization.

**Practical implementation.** Rather than searching over a grid of $K$ values (which would require multiple full MCMC runs), $K$ can be treated as a Bayesian hyperparameter with a log-uniform prior $p(K) \propto 1/K$ over some range $[K_{\min}, K_{\max}]$ and updated adaptively during the warm-up phase by monitoring $\mathcal{L}_{\text{OOB}}(K, t)$ at each iteration. A gradient-free line search (e.g., golden section) on the cumulative OOB loss over the burn-in period is sufficient in practice.

---

## 4. Analysis of Alternative Normalization Schemes


To understand what is lost by normalizing the bootstrap weights, note that if weights are scaled to sum to a target quantity $M$, the weighted estimator is scale-invariant so **Term 1 remains fixed at $\Sigma_{\text{sandwich}}$**. However, Term 2 scales inversely with $M$. To see this, observe that the average leaf weight sum is $E_*[W_l^*] = M \cdot n_l / n$, so:

$$
E_{\mathbf{w}^*}\bigl[\text{Var}(\theta \mid \mathbf{w}^*)\bigr] \approx \frac{\sigma^2}{M n_l / n} \sim \frac{n}{M} V_{\text{model}}
$$

where $V_{\text{model}} = \sigma^2 n / (n_l \cdot n) = \sigma^2 / n_l$ is the standard unweighted posterior variance. This yields the following consequences:

### Case A: Normalizing to Sample Size ($M = n$)

If weights are normalized to sum to $n$ (i.e., $\tilde{w}_i = w_i^* \cdot n / \hat{N}^*$), then $\frac{n}{M} = 1$ and:

$$
\text{Var}(\theta) \approx \Sigma_{\text{sandwich}} + V_{\text{model}}
$$

Since $\Sigma_{\text{sandwich}} = \text{Deff} \cdot V_{\text{model}}$ and $\text{Deff} > 1$ in complex designs, the over-coverage is modest but systematic:

$$
\frac{\text{Var}(\theta)}{\Sigma_{\text{sandwich}}} = 1 + \frac{1}{\text{Deff}} \in \left(1,\, 2\right)
$$

### Case B: Normalizing to Effective Sample Size ($M = n_{\text{eff}} = n / \text{Deff}$)

If weights are normalized to sum to $n_{\text{eff}}$ (the recommendation sometimes made to *inflate* posterior variance to match design uncertainty in fixed-weight models), then $\frac{n}{M} = \text{Deff}$ and:

$$
E_{\mathbf{w}^*}\bigl[\text{Var}(\theta \mid \mathbf{w}^*)\bigr] \approx \text{Deff} \cdot V_{\text{model}} \approx \Sigma_{\text{sandwich}}
$$

Substituting into the Law of Total Variance:

$$
\text{Var}(\theta) \approx \Sigma_{\text{sandwich}} + \Sigma_{\text{sandwich}} = 2\,\Sigma_{\text{sandwich}}
$$

The total MCMC variance is **exactly double** the target, inflating credible intervals by a factor of $\sqrt{2} \approx 1.41$. This is the worst possible normalization choice under the one-step approach.

---

## 5. Contrast with the Two-Step Approach of Savitsky & Williams

The distinction between the one-step and two-step approaches reflects a fundamental inversion of the role of Term 1 and Term 2:

| | **Savitsky & Toth (2016)** | **Williams & Savitsky (2021)** | **One-Step (Proposed)** |
|:---|:---:|:---:|:---:|
| Bootstrap resampling at each MCMC step? | ✗ | ✗ | ✓ |
| Weight normalization | $\sum \tilde{w}_i = n$ | $\sum \tilde{w}_i = n$ | None (raw, or $\times K$) |
| Term 1 | $0$ | $0$ | $\Sigma_{\text{sandwich}}$ |
| Term 2 | $\approx V_{\text{model}} \ll \Sigma_{\text{sandwich}}$ | $\approx V_{\text{model}}$ | $\approx 0$ |
| Net result | **Under-coverage** | **Under-coverage** (corrected post-hoc) | **Exact** |

**Savitsky & Toth (2016)** normalize weights to sum to $n$ so that the model-based posterior variance $V_{\text{model}}$ is on the scale of $1/n$. Because $H_0 \neq J_0$ under informative sampling, the resulting posterior **undercovers** — the credible intervals are too narrow by a factor proportional to $\text{Deff}$.

**Williams & Savitsky (2021)** keep the same normalization but apply a *post-hoc* affine transformation to the MCMC draws to rescale Term 2 up to $\Sigma_{\text{sandwich}}$. This requires explicit computation of the sandwich matrices $H_0$ and $J_0$, which is analytically intractable for non-parametric models like BART, and can distort non-Gaussian posterior shapes.

In the **one-step Random Weight MCMC**, Term 1 is non-zero and already equals $\Sigma_{\text{sandwich}}$ by the Kim-Rao CLT for bootstrap score functions. Term 2 must therefore be collapsed to zero to avoid over-coverage. This is achieved by using raw (or $K$-scaled) bootstrap weights, which naturally place the leaf weight sum on a scale that renders the model-based posterior variance negligible.
