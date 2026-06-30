# Non-Parametric Bayesian Inference under Complex Survey Designs via Random Weight MCMC: Applications to DPMM and BART

## Introduction and Superpopulation Framework

In finite population inference under complex survey designs, classical pseudo-Bayesian methods struggle with covariance matrix estimation due to informative sampling. We formalize this problem using a superpopulation framework. Let $\mathcal{M}$ denote the superpopulation model governing the data-generating mechanism. The finite population $U = \{1, \dots, N\}$ is a realization from $\mathcal{M}$. A complex sample $S \subset U$ of size $n$ is subsequently drawn from the finite population according to a probabilistic sampling design $p(S \mid \mathbf{y}_U, \mathbf{z}_U)$, where $\mathbf{z}_U$ denotes design variables. 

To draw inference on the superpopulation parameter $\theta$, the pseudo-likelihood approach evaluates the posterior using fixed survey weights $w_i \propto 1/\pi_i$:
$$p_{\text{pseudo}}(\theta \mid \mathbf{y}_S, \mathbf{w}) \propto p(\theta) \prod_{i \in S} p(y_i \mid \theta_i)^{w_i}$$

However, as established by Savitsky and Toth (2016), while this fixed-weight pseudo-posterior is consistent for $\theta$, it fails Bartlett's second identity ($H_0 \neq J_0$). It systematically underestimates the posterior variance and warps the shape of the parameter covariance matrix because it treats the survey weights $\mathbf{w}$ as fixed, precise information, thereby ignoring the complex sampling design variance. To resolve this, Williams and Savitsky (2021) proposed a two-step post-hoc affine transformation. However, this two-step approach requires explicit computation of the Hessian and sandwich matrices and can distort non-normal posteriors. The Random Weight MCMC proposed here corrects this by replacing the two-step adjustment with a unified, one-step simulation-based variance resolution.

This proposal introduces the **Random Weight MCMC**, a novel unified methodological framework that structurally bypasses the two-step problem by natively integrating the unified bootstrap weight mechanism of Kim, Rao, and Wang (2024) directly into the transition steps of non-parametric Bayesian algorithms. We demonstrate this framework across two major applications: unsupervised density estimation via Dirichlet Process Mixture Models (DPMM) and supervised predictive modeling via Bayesian Additive Regression Trees (BART).

---

## PART I: The Random Weight MCMC Core Methodology

### 1. Inner-Loop Bootstrap Weight Resampling
Instead of relying on fixed survey weights, we target a randomized pseudo-posterior by evaluating the transition probabilities using dynamically fluctuating weights:
$$p_{\text{pseudo}}(\mathbf{y}_S \mid \theta, \mathbf{w}^*) = \prod_{i \in S} p(y_i \mid \theta_i)^{w_i^*}$$

At each MCMC iteration, we draw a fresh set of bootstrap weights $\mathbf{w}^* = (w_1^*, \dots, w_n^*)$ from the Kim et al. (2024) resampling distribution. These weights are constructed as $w_i^* = r_i^* w_i$, where the rescaling factors $r_i^*$ are drawn via a multinomial distribution over the primary sampling units within each stratum to replicate the exact design variance. By redrawing $w_i^*$ at every transition step, the score function is stochastically perturbed, natively injecting the complex survey design uncertainty into the Markov chain.

**Computational Sparsity:** A natural property of the bootstrap multinomial resampling is that approximately $1/e \approx 37\%$ of the primary sampling units will receive a rescaling factor of $r_i^* = 0$ at any given MCMC iteration. This "dropout" provides a massive computational advantage: the pseudo-likelihood evaluation inside the MCMC inner loop can entirely skip ~37% of the dataset at every transition step, significantly reducing the floating-point operations required for large datasets like HCHS/SOL while natively regularizing the model against overfitting to specific observations.

### 2. Theoretical Justification: Point Consistency and One-Step Variance Resolution
The frequentist validity of the Random Weight MCMC is formally established by bridging the Horvitz-Thompson uniform empirical process theory of Han and Wellner (2021), the pseudo-posterior consistency of Savitsky and Toth (2016), and the bootstrap central limit theory of Kim, Rao, and Wang (2024). 

#### Non-parametric Consistency
Unlike traditional parametric models, non-parametric Bayesian learners like DPMMs and BART operate over infinite-dimensional or highly complex function spaces. We justify the design-consistency of our point estimators by modeling our tree structures and mixture components as functions belonging to a VC-subgraph class $\mathcal{F}$ with bounded bracketing entropy. By applying the uniform Glivenko-Cantelli and Donsker theorems for Horvitz-Thompson empirical processes (Han and Wellner, 2021), the weighted empirical measure:
$$\mathbb{P}_N^{\pi} = \frac{1}{N} \sum_{i \in S} w_i \delta(X_i)$$
converges uniformly over $\mathcal{F}$, ensuring that our estimated conditional expectations converge in probability to the true population conditional expectations under the design-weighted $L_2$ norm:
$$\int \left(\hat{m}(X) - m(X)\right)^2 dF_X(X) \xrightarrow{p} 0$$

#### Correcting the BvM Failure via One-Step Variance Resolution
A central challenge in pseudo-Bayesian inference is the failure of the Bernstein-von Mises (BvM) theorem under complex sampling. As demonstrated by Williams and Savitsky (2021), the pseudo-posterior variance from a fixed-weight likelihood concentrates around the inverse Fisher Information $H_0^{-1}$. However, the true design-consistent variance of the M-estimator is the sandwich form:
$$\Sigma_{\text{sandwich}} = H_0^{-1} J_0^{\pi} H_0^{-1}$$
where $J_0^{\pi}$ is the design variance of the score function. Because $H_0 \neq J_0^{\pi}$ under complex designs, fixed-weight pseudo-posteriors exhibit severe undercoverage.

To resolve this, Williams and Savitsky (2021) proposed a two-step post-hoc affine transformation. We replace this with a **unified one-step variance resolution** by drawing fresh bootstrap weights $\mathbf{w}^{*(t)}$ from the Kim-24 stratified cluster bootstrap at each MCMC step. Let $\mathbb{P}_N^{*(b)}$ denote the bootstrap empirical measure. For any parameter $\theta$ estimated by the score equations, the bootstrap variance of the score function conditional on the sample satisfies:
$$V_* \left[ \sqrt{N} \mathbb{P}_N^{*(b)} S(\theta_0) \right] = \hat{V} \left[ \sqrt{N} \mathbb{P}_N^{\pi} S(\theta_0) \;\middle|\; \mathcal{F}_N \right] \xrightarrow{p} J_0^{\pi}$$

Applying a first-order Taylor expansion to the bootstrap score equation $\mathbb{P}_N^{*(b)} S\left(\hat{\theta}^{*(b)}\right) = 0$ yields:
$$\sqrt{N} \left( \hat{\theta}^{*(b)} - \hat{\theta}_{\pi} \right) = \left[ \mathbb{I}_N^{*(b)}(\bar{\theta}) \right]^{-1} \sqrt{N} \mathbb{P}_N^{*(b)} S\left(\hat{\theta}_{\pi}\right)$$
Taking the conditional bootstrap variance on both sides:
$$\text{Var}_* \left[ \sqrt{N} \left( \hat{\theta}^{*(b)} - \hat{\theta}_{\pi} \right) \right] = \left[ \mathbb{I}_N^{*(b)}(\bar{\theta}) \right]^{-1} \left( \text{Var}_* \left[ \sqrt{N} \mathbb{P}_N^{*(b)} S\left(\hat{\theta}_{\pi}\right) \right] \right) \left[ \mathbb{I}_N^{*(b)}(\bar{\theta}) \right]^{-T} + o_p(1)$$
Since the bootstrap Hessian $\mathbb{I}_N^{*(b)}(\bar{\theta}) \xrightarrow{p} H_0$ and the score variance converges to $J_0^{\pi}$, the variance of the Random Weight MCMC draws automatically converges to:
$$\text{Var}_* \left[ \sqrt{N} \left( \hat{\theta}^{*(b)} - \hat{\theta}_{\pi} \right) \right] \xrightarrow{p} H_0^{-1} J_0^{\pi} H_0^{-1}$$
This one-step perturbation bypasses the need for post-processing matrix calculations, preserves the natural boundaries of the parameter space, and naturally generalizes to complex, variable-dimension Bayesian models like DPMMs and BART.

### 3. Extending Kim-24 to Non-Concave Likelihoods
A notable challenge in applying the Kim-24 bootstrap to non-parametric models (like DPMMs or BART) is the requirement of global log-likelihood concavity. However, the BvM synthesis holds through two mechanisms:
1. **Conditional Concavity:** The Kim-24 weights are injected *inside* the Gibbs sampler, which operates conditionally. Conditionally on the latent structures (cluster assignments in DPMM, or tree structures in BART), the likelihoods for the continuous parameters reduce to standard exponential families, ensuring strict conditional concavity at every transition step.
2. **Local Asymptotic Normality (LAN):** Modern Bayesian asymptotics rely only on LAN (Kleijn and van der Vaart, 2012). Because the pseudo-posterior contracts into a shrinking neighborhood around the true data-generating parameter, the log-likelihood surface is locally quadratic (strictly concave) in the stationary region.

### 4. Robustness via Generalized MMD-Bayes
Standard pseudo-Bayesian consistency proofs (e.g., Savitsky and Toth, 2016) implicitly rely on minimizing the Kullback-Leibler (KL) divergence between the true superpopulation and the model family. However, KL-divergence is highly sensitive to tail behavior; under severe model misspecification or data contamination, KL-based posteriors fail to contract meaningfully. To circumvent this, the Random Weight MCMC can be extended into a **Generalized Bayes** framework by abandoning the strict pseudo-likelihood in favor of an integral probability metric loss function. 

Specifically, we employ the Maximum Mean Discrepancy (MMD) formulated within a Reproducing Kernel Hilbert Space (RKHS) (Chérief-Abdellatif and Alquier, 2020). By replacing the standard empirical measure with the Kim-24 randomized bootstrap empirical measure $\mathbb{P}_N^{\pi*} = \frac{1}{N} \sum_{i \in S} w_i^* \delta_{y_i}$, we construct a survey-weighted MMD pseudo-loss:
$$ p_{\text{robust}}(\theta \mid \mathbf{y}_S, \mathbf{w}^*) \propto \exp\left( -\lambda \cdot \text{MMD}^2(\mathbb{P}_N^{\pi*}, P_\theta) \right) \pi(\theta) $$
This Generalized MMD-Bayes formulation achieves a theoretical synthesis: the kernel mean embeddings confer strict asymptotic robustness against model misspecification, while the dynamic injection of the Kim-24 bootstrap weights ensures the stationary distribution natively captures the target frequentist sandwich covariance $\Sigma = H_0^{-1} J_0 H_0^{-1}$ of the complex survey design.

**Practical Considerations: Tuning MMD Hyperparameters**
A known computational hurdle of the Generalized Bayes MMD framework is the reliance on un-estimable hyperparameters: the temperature scaling $\lambda$ and the RKHS kernel bandwidth. Because we have abandoned the true generative likelihood, standard Bayesian model selection via marginal likelihood is invalid. To resolve this, the kernel bandwidth is typically anchored using the median heuristic. Meanwhile, the temperature parameter $\lambda$, which dictates the spread of the pseudo-posterior, must be carefully tuned. 

Brilliantly, the Random Weight MCMC offers a "free" computational solution to tune $\lambda$ by exploiting the computational sparsity of the bootstrap. Because approximately $37\%$ of the data receives a weight of $w_i^* = 0$ at any given MCMC iteration, these unselected units natively form an **Out-Of-Bag (OOB) validation set**. Rather than relying on computationally heavy K-fold cross-validation or PSIS-LOO, $\lambda$ can be tuned by evaluating the predictive MMD loss of the current MCMC state $\theta^{(t)}$ exclusively on the ~37% OOB samples at each iteration. Selecting the $\lambda$ that minimizes this cumulative OOB loss ensures the Generalized posterior achieves optimal out-of-sample generalization without any additional computational overhead.

---

## PART II: Application 1 - CKM Phenotyping and Synthetic Data Generation via DPMM

To perform multivariate joint modeling and generate synthetic finite populations $U_{syn}$, we specify a full DPMM as our superpopulation model $\mathcal{M}$:
$$G \sim DP(\alpha, G_0)$$
$$\theta_i \mid G \sim G$$
$$y_i \mid \theta_i \sim F(\theta_i)$$

Setting $\alpha > 0$ allows the superpopulation model to learn latent cluster structure and generate diverse, continuously distributed synthetic populations that extend beyond the support of the observed sample $S$ (overcoming the limitations of the Pólya posterior where $\alpha \to 0$).

### Random Weight DPMM Inference
1. **Cluster Assignment Update:** During the Gibbs step for cluster assignments $z_i$, calculate the aggregate bootstrap weight for each existing DP cluster $k$: $w_k^* = \sum_{i \in S : z_i = k} w_i^*$. The predictive probabilities are evaluated using these dynamically fluctuating aggregate weights.
2. **Synthetic Population Generation:** For the unobserved population units $m = N - n$, draw sequentially from the posterior predictive distribution, re-weighting the allocation probabilities $n_k$ by the aggregate bootstrap weights $w_k^*$.

### Empirical Application: CKM Phenotyping in HCHS/SOL
The primary clinical application of this framework is identifying Cardiovascular-Kidney-Metabolic (CKM) phenotypes within the Hispanic Community Health Study / Study of Latinos (HCHS/SOL). Phenotyping requires joint multivariate modeling of both continuous markers (e.g., blood pressure, HbA1c, eGFR, BMI) and categorical markers (e.g., smoking status, medication adherence).

By specifying the base measure $G_0$ as a multivariate joint distribution (a mixture of multivariate Normals and multinomials), the DPMM organically learns the optimal number of latent CKM phenotypes present in the population without requiring *a priori* class specification like standard Latent Class Analysis (LCA) or K-means clustering. 

Crucially, HCHS/SOL employs a highly complex sampling design with deep geographic stratification and multi-stage clustering. A standard Bayesian DPMM would yield biased phenotype prevalences and artificially narrow credible intervals. The **Random Weight DPMM** resolves this: the dynamically fluctuating Kim-24 bootstrap weights ensure that the discovered CKM phenotypes accurately reflect the true target population of US Hispanics/Latinos, while natively inflating the posterior credible intervals of the phenotype profiles to account for intra-household and geographic correlations.

---

## PART III: Application 2 - Predictive Modeling via BART

Bayesian Additive Regression Trees (BART) model the conditional mean using a sum-of-trees: $f(x) = \sum_{j=1}^m g(x; T_j, M_j)$, where $T_j$ denotes the tree structure and $M_j$ denotes the terminal node parameters. Fitting BART directly to survey data typically ignores design variance, underestimating predictive credible intervals. 

### Random Weight BART Inference
We apply the Random Weight MCMC to the standard Bayesian backfitting algorithm:
1. **Tree Structure Proposals (Metropolis-Hastings):** When proposing a change to a tree structure (grow, prune, change), evaluate the marginal pseudo-likelihood for the MH acceptance ratio using the dynamically fluctuating weights $w_i^*$.
2. **Terminal Node Updates (Gibbs Step):** When drawing the terminal node parameters $\mu_{jl}$, evaluate the sufficient statistics (the sum of residuals and sum of squared residuals) using the **$w_i^*$-weighted average**. Because the weights dynamically inflate and deflate according to the survey design variance, the terminal nodes natively yield wider, correctly-calibrated credible intervals, and the tree structures organically resist overfitting to unrepresentative survey clusters.

### Empirical Application: Predicting Longitudinal CKM Outcomes in HCHS/SOL
While the DPMM serves to discover unsupervised CKM phenotypes cross-sectionally, the **Random Weight BART** serves as a powerful supervised engine to predict longitudinal outcomes—such as the incidence of major adverse cardiovascular events (MACE), kidney failure, or all-cause mortality—among the Hispanic/Latino population. 

BART is uniquely suited for this because of its non-parametric sum-of-trees architecture, which automatically captures complex, high-order, non-linear interactions between sociodemographic variables, baseline CKM markers, and the previously discovered latent phenotypes. However, predicting longitudinal outcomes in HCHS/SOL requires strict adherence to the survey design; ignoring the stratification and clustering leads to overconfident predictions and severe under-coverage of the predictive credible intervals.

By employing the Random Weight BART, the predictive modeling natively incorporates the complex survey design. The dynamically fluctuating Kim-24 bootstrap weights $w_i^*$ ensure that the regression trees do not overfit to the idiosyncratic noise of oversampled geographic clusters or households. Consequently, the resulting posterior predictive distributions provide robust, population-representative risk stratification for longitudinal CKM outcomes, complete with valid frequentist coverage for the predictive credible intervals.

---

## References
- **Breslow, N. E., & Wellner, J. A. (2007).** Weighted Likelihood for Semiparametric Models and Two-phase Stratified Samples. *Scandinavian Journal of Statistics*.
- **Han, Q., & Wellner, J. A. (2021).** Complex Sampling Designs: Uniform Limit Theorems and Applications. *Annals of Statistics*.
- **Kim, J. K., Rao, J. N. K., & Wang, Z. (2024).** Hypotheses Testing from Complex Survey Data Using Bootstrap Weights: A Unified Approach. *Journal of the American Statistical Association*.
- **Savitsky, T. D., & Toth, D. (2016).** Bayesian estimation under informative sampling. *Electronic Journal of Statistics*.
- **Williams, M. R., & Savitsky, T. D. (2021).** Uncertainty Estimation for Pseudo-Bayesian Inference Under Complex Sampling. *International Statistical Review*.
