## Supplemental Material for "Hypotheses Testing from Complex Survey Data Using Bootstrap Weights: A Unified Approach"

Without loss of generality, we use "\*" to denote the quantities for the bootstrap sample, use y1, . . . , y<sup>N</sup> to denote the elements of the finite population and use Y1, . . . , Y<sup>N</sup> to denote the corresponding random variables with respect to the super-population model. Besides, we use the same notation, "∥·∥", to denote the operator norm for a matrix and Euclidean norm for a vector. Specifically, for a vector x, ∥x∥ = (x <sup>T</sup>x) 1/2 , and for an m × m matrix M, ∥M∥ = supx∈Rm,∥x∥=1∥Mx∥.

#### S1 Comments on regularity conditions

In Condition C1, the concavity of l(θ; x, y) and the pointwise convergence of the design variance V {lw(θ) | F<sup>N</sup> } in Condition C2 are sufficient for the uniform convergence supθ∈A|lw(θ)− E{l(θ; X, Y )}| → 0 in probability with respect to the super-population model and the sampling mechanism over any compact set A ⊂ Θ (Andersen and Gill, 1982, Theorem II.1), and such a uniform convergence of lw(θ) guarantees the consistency of the point estimator θˆ as well as its bootstrap counterpart in Lemmas S1–S2 below. The twice continuous differentiability condition of l(θ; x, y) guarantees that the weighted score function Sˆ <sup>w</sup>(θ) is continuously differentiable with respect to θ, and such a condition applies to many classical models; refer to Assumption 2 of Yuan and Jennrich (1998) and Theorem 6.1 of Rubin-Bleuer and Schiopu-Kratina (2005) for a similar condition. However, the twice continuous differentiability condition may fail for some cases, such as l(θ; y) = |y −θ|, and it is not in our scope to investigate these cases; see Theorem 5.21 of van der Vaart (2000) for detailed discussion. If the finite population is independently generated from a super-population model with density function f(y; θ0) and if we are interested in θ0, then l(θ; y) = log f(y; θ). The unique maximizer of E{l(θ; Y )} is equal to the true parameter θ<sup>0</sup> when the model is identifiable and correctly specified.

The design-unbiased variance estimator of lw(θ) is used to guarantee the consistency of the bootstrap estimator θˆ ∗ , and it can be obtained under regularity conditions for general sampling designs; see Fuller (2009, Chapte 1) for details. Specifically, the sampling weights can be random with respect to the super-population model; refer to Section S4 for details. Condition C3 involves convergence property of the design variance, and it is required to get the central limit theorems in Condition C4 and Condition C7; see Theorem 5.1 of Rubin-Bleuer and Schiopu-Kratina (2005) and the very bottom of page 256 of Dumitrescu et al. (2020) for details. By "design variance V {Sˆ <sup>w</sup>(θ) | F<sup>N</sup> }", we mean the variance of the weighted score function Sˆ <sup>w</sup>(θ) with respect to the sampling mechanism conditional on the finite population. Specifically, since θˆ is random, we assume a uniform central limit theorem for θ ∈ B in Condition C7, and this is the reason why we need ΣS(θ) to be invertible for θ ∈ B. Condition C4 assumes the central limit theorem for Sˆ <sup>w</sup>(θ0) with respect to the super-population model and the sampling mechanism, and it is used to derive the asymptotic distribution of θˆ. In Condition C4, we implicitly assume that the variance of S<sup>N</sup> (θ0) = N <sup>−</sup><sup>1</sup> P<sup>N</sup> <sup>i</sup>=1 S(θ0; x<sup>i</sup> , yi) is negligible compared with the design V {Sˆ <sup>w</sup>(θ0) | F<sup>N</sup> } conditional on the finite population asymptotically, and such an assumption holds when the sample size is negligible compared with the finite population size under regularity conditions for general sampling designs. For example, if the finite population is independently generated from a super-population model, the variance of S<sup>N</sup> (θ0) is V {S<sup>N</sup> (θ0)} = N <sup>−</sup><sup>1</sup>Σ<sup>0</sup> with respect to the super-population model, where Σ<sup>0</sup> = V {S(θ0; X, Y )} is assumed to be finite. Under simple random sampling without replacement, the design variance is V {Sˆ <sup>w</sup>(θ0) | F<sup>N</sup> } = n −1 (1 − nN <sup>−</sup><sup>1</sup> )Σˆ <sup>0</sup>, where Σˆ <sup>0</sup> = (N − 1)<sup>−</sup><sup>1</sup> P<sup>N</sup> <sup>i</sup>=1{S(θ0; x<sup>i</sup> , yi) − S<sup>N</sup> (θ0)} ⊗2 is a consistent estimator of Σ<sup>0</sup> with respect to the super-population model, and A<sup>⊗</sup><sup>2</sup> = AA<sup>T</sup> . Thus, if nN <sup>−</sup><sup>1</sup> → 0, V {S<sup>N</sup> (θ0)} is asymptotically negligible compared with V {Sˆ <sup>w</sup>(θ0) | F<sup>N</sup> } under simple random sampling without replacement. In Condition C5, the uniform convergence of ˆIw(θ) is used to derive the limiting distributions of θˆ by the Slutsky's theorem (Athreya and Lahiri, 2006, Theorem 9.1.6); refer to Assumption 2 of Yuan and Jennrich (1998) and Theorem 6.1 of Rubin-Bleuer and Schiopu-Kratina (2005) for a similar assumption.

The first part of Condition C6 is a general condition for the bootstrap weights, which are determined only by a specified sampling design. The bootstrap weights also guarantee E∗{l ∗ <sup>w</sup>(θ)} = lw(θ) and V∗{l ∗ <sup>w</sup>(θ)} = Vˆ {lw(θ) | F<sup>N</sup> } conditional on the realized sample A for θ ∈ Θ. Besides, the non-negative condition on the bootstrap weights guarantees that l ∗ <sup>w</sup>(θ) is concave, which is used to show the consistency of the bootstrap estimator <sup>θ</sup><sup>ˆ</sup> ∗ . The second part of Condition C6 guarantees a design-consistent variance estimator, and it is satisfied under regularity conditions for commonly used sampling designs, including simple random sampling, probability proportional to size sampling, and stratified two-stage sampling. Condition C7 is the bootstrap central limit theorem for the bootstrap weighted score function Sˆ ∗ <sup>w</sup>(θ) for θ ∈ B. Different from Condition C4, the uniform central limit theorem in Condition C7 guarantees

$$
n^{1/2}\hat{\boldsymbol{S}}_{w}^{*}(\hat{\boldsymbol{\theta}}) \mid A \longrightarrow N(\mathbf{0}, \boldsymbol{\Sigma}_{S}(\boldsymbol{\theta}_{0})) \tag{S.1}
$$

in distribution with respect to the bootstrap sampling distribution conditional on the realized sample A, and it is a building block to prove Theorem 2. Condition C8 is a bootstrap version of Condition C5, and it is used to derive the limiting distribution of θˆ ∗ . To get the asymptotic distribution on the right-hand side of (S.1), we have used the convergence results in Condition C3 and Condition C6.

#### S2 Proof of Theorem 1

First, we show the consistency of the point estimator θˆ.

Lemma S1. Under Conditions C1–C2, the estimator θˆ is a consistent estimator of θ<sup>0</sup> with respect to the super-population model and the sampling mechanism.

Proof of Lemma S1. Since E{lw(θ) | F<sup>N</sup> } = l<sup>N</sup> (θ), it follows from Condition C2 that

$$
l_w(\boldsymbol{\theta}) - l_N(\boldsymbol{\theta}) \to 0 \tag{S.2}
$$

in probability with respect to the super-population model and the sampling mechanism as n → ∞ for θ ∈ Θ. Since lw(θ) is concave in Θ according to Condition C1, Lemma S1 can be proved by Theorem 2.1 and Theorem 2.7 of Engle and McFadden (1994, Chapter 36).

Proof of Theorem 1. We can adopt a similar proof for Theorem 4 of Yuan and Jennrich (1998), so its proof is omitted.

#### S3 Proof of Theorem 2

First, we show the consistency of the bootstrap estimator θˆ ∗ .

Lemma S2. Under Conditions C1–C2 and C6, the bootstrap point estimator θˆ ∗ is consistent for θ0.

Proof of Lemma S2. By Condition C1 and Condition C6, l ∗ <sup>w</sup>(θ) = N <sup>−</sup><sup>1</sup> P <sup>i</sup>∈<sup>A</sup> w ∗ i l(θ; yi) is concave and twice continuously differentiable. By Condition C2, E[Vˆ {lw(θ) | F<sup>N</sup> } | F<sup>N</sup> ] = V {lw(θ) | F<sup>N</sup> }, so E[Vˆ {lw(θ) | F<sup>N</sup> } | F<sup>N</sup> ] → 0 in probability with respect to the superpopulation model. By the Markov's inequality (Athreya and Lahiri, 2006, Theorem 3.1.1), we conclude that Vˆ {lw(θ) | F<sup>N</sup> } → 0 in probability with respect to the super-population model and the sampling mechanism by noting the fact that Vˆ {lw(θ) | F<sup>N</sup> } is non-negative. By Condition C2 and Condition C6, we have V∗{l ∗ <sup>w</sup>(θ)} → 0 in probability with respect to the bootstrap sampling distribution conditional on the realized sample A, and it leads to

$$
l_w^*(\boldsymbol{\theta}) - E_*\{l_w^*(\boldsymbol{\theta})\} = l_w^*(\boldsymbol{\theta}) - l_w(\boldsymbol{\theta}) \to 0
$$
\n(S.3)

in probability conditional on the realized sample A. By Condition C2 and (S.3), we have

$$
l_w^*(\boldsymbol{\theta}) - E\{l(\boldsymbol{\theta}; \boldsymbol{X}, Y)\} \to 0
$$
\n(S.4)

in probability with respect to the super-population model, the sampling mechanism as well as the bootstrap sampling distribution for θ ∈ Θ. Thus, by (S.4) and the fact that l ∗ <sup>w</sup>(θ) is concave, the convergence in (S.4) is uniform over any bounded and closed subset of Θ (Andersen and Gill, 1982, Theorem II.1). Thus, Lemma S2 can be proved by Theorem 2.1 of Engle and McFadden (1994, Chapter 36).

Proof of Theorem 2. By Lemmas S1–S2, we conclude that

$$
\hat{\boldsymbol{\theta}}^* - \hat{\boldsymbol{\theta}} \to 0 \tag{S.5}
$$

in probability with respect to the super-population model, the sampling mechanism as well as the bootstrap sampling distribution as n → ∞.

By the mean value theorem, we have

$$
\hat{\boldsymbol{S}}_{w}^{*}(\hat{\boldsymbol{\theta}}^{*}) = \hat{\boldsymbol{S}}_{w}^{*}(\hat{\boldsymbol{\theta}}) - \hat{\boldsymbol{I}}_{w}^{*}(\tilde{\boldsymbol{\theta}}^{*})(\hat{\boldsymbol{\theta}}^{*} - \hat{\boldsymbol{\theta}}),
$$
\n(S.6)

where ˆI ∗ <sup>w</sup>(θ) = <sup>−</sup>∂S<sup>ˆ</sup> ∗ <sup>w</sup>(θ)/∂θ T , and θ˜ ∗ lies on the segment joining θˆ and θˆ ∗ .

Now, if we can show that

$$
\hat{\boldsymbol{I}}_w^*(\tilde{\boldsymbol{\theta}}^*) \to \boldsymbol{\mathcal{I}}(\boldsymbol{\theta}_0) \tag{S.7}
$$

in probability with respect to the super-population model, the sampling mechanism as well as the bootstrap sampling distribution, then, by (S.6)–(S.7), we have

$$
\hat{\boldsymbol{S}}_{w}^{*}(\hat{\boldsymbol{\theta}}^{*}) = \hat{\boldsymbol{S}}_{w}^{*}(\hat{\boldsymbol{\theta}}) - \mathcal{I}(\boldsymbol{\theta}_{0})(\hat{\boldsymbol{\theta}} - \hat{\boldsymbol{\theta}}^{*}) + o_{p}(\|\hat{\boldsymbol{\theta}} - \hat{\boldsymbol{\theta}}^{*}\|),
$$
\n(S.8)

where the op-term is with respect to the super-population model, the sampling mechanism

as well as the bootstrap sampling distribution. Since Sˆ ∗ <sup>w</sup>(θ<sup>ˆ</sup> ∗ ) = 0, by (S.5) and (S.8), we have

$$
\hat{\boldsymbol{\theta}} - \hat{\boldsymbol{\theta}}^* = \boldsymbol{\mathcal{I}}^{-1}(\boldsymbol{\theta}_0) \hat{\boldsymbol{S}}_w^*(\hat{\boldsymbol{\theta}}) + o_p(1),
$$
\n(S.9)

where the op-term is with respect to the super-population model, the sampling mechanism as well as the bootstrap sampling distribution. Thus, by Conditions C4–C7, (S.9) and the fact that Sˆ <sup>w</sup>(θˆ) = 0, we have proved Theorem 2.

It remains to show (S.7). Since p is fixed, it is sufficient to show (S.7) for p = 1. For p = 1, we have

$$
|\hat{I}^*_w(\tilde{\theta}^*) - \mathcal{I}(\theta_0)| \le |\hat{I}^*_w(\tilde{\theta}^*) - \hat{I}_w(\tilde{\theta}^*)| + |\hat{I}_w(\tilde{\theta}^*) - \mathcal{I}(\tilde{\theta}^*)| + |\mathcal{I}(\tilde{\theta}^*) - \mathcal{I}(\theta_0)|. \tag{S.10}
$$

By Lemma S1, (S.5) and the fact that θ<sup>0</sup> is an interior point of B from Condition C3, we can show that the probability for ˜θ <sup>∗</sup> ∈ B converges to 1 with respect to the super-population model, the sampling mechanism as well as the bootstrap sampling distribution. Thus, by Condition C8, we have

$$
\hat{I}^*_{w}(\tilde{\theta}^*) - \hat{I}_{w}(\tilde{\theta}^*) \to 0
$$
\n(S.11)

in probability with respect to the bootstrap sampling distribution conditional on the realized sample A. Similarly, by Condition C5, we can also claim that

$$
\hat{I}_w(\tilde{\theta}^*) - \mathcal{I}(\tilde{\theta}^*) \to 0 \tag{S.12}
$$

in probability with respect to the super-population model and the sampling mechanism. By Lemma S1, (S.5) and the fact that ˜θ ∗ lies on the segment joining ˆθ and ˆθ ∗ , we conclude that ˜θ <sup>∗</sup> − θ<sup>0</sup> → 0 in probability with respect to the super-population model, the sampling mechanism as well as the bootstrap sampling distribution. Thus, by Condition C1, we can show that I(θ) is continuous on B, so

$$
\mathcal{I}(\tilde{\theta}^*) - \mathcal{I}(\theta_0) \to 0 \tag{S.13}
$$

in probability with respect to the super-population model, the sampling mechanism as well as the bootstrap sampling distribution. By (S.10)–(S.13), we have validated (S.7).

### S4 Stratified Multi-Stage Sampling

Krewski and Rao (1981) and Shao and Tu (1995, Chapter 6) discussed design limiting properties for stratified multi-stage sampling with clusters sampled with replacement conditional on the finite population, and Rubin-Bleuer and Schiopu-Kratina (2005) considered a superpopulation model for a stratified two-stage sampling design. In this section, we mainly extend the discussion of Krewski and Rao (1981) to a model-based framework and discuss the theoretical properties of parameters obtained by certain estimating equations; also see Example 3.1 and Example 6.1 of Rubin-Bleuer and Schiopu-Kratina (2005) for details.

We first consider a super-population model similar to Example 3.1 of Rubin-Bleuer and Schiopu-Kratina (2005) as Example 1 under stratified multi-stage sampling. Assume that the cluster size measures Mhi are generated from a super-population model for h = 1, . . . , H and i = 1, . . . , Nh, where Mhi is the size of the (hi)th cluster, containing elements in the second and sub-sequential stages, and N<sup>h</sup> is the number of the clusters (primary sampling units, PSU's) with respect to the first sampling design. Under stratified multi-stage sampling, let N = P<sup>H</sup> <sup>h</sup>=1 N<sup>h</sup> be the total number of clusters in the finite population. Conditional on the generated cluster sizes, we generate the finite population F<sup>N</sup> from a super-population model. We implicitly assume that H → ∞ and the number of stages are fixed for the sampling design. Besides, the population quantities are implicitly indexed by H. The following specific conditions are assumed for stratified multi-stage sampling; see Example 1

of the paper for details about the notations.

- A1. The parameter space Θ is an open and convex set containing θ<sup>0</sup> as an interior point. For θ ∈ Θ, l(θ; x, y) is concave and twice-continuously differentiable with respect to θ. There exists a non-stochastic function L(θ), such that l<sup>N</sup> (θ) → L(θ) in probability with respect to the super-population model for θ ∈ Θ. Besides, L(θ) is uniquely maximized at θ0.
- A2. Condition C2 holds.
- A3. max1≤h≤<sup>H</sup> n<sup>h</sup> = O(1) and n<sup>h</sup> ≥ 2 for h = 1, 2, . . ..
- A4. max1≤h≤<sup>H</sup> W<sup>h</sup> = Op(H<sup>−</sup><sup>1</sup> ) with respect to the super-population model.
- A5. There exists a δ > 0, such that P<sup>H</sup> <sup>h</sup>=1 <sup>W</sup>hE{∥S˜ w,h1(θ) − Sh(θ)∥ 2+δ | F<sup>N</sup> } = Op(1) with respect to the super-population model uniformly for θ ∈ B, where B is a compact set containing θ<sup>0</sup> as an interior point, S˜ w,h1(θ) = da(h1)p −1 a(h1)S<sup>ˆ</sup> w,a(h1)(θ), Sh(θ) = M<sup>−</sup><sup>1</sup> h P<sup>N</sup><sup>h</sup> <sup>i</sup>=1 MhiShi(θ) and Shi(θ) are the mean of the score functions of the h-th stratum and (hi)-th cluster, respectively.
- A6. supθ∈B n P<sup>H</sup> <sup>h</sup>=1 W<sup>2</sup> <sup>h</sup><sup>V</sup> {S˜ w,a(h1)(θ) | F<sup>N</sup> }/n<sup>h</sup> − ΣS(θ) → 0 in probability with respect to the super-population model, where ΣS(θ) is non-stochastic and positive definitive for θ ∈ B, and n = P<sup>H</sup> <sup>h</sup>=1 n<sup>h</sup> is the number of sampled clusters.
- A7. supθ∈B V {S(θ)} = o(n −1 ) with respect to the super-population model, where S(θ) = E{Sˆ <sup>w</sup>(θ) | F<sup>N</sup> }.
- A8. There exist two constants 0 < C<sup>5</sup> < C<sup>6</sup> < ∞, such that C<sup>5</sup> < Mhphi < C<sup>6</sup> almost surely with respect to the super-population model for h = 1, . . . , H and i = 1, . . . , Nh.
- A9. max1≤h≤<sup>H</sup> supθ∈B E{∥ˆIw,a(h1)(θ)∥ 2 | F<sup>N</sup> } = Op(1) with respect to the super-population model, where ˆIw,hi(θ) = ∂Sˆ w,hi(θ)/∂θ T .

A10. There exists a non-stochastic function I(θ) such that supθ∈B∥I(θ) − I(θ)∥ → 0 in probability with respect to the super-population model and the sampling mechanism, where I(θ) = E{ ˆIw(θ) | F<sup>N</sup> }, and I(θ0) is invertible.

Condition A1 is a counterpart of Condition C1. Since cluster-level random effects are allowed when generating finite populations under stratified multistage sampling, the elements in the finite population are no longer realizations of IID random variables. Thus, instead of E{l(θ; X, Y )}, we assume the existence of a non-stochastic function L(θ) in Condition A1. By Conditions A1–A2, we can validate the consistency of θˆ and θˆ ∗ under stratified multistage sampling; see Lemma S4 for details. Although we may have more than two sampling stages, we mainly focus on the PSU's in the first stage of sampling. Condition A3 shows that we focus on the stratified multi-stage sampling where the number of selected PSUs is bounded within each stratum, but we assume that the number of strata diverges; see Krewski and Rao (1981) for details. Condition A4 rules out the existence of extremely large or small strata. Condition A5 is a typical Lyapounov-type condition, and it is used to validate the central limit theorem in Condition C4. Condition A6 shows the convergence rate of the variance of the estimation equation. Condition A7 shows that V {S(θ)} is negligible compared with the design variance of Sˆ <sup>w</sup>(θ) conditional on the finite population. Conditions A6–A7 are mainly used to verify Condition C4; see Section S1 for more detailed discussion. Condition A8 is used to regulate the selection probability of each PSU. Conditions A9–A10 are used to show the convergence property of the negative Hessian matrix.

Remark S1. As mentioned above, we do not assume that the finite population is a random sample from a super-population model in this section. To make this point clear, consider Sˆ <sup>w</sup>(θ0) = Sˆ <sup>w</sup>(θ0) − S(θ0) + S(θ0) − E{Sˆ <sup>w</sup>(θ0)}. Conditional on a finite population, we can show that {S˜ w,hi(θ) : i = 1, . . . , nh} are IID for each h = 1, . . . , H, and S˜ w,h1<sup>i</sup>(θ) are independent of S˜ w,h2<sup>j</sup> (θ) for h<sup>1</sup> ̸= h<sup>2</sup> under stratified multi-stage sampling with clusters selected with replacement; we refer to a similar discussion on page 1012 of Krewski and Rao (1981). It may be noted that such independence is only due to the sampling design given the finite population. Thus, under Condition A1–A6, we can obtain a similar result as Theorem 3.1 of Krewski and Rao (1981) for the estimating function in probability. Under Condition A7, we can show that S(θ0) − E{Sˆ <sup>w</sup>(θ0)} is asymptotically negligible compared with Sˆ <sup>w</sup>(θ0)−S(θ0). Thus, we can establish the central limit theorem in Condition C4 under a more generate setup for the super-population model, and this is also the intuition to prove Lemma S5 below.

To make the discussion more clear, we consider a stratified two-stage sampling design similar as Example 6.1 of Rubin-Bleuer and Schiopu-Kratina (2005). Let the finite population be F<sup>N</sup> = {yhij : h = 1, . . . , H;i = 1, . . . , Nh; j = 1, . . . , Mhi}. For h = 1, . . . , H and i = 1, . . . , Nh, {yhij : j = 1, . . . , Mhi} are IID conditional on the cluster; also see Rabe-Hesketh and Skrondal (2006) for a similar model assumption under multilevel modeling. Let l(θ; y) be a concave and twice-differentiable objective function. Then, l<sup>N</sup> (θ) = M<sup>−</sup><sup>1</sup> P<sup>H</sup> h=1 P<sup>N</sup><sup>h</sup> i=1 P<sup>M</sup>hi <sup>j</sup>=1 l(θ; yhij ). Due to the independence among sample clusters, Condition A1 can be satisfied under regularity conditions on the super-population model. Conditions A5–A6 are used to show the central limit theorem of the weighted score function; refer to Condition 5 of Rubin-Bleuer and Schiopu-Kratina (2005) and Lemma 3.1 of Krewski and Rao (1981) for details. By a similar argument of Example 6.1 of Rubin-Bleuer and Schiopu-Kratina (2005), we can show that V {S<sup>N</sup> (θ)} = O(N <sup>−</sup><sup>1</sup> ) for θ ∈ Θ. Thus, if n = o(N), Condition A7 can be guaranteed under certain moment conditions. Recall that N is the total number of population clusters, and n is the number of sampled PSUs. The discussion for Conditions A9–A10 is the same as the preceding paragraphs.

Before validating the regular conditions, we need the following result.

Lemma S3. Consider a sequence of m×m real-valued and symmetric matrices, M1,M2, . . .. If for any a ∈ R <sup>m</sup> satisfying ∥a∥ = 1,

$$
\mathbf{a}^{\mathrm{T}}\mathbf{M}_n\mathbf{a}\to 0, \quad n\to\infty,\tag{S.14}
$$

then we have

$$
||M_n|| \to 0, \quad n \to \infty.
$$

Proof of Lemma S3. By the definition of the operator norm at the very beginning of the Supplementary Material, we can conclude that

$$
\|\mathbf{M}\| = \max\{|\lambda| : \lambda \text{ is an eigenvalue of } \mathbf{M}\}\
$$
\n(S.15)

for any real-valued and symmetric matrix M; see Section 5.6 of Horn and Johnson (2012) for details. By (S.14), all eigenvalues of M<sup>n</sup> converges to 0 as n → ∞; refer to the Rayleigh-Ritz Theorem (Horn and Johnson, 2012, Section 4.2) for details. Thus, we have proved Lemma S3 by (S.15).

Lemma S4. Under stratified multi-stage sampling and Conditions A1–A2, the point estimator θˆ is consistent of θ<sup>0</sup> with respect to the super-population model and the sampling mechanism, and the bootstrap point estimator θˆ ∗ is consistent for θ<sup>0</sup> as well.

The proof of Lemma S4 is essentially the same as Lemmas S1–S2 by noting the fact that the bootstrap variance is unbiased for the variance estimator of the original sample; see (S.23) in Lemma S7 below for details. Thus, its proof is omitted.

Lemma S5. Under stratified multi-stage sampling and Condition A1–A7, Condition C4 holds.

Proof of Lemma S5. By Conditions A3–A5, Krewski and Rao (1981) established the following central limit theorem:

$$
n^{1/2}\{\hat{\mathbf{S}}_{w}(\boldsymbol{\theta}_{0}) - \mathbf{S}(\boldsymbol{\theta}_{0})\} \mid \mathcal{F}_{N} \to N(\mathbf{0}, \Sigma_{S}(\boldsymbol{\theta}_{0}))
$$
\n(S.16)

in distribution with respect to the sampling mechanism in probability regarding the superpopulation model, where S(θ) = E{Sˆ <sup>w</sup>(θ) | F<sup>N</sup> }.

Besides, by Condition A7, we have

$$
\mathbf{S}(\boldsymbol{\theta}_0) = o_p(n^{-1}) \tag{S.17}
$$

with respect to the super-population model, since E{S(θ0)} = 0 by Condition C1.

By Condition A6 and (S.16)–(S.17), Theorem 5.1 of Rubin-Bleuer and Schiopu-Kratina (2005) validates Lemma S5.

Lemma S6. Under stratified multi-stage sampling, Conditions A3–A4 and Conditions A9– A10, Condition C5 holds.

Proof of Lemma S6. By Condition A10, it is sufficient to show

$$
\sup_{\theta \in \mathcal{B}} \|\hat{\boldsymbol{I}}_{w}(\boldsymbol{\theta}) - \boldsymbol{I}(\boldsymbol{\theta})\| \to 0 \tag{S.18}
$$

in probability with respect to the super-population and the sampling mechanism. Since

$$
\mathbf{I}(\boldsymbol{\theta}) = E\{\hat{\boldsymbol{I}}_{w}(\boldsymbol{\theta}) \mid \mathcal{F}_{N}\},\tag{S.19}
$$

by Lemma S3, it remains to show

$$
\sup_{\theta \in \mathcal{B}} V\{\mathbf{a}^{\mathrm{T}}\hat{\boldsymbol{I}}_{w}(\theta)\mathbf{a} \mid \mathcal{F}_{N}\} \to 0
$$
\n(S.20)

in probability with respect to the super-population model for a ∈ R p satisfying ∥a∥ = 1.

Consider

$$
V\{\boldsymbol{a}^{\mathrm{T}}\hat{\boldsymbol{I}}_{w}(\boldsymbol{\theta})\boldsymbol{a} \mid \mathcal{F}_{N}\} = \sum_{h=1}^{H} W_{h}^{2} m_{h}^{-1} V\{\boldsymbol{a}^{\mathrm{T}}\hat{\boldsymbol{I}}_{w,a(h1)}(\boldsymbol{\theta})\boldsymbol{a} \mid \mathcal{F}_{N}\}
$$
  
\n
$$
= O(H^{-1}) \sum_{h=1}^{H} W_{h} m_{h}^{-1} V\{\boldsymbol{a}^{\mathrm{T}}\hat{\boldsymbol{I}}_{w,a(h1)}(\boldsymbol{\theta})\boldsymbol{a} \mid \mathcal{F}_{N}\}
$$
  
\n
$$
= O_{p}(H^{-1})
$$
\n(S.21)

with respect to the super-population model uniformly over B and a ∈ R p satisfying ∥a∥ = 1, where the second equality holds by Condition A4, and the last equality holds by Conditions A3–A4 and Condition A9. Thus, we have shown (S.20) by (S.21) and proved Lemma S6.

Lemma S7. The proposed bootstrap satisfies the first part of Condition C6, involving the bootstrap expectation and variance.

Proof of Lemma S7. Recall that Sˆ <sup>w</sup>(θ) = P<sup>H</sup> <sup>h</sup>=1 <sup>W</sup>hS<sup>ˆ</sup> w,h(θ), where W<sup>h</sup> = Mh/M, M = P<sup>H</sup> <sup>h</sup>=1 Mh, M<sup>h</sup> = P<sup>N</sup><sup>h</sup> <sup>i</sup>=1 <sup>M</sup>hi, <sup>M</sup>hi is the size of the (hi)-th cluster, <sup>S</sup><sup>ˆ</sup> w,h(θ) = n −1 h P<sup>n</sup><sup>h</sup> <sup>i</sup>=1 <sup>S</sup>˜ w,hi(θ), S˜ w,hi(θ) = da(hi)p −1 a(hi)S<sup>ˆ</sup> w,a(hi)(θ), da(hi) = Ma(hi)/Mh, and a(hi) is the index of the i-th selected cluster in the h-th stratum. From Example 1, the bootstrap replicate of Sˆ <sup>w</sup>(θ) is Sˆ ∗ <sup>w</sup>(θ) = P<sup>H</sup> <sup>h</sup>=1 <sup>W</sup>hS<sup>ˆ</sup> ∗ w,h(θ), where <sup>S</sup><sup>ˆ</sup> ∗ w,h(θ) = n −1 h P<sup>n</sup><sup>h</sup> <sup>i</sup>=1 r ∗ hiS˜ w,hi(θ), r ∗ hi = khn ∗ hi, k<sup>h</sup> = nh/(n<sup>h</sup> − 1), and n ∗ <sup>h</sup> = (n ∗ h1 , . . . , n<sup>∗</sup> hn<sup>h</sup> ) T is generated by a multinomial distribution with n<sup>h</sup> − 1 trials and a success probability vector n −1 h (1, . . . , 1)<sup>T</sup> of length nh. Since E∗(r ∗ hi) = 1 for h = 1, . . . , H and i = 1, . . . , mh, we have

$$
E_*\{\hat{\boldsymbol{S}}_w^*(\boldsymbol{\theta})\} = \hat{\boldsymbol{S}}_w(\boldsymbol{\theta}).
$$
\n(S.22)

By properties of a multinomial distribution, we have

$$
\mathbf{V}_{*}\{\hat{\boldsymbol{S}}_{w}^{*}(\boldsymbol{\theta})\} = \sum_{h=1}^{H} W_{h}^{2} \{n_{h}(n_{h}-1)\}^{-1} \tilde{\boldsymbol{S}}_{h,vec}(\boldsymbol{\theta}) (\boldsymbol{I}_{n_{h}} - \mathcal{P}_{1,n_{h}}) \tilde{\boldsymbol{S}}_{h,vec}^{T}(\boldsymbol{\theta})
$$
  
=  $\hat{V} \{\hat{\boldsymbol{S}}_{w}(\boldsymbol{\theta}) \mid \mathcal{F}_{N}\},$  (S.23)

where P1,n<sup>h</sup> = 1<sup>n</sup><sup>h</sup> (1 T n<sup>h</sup> 1<sup>n</sup><sup>h</sup> ) <sup>−</sup><sup>1</sup>1 T n<sup>h</sup> is the projection matrix to the linear space spanned by 1<sup>n</sup><sup>h</sup> and S˜ h,vec(θ) = (da(h1)p −1 a(h1)S<sup>ˆ</sup> w,a(h1)(θ), . . . , da(hnh)p −1 <sup>a</sup>(hnh)S<sup>ˆ</sup> w,a(hnh)(θ)).

By (S.22)–(S.23), we have completed the proof of Lemma S7.

The second part of Condition C6, involving design variance estimator, can be validated in a similar manner as Theorem 3.2 of Krewski and Rao (1981) based on Conditions A3–A5, so we omit the proof here.

Lemma S8. Under stratified multi-stage sampling and Conditions A3–A6, Condition C7 holds.

Proof of Lemma S8. We prove Lemma S8 by validating the Lyapounov condition for Sˆ ∗ <sup>w</sup>(θ) with respect to the super-population model, the sampling mechanism and the bootstrap sampling distribution. Notice that Sˆ ∗ w,h(θ) = n −1 h kh P<sup>n</sup>h−<sup>1</sup> <sup>i</sup>=1 <sup>S</sup>˜ ∗ w,hi(θ), where {S˜ ∗ w,h1 (θ), . . . ,S˜ ∗ w,h(nh−1)(θ)} is a random sample of size nh−1 from {S˜ w,h1(θ), . . . ,S˜ w,hn<sup>h</sup> (θ)}. Thus, we have E∗{S˜ ∗ w,h1 (θ)} = Sˆ w,h(θ), where Sˆ w,h(θ) = n −1 h P<sup>n</sup><sup>h</sup> <sup>i</sup>=1 <sup>S</sup>˜ w,hi(θ). For δ > 0 in Condition A5, consider

$$
\sum_{h=1}^{H} W_h^{2+\delta} E_* \{ \|\tilde{\boldsymbol{S}}_{w,h1}^* (\boldsymbol{\theta}) - \hat{\boldsymbol{S}}_{w,h} (\boldsymbol{\theta}) \|^{2+\delta} \}
$$
\n
$$
= \sum_{h=1}^{H} W_h^{2+\delta} \frac{1}{n_h} \sum_{i=1}^{n_h} \|\tilde{\boldsymbol{S}}_{w,hi} (\boldsymbol{\theta}) - \hat{\boldsymbol{S}}_{w,h} (\boldsymbol{\theta}) \|^{2+\delta}
$$
\n
$$
\leq C(\delta) \left\{ \sum_{h=1}^{H} W_h^{2+\delta} \frac{1}{n_h} \sum_{i=1}^{n_h} \|\tilde{\boldsymbol{S}}_{w,hi} (\boldsymbol{\theta}) - \boldsymbol{S}_h (\boldsymbol{\theta}) \|^{2+\delta} + \sum_{h=1}^{H} W_h^{2+\delta} \|\hat{\boldsymbol{S}}_{w,h} (\boldsymbol{\theta}) - \boldsymbol{S}_h (\boldsymbol{\theta}) \|^{2+\delta} \right\}
$$
\n(S.24)

by the triangular inequality and the H¨older's inequality, where C(δ) is a constant determined by δ.

Consider the first term of (S.24),

$$
E\left\{\sum_{h=1}^{H}W_h^{2+\delta}\frac{1}{n_h}\sum_{i=1}^{n_h} \|\tilde{\boldsymbol{S}}_{w,hi}(\boldsymbol{\theta})-\boldsymbol{S}_h(\boldsymbol{\theta})\|^{2+\delta} \|\mathcal{F}_N\right\} = \sum_{h=1}^{H}W_h^{2+\delta}E\{\|\tilde{\boldsymbol{S}}_{w,hi}(\boldsymbol{\theta})-\boldsymbol{S}_h(\boldsymbol{\theta})\|^{2+\delta} \|\mathcal{F}_N\}
$$
  
=  $O_p(H^{-1-\delta})$  (S.25)

with respect to the super-population model and the sampling mechanism, where the first equality holds by the fact that we use probability proportional to size sampling with replacement in the first stage sampling, and the second equality holds by Conditions A4–A5.

Before considering the second term of (S.24), we have

$$
\begin{array}{rcl} \|\hat{\bm{S}}_{w,h}(\bm{\theta})-\bm{S}_h(\bm{\theta})\|^{2+\delta} & = & \left\|\frac{1}{n_h}\sum_{i=1}^{n_h}\{\tilde{\bm{S}}_{w,hi}(\bm{\theta})-\bm{S}_h(\bm{\theta})\}\right\| \\ \\ & \leq & \left\{\frac{1}{n_h}\sum_{i=1}^{n_h}\left\|\tilde{\bm{S}}_{w,hi}(\bm{\theta})-\bm{S}_h(\bm{\theta})\right\|\right\}^{2+\delta} \\ \\ & \leq & \frac{1}{n_h}\sum_{i=1}^{n_h}\left\|\tilde{\bm{S}}_{w,hi}(\bm{\theta})-\bm{S}_h(\bm{\theta})\right\|^{2+\delta}, \end{array}
$$

where the second inequality holds by the triangular inequality for norms, and the third inequality holds by the H¨older's inequality. Thus, we can use a similar argument to (S.25) to show that

$$
E\left\{\sum_{h=1}^{H} W_h^{2+\delta} \|\hat{\boldsymbol{S}}_{w,h}(\boldsymbol{\theta}) - \boldsymbol{S}_h(\boldsymbol{\theta})\|^{2+\delta} \mid \mathcal{F}_N\right\} = O_p(H^{-1-\delta})
$$
(S.26)

with respect to the super-population model and the sampling mechanism. Thus, by (S.24)– (S.26) and the Markov's inequality, we can show that

$$
\sum_{h=1}^{H} W_h^{2+\delta} E_* \{ \|\tilde{\boldsymbol{S}}_{h1}^*(\boldsymbol{\theta}) - \hat{\boldsymbol{S}}_{w,h}(\boldsymbol{\theta})\|^{2+\delta} \} = O_p(H^{-1-\delta})
$$
\n(S.27)

with respect to the super-population model and the sampling mechanism. On the other hand, by Lemma S7, Condition A3 and Condition A6, we have shown that

$$
\boldsymbol{V}_* \{\hat{\boldsymbol{S}}_w^*(\boldsymbol{\theta})\} \asymp H^{-1} \tag{S.28}
$$

in probability with respect to the super-population model and the sampling mechanism. Thus, for any a ∈ R <sup>p</sup> with ∥a∥ = 1, we have

$$
[V_*\{\boldsymbol{a}^{\mathrm{T}}\hat{\boldsymbol{S}}_w^*(\boldsymbol{\theta})\}]^{-(2+\delta)/2} \sum_{h=1}^H W_h^{2+\delta} E_*\{\|\tilde{\boldsymbol{S}}_h^*(\boldsymbol{\theta}) - \hat{\boldsymbol{S}}_{w,h}(\boldsymbol{\theta})\|^{2+\delta}\} = O_p(H^{-\delta/2}) = o_p(1) \quad (S.29)
$$

with respect to the super-population model and the sampling mechanism, where the last

equality holds by the fact that H → ∞. Thus, by (S.29) and the Cram´er-Wold device (Athreya and Lahiri, 2006, Theorem 10.4.5), we have completed the proof of Lemma S8.

Lemma S9. Under stratified multi-stage sampling and Conditions A3–A4 and Condition A9, Condition C8 holds.

Proof of Lemma S9. Recall that ˆIw(θ) = P<sup>H</sup> <sup>h</sup>=1 W<sup>h</sup> ˆIw,h(θ), where ˆIw,h(θ) = n −1 h P<sup>n</sup><sup>h</sup> i=1 ˜Iw,hi(θ), ˜Iw,hi(θ) = da(hi)p −1 a(hi) ˆIw,a(hi)(θ). The bootstrap replicate of ˆIw(θ) is ˆI ∗ <sup>w</sup>(θ) = P<sup>H</sup> <sup>h</sup>=1 W<sup>h</sup> ˆI ∗ w,h(θ), where ˆI ∗ w,h(θ) = n −1 h P<sup>n</sup><sup>h</sup> <sup>i</sup>=1 r ∗ hi ˜Iw,hi(θ).

First, we can use a similar argument to (S.22) to show that

$$
E_*\{\hat{\boldsymbol{I}}_w^*(\boldsymbol{\theta})\} = \hat{\boldsymbol{I}}_w(\boldsymbol{\theta}).
$$
\n(S.30)

For a ∈ R <sup>p</sup> with ∥a∥ = 1, denote Za,hi(θ) = a <sup>T</sup> ˜Iw,hi(θ)a. Besides, we can show that

$$
\mathbf{V}_{*}\{\boldsymbol{a}^{T}\hat{\boldsymbol{I}}_{w}^{*}(\boldsymbol{\theta})\boldsymbol{a}\} = \sum_{h=1}^{H} W_{h}^{2}\{n_{h}(n_{h}-1)\}^{-1}\tilde{\boldsymbol{Z}}_{h,vec}(\boldsymbol{\theta})(\boldsymbol{I}_{n_{h}} - \mathcal{P}_{1,n_{h}})\tilde{\boldsymbol{Z}}_{h,vec}^{T}(\boldsymbol{\theta})
$$
\n
$$
\leq \sum_{h=1}^{H} W_{h}^{2}\{n_{h}(n_{h}-1)\}^{-1}\sum_{i=1}^{n_{h}} Z_{\boldsymbol{a},hi}^{2}(\boldsymbol{\theta}), \qquad (S.31)
$$

where Z˜ h,vec(θ) = (Za,h1(θ), . . . , Za,hn<sup>h</sup> (θ)), and the second inequality holds since P1,n<sup>h</sup> is non-negative definitive.

Consider

$$
E\left[\sum_{h=1}^{H} W_h^2\{n_h(n_h-1)\}^{-1} \sum_{i=1}^{n_h} Z_{a,hi}^2(\boldsymbol{\theta}) \mid \mathcal{F}_N\right] \leq \sum_{h=1}^{H} W_h^2 E\{Z_{a,hi}^2(\boldsymbol{\theta}) \mid \mathcal{F}_N\}
$$
  
=  $O_p(H^{-1})$  (S.32)

with respect to the super-population model uniformly over B, where the first inequality holds by Condition A3, and the second inequality holds by Condition A4 and Condition A9.

By (S.31)–(S.32) to show

$$
\boldsymbol{V}_*\{\boldsymbol{a}^{\mathrm{T}}\hat{\boldsymbol{I}}_w^*(\boldsymbol{\theta})\boldsymbol{a}\} = o_p(1)
$$

with respect to the super-population model and the sampling mechanism uniformly for θ ∈ B by Markov's inequality. Thus, by Lemma S3, (S.30) and (S.33), we have completed the proof of Lemma S9.

### S5 Poisson sampling

Poisson sampling is often applied when administrative files are available, and it can also simplify the variance estimator under some complex sampling designs; see Beaumont and Patak (2012) for details. Under Poisson sampling, the sample selection is made from independent Bernoulli distributions with parameter π<sup>i</sup> for i = 1, . . . , N. Specifically, let I<sup>i</sup> , following Bernoulli(πi) distribution, be the sampling indicator of y<sup>i</sup> , and I<sup>i</sup> = 1 if and only if yi is selected. The rescaling factor is r ∗ <sup>i</sup> = 1+m<sup>∗</sup> <sup>i</sup> −π<sup>i</sup> , where m<sup>∗</sup> <sup>i</sup> ∼ Bernoulli(πi). The bootstrap replicate of Sˆ <sup>w</sup>(θ) is Sˆ ∗ <sup>w</sup>(θ) = N <sup>−</sup><sup>1</sup> P i∈A r ∗ <sup>i</sup> wiS(θ; yi) with w<sup>i</sup> = π −1 i , and it satisfies E∗{Sˆ ∗ <sup>w</sup>(θ)} <sup>=</sup> <sup>S</sup><sup>ˆ</sup> <sup>w</sup>(θ) and V∗{Sˆ ∗ <sup>w</sup>(θ)} = N <sup>−</sup><sup>2</sup> P i∈A π −1 i (1 − πi)S(θ; yi) ⊗2 , where x <sup>⊗</sup><sup>2</sup> = xx<sup>T</sup> for a vector x.

## S6 Probability proportional to size sampling with replacement

Probability proportional to size sampling is widely used as a building block for stratified multi-stage sampling. Under probability proportional to size sampling with replacement, a sample of size n is independently generated from F<sup>N</sup> with replacement, and the selection probability of y<sup>i</sup> is p<sup>i</sup> with P<sup>N</sup> <sup>i</sup>=1 <sup>p</sup><sup>i</sup> = 1. Then, we have <sup>S</sup><sup>ˆ</sup> <sup>w</sup>(θ) = (Nn) <sup>−</sup><sup>1</sup> P<sup>n</sup> <sup>i</sup>=1 p −1 <sup>a</sup>(i)S(θ; ya(i)), where a(i) is the index of the element selected at the i-th draw. A design-unbiased variance estimator is Vˆ {Sˆ <sup>w</sup>(θ)} = {N<sup>2</sup>n(n − 1)} <sup>−</sup><sup>1</sup>Svec(θ)(I<sup>n</sup> − P1,n)S T vec(θ), where Svec(θ) = (p −1 <sup>a</sup>(1)S(θ; ya(1)), . . . , p<sup>−</sup><sup>1</sup> <sup>a</sup>(n)S(θ; ya(n))) is a p × n matrix, I<sup>n</sup> is the n × n identity matrix, and P1,n = 1n(1 T <sup>n</sup>1n) <sup>−</sup><sup>1</sup>1 T n is the projection matrix to the linear space spanned by 1n.

Then, the bootstrap replicate of Sˆ <sup>w</sup>(θ) is Sˆ ∗ <sup>w</sup>(θ) = (Nn) <sup>−</sup><sup>1</sup> P<sup>n</sup> <sup>i</sup>=1 r ∗ i p −1 <sup>a</sup>(i)S(θ; ya(i)) with the rescaling factor is r ∗ <sup>i</sup> = knm<sup>∗</sup> i , where m<sup>∗</sup> = (m<sup>∗</sup> 1 , . . . , m<sup>∗</sup> n ) T is generated using a multinomial distribution with n − 1 trials and a success probability vector n <sup>−</sup><sup>1</sup>1n, and k<sup>n</sup> = n/(n − 1). Using the property of the multinomial distribution, we can verify that E∗{Sˆ ∗ <sup>w</sup>(θ)} <sup>=</sup> <sup>S</sup><sup>ˆ</sup> <sup>w</sup>(θ) and V <sup>∗</sup>{Sˆ ∗ <sup>w</sup>(θ)} <sup>=</sup> <sup>V</sup><sup>ˆ</sup> {S<sup>ˆ</sup> <sup>w</sup>(θ) | F<sup>N</sup> }. Under generality conditions, the proposed bootstrap method can be validated in a similar manner as in Section S4.

#### S7 Stratified Simple Random Sampling

In this section, we propose a bootstrap method for stratified random sampling without replacement, which is commonly used and has been investigated by Bickel and Freedman (1984). Specifically, denote N<sup>h</sup> and n<sup>h</sup> to be the stratum population size and sample size of the hth stratum, respectively. Let N = P<sup>H</sup> <sup>h</sup> N<sup>h</sup> and n = P<sup>H</sup> <sup>h</sup>=1 n<sup>h</sup> be the population size and sample size, where H is the number of strata. Simple random sampling without replacement is applied within each stratum, and we assume that the first n<sup>h</sup> elements are sampled within the hth stratum, for simplicity. Let H, N<sup>h</sup> and n<sup>h</sup> depend on an index ν such that n(ν) = n1(ν) + · · · + nH(ν) → ∞ as ν → ∞, and we omit the index ν without loss of generality. Bickel and Freedman (1984) assumed Lindeberg conditions to get the limiting results.

Example 3. Under stratified simple random sampling without replacement, the estimator θˆ is obtained by solving

$$
\hat{\boldsymbol{S}}_w(\boldsymbol{\theta}) = \sum_{h=1}^H W_h N_h^{-1} \sum_{i=1}^{n_h} w_{hi} \boldsymbol{S}(\boldsymbol{\theta}; \boldsymbol{x}_{hi}, y_{hi}) = \sum_{h=1}^H W_h n_h^{-1} \sum_{i=1}^{n_h} \boldsymbol{S}(\boldsymbol{\theta}; \boldsymbol{x}_{hi}, y_{hi}) = \boldsymbol{0},
$$

where W<sup>h</sup> = NhN <sup>−</sup><sup>1</sup> , whi = Nhn −1 h . Then, we can show that

$$
\hat{\boldsymbol{V}}\{\hat{\boldsymbol{S}}_w(\boldsymbol{\theta})\mid\mathcal{F}_N\}=\sum_{h=1}^HW_h^2n_h^{-1}(1-n_hN_h^{-1})(n_h-1)^{-1}\sum_{i=1}^{n_h}\{\boldsymbol{S}(\boldsymbol{\theta};\boldsymbol{x}_{hi},y_{hi})-\bar{\boldsymbol{S}}_h(\boldsymbol{\theta})\}^{\otimes2},
$$

where S¯ <sup>h</sup>(θ) = n −1 h P<sup>n</sup><sup>h</sup> <sup>i</sup>=1 S(θ; xhi, yhi).

The rescaling factor of the bootstrap method is r ∗ hi = 1 + {m<sup>∗</sup> hi − (n<sup>h</sup> − 1)n −1 h }kh, where m<sup>∗</sup> <sup>h</sup> = (m<sup>∗</sup> h1 , . . . , m<sup>∗</sup> hn<sup>h</sup> ) T is generated using a multinomial distribution with n<sup>h</sup> − 1 trials and a success probability vector n −1 h (1, . . . , 1)<sup>T</sup> of length nh, and k 2 <sup>h</sup> = n 2 h (1 − nhN −1 h )(n<sup>h</sup> − 1)<sup>−</sup><sup>2</sup> .

The proposed bootstrap method is essentially the same as the one considered in Section 4 of Rao and Wu (1988) when estimating the population mean.

#### S8 Design effect

#### S8.1 Tests for a scalar parameter

For the first simulation and the additional simulation study in Section S9 below, we are interested in testing hypotheses with respect to a single regression parameter θ2. The design effect for estimating θ<sup>2</sup> is

$$
\text{deff}(\theta_2) = V_p(\hat{\theta}_{2,p})/V_{srs}(\hat{\theta}_{2,srs}),
$$

where ˆθ2,p and ˆθ2,srs are the point estimators of θ<sup>2</sup> under the current sampling design and simple random sampling with replacement with the same sample size, respectively, and Vp( ˆθ2,p) and Vsrs( ˆθ2,srs) are their variances with respect to the specific sampling design as well as the super-population model. In this paper, we use 50 000 Monte Carlo simulations to estimate Vp( ˆθ2,p) and Vsrs( ˆθ2,srs).

#### S8.2 Tests of independence in a two-way table

In the simulation study, we consider a 3 × 3 table. For i = 1, 2 and j = 1, 2, let hij (p) = pij − pi+p+<sup>j</sup> , where p = (p11, p12, p13, . . . , p31, p32) T , pi<sup>+</sup> = P<sup>3</sup> <sup>j</sup>=1 pij for i = 1, 2 and p+<sup>j</sup> = P<sup>3</sup> <sup>i</sup>=1 pij for j = 1, 2. Then, denote h(p) = (h11(p), h12(p), h21(p), h22(p))<sup>T</sup> , and let V <sup>h</sup> be the variance of h(pˆ) with respect to the sampling design as well as the super-population model, where pˆ is a design-unbiased estimator of p. Then, the asymptotic distribution of the test statistic X<sup>2</sup> I in (18) can be approximated by P<sup>4</sup> <sup>k</sup>=1 δkWk, where δ1, . . . , δ<sup>4</sup> are the eigenvalues of (P −1 <sup>r</sup> ⊗P −1 c )V <sup>h</sup>, P <sup>r</sup> = diag(p<sup>r</sup> )−prp T r , P <sup>c</sup> = diag(p<sup>c</sup> )−pcp T c , p<sup>r</sup> = (p1+, p2+) T , p<sup>c</sup> = (p+1, p+2) T , and W1, . . . , W<sup>4</sup> are independent χ 2 (1) random variables.

We use 50 000 Monte Carlo simulations to estimate V <sup>h</sup>, and it can be shown that all the estimated eigenvalues are greater than 1.3 under both scenarios. Thus, the asymptotic distribution of the test statistic X<sup>2</sup> I in (18) is more dispersed than a χ 2 (4) distribution, which is the reference distribution for the naive Pearson method.

# S9 Additional Simulation Study: Stratified Two-Stage Sampling

In this section, we consider a stratified two-stage sampling design to study the performance of the proposed bootstrap testing method. A finite population F<sup>N</sup> = {(xhij , yhij ) : h = 1, . . . , H;i = 1, . . . , Mh; j = 1, . . . , Mhi} is generated based on the following steps, where H is the number of strata, M<sup>h</sup> is the number of clusters in the h-th stratum, and Mhi is the size of the (hi)-th cluster. Let a<sup>h</sup> ∼ Ex(1) be the stratum effect, bhi ∼ Ex(1) be the cluster effect, M<sup>h</sup> | a<sup>h</sup> ∼ 5Po(ah) + 20, Mhi ∼ 5Po(a<sup>h</sup> + bhi) + 30, xhij ∼ U(0, 10), and yhij | (ah, bhi, xhij ) ∼ N(µhij , 1), where Ex(λ) is an exponential distribution with rate parameter λ, Po(λ) is a Poisson distribution with parameter λ, µij = θ<sup>1</sup> + θ2xhij + 0.3a<sup>h</sup> + 0.3bhi, and (θ1, θ2) = (0, 1). From the finite population, we use a stratified two-stage sampling design to obtain a sample. The first-stage sampling design is probability proportional to size sampling with replacement, where the selection probability is proportional to the cluster size, and simple random sampling without replacement is applied at the second stage. We consider two scenarios for the number of strata: H = 20 and H = 50. For each scenario, (m1, m2) = (5, 10), where m<sup>1</sup> and m<sup>2</sup> are the sample sizes for the first-stage and second-stage sampling.

We test H<sup>0</sup> : θ<sup>2</sup> = θ (0) <sup>2</sup> with α = 0.05 nominal significance level and consider two different values for θ (0) 2 : 1 and 1.01. Since the likelihood function involves intractable integral forms, we do not consider the likelihood-ratio test and only compare the bootstrap quasi-score method with the naive quasi-score method.

For each scenario, we generate 1 000 Monte Carlo samples, and we consider M = 500 and M = 1 000 iterations for the bootstrap quasi-score method. Table S1 summarizes the simulation results. For both scenarios, the naive quasi-score method has a larger type I error rate under H0. The design effects for estimating θ<sup>2</sup> are larger than 1.2 under both scenarios, indicating that the distribution of the test statistic is more "dispersed" than the χ 2 (1) distribution. In contrast, the type I error rates of the bootstrap quasi-score method are close to the nominal significance level when the null hypothesis is true for both scenarios. Besides, the power of the proposed bootstrap testing methods is reasonable compared with the naive quasi-score method, which has an inflated type error rate. In addition, we get similar test results regardless of the number of bootstrap repetitions.

Table S1: Test power for the hypothesis test H<sup>0</sup> : θ<sup>2</sup> = θ (0) <sup>2</sup> based on 1 000 Monte Carlo simulations, the true parameter is θ<sup>2</sup> = 1.00, and the nominal significance level is 0.05. The numbers of bootstrap repetitions are set to be M = 500 and M = 1 000.

NOTE: NQS, naive quasi-score method; BQS I: bootstrap quasi-score method with M = 500 bootstrap replications; BQS II: bootstrap quasi-score method with M = 1 000 bootstrap replications.

#### References

Andersen, P. K. and Gill, R. D. (1982), "Cox's regression model for counting processes: a large sample study," Annals of Statistics, 10, 1100–1120.

- Athreya, K. B. and Lahiri, S. N. (2006), Measure Theory and Probability Theory, New York: Springer.
- Beaumont, J.-F. and Patak, Z. (2012), "On the generalized bootstrap for sample surveys with special attention to Poisson sampling," International Statistical Review, 80, 127–148.
- Bickel, P. J. and Freedman, D. A. (1984), "Asymptotic normality and the bootstrap in stratified sampling," Annals of Statistics, 12, 470–482.
- Dumitrescu, L., Qian, W., and Rao, J. N. K. (2020), "Inference for longitudinal data from complex sampling surveys: An approach based on quadratic inference functions," Scandinavian Journal of Statistics, 48, 1–29.
- Engle, R. F. and McFadden, D. L. (1994), Handbook of Econometrics: Volume IV, Amsterdam.
- Fuller, W. A. (2009), Sampling Statistics, New Jersey: Wiley.
- Horn, R. A. and Johnson, C. R. (2012), Matrix Analysis, New York: Cambridge University Press.
- Krewski, D. and Rao, J. N. K. (1981), "Inference from stratified samples: properties of the linearization, jackknife and balanced repeated replication methods," Annals of Statistics, 9, 1010–1019.
- Rabe-Hesketh, S. and Skrondal, A. (2006), "Multilevel modelling of complex survey data," Journal of the Royal Statistical Society: Series A (Statistics in Society), 169, 805–827.
- Rao, J. N. K. and Wu, C. F. J. (1988), "Resampling inference with complex survey data," Journal of the American Statistical Association, 83, 231–241.
- Rubin-Bleuer, S. and Schiopu-Kratina, I. (2005), "On the two-phase framework for joint model and design-based inference," Annals of Statistics, 33, 2789–2810.

Shao, J. and Tu, D. (1995), The Jackknife and Bootstrap, New York: Springer.

van der Vaart, A. W. (2000), Asymptotic Statistics, New York: Cambridge University Press.

Yuan, K.-H. and Jennrich, R. I. (1998), "Asymptotics of estimating equations under natural conditions," Journal of Multivariate Analysis, 65, 245–260.