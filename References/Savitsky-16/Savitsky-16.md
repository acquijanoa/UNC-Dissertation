**[Electronic Journal of Statistics](http://projecteuclid.org/ejs)** Vol. 10 (2016) 1677–1708 ISSN: 1935-7524 DOI: [10.1214/16-EJS1153](http://dx.doi.org/10.1214/16-EJS1153)

# **Bayesian estimation under informative sampling**

### **Terrance D. Savitsky and Daniell Toth**

U.S. Bureau of Labor Statistics, 2 Massachusetts Ave. N.E Washington, D.C. 20212 USA e-mail: [Savitsky.Terrance@bls.gov](mailto:Savitsky.Terrance@bls.gov) [Toth.Daniell@bls.gov](mailto:Toth.Daniell@bls.gov)

**Abstract:** Bayesian analysis is increasingly popular for use in social science and other application areas where the data are observations from an informative sample. An informative sampling design leads to inclusion probabilities that are correlated with the response variable of interest. Model inference performed on the observed sample taken from the population will be biased for the population generative model under informative sampling since the balance of information in the sample data is different from that for the population. Typical approaches to account for an informative sampling design under Bayesian estimation are often difficult to implement because they require re-parameterization of the hypothesized generating model, or focus on design, rather than model-based, inference. We propose to construct a pseudo-posterior distribution that utilizes sampling weights based on the marginal inclusion probabilities to exponentiate the likelihood contribution of each sampled unit, which weights the information in the sample back to the population. Our approach provides a nearly automated estimation procedure applicable to any model specified by the data analyst for the population and retains the population model parameterization and posterior sampling geometry. We construct conditions on known marginal and pairwise inclusion probabilities that define a class of sampling designs where L<sup>1</sup> consistency of the pseudo posterior is guaranteed. We demonstrate our method on an application concerning the Bureau of Labor Statistics Job Openings and Labor Turnover Survey.

**Keywords and phrases:** Survey sampling, Gaussian process, Dirichlet process, Bayesian hierarchical models, Latent models, Markov Chain Monte Carlo.

Received July 2015.

## **1. Introduction**

Bayesian formulations are increasingly popular for modeling hypothesized distributions with complicated dependence structures. Their popularity stems from the ease of capturing this dependence by employing models with random effects parameters with a hierarchical construction that regulates the borrowing of information for estimation. Latent parameters are often used in the model to permit flexibility in the estimation of the dependencies among the observations (Dunson, [2010](#page-30-0)). In social science applications, utilization of latent parameters may be useful for making inference about intrinsic belief states of people from their observed actions(see for example, Savitsky & Dalal [\(2013\)](#page-31-0)) Other application areas in which latent parameters may be employed include, engineering and natural science, which use them to parameterize elements of an evolving process.

Data used in these type of applications are often acquired through a complex sample design, resulting in probabilities of inclusion that are associated with the variable of interest. This association could result in an observed data set consisting of units that are not independent and identically distributed. A sampling design that produces a correlation between selection probabilities and observed values is referred to as informative. Failure to account for this dependence caused by the sampling design could bias estimation of parameters that index the joint distribution hypothesized to have generated the population (Holt et al., [1980\)](#page-30-1).

### <span id="page-1-0"></span>*1.1. Examples*

We next outline some examples of survey instruments that employ informative sampling designs and associated inferential goals for models estimated on observed samples realized from these surveys.

Example 1: The Survey of Occupational Illnesses and Injuries (SOII) is administered to U.S. business establishments by the U.S. Bureau of Labor Statistics (BLS), in partnership with individual states, in order to capture workplace induced injuries and illnesses. A stratified sampling design is used where strata are indexed by state-industry-size-injury rate. Strata containing establishments that historically express higher injury rates are assigned higher sample inclusion probabilities. The resulting sample will contain a larger proportion of establishments that express higher injury rates than the population, as a whole. States desire to perform regression modeling with variable selection to discover the root causes that predict illnesses and injuries among the population of establishments, estimated from the observed sample. The model-estimated coefficients from the sample will be biased absent correction for over-representation of establishments that tend to express relatively high injury rates.

Example 2: The Current Establishment Statistics (CES) is a BLS survey of U.S. business establishments that collects employment count data across states and industries under a stratified sampling design with strata indexed by the number of employees in each establishment. Strata containing relatively larger establishments are assigned higher inclusion probabilities than those which hold establishments with relatively fewer employees. The distribution of employment totals in the observed sample of establishments will be skewed towards relatively larger values as compared to the population of establishments. An important area of modeling inference is to understand industry-indexed differences in monthly employment trends and correlations among industries in the population. We would use a mixed effects model, parameterized with random effects indexed by industry and month. Estimation of the population distribution under our model from the observed sample will be biased absent some correction for the skewness in the sample towards larger-sized establishments.

Example 3: BLS collects establishment-indexed employment totals in both the Quarterly Census of Employment and Wages (QCEW) and the CES survey. CES survey participants also provide submissions to the QCEW, such that their reported monthly employment totals for an overlapping time period of interest should be equal between the two instruments, but they are not for approximately 10000 establishments, indicating one or more employment count submission errors for those respondents. A response variable of interest, termed the "error time series", was created by taking the absolute value of the difference in reported employment totals among the 10000 establishments for each month over a 12 month period. A "response analysis survey" (RAS) of approximately 2000 establishments was taken from this population with the goal to understand the process drivers for committing errors so that BLS may target resources to establishments that mitigate them. The modeling focus is to identify probabilistic clusters of establishments with similar error patterns over the 12 month period and to examine the process by which establishments in each cluster construct their data submissions to BLS. The RAS survey design stratified the population of 10000 establishments based on phenomena of interest expressed in portions of each time series; for example, a big jump in the reported difference at yearend may indicate establishments who count checks that include regular pay and bonuses for each employee, instead of counting employees. Higher inclusion probabilities were assigned to those strata expressing phenomena of relatively greater interest to BLS researchers. Modeling the number of and memberships in probabilistic clusters of error patterns expressed in the population from the RAS sample may be biased because the proportions of error patterns expressed in the sample are designed to be different from the population.

Example 4: The Current Expenditure (CE) survey is administered to U.S. households by BLS for the purpose of determining the amount of spending for a broad collection of goods and service categories and it serves as the main source used to construct the basket of goods later used to formulate the Consumer Price Index. The CE employs a multi-stage sampling design that draws clusters of core-based statistical areas (CBSAs), such as metropolitan and micropolitan areas, from which Census blocks and, ultimately, households are sampled. Economists desire to model the propensity or probability of purchase for a variety of goods and services. The balance of sampled clusters may not be reflective of those in the population; for example, if particularly high income ares are included in the sample. So inference on purchase propensities for the population made from the observed sample will be biased absent correction for the informative sampling design.

Example 5: BLS administers the Job Openings and Labor Turnover survey (JOLTS) to business establishments with the focus to measure labor market dynamics by reporting the number of job openings, hires and separations, which is a leading indicator for employment trends. The sampling design assigns larger inclusion probabilities to establishments with relatively more employees because larger establishments drive the variance in the reported statistics. Our modeling goals are to understand differences in labor force dynamics based on employment ownership (e.g., private, public) and region as part of imputing missing values with respect to the population generating distribution. As with the CES sampling design, however, our sample will tend to over-represent relatively largersized establishments, so that inference and imputation using the sample will be biased for the population. We develop a multivariate count data population generating model in Section [4,](#page-11-0) where we illustrate the resulting estimation bias from failure to account for the correlations between assigned inclusion probabilities and the response variables of interest for our sample.

The target audience for this article are data analysts who wish to perform some distributional inference using data obtained from an informative sample design on a population using a model they specify, p (yi|*λ*), *λ* ∈ Λ, for density, p. We discuss, in the next section, how the limited literature on this topic does not adequately provide a general method for making distributional inference on a population while adjusting for the unequal probabilities of selection.

In this article, we propose an approach that replaces the likelihood with the "pseudo" likelihood (Chambers & Skinner, [2003\)](#page-30-2), p (yi|δ<sup>i</sup> = 1, *λ*) w<sup>i</sup> , using sampling weight, w<sup>i</sup> ∝ 1/πi. This re-weights the likelihood contribution for each observed unit with intent to re-balance the information in the observed sample to approximate the balance of information in the target finite population; correcting for the informativeness. We show that the proposed method for Bayesian estimation on complex sample data allows for asymptotically consistent inference on any population-generating model specified by the data analyst.

Additionally, this method does not require information about the complex design, other than the probabilities of selection, or about the full population, other than the observed data. We believe this makes the method applicable to more situations. Indeed, it is often the case that the data analyst does not have access to the full design information or auxiliary variables on the population, z1,...,z<sup>N</sup> , used to assign the probabilities of selection π1,...,π<sup>N</sup> . However, it is common for the probabilities of selection for the units in the sample, π1,...,πn, to be provided with the observed sample data.

### *1.2. Review of Methods to Account for Informative Sampling*

One current approach is to account for the informativeness by parameterizing the sampling design into the model (Little, [2004](#page-30-3)). Parameterizing even a simple informative design is often difficult to accomplish and may disrupt desired inference by requiring a change to the underlying population model parameterization. The analyst in Example 3, above, desires to perform inference on an a priori unknown clustering of sampled units with their population model for data acquired under a stratified sampling design. Specifying random effects to be indexed by strata will likely conflict with the identification and composition of inferred clusters. Further, the data analyst may not have access to the sampling design, but only indirect information in form of sampling weights. Lastly, the analyst is sometimes required to impute the unobserved units in the finite population, which may be computationally infeasible.

Another approach incorporates the sampling weights into inference about the population, as is our intent, but requires a particular form for the likelihood that

does not allow the analyst to impose their own population model formulation of inferential interest. For example, Dong et al. [\(2014\)](#page-30-4) specifies an empirical likelihood, while Kunihama et al. [\(2014\)](#page-30-5) constructs a non-parametric mixture for the likelihood and Rao & Wu [\(2010](#page-30-6)) uses a sampling-weighted (pseudo) empirical likelihood. All of these approaches impose Dirichlet distribution priors for the mixture components with hyperparameters specified as a function of the first-order sampling weights. Si et al. [\(2015](#page-31-1)) regress the response variable on a Gaussian process function of the weights for sampling designs where subgroups of sampled units have equal weights (e.g., a stratified sampling design). These approaches are designed for inference about simple mean and total statistics, rather than inference for parameters that characterize an analyst-specified population model that is the focus for our proposed method.

One method that uses a plug-in estimator, as do we in our method, is to construct a joint likelihood of the population distribution and sample inclusion in a simple logistic regression model (Malec et al., [1999](#page-30-7)). This allows one to analytically marginalize over the parameters indexed by the non-sampled units. This approach is limited in application to a class of simple population models that permit analytic integration and may not be applied to more general classes of Bayesian models for the population that we envision in development of our approach.

Perhaps the most general Bayesian approach constructs models to co-estimate parameters for conditional expectations of inclusion probabilities jointly with the population-generating model parameters at each level of a hierarchical construction (Pfeffermann et al., [2006\)](#page-30-8). This formulation is fully Bayesian such that it accounts for all sources of uncertainty in population generation and inclusion of units, but requires a custom implementation of an MCMC sampler for each specified population model, such as their simple two-level linear regression model. The implementations may increase the complexity of the specified model and reduce the quality of posterior mixing in the MCMC, so that they are suitable for relatively simple population probability models.

The method we propose is intended to allow Bayesian inference from any population model that may be specified by the the data analyst under an informative sampling design, unlike the alternative methods. It provides asymptotically unbiased estimation using only the distribution for the observed sample units and normalized H´ajek-like sampling weights. The "plug-in" type method accounts for the informative sampling design by raising the likelihood contribution of each sampled observation to the power of their associated sampling weight. The implementation of the plug-in procedure for Bayesian estimation multiplies the sampling weight into each full conditional log-posterior density. This can then sampled in the typical sequential scan MCMC.

Unlike these other methods that are prominent in the literature, this method: 1) does not impose a population model (implicitly or explicitly), unlike the most recently-developed methods (Dong et al., [2014](#page-30-4); Kunihama et al., [2014;](#page-30-5) Rao & Wu, [2010](#page-30-6); Si et al., [2015\)](#page-31-1); 2) requires only the sampling weights and does not require parameterizing the sampling design unlike Little [\(2004](#page-30-3)); 3) does not require a customized MCMC sampling procedure unlike Pfeffermann et al. [\(2006](#page-30-8)), so can be done automatically; 4) does not require imputing the non-sampled units in the finite population. Our data application and estimation model in the sequel are intended to be representative of common problems for Bayesian inference, and the application data are not readily estimated with these other methods that account for informative sampling.

We formulate the pseudo-posterior density as sampling weight-adjusted plugin from which we conduct model inference about the population under a dependent, informative sampling design in Section [2.](#page-5-0) Conditions are constructed that guarantee a frequentist L<sup>1</sup> contraction of the pseudo posterior distribution on the true generating distribution in Section [3.](#page-6-0) We make an application of the pseudo posterior estimator to construct a regression model for count data using a dataset of monthly job hires and separations collected by the U.S. Bureau of Labor Statistics in Section [4.](#page-11-0) We reveal large differences for parameter estimates between incorporation versus ignoring the sampling weights. This section also includes a simulation study that compares the pseudo posterior estimated on the observed sample to the posterior estimated on the entire finite population. The paper concludes with a discussion in Section [5.](#page-19-0) The proofs for the main result, along with two enabling results are contained in an Appendix.

## <span id="page-5-0"></span>**2. Method to account for informative sampling**

We begin by constructing the pseudo likelihood and associated pseudo posterior density under any analyst-specified prior formulation on the model, *λ* ∈ Λ.

### *2.1. Pseudo posterior*

Suppose there exists a Lebesgue measurable population-generating density, π (y|*λ*), indexed by parameters, *λ* ∈ Λ. Let δ<sup>i</sup> ∈ {0, 1} denote the sample inclusion indicator for units i = 1,...,N from the population under sampling without replacement. The density for the observed sample is denoted by, π (yo|*λ*) = π (y|δ<sup>i</sup> = 1,*λ*), where "o" indicates "observed".

The plug-in estimator for posterior density under the analyst-specified model for *λ* ∈ Λ is

$$
\hat{\pi}(\lambda|\mathbf{y}_o, \tilde{\mathbf{w}}) \propto \left[\prod_{i=1}^n p(y_{o,i}|\lambda)^{\tilde{w}_i}\right] \pi(\lambda), \tag{1}
$$

where <sup>n</sup> <sup>i</sup>=1 p (yo,i|*λ*) <sup>w</sup>˜<sup>i</sup> denotes the pseudo likelihood for observed sample responses, **y**o. The joint prior density on model space assigned by the analyst is denoted by π (*λ*). This pseudo likelihood employs sampling weights, {w˜<sup>i</sup> ∝ 1/πi}, constructed to be inversely proportional to unit inclusion probabilities. Each sampling weight assigns the relative importance of the likelihood contribution for each sample observation to approximate the likelihood for the population. We use ˆπ to denote the noisy approximation to posterior distribution, π, and we make note that the approximation is based on the data, **y**<sup>o</sup> , and sampling weights, {**w**˜ }, confined to those units included in the sample, S.

The total estimated posterior variance is regulated by the sum of the sampling weights. We define unnormalized weights, {w<sup>i</sup> = 1/πi}, and subsequently normalize them, ˜w<sup>i</sup> = w<sup>i</sup> wi n , i = 1,...,n, to sum to the sample size, n, the asymptotic units of information in the sample. Incorporation of the sampling weights to formulate the pseudo posterior estimator is expected to increase the estimated parameter posterior variances relative to the (unweighted) posterior estimated on a simple random (non-informative) sample because the weights encode the uncertainty with which samples represent the finite population under repeated sampling. This increase in estimated posterior variance may be partly or wholly offset to the extent that the informative sampling design is more efficient than simple random sampling; for example, a stratified sampling design that takes simple random samples within each stratum may produce samples that provide better coverage of the population. Although our method utilizes the weights as a "plug-in", rather than imposing a prior, Pfeffermann & Sverchkov [\(2009\)](#page-30-9) use Bayes rule to demonstrate one may replace the weights with their conditional expectation given the observed response to correct for informative sampling. Replacing the raw weights with their conditional expectation given the observed response may serve to reduce the total variation attributed to weighting (and the resulting posterior uncertainty) in the case where the actual sampled observations express information in different proportions than intended in the sampling design. Even though the conditional distribution of the weights given the response is generally different for the observed sample than for the population, nevertheless their conditional expectations are equal.

### <span id="page-6-0"></span>**3. Pseudo posterior consistency**

We formulate a pseudo posterior distribution in this section and specify conditions under which it contracts on the true generating distribution in L1. Let ν ∈ Z<sup>+</sup> index a sequence of finite populations, {Uν}<sup>ν</sup>=1,...,N<sup>ν</sup> , each of size, |Uν| = Nν, such that N<sup>ν</sup> < Nν- , for ν<ν-, so that the finite population size grows as ν increases. Suppose that **X**ν,1,..., **X**ν,N<sup>ν</sup> are independently distributed according to some unknown distribution P, (with density, p) defined on the sample space, (X , A). If Π is a prior distribution on the model space, (P, C) to which P is known to belong, then the posterior distribution is given by

<span id="page-6-1"></span>
$$
\Pi\left(B|\mathbf{X}_{1},\ldots,\mathbf{X}_{N_{\nu}}\right) = \frac{\int_{P\in B} \prod_{i=1}^{N_{\nu}} \frac{p}{p_{0}}(\mathbf{X}_{i}) d\Pi(P)}{\int_{P\in\mathcal{P}} \prod_{i=1}^{N_{\nu}} \frac{p}{p_{0}}(\mathbf{X}_{i}) d\Pi(P)},\tag{2}
$$

for any B ∈ C, where we refer to {**X**ν,i}<sup>i</sup>=1,...,N<sup>ν</sup> as {**X**i}<sup>i</sup>=1,...,N<sup>ν</sup> for readability when the context is clear.

Ghosal & van der Vaart [\(2007](#page-30-10)) study the rate at which this posterior distribution converges to the assumed true (and fixed) generating distribution P0. They prove, under certain conditions on the model space, P, and the prior distribution, Π, that in P0-probability, the posterior distribution concentrates on an arbitrarily small neighborhood of P<sup>0</sup> as N<sup>ν</sup> ↑ ∞.

The observed data on which we focus is not the entire finite population, **X**1,..., **X**N<sup>ν</sup> , but rather a sample, **X**1,..., **X**n<sup>ν</sup> , with n<sup>ν</sup> ≤ Nν, drawn under a sampling design distribution applied to the finite population under which each unit, i ∈ (1,...,Nν), is assigned a probability of inclusion in the sample. These unit inclusion probabilities are constructed to depend on the realized finite population values, **X**1,..., **X**N<sup>ν</sup> , at each ν.

### *3.1. Pseudo posterior distribution*

A sampling design is defined by placing a known distribution on a vector of inclusion indicators, *δ*<sup>ν</sup> = (δν1,...,δνN<sup>ν</sup> ), linked to the units comprising the population, Uν. The sampling distribution is subsequently used to take an observed random sample of size n<sup>ν</sup> ≤ Nν. Our conditions needed for the main result employ known marginal unit inclusion probabilities, πνi = Pr{δνi = 1} for all i ∈ U<sup>ν</sup> and the second-order pairwise probabilities, πνij = Pr{δνi = 1∩δνj = 1} for i, j ∈ Uν, which are obtained from the joint distribution over (δ<sup>ν</sup>1,...,δνN<sup>ν</sup> ). The dependence among unit inclusions in the sample contrasts with the usual iid draws from P. We denote the sampling distribution by Pν.

Under informative sampling, the marginal inclusion probabilities, πνi = P{δνi = 1}, i ∈ (1,...,Nν), are formulated to depend on the finite population data values, **X**<sup>N</sup><sup>ν</sup> = (**X**1,..., **X**<sup>N</sup><sup>ν</sup> ). Since the resulting balance of information would be different in the sample, the posterior distribution for (**X**1δ<sup>ν</sup>1,..., **X**<sup>N</sup><sup>ν</sup> δνN<sup>ν</sup> ), that we employ for inference about P0, is not equal to that of Equation [2.](#page-6-1)

Our task is to perform inference about the population generating distribution, P0, using the observed data taken under an informative sampling design. We account for informative sampling by "undoing" the sampling design with the weighted estimator,

$$
p^{\pi} \left( \mathbf{X}_{i} \delta_{\nu i} \right) := p \left( \mathbf{X}_{i} \right)^{\delta_{\nu i} / \pi_{\nu i}}, \ i \in U_{\nu}, \tag{3}
$$

which weights each density contribution, p(**X**i), by the inverse of its marginal inclusion probability. This construction re-weights the likelihood contributions defined on those units randomly-selected for inclusion in the observed sample ({i ∈ U<sup>ν</sup> : δνi = 1}) to approximate the balance of information in Uν. This approximation for the population likelihood produces the associated pseudo posterior,

<span id="page-7-0"></span>
$$
\Pi^{\pi}\left(B|\mathbf{X}_{1}\delta_{\nu1},\ldots,\mathbf{X}_{N_{\nu}}\delta_{\nu N_{\nu}}\right) = \frac{\int_{P\in B} \prod_{i=1}^{N_{\nu}} \frac{p^{\pi}}{p_0^{\pi}}(\mathbf{X}_{i}\delta_{\nu i}) d\Pi(P)}{\int_{P\in\mathcal{P}} \prod_{i=1}^{N_{\nu}} \frac{p^{\pi}}{p_0^{\pi}}(\mathbf{X}_{i}\delta_{\nu i}) d\Pi(P)},\tag{4}
$$

that we use to achieve our required conditions for the rate of contraction of the pseudo posterior distribution on P0. We recall that both P and *δ*<sup>ν</sup> are random variables defined on the space of measures and possible samples, respectively. Additional conditions are later formulated for the distribution over

samples, Pν, drawn under the known sampling design, to achieve contraction of the pseudo posterior on P0. We assume measurability for the sets on which we compute prior, posterior and pseudo posterior probabilities on the joint product space, X ×P. For brevity, we use the superscript, π, to denote the dependence on the known sampling probabilities, {πνi}i=1,...,N<sup>ν</sup> ; for example, Π<sup>π</sup> (B|**X**1δν1,..., **X**N<sup>ν</sup> δνN<sup>ν</sup> ) := Π (B|(**X**1δν1,..., **X**N<sup>ν</sup> δνN<sup>ν</sup> ),(πν1,...,πνN<sup>ν</sup> )).

Our main result is achieved in the limit as ν ↑ ∞, under the countable set of successively larger-sized populations, {Uν}ν∈Z<sup>+</sup> . We define the associated rate of convergence notation, O(bν), to denote limν↑∞ O(b<sup>ν</sup> ) <sup>b</sup><sup>ν</sup> = 0.

### <span id="page-8-1"></span>*3.2. Empirical process functionals*

We employ the empirical distribution approximation for the joint distribution over population generation and the draw of an informative sample that produces our observed data to formulate our results. Our empirical distribution construction follows Breslow & Wellner [\(2007\)](#page-29-0) and incorporates inverse inclusion probability weights, {1/πνi}<sup>i</sup>=1,...,N<sup>ν</sup> , to account for the informative sampling design,

$$
\mathbb{P}_{N_{\nu}}^{\pi} = \frac{1}{N_{\nu}} \sum_{i=1}^{N_{\nu}} \frac{\delta_{\nu i}}{\pi_{\nu i}} \delta\left(\mathbf{X}_{i}\right),\tag{5}
$$

where δ (**X**i) denotes the Dirac delta function, with probability mass 1 on **X**<sup>i</sup> and we recall that N<sup>ν</sup> = |Uν| denotes the size of of the finite population. This construction contrasts with the usual empirical distribution, P<sup>N</sup><sup>ν</sup> = <sup>1</sup> N<sup>v</sup> <sup>N</sup><sup>ν</sup> <sup>i</sup>=1 δ (**X**i), used to approximate P ∈ P, the distribution hypothesized to generate the finite population, Uν.

We follow the notational convention of Ghosal et al. [\(2000](#page-30-11)) and define the associated expectation functionals with respect to these empirical distributions by P<sup>π</sup> <sup>N</sup><sup>ν</sup> f = <sup>1</sup> N<sup>ν</sup> <sup>N</sup><sup>ν</sup> i=1 δνi <sup>π</sup>νi <sup>f</sup> (**X**i). Similarly, <sup>P</sup><sup>N</sup><sup>ν</sup> <sup>f</sup> <sup>=</sup> <sup>1</sup> N<sup>ν</sup> <sup>N</sup><sup>ν</sup> <sup>i</sup>=1 f (**X**i). Lastly, we use the associated centered empirical processes, G<sup>π</sup> <sup>N</sup><sup>ν</sup> <sup>=</sup> <sup>√</sup>N<sup>ν</sup> P<sup>π</sup> <sup>N</sup><sup>ν</sup> − P<sup>0</sup> and <sup>G</sup><sup>N</sup><sup>ν</sup> <sup>=</sup> <sup>√</sup>N<sup>ν</sup> (P<sup>N</sup><sup>ν</sup> <sup>−</sup> <sup>P</sup>0).

The sampling-weighted, (average) pseudo Hellinger distance between distributions, P1, P<sup>2</sup> ∈ P, dπ,<sup>2</sup> <sup>N</sup><sup>ν</sup> (p1, p2) = <sup>1</sup> N<sup>ν</sup> <sup>N</sup><sup>ν</sup> i=1 δνi <sup>π</sup>νi <sup>d</sup><sup>2</sup> (p1(**X**i), p2(**X**i)), where <sup>d</sup> (p1, p2) = <sup>√</sup>p<sup>1</sup> <sup>−</sup> <sup>√</sup>p<sup>2</sup> <sup>2</sup> dμ <sup>1</sup> 2 (for dominating measure, μ). We need this empirical average distance metric because the observed (sample) data drawn from the finite population under P<sup>ν</sup> are no longer independent. The implication is that our result apply to finite populations generated as inid from which informative samples are taken. The associated non-sampling Hellinger distance is specified with, d<sup>2</sup> <sup>N</sup><sup>ν</sup> (p1, p2) = <sup>1</sup> N<sup>ν</sup> <sup>N</sup><sup>ν</sup> <sup>i</sup>=1 d<sup>2</sup> (p1(**X**i), p2(**X**i)).

# <span id="page-8-0"></span>*3.3. Main result*

We proceed to construct associated conditions and a theorem that contain our main result on the consistency of the pseudo posterior distribution under a class of informative sampling designs at the true generating distribution, P0. Our approach extends the main in-probability convergence result of Ghosal & van der Vaart [\(2007\)](#page-30-10) by adding new conditions that restrict the distribution of the informative sampling design. Suppose we have a sequence, ξN<sup>ν</sup> ↓ 0 and Nνξ<sup>2</sup> <sup>N</sup><sup>ν</sup> ↑ ∞ and nνξ<sup>2</sup> <sup>N</sup><sup>ν</sup> ↑ ∞ as <sup>ν</sup> <sup>∈</sup> <sup>Z</sup><sup>+</sup> ↑ ∞ and any constant, C > 0,

**(A1)** (Local entropy condition – Size of model)

$$
\sup_{\xi > \xi_{N_{\nu}}} \log N(\xi/36, \{ P \in \mathcal{P}_{N_{\nu}} : d_{N_{\nu}}(P, P_0) < \xi \}, d_{N_{\nu}}) \le N_{\nu} \xi_{N_{\nu}}^2,
$$

**(A2)** (Size of space)

$$
\Pi(\mathcal{P}\backslash\mathcal{P}_{N_{\nu}})\leq \exp\left(-N_{\nu}\xi_{N_{\nu}}^{2}\left(2(1+2C)\right)\right)
$$

**(A3)** (Prior mass covering the truth)

$$
\Pi\left(P: -P_0 \log \frac{p}{p_0} \le \xi_{N_\nu}^2 \cap P_0 \left[ \log \frac{p}{p_0} \right]^2 \le \xi_{N_\nu}^2 \right) \ge \exp\left(-N_\nu \xi_{N_\nu}^2 C\right)
$$

**(A4)** (Non-zero Inclusion Probabilities)

$$
\sup_{\nu} \left[ \frac{1}{\min_{i \in U_{\nu}} \pi_{\nu i}} \right] \leq \gamma, \text{ with } P_0\text{-probability 1.}
$$

**(A5)** (Asymptotic Independence Condition)

$$
\limsup_{\nu \uparrow \infty} \max_{i \neq j \in U_{\nu}} \left| \frac{\pi_{\nu ij}}{\pi_{\nu i} \pi_{\nu j}} - 1 \right| = \mathcal{O}(N_{\nu}^{-1}), \text{ with } P_0\text{-probability 1}
$$

such that for some constant, C<sup>3</sup> > 0,

$$
N_{\nu} \sup_{\nu} \max_{i \neq j \in U_{\nu}} \left[ \frac{\pi_{\nu ij}}{\pi_{\nu i} \pi_{\nu j}} \right] \le C_3, \text{ for } N_{\nu} \text{ sufficiently large.}
$$

**(A6)** (Constant Sampling fraction) For some constant, f ∈ (0, 1), that we term the "sampling fraction",

$$
\limsup_{\nu} \left| \frac{n_{\nu}}{N_{\nu}} - f \right| = \mathcal{O}(1), \text{ with } P_0\text{-probability 1.}
$$

Condition [\(A1\)](#page-8-0) denotes the logarithm of the covering number, defined as the minimum number of balls of radius ξ/36 needed to cover

{P ∈ P<sup>N</sup><sup>ν</sup> : d<sup>N</sup><sup>ν</sup> (P, P0) < ξ} under distance metric, d<sup>N</sup><sup>ν</sup> . This condition restricts the growth in the size of the model space, or as noted by Ghosal et al. [\(2000\)](#page-30-11), the space, P<sup>N</sup><sup>ν</sup> , must be not too big in order that the condition specifies an optimal convergence rate (Wong & Shen, [1995\)](#page-31-2). This condition guarantees the existence of test statistics, φ<sup>n</sup><sup>ν</sup> (**X**1δ<sup>ν</sup>1,..., **X**<sup>N</sup><sup>ν</sup> δνN<sup>ν</sup> ) ∈ (0, 1), needed for enabling Lemma [B.1,](#page-24-0) stated in the Appendix, that bounds the expectation of the

pseudo posterior mass assigned on the set {P ∈ PN<sup>ν</sup> : dn<sup>ν</sup> (P, P0) ≥ ξN<sup>ν</sup> }. Condition [\(A3\)](#page-8-0) ensures the prior, Π, assigns mass to convex balls in the vicinity of P0. Conditions [3.3](#page-8-0) and [3.3,](#page-8-0) together, define the minimum value of ξN<sup>ν</sup> , where if these conditions are satisfied for some ξN<sup>ν</sup> , then they are also satisfied for any ξ>ξN<sup>ν</sup> . Condition [\(A2\)](#page-8-0) allows, but restricts, the prior mass placed on the uncountable portion of the model space, such that we may direct our inference to an approximating sieve, PN<sup>ν</sup> . This sequence of spaces "trims" away a portion of the space that is not entropy bounded (in condition [\(A1\)\)](#page-8-0). In practice, trimming the space may usually be performed to ensure the entropy bound.

The next three new conditions impose restrictions on the sampling design and associated known distribution, Pν, used to draw the observed sample data that, together, define a class of allowable sampling designs on which the contraction result for the pseudo posterior is guaranteed. Condition [\(A4\)](#page-8-0) requires the sampling design to assign a positive probability for inclusion of every unit in the population because the restriction bounds the sampling inclusion probabilities away from 0. Since the maximum inclusion probability is 1, the bound, γ ≥ 1. No portion of the population may be systematically excluded, which would prevent a sample of any size from containing information about the population from which the sample is taken. Condition [\(A5\)](#page-8-0) restricts the result to sampling designs where the dependence among lowest-level sampled units attenuates to 0 as ν ↑ ∞; for example, a two-stage sampling design of clusters within strata would meet this condition if the number of population units nested within each cluster (from which the sample is drawn) increases in the limit of ν. Such would be the case in a survey of households within each cluster if the cluster domains are geographically defined and would grow in area as ν increases. In this case of increasing cluster area, the dependence among the inclusion of any two households in a given cluster would decline as the number of households increases with the size of the area defined for that cluster. Condition [\(A6\)](#page-8-0) ensures that the observed sample size, nν, limits to ∞ along with the size of the partially-observed finite population, Nν.

<span id="page-10-1"></span>**Theorem 3.1.** Suppose conditions [\(A1\)–\(A6\)](#page-8-0) hold. Then for sets P<sup>N</sup><sup>ν</sup> ⊂ P, constants, K > 0, and M sufficiently large,

<span id="page-10-0"></span>
$$
\mathbb{E}_{P_0, P_\nu} \Pi^{\pi} \left( P : d_{N_\nu}^{\pi} \left( P, P_0 \right) \ge M \xi_{N_\nu} | \mathbf{X}_1 \delta_{\nu 1}, \dots, \mathbf{X}_{N_\nu} \delta_{\nu N_\nu} \right) \le
$$
  

$$
\frac{16\gamma^2 \left[ \gamma + C_3 \right]}{\left( K f + 1 - 2\gamma \right)^2 N_\nu \xi_{N_\nu}^2} + 5\gamma \exp\left( -\frac{K n_\nu \xi_{N_\nu}^2}{2\gamma} \right),
$$
 (6)

which tends to 0 as (nν, Nν) ↑ ∞.

We note that the rate of convergence is injured for a sampling distribution, Pν, that assigns relatively low inclusion probabilities to some units in the finite population such that γ will be relatively larger. Samples drawn under a design that expresses a large variability in the sampling weights will express more dispersion in their information similarity to the underlying finite population. Similarly, the larger the dependence among the finite population unit inclusions induced by Pν, the higher will be C<sup>3</sup> and the slower will be the rate of contraction.

The separability of the conditions on P and Π (P), on the one hand, from those on the sampling design distribution, Pν, on the other hand, coupled with the sequential process of taking the observed sample from the finite population reveal that the pseudo posterior, defined on the partially-observed sample from a population, contracts on P<sup>0</sup> through converging to the posterior distribution defined on each fully-observed population. We demonstrate this property of the pseudo posterior in a simulation study conducted in Section [4.1.](#page-16-0) By contrast, if the posterior distribution, defined on each fully-observed finite population, fails to meet conditions [3.3,](#page-8-0) [3.3](#page-8-0) and [3.3](#page-8-0) for the main result from Equation [6,](#page-10-0) such that it fails to contract on P0, then the associated pseudo posterior cannot contract on P0, even if the sampling design satisfies conditions [\(A4\),](#page-8-0) [\(A5\)](#page-8-0) and [\(A6\).](#page-8-0)

The proof generally follows that of Ghosal et al. [\(2000\)](#page-30-11) with substantial modification to account for informative sampling. The L<sup>1</sup> rate of contraction of the pseudo posterior distribution with respect to the joint distribution for population generation and the taking of informative samples is derived. Our approach includes two unique enabling results. Please see Appendix sections [A](#page-21-0) and [B](#page-24-1) for details.

## <span id="page-11-0"></span>**4. Application**

We construct a model for count data and perform inference on survey responses collected by the Job Openings and Labor Turnover Survey (JOLTS), introduced in Example 5 of Section [1.1,](#page-1-0) which is administered by BLS on a monthly basis to a randomly-selected sample from a frame composed of non-agricultural U.S. private (business) and public establishments. JOLTS focuses on the demand side of U.S. labor force dynamics and measures job hires, separations (e.g. quits, layoffs and discharges) and openings. The JOLTS sampling design assigns inclusion probabilities (under sampling without replacement) to establishments to be proportional to the number of employees for each establishment (as obtained from the Quarterly Census of Employment and Wages (QCEW)). This design is informative in that the number of employees for an establishment will generally be correlated with the number of hires, separations and openings. We perform our modeling analysis on a May, 2012 data set of n = 8595 responding establishments.

We begin by specifying a finite population regression probability model from which we formulate the sampling-weighted pseudo posterior joint distribution that we use to make inference on model parameters from the population generating distribution with only the observed sample of a finite population. We demonstrate that failing to incorporate sampling weights (e.g. by estimating the posterior distribution defined for the finite population on the observed sample) produces large differences in estimates of parameters.

Our regression model defines a multivariate response as the number of job hires (Hires) for the first response variable and total separations (Seps) as the

second response variable. We construct a single multivariate model (as contrasted with the specification of two univariate models) because these variables of interest tend to be highly correlated such that we expect the regression parameters to express dependence; for example, these two variables are correlated at 60% in our May 2012 dataset.

We formulate a model for count data that accommodates the high degree of over-dispersion expressed in our establishment-indexed multivariate responses due to the large employment size differences across the establishments. Were we working with domain-indexed (e.g., by state or county) responses, we may consider to use a Gaussian approximation for the count data likelihood, but such is not appropriate for us due to the presence of many small-sized establishments. The modeling of count data outcomes is very typical for the analysis of BLS survey data for establishments focused on (un)employment.

We specify the following count data model for the population,

$$
y_{id} \stackrel{\text{ind}}{\sim} \text{Pois} \left( \exp \left( \psi_{id} \right) \right) \tag{7}
$$

$$
\Psi \sim \mathbf{X}^{N \times D} \mathbf{B}^{N \times D} + \mathcal{N}_{N \times D} \left( \mathbb{I}_N, \Lambda^{-1} \right)
$$
\n(8)

<span id="page-12-0"></span>
$$
\mathbf{B} \sim \mathbf{0} + \mathcal{N}_{P \times D} \left( \mathbf{M}^{-1}, \left[ \tau_B \Lambda \right]^{-1} \right) \tag{9}
$$

<span id="page-12-1"></span>
$$
\Lambda \sim \mathcal{W}_D\left( (D+1), \mathbf{I}_D \right) \tag{10}
$$

$$
\tau_B \sim \mathcal{G}(1,1) \tag{11}
$$

$$
\mathbf{M} \sim \mathcal{W}_P\left((P+1), \mathbf{I}_P\right),\tag{12}
$$

where i = 1,...,N indexes the number of establishments and d = 1,...,D indexes the number of dimensions for the multivariate response, **Y**. The N × D log-mean, Ψ = <sup>D</sup>×<sup>1</sup> *ψ*-<sup>1</sup> ,..., *ψ*-N , may be viewed as a latent response whose

columns index the number of job hires (Hires) and total separations (Seps) under our JOLTS application, so that D = 2. The number of predictors in the design matrix, **X**, is denoted by P and **B** are the unknown matrix of population coefficients that serve as the focus for our inference. Our model is formulated as a multivariate Poisson-lognormal model, under which the Gaussian prior of Equation [8](#page-12-0) for the logarithm of the Poisson mean allows for over-dispersion (of different degrees) in each of the D dimensions. The priors in Equation [8](#page-12-0) and Equation [9](#page-12-1) are formulated in matrix variate (or, more generally, tensor product) Gaussian distributions using the notation of Dawid [\(1981](#page-30-12)); for example, the prior for the P × D matrix of coefficients, **B**, assigns the P × D mean **0** for a Gaussian distribution that employs a separable covariance structure where the P × P, **M**, denotes the precision matrix for the columns of **B**, and the D×D, τBΛ, denotes the precision matrix for the rows. This prior formulation is the equivalent of assigning a P D dimensional Gaussian prior to a vectorization of **B** accomplished by stacking its columns with P D × P D precision matrix, **M**⊗(τBΛ). (See Hoff [\(2011](#page-30-13)) for more background). Precision matrices, (**M**,Λ), each receive Wishart priors with hyperparameter values that impose uniform marginal prior distributions on the correlations (Barnard et al., [2000\)](#page-29-1).

We regress the multivariate latent response, Ψ, on predictors representing the logarithm of the overall establishment-indexed number of employees (Emp), obtained from the QCEW, the logarithm of the number of job openings (Open), region (Northeast, South, West, Midwest (Midw)) and ownership type (Private, Federal Government, State Government (State), Local Government (Local)). We convert region and ownership type to binary indicators and leave out the Northeast region and Federal Government ownership to provide the baseline of a full-column rank predictor matrix. We summarize our regression model on the logarithm scale by: (ψHires, ψSeps) ∼ 1 + West + Midw + South + State + Local + Private + log(Emp) + log(Opens) + error, where 1 denotes an intercept (Int).

Our population model is hypothesized to generate the finite population of the U.S. non-agricultural establishments, from which we have taken a sample of size n = 8595 for May, 2012 as our observations. For ease of reading, we will continue to use **Y** and **X**, to next define the associated pseudo posterior, though each possesses n<N rows representing the sampled observations, in this context.

The population model likelihood contribution for establishment, i, on dimension, d, is formed with the integration,

$$
p(y_{id}|\mathbf{x}_i, \mathbf{B}, \Lambda) = \int_{\mathbb{R}} p(y_{id}|\psi_{id}) \times p(\psi_{id}|\mathbf{x}_i, \mathbf{B}, \Lambda) d\psi_{id}, \qquad (13)
$$

where sampling weight, w<sup>i</sup> = 1/π<sup>i</sup> and ˜w<sup>i</sup> = n × wi/ <sup>n</sup> <sup>i</sup>=1 wi, such that the adjusted weights sum to n, the asymptotic amount of information contained in the sample (under a sampling design that obeys condition [\(A5\)\)](#page-8-0). This integrated likelihood induces the following pseudo likelihood,

$$
p^{\pi} (y_{id} | \mathbf{x}_i, \mathbf{B}, \Lambda) = \left[ \int_{\mathbb{R}} p (y_{id} | \psi_{id}) \times p (\psi_{id} | \mathbf{x}_i, \mathbf{B}, \Lambda) d\psi_{id} \right]^{\tilde{w}_i}, \quad (14)
$$

which is analytically intractable, so we perform the integration, numerically, in our MCMC using the prior for each ψid exponentiated by the normalized sampling weight, ˜wi, which we use to construct its pseudo posterior distribution. Using Bayes rule we present the logarithm of the pseudo posteriors for the latent set of D × 1 log-mean parameters, {*ψ*i}, (which are a posteriori independent over i = 1,...,n), with,

$$
\log p^{\pi} (\psi_i | \mathbf{y}_i, \mathbf{x}_i, \mathbf{B}, \Lambda) \propto
$$
\n
$$
\log \left\{ \left[ \prod_{d=1}^D \exp (\psi_{id})^{y_{id}} \exp (-\exp (\psi_{id})) \right]^{ \tilde{w}_i} \times \left[ \mathcal{N}_D \left( \psi_i | \mathbf{x}_i' \mathbf{B}, \Lambda^{-1} \right) \right]^{ \tilde{w}_i} \right\}
$$
\n(15b)

Bayesian estimation under informative sampling 1691

$$
\propto \tilde{w}_i \sum_{d=1}^{D} \left[ y_{id} \psi_{id} - \exp \left( \psi_{id} \right) \right] - \frac{1}{2} \left( \boldsymbol{\psi}_i - \mathbf{x}_i' \mathbf{B} \right)' \tilde{w}_i \Lambda \left( \boldsymbol{\psi}_i - \mathbf{x}_i' \mathbf{B} \right), \tag{15c}
$$

where we note in the second expression in Equation [15c](#page-14-0) that the sampling weights influence the prior precision for each *ψ*i, such that a higher-weighted observation will exert relatively more influence on posterior inference because this observation is relatively more representative of the population. We take samples from the pseudo posterior distribution specified Equation [15c](#page-14-0) in our MCMC using the elliptical slice sampler of Murray et al. [\(2010](#page-30-14)), where we draw *ψ*i ind∼ N<sup>D</sup> **x** i**B**,( ˜wiΛ)−<sup>1</sup> and formulate a proposal as a convex combination (parameterized on an ellipse) of this draw from the prior and the value selected on the previous iteration of the MCMC. We evaluate each proposal using the weighted likelihood in the first expression of Equation [15c.](#page-14-0)

We next illustrate the construction of the pseudo posterior distribution for the P ×D matrix of regression coefficients, **B**, (which by D-separation is independent of the observations, ((yid), given (ψid)),

$$
p^{\pi}(\mathbf{B}|\mathbf{Y}, \mathbf{X}, \Psi, \Lambda, \mathbf{M}, \tau_B) \propto \left[ \prod_{i=1}^{n} \mathcal{N}_{n \times D} \left( \psi_i | \mathbf{B}' \mathbf{x}_i, \mathbb{I}_n, \Lambda^{-1} \right)^{\tilde{w}_i} \right] \times \mathcal{N}_{P \times D} \left( \mathbf{B} | \mathbf{M}^{-1}, (\tau_B \Lambda)^{-1} \right)
$$
(16a)

log p<sup>π</sup> (**B**|**Y**, **X**, Ψ,Λ,**M**, τB) ∝

<span id="page-14-1"></span>
$$
\sum_{i=1}^{n} \left[ \frac{\tilde{w}_i}{2} \log |\Lambda| - \frac{\tilde{w}_i}{2} \left( \boldsymbol{\psi}_i - \mathbf{B}' \mathbf{x}_i \right)' \Lambda \left( \boldsymbol{\psi}_i - \mathbf{B}' \mathbf{x}_i \right) \right]
$$
  
+  $\log \mathcal{N}_{P \times D} \left( \mathbf{B} |\mathbf{M}^{-1}, (\tau_B \Lambda)^{-1} \right).$  (16b)

In a Bayesian setting, the sum of the weights (n = <sup>n</sup> <sup>i</sup>=1 w˜i) impacts the estimated posterior variance as we observe in Equation [16b.](#page-14-1) We see that weights scale the quadratic product of the Gaussian kernel in Equation [16b](#page-14-1) such that we may accomplish the same result using the matrix variate formation to define the pseudo likelihood, N<sup>n</sup>×<sup>D</sup> Ψ − **XB**|**W˜** ,Λ<sup>−</sup><sup>1</sup> , where **W˜** = diag ( ˜w1,..., w˜n), the weights for the sampled observations, from which we compute the following conjugate conditional pseudo posterior distribution defined on the n observations,

<span id="page-14-2"></span>
$$
p^{\pi}(\mathbf{B}|\mathbf{Y}, \mathbf{X}, \Psi, \Lambda, \mathbf{M}, \tau_B) = \mathbf{h}_B^{\pi} + \mathcal{N}_{P \times D}(\mathbf{B} | (\phi_B^{\pi})^{-1}, \Lambda^{-1}), \tag{17}
$$

where *φ*<sup>π</sup> <sup>B</sup> = **X**-**WX˜** + τB**M** and **h**<sup>π</sup> <sup>B</sup> = (*φ*<sup>π</sup> B)<sup>−</sup><sup>1</sup>**X**-**W˜** Ψ.

Under employment of a simpler continuous response framework, the conditional posterior for **B** retains the same form as Equation [17,](#page-14-2) except the latent response on the logarithm scale, Ψ, would be replaced by the observed data, **Y**. Intuitively, we note using a sampling-weighted pseudo prior for the latent response, Ψ, for sampling coefficients, **B**, is analogous to using the samplingweighted likelihood in the case of an observed, continuous response, **Y**.

<span id="page-14-0"></span>

![](_page_15_Figure_1.jpeg)

<span id="page-15-0"></span>Fig 1. Comparison of posterior densities for the each coefficient in the (<sup>P</sup> = 9) <sup>×</sup> (<sup>D</sup> = 2) coefficient matrix, **B**, within 95% credible intervals, based on inclusion sampling weights in a pseudo posterior (the left-hand plot in each panel) and exclusion of the sampling weights using the posterior distribution defined for the population (in the right-hand plot). Each plot panel is labeled by "predictor,response" for the two included response variables, "Hires", and "Seps" (total separations).

Each plot panel in Figure [1](#page-15-0) compares estimated posterior distributions for a coefficient in **B** (within 95% credible intervals), labeled by "predictor, dimension (of the multivariate response)", when applied to the May, 2012 JOLTS dataset between two estimation models: 1. The left-hand plot in each panel employs the sampling weights to estimate the pseudo posterior for **B**, induced by the pseudo posterior for the latent response in Equation [15c;](#page-14-0) 2. The right-hand plot estimates the coefficients using the posterior distribution defined on the finite population, which may be achieved by replacing **W˜** by the identity matrix to equally weight establishments. Equal weighting of establishments assumes that the sample represents the same balance of information as the population from which it was drawn, which is not the case under an informative sampling design. Comparing estimation results from the pseudo posterior and population posterior distributions provides one method to assess the sensitivity of estimated parameter distributions to the sampling design.

We observe that the estimated results are quite different in both location and variation between estimations performed under the pseudo posterior and population posterior distributions, indicating a high degree of informativeness in the sampling design. The 95% credible intervals for the coefficients of the continuous predictors – (the log of) job openings (Opens) and employment (Emp) – don't even overlap on both the number of hires (Hires) and separations (Seps)

<span id="page-16-1"></span>Characteristics of single stage, fixed size pps sampling design used in simulation study. nν denotes the sample size. CUs denotes the number of certainty units (with inclusion probabilities equal to 1). πν denotes the inclusion probabilities (proportional to square root of JOLTS employment), CV(πν) denotes the coefficient of variation of πν, Cor(yhires,πν) denotes correlation of the number of hires and πν and Cor(ySeps,πν) denotes the correlation of the number of separations and πν.

responses. The coefficient for the State ownership predictor and the number of hires response is bounded away from 0 when estimated under the (unweighted) population posterior, but is centered on 0 under the sampling-weighted, pseudo posterior. The coefficient posterior variances estimated on the observed sample under the population posterior are understated because they don't reflect the uncertainty with which the information in the sample expresses that in the population (which is captured through the sampling weights).

# <span id="page-16-0"></span>*4.1. Simulation Study*

We implement a simulation study to compare the marginal pseudo posterior distributions to the (unweighted) population posterior distributions for the regression coefficients, where both are estimated on the observed sample drawn under an informative sampling design. For this study we use the N = 8595 observations from the JOLTS May, 2012 data as our population. We take 100 Monte Carlo samples of size **n**<sup>ν</sup> = (500, 1000, 1500, 2500) establishments using an informative single-stage sample design with unequal inclusion probabilities based on the proportional to size sample used for the real JOLTS survey. Characteristics of the the sampling design, used for this study, at each sample size are presented in Table [1.](#page-16-1)

This sampling design will induce distributions of the observed samples that will be different from those for the population. The designed correlation between the response and inclusion probabilities will produce observed samples with values skewed towards higher numbers of hires and separations than in the population. Figure [2](#page-17-0) demonstrates this difference between the distributions for realized samples under the informative sampling design compared to those for the finite population. The left-most box plot in each of the two panels displays the population distribution for a response value. A single sample is drawn under a sequence of increasing sample sizes for illustration. The next set of box plots displays the resulting distributions for the response values in each sample with size increasing from left-to-right. The left-hand plot panel displays the distributions for the Hires response, while the right-hand panel displays those for the Seps (separations) response variable.

1694 T. D. Savitsky and D. Toth

![](_page_17_Figure_1.jpeg)

<span id="page-17-0"></span>Fig 2. Distributions of response values for population compared to informative samples. The left-most box plot in each of the two plot panels contains the distribution for the JOLTS sample that we use as our "population" in the simulation study. The next set of box plots show the distribution for the response values for increasing sample sizes (from left-to-right) for each sample drawn under our single stage proportion-to-size design. The left-hand plot panel displays the Hires response variable and the right-hand panel displays the Seps (separations) response variable.

Pseudo posterior and population posterior distributions are estimated on each Monte Carlo sample at each sample size in **n**ν. Figure [3](#page-18-0) compares estimation of the posterior distribution from the fully-observed population (left-hand box plot) to estimation using the pseudo posterior from sample observations taken under the proportional-to-size sampling design. The third box plot in each panel shows the estimation of the posterior distribution estimated on the same sample ignoring the informative sampling design. The last box plot in each panel displays the estimates of the posterior distribution from a simple random sample of the same size, where no correction for the sampling design is required, as a gold standard against which to measure the performance of the pseudo posterior distribution. We estimate the distributions on each of the 100 Monte Carlo draws for each sample size and concatenate the results such that they incorporate both the variation of population generation and repeated sampling from that population. The sample sizes, nν, increase from left-to-right across the plot panels. The top set of plot panels display the posterior distributions of the regression coefficient for the employment predictor (Emp) and the hires response (Hires), while the bottom set of panels display the coefficient distributions for the employment predictor (Emp) and the total separations response (Seps).

Scanning from left-to-right in each row across the increasing sample sizes, we readily note a consistent difference in the estimated posterior mean, as expected, between the population model estimated on the samples without ad-

![](_page_18_Figure_1.jpeg)

<span id="page-18-0"></span>Fig 3. Comparison of posterior densities for 2 coefficients, Employment-Hires (top row of plot panels) and Employment-Separation (bottom row of plot panels) in **B**, within 95% credible intervals, between estimation on the population (left-hand plot in each panel), estimations from informative samples data taken from that population, which include sampling weights in a pseudo posterior (the second plot from the left in each panel) and exclusion of the sampling weights using the population posterior distribution (the third plot from the left) under a simulation study. The right-most plot presents the posterior density estimated from a simple random sample of the same size for comparison. The simulation study uses the May, 2012 JOLTS sample as the "population" and generates 500 informative samples for a range of sample sizes (of 500, 1000, 1500, 2500, from left-to-right) under a sampling without replacement design with inclusion probabilities set proportionally to the square root of employment levels. A separate estimation is performed on each Monte Carlo sample and the draws from estimated distributions are concatenated over the samples.

justment for the informative sampling design as compared to the mean of the posterior distribution estimated on the entire finite population. The application of the pseudo posterior model, however, produces much less difference (relative to estimation on the fully observed population), though the difference between the estimated pseudo-posterior and the population posterior is yet notably more than that for the simple random sampling result (estimated on samples of the same size as the informative sample). The estimated difference for the pseudoposterior converges to 0, however, as the sample size increases. The posterior variance for the estimated posterior under simple random sampling remains larger than that for the pseudo posterior estimated on the informative sample because our proportion-to-size sampling design over-samples the highest variance units, which provides better capture of information in the population (which is why this design is used). So, in this case, the improved capture of information in the finite population provided by our sampling design more than overcomes the added variation induced by estimation with the sampling weights. In summary, this simulation study demonstrates the contraction of the pseudo posterior distribution estimated on the sample onto the posterior distribution estimated on a fully-observed finite population.

1696 T. D. Savitsky and D. Toth

We were able to directly perform posterior inference about the population using only quantities available for the observed sample under the pseudo posterior full conditional distributions outlined in Section [4.](#page-11-0) By contrast, Little [\(2004](#page-30-3)) offer no modeling approach that parameterizes a proportion-to-size sampling design because they note that each unit is in its own group under the stratum-indexed construction they generally suggest. A typical naive approach, however, is to simply include the sampling weights or a variable highly corrected with them as a predictor only for observed units with no imputation of the non-sampled units. This is precisely the construction of the alternative that estimates the population posterior distribution on the informative sample, which is shown as the third box plot in each plot panel of Figure [3,](#page-18-0) because the employment variable, Emp, is included as a predictor and is highly correlated with the sampling weights. This option includes Emp only for the sampled units as does the model for the pseudo posterior. The reason for biased inference, even when including a predictor that is highly correlated with sampling weights, is because the distribution for the sampled data conditioned on the sampling weights is not generally equal to the distribution for the population conditioned on the sampling weights by Bayes rule (Pfeffermann & Sverchkov, [2009\)](#page-30-9).

The JOLTS respondent-level data from which samples were drawn for our Monte Carlo simulation study may not be publicly released due to restrictions that protect confidentiality of survey participants. A Monte Carlo simulation study using our pseudo posterior plug-in method may, however, be generated under a Bayesian nonparametric model for functional data that is available from the growfunctions package for R (Savitsky, [2015\)](#page-31-3). The package includes a synthetic data generator whose use is illustrated in Savitsky [\(2014](#page-31-4)), along with a Monte Carlo simulator that compares parameter estimates when correcting versus ignoring an informative sampling design. Although the functional data model is more complicated than our count data application, the Monte Carlo simulation function available in growfunctions produces a figure that demonstrates results very similar to those displayed in Figure [3.](#page-18-0)

### <span id="page-19-0"></span>**5. Discussion**

A variety of broadly applicable approaches are available for incorporating sampling weights into weighted maximum likelihood estimation procedures (Pfeffermann & Sverchkov, [2009\)](#page-30-9) to account for an informative sampling design. Defining easily adaptable algorithms that account for sampling design informativeness under any model for the population specified by the data analyst has proved more challenging for estimating Bayesian probability models. Solutions have focused on parameterizing the sampling design or co-estimating the conditional expectation of inclusion, along with the population-generating model. While these approaches allow estimation using the sampled observations, the implementations typically require a high degree of customization to each population model. We take a different approach that constructs a sample-weighted pseudoposterior to account for an informative sampling design in our plug-in method that is readily accommodated to any Bayesian population probability model.

We demonstrated the applicability of the plug-in method to a poisson – lognormal model for count data. We showed that the plug-in method reduces estimation bias and posterior estimation includes the uncertainty with which the sample reflects the population on these covariance parameters.

The plug-in method is as easily-implemented and broadly applicable as those methods used for likelihood based optimization. We illustrated in Section [4](#page-11-0) that the full conditional posterior distributions defined for the population generating model are easily updated by multiplying the log-likelihoods for ({yid}, ψid), by the sampling weights, {w˜i}, without changing the constructions for full conditional posterior distributions. The same concerns that apply in the use of sampling weights under likelihood optimization also apply for Bayesian estimation. The quality of posterior estimation is highly dependent on the population representativeness of the realized sample. Sampling weights may be adjusted based on the composition of the realized sample through estimating the conditional expectation of the weights, given the response values for the observed units. Regressing the weights on the response variables using the observed data and replacing the raw weights with their conditional expectation, known as "weight smoothing", would be expected to reduce the posterior variances for estimated parameters to the extent that the weights express variance unrelated to the response. While the conditional distribution for the sampling weights given the response under the realized sample is not generally expected to be the same as that for the finite population, their expectations are equal (Pfeffermann & Sverchkov, [2009](#page-30-9)). We explored such smoothing for our sampling weights for the JOLTS application, but there was little reduction in variance, so we employed the published weights for simplicity.

Even after adjustment, if the composition of the realized sample unevenly reflects information in the population, the weights would express a high variation. Approaches that calibrate the sampling weights to actual population totals, where known, may improve the quality of estimation produced from the plug-in method. BLS performs a calibration adjustment of the JOLTS sampling weights such that the weighted difference of hires and separations reported in the sample ties to the monthly total employment change from the CES survey. (The CES survey is introduced and discussed in Section [1.1\)](#page-1-0). This step adjusts the sampling weights computed from inclusion probabilities under the JOLTS sampling design based on the actual sample achieved in each month.

One may, alternatively, implement a more fully Bayesian approach that parameterizes a joint model for the response and sampling weights, specific to a given population generating model, as a method that smoothes the weights in the presence of the response values. Doing so, however, requires imputation of weights and response values for non-sampled units, which can be computationally expensive for a survey that samples from the entire U.S. population of business establishments, as does JOLTS.

Lastly, we construct conditions which, together, define a class of sampling designs under which L<sup>1</sup> consistency of the pseudo posterior is guaranteed. One of these conditions requires that the pairwise sample inclusion dependencies asymptotically decrease to 0. While there are many sampling designs, in prac1698 T. D. Savitsky and D. Toth

tice, which are members of this class, including the proportion-to-size sampling design used for our JOLTS application, there are some designs which are not – such as a cluster sampling design where the number of clusters grows, but the number of units in each cluster remains relatively fixed. A direction for future research will be to widen the class of allowable designs by incorporating second order (or pairwise) inclusion probabilities for inference, though doing so will introduce some practical restrictions on the specifications for the population generating model.

### <span id="page-21-0"></span>**Appendix A: Proof of Theorem [3.3](#page-8-0)**

Proof. Condition [\(A1\)](#page-8-0) establishes the existence of test statistics, φn<sup>ν</sup> (X1δν1,...,XN<sup>ν</sup> δνN<sup>ν</sup> ) ∈ (0, 1) used to achieve the following result,

<span id="page-21-2"></span><span id="page-21-1"></span>
$$
\mathbb{E}_{P_0, P_\nu} \phi_{n_\nu}
$$
\n
$$
\leq \exp\left(n_\nu \xi_{N_\nu}^2\right) \cdot \frac{\exp\left(-K n_\nu M^2 \xi_{n_\nu}^2\right)}{1 - \exp\left(-K n_\nu M^2 \xi_{N_\nu}^2\right)}
$$
\n
$$
\leq 2 \exp\left(-K n_\nu \xi_{N_\nu}^2\right),\tag{18}
$$

in Lemmas 2 and 9 of Ghosal & van der Vaart [\(2007\)](#page-30-10) by setting ξ = Mξ<sup>N</sup><sup>ν</sup> , and by choosing constant M > 0 sufficiently large, such that KM<sup>2</sup> − 1 > K.

We will bound the expectation (under (P0, Pν), jointly) of the mass assigned by pseudo posterior distribution for those P at some minimum distance from P0,

$$
\Pi^{\pi} (P \in \mathcal{P} : d_{N_{\nu}}^{\pi} (P, P_0) \ge M \xi_{N_{\nu}} | X_1 \delta_{\nu 1}, \dots, X_{N_{\nu}} \delta_{\nu N_{\nu}}) =   
\Pi^{\pi} (P \in \mathcal{P} : d_{N_{\nu}}^{\pi} (P, P_0) \ge M \xi_{N_{\nu}} | X_1 \delta_{\nu 1}, \dots, X_{N_{\nu}} \delta_{\nu N_{\nu}}) (\phi_{n_{\nu}} + 1 - \phi_{n_{\nu}}).
$$
(19)

Equation [18](#page-21-1) establishes the bound,

$$
\mathbb{E}_{P_0, P_{\nu}} \Pi^{\pi} \left( P \in \mathcal{P} : d_{N_{\nu}}^{\pi} \left( P, P_0 \right) \ge M \xi_{N_{\nu}} \middle| X_1 \delta_{\nu 1}, \dots, X_{N_{\nu}} \delta_{N_{\nu}} \right) \phi_{n_{\nu}} \le \mathbb{E}_{P_0} \phi_{n_{\nu}} \le
$$
  
2 
$$
2 \exp \left( -K n_{\nu} \xi_{N_{\nu}}^2 \right), \quad (20)
$$

since the pseudo posterior mass is bounded from above by 1. We next enumerate the pseudo posterior distribution for the second term of Equation [19,](#page-21-2)

<span id="page-21-3"></span>
$$
\Pi^{\pi} \left( P \in \mathcal{P} : d_{N_{\nu}}^{\pi} \left( P, P_0 \right) \ge M \xi_{N_{\nu}} \left| \mathbf{X}_{N_{\nu}} \delta_{N_{\nu}} \right) \left( 1 - \phi_{n_{\nu}} \right) =
$$
\n
$$
\frac{\int_{P \in \mathcal{P} : d_{N_{\nu}}^{\pi} \left( P, P_0 \right) \ge M \xi_{N_{\nu}} \prod_{i=1}^{N_{\nu}} \frac{p^{\pi}}{p_0^{\pi}} \left( X_i \delta_{\nu i} \right) d\Pi \left( P \right) \left( 1 - \phi_{n_{\nu}} \right)}{\int_{P \in \mathcal{P}} \prod_{i=1}^{N_{\nu}} \frac{p^{\pi}}{p_0^{\pi}} \left( X_i \delta_{\nu i} \right) d\Pi \left( P \right)}.
$$
\n
$$
(21)
$$

We may bound the denominator of Equation [21](#page-21-3) from below, in probability. Define the event,

$$
B_{N_{\nu}} = \left\{ P : -P_0 \log \left( \frac{p}{p_0} \right) \le \xi_{N_{\nu}}^2, P_0 \left( \log \frac{p}{p_0} \right)^2 \le \xi_{N_{\nu}}^2 \right\}
$$

We have from Lemma [B.2,](#page-27-0)

$$
\Pr\left\{\int\limits_{P\in\mathcal{P}}\prod_{i=1}^{N_{\nu}}\frac{p^{\pi}}{p_0^{\pi}}\left(X_i\delta_{\nu i}\right)d\Pi\left(P\right)\geq \exp\left[-(1+C)N_{\nu}\xi^2\right]\right\}\geq 1-\frac{\gamma+C_3}{C^2N_{\nu}\xi^2},
$$

for every P ∈ BN<sup>ν</sup> and any C > 0, γ > 1, where γ may be set closer to 1 for sampling designs that define a low gradient for inclusion probabilities, {πνi}. The constant, C<sup>3</sup> > 0, and will be close to 1 for sufficiently large ν. Condition [\(A3\)](#page-8-0) restricts the prior on B<sup>N</sup><sup>ν</sup> ,

$$
\Pi(B_{N_{\nu}}) \ge \exp\left(-N_{\nu}\xi_{N_{\nu}}^2 C\right)
$$

.

Then with probability at least 1 <sup>−</sup> <sup>16</sup>γ2[γ+C3] (KM2f−2γ)2N<sup>ν</sup> <sup>ξ</sup><sup>2</sup> ,

$$
\int_{P \in \mathcal{P}} \prod_{i=1}^{N_{\nu}} \frac{p^{\pi}}{p_0^{\pi}} \left( X_i \delta_{\nu i} \right) d\Pi \left( P \right) \ge \exp \left[ - (1 + C) N_{\nu} \xi^2 \right] \Pi \left( B_{N_{\nu}} \right)
$$
\n
$$
\ge \exp \left( - (1 + 2C) N_{\nu} \xi^2 \right)
$$
\n
$$
\ge \exp \left( - \frac{K M^2 n_{\nu} \xi_{N_{\nu}}^2}{2\gamma} \right),
$$

where we set 1 + 2C = KM2<sup>f</sup> <sup>2</sup><sup>γ</sup> , where we use condition [\(A6\)](#page-8-0) to replace f × N<sup>ν</sup> with n<sup>ν</sup> for ν sufficiently large.

Denote this event by,

$$
A_{N_{\nu}}^{\pi} = \left\{ \int\limits_{P \in \mathcal{P}} \prod_{i=1}^{N_{\nu}} \frac{p^{\pi}}{p_0^{\pi}} \left( X_i \delta_{\nu i} \right) d\Pi \left( P \right) \ge \exp \left( -\frac{K M^2 n_{\nu} \xi_{N_{\nu}}^2}{2\gamma} \right) \right\},\qquad(22)
$$

such that,

$$
\Pi^{\pi} \left( P \in \mathcal{P} : d_{N_{\nu}}^{\pi} \left( P, P_0 \right) \ge M \xi_{N_{\nu}} \middle| \mathbf{X}_{N_{\nu}} \delta_{N_{\nu}} \right) (1 - \phi_{n_{\nu}}) =
$$
\n
$$
\left[ \underbrace{\int \int \int \int \int \int \int \int \int \int \int \int \int \int \int \int \int \int \int \
$$

1700 T. D. Savitsky and D. Toth

$$
\leq \mathbb{I}\left(\left[A_{N_{\nu}}^{\pi}\right]^{c}\right) + \mathbb{I}\left(A_{N_{\nu}}^{\pi}\right) \times \left[\exp\left(\frac{KM^{2}n_{\nu}\xi_{N_{\nu}}^{2}}{2\gamma}\right)\Pi\left(\mathcal{P}\backslash\mathcal{P}_{N_{\nu}}\right) \right]   
+ \exp\left(\frac{KM^{2}n_{\nu}\xi_{N_{\nu}}^{2}}{2\gamma}\right) \int\limits_{\left\{P\in\mathcal{P}_{N_{\nu}}:d_{N_{\nu}}^{\pi}(P,P_{0})\geq M\xi_{N_{\nu}}\right\}} \prod_{i=1}^{N_{\nu}} \frac{p^{\pi}}{p_{0}^{\pi}}\left(X_{i}\delta_{\nu i}\right) d\Pi\left(P\right)\left(1-\phi_{n_{\nu}}\right)\right]
$$

Taking the expectation of both sides with respect to the joint distribution, (P0, Pν),

$$
\mathbb{E}_{P_0, P_{\nu}} \Pi^{\pi} \left( P \in \mathcal{P} : d_{N_{\nu}}^{\pi} (P, P_0) \ge M \xi_{N_{\nu}} \left| \mathbf{X}_{N_{\nu}} \delta_{N_{\nu}} \right) (1 - \phi_{n_{\nu}}) \right.   
\n\le P \left( \left[ A_{N_{\nu}}^{\pi} \right]^c \right) + \exp \left( \frac{K M^2 n_{\nu} \xi_{N_{\nu}}^2}{2\gamma} \right) \Pi \left( \mathcal{P} \backslash \mathcal{P}_{N_{\nu}} \right)   
\n+ \exp \left( \frac{K M^2 n_{\nu} \xi_{N_{\nu}}^2}{2\gamma} \right) \times   
\n\mathbb{E}_{P_0, P_{\nu}} \int_{\{P \in \mathcal{P}_{N_{\nu}} : d_{N_{\nu}}^{\pi} (P, P_0) \ge M \xi_{N_{\nu}}\}} \prod_{i=1}^{N_{\nu}} \frac{p^{\pi}}{p_0^{\pi}} \left( X_i \delta_{\nu i} \right) d\Pi \left( P \right) (1 - \phi_{n_{\nu}})   
\n\le \frac{(i)}{(K M^2 f - 2\gamma)^2 n_{\nu} \xi_{N_{\nu}}^2} + \exp \left( -\frac{K M^2 n_{\nu} \xi_{N_{\nu}}^2}{2\gamma} \right)   
\n+ \exp \left( \frac{K M^2 n_{\nu} \xi_{N_{\nu}}^2}{2\gamma} \right) \times   
\n\mathbb{E}_{P_0, P_{\nu}} \int_{\{P \in \mathcal{P}_{N_{\nu}} : d_{N_{\nu}}^{\pi} (P, P_0) \ge M \xi_{N_{\nu}}\}} \prod_{i=1}^{N_{\nu}} \frac{p^{\pi}}{p_0^{\pi}} \left( X_i \delta_{\nu i} \right) d\Pi \left( P \right) (1 - \phi_{n_{\nu}}), \tag{23}
$$

where in (i) we have used condition [3.3](#page-8-0) that bounds from above Π (P\P<sup>N</sup><sup>ν</sup> ), the prior mass assigned on the portion of the model space that lies outside the sieve, and have plugged in for constant, C.

By conditions [\(A1\),](#page-8-0) [\(A4\)](#page-8-0) and Lemma [B.1,](#page-24-0)

<span id="page-23-0"></span>
$$
\mathbb{E}_{P_0, P_\nu} \int_{\{P \in \mathcal{P}_{N_\nu} : d_{N_\nu}^{\pi}(P, P_0) \ge M \xi_{N_\nu}\}} \prod_{i=1}^{N_\nu} \frac{p^{\pi}}{p_0^{\pi}} \left(X_i \delta_{\nu i}\right) d\Pi\left(P\right) \left(1 - \phi_{N_\nu}\right)
$$
  

$$
\le 2\gamma \exp\left(\frac{-KM^2 n_\nu \xi_{N_\nu}^2}{\gamma}\right)
$$

Returning to the expectation in Equation [23,](#page-23-0)

$$
\mathbb{E}_{P_0, P_\nu} \Pi^{\pi} \left( P \in \mathcal{P} : d_{N_\nu}^{\pi} (P, P_0) \ge M \xi_{N_\nu} \middle| \mathbf{X}_{N_\nu} \delta_{N_\nu} \right) (1 - \phi_{n_\nu})
$$
\n
$$
\le \frac{16\gamma^2 \left[ \gamma + C_3 \right]}{\left( K M^2 - 2\gamma \right)^2 N_\nu \xi_{N_\nu}^2} + \exp\left( -\frac{K M^2 n_\nu \xi_{N_\nu}^2}{2\gamma} \right)
$$
\n
$$
+ \exp\left( \frac{K M^2 n_\nu \xi_{N_\nu}^2}{2\gamma} \right) \times 2\gamma \exp\left( -\frac{K M^2 n_\nu \xi_{N_\nu}^2}{\gamma} \right)
$$

Bayesian estimation under informative sampling 1701

$$
\stackrel{(i)}{\leq} \frac{16\gamma^2 \left[\gamma + C_3\right]}{\left(Kf - 2\gamma\right)^2 N_\nu \xi_{N_\nu}^2} + 3\gamma \exp\left(-\frac{K M^2 n_\nu \xi_{N_\nu}^2}{2\gamma}\right),\tag{24}
$$

where in (i) we use our earlier stated bound, KM<sup>2</sup> − 1 > K → KM<sup>2</sup> > K + 1. Bringing all the pieces together,

$$
\mathbb{E}_{P_0, P_{\nu}} \Pi^{\pi} \left( P \in \mathcal{P} : d_{N_{\nu}}^{\pi} \left( P, P_0 \right) \ge M \xi_{N_{\nu}} \middle| X_1 \delta_{\nu 1}, \dots, X_{N_{\nu}} \delta_{N_{\nu}} \right)   
\n\le 2 \exp \left( -Kn_{\nu} \xi_{N_{\nu}}^2 \right) + \frac{16\gamma^2 \left[ \gamma + C_3 \right]}{\left( Kf - 2\gamma \right)^2 N_{\nu} \xi_{N_{\nu}}^2} + 3\gamma \exp \left( -\frac{KM^2 n_{\nu} \xi_{N_{\nu}}^2}{2\gamma} \right)   
\n\le \frac{16\gamma^2 \left[ \gamma + C_3 \right]}{\left( Kf - 2\gamma \right)^2 N_{\nu} \xi_{N_{\nu}}^2} + 5\gamma \exp \left( -\frac{Kn_{\nu} \xi_{N_{\nu}}^2}{2\gamma} \right)
$$
\n(25)

where γ ≥ 1 and C<sup>3</sup> > 0. The right-hand side of Equation [25](#page-24-2) tends to 0 (as ν ↑ ∞) in P<sup>0</sup> probability. This concludes the proof.

### <span id="page-24-1"></span>**Appendix B: Enabling Lemmas**

We next construct two enabling results needed to prove Theorem [3.1](#page-10-1) to account informative sampling under [\(A4\),](#page-8-0) [\(A5\)](#page-8-0) and [\(A6\).](#page-8-0) The first enabling result, Lemma [B.1,](#page-24-0) extends the applicability of Ghosal & van der Vaart [\(2007\)](#page-30-10) – Lemmas 2 and 9 for inid models to informative sampling without replacement. This result is used to bound from above the numerator for the expectation with respect to the joint distribution for population generation and the taking of the informative sample, (P0, Pν), of the pseudo posterior distribution in Equation [4](#page-7-0) on the restricted set of measures, {P ∈ B}, where B = {P ∈ P : d<sup>N</sup><sup>ν</sup> (P, P0) > δξ<sup>N</sup><sup>ν</sup> }, (for any δ > 0). The restricted set includes those P that are at some minimum distance, δξ<sup>N</sup><sup>ν</sup> , from P<sup>0</sup> under pseudo Hellinger metric, d<sup>π</sup> <sup>N</sup><sup>ν</sup> . The second result, Lemma [B.2,](#page-27-0) extends Lemma 8.1 of Ghosal et al. [\(2000\)](#page-30-11) to bound the probability of the denominator of Equation [4](#page-7-0) with respect to (P0, Pν), from below.

<span id="page-24-0"></span>**Lemma B.1.** Suppose conditions [\(A1\)](#page-8-0) and [\(A4\)](#page-8-0) hold. Then for every ξ>ξ<sup>N</sup><sup>ν</sup> , a constant, K > 0, and any constant, δ > 0,

<span id="page-24-4"></span><span id="page-24-3"></span>
$$
\mathbb{E}_{P_0, P_{\nu}} \left[ \int_{P \in \mathcal{P} \setminus \mathcal{P}_{N_{\nu}}} \prod_{i=1}^{N_{\nu}} \frac{p^{\pi}}{p_0^{\pi}} \left( \mathbf{X}_i \delta_{\nu i} \right) d\Pi \left( P \right) \left( 1 - \phi_{n_{\nu}} \right) \right] \leq \Pi \left( \mathcal{P} \setminus \mathcal{P}_{N_{\nu}} \right) \tag{26}
$$
\n
$$
\mathbb{E}_{P_0, P_{\nu}} \left[ \int_{P \in \mathcal{P}_{N_{\nu}} : d_{N_{\nu}}^{\pi} \left( P, P_0 \right) > \delta_{\xi}} \prod_{i=1}^{N_{\nu}} \frac{p^{\pi}}{p_0^{\pi}} \left( \mathbf{X}_i \delta_{\nu i} \right) d\Pi \left( P \right) \left( 1 - \phi_{n_{\nu}} \right) \right] \leq 2\gamma \exp \left( \frac{-K n_{\nu} \delta^2 \xi^2}{\gamma} \right). \tag{27}
$$

<span id="page-24-2"></span>

The constant multiplier, γ ≥ 1, defined in condition [\(A4\),](#page-8-0) restricts the distribution of the sampling design by bounding all marginal inclusion probabilities for population units away from 0. As with the main result, the upper bound is injured by γ.

Proof. We proceed constructively to simplify the form of the expectations on the left-hand side of both Equations [26](#page-24-3) and [27](#page-24-4) and follow with an application of Lemma 2 (and result 2.2) and Lemma 9 of Ghosal & van der Vaart [\(2007\)](#page-30-10), which is used to establish the right-hand bound of Equation [27](#page-24-4) (based on the existence of tests, φn<sup>ν</sup> ).

Fixing ν, we index units that comprise the population with, U<sup>ν</sup> = {1,...,Nν}. Next, draw a single observed sample of n<sup>ν</sup> units from Uν, indexed by subsequence,

{i ∈ U<sup>ν</sup> : δνi = 1, = 1,...,nν}. Without loss of generality, we simplify notation to follow by indexing the observed sample, sequentially, with = 1,...,nν.

We next decompose the expectation under the joint distribution with respect to population generation, P0, and the drawing of a sample, Pν,

Suppose we draw P from some set B ⊂ P. By Fubini,

$$
\mathbb{E}_{P_0, P_{\nu}} \left[ \int_{P \in B} \prod_{i=1}^{N_{\nu}} \frac{p^{\pi}}{p_0^{\pi}} \left( \mathbf{X}_i \delta_{\nu i} \right) d\Pi \left( P \right) \left( 1 - \phi_{n_{\nu}} \right) \right]
$$
\n
$$
\leq \int_{P \in B} \left[ \mathbb{E}_{P_0, P_{\nu}} \prod_{i=1}^{N_{\nu}} \frac{p^{\pi}}{p_0^{\pi}} \left( \mathbf{X}_i \delta_{\nu i} \right) \left( 1 - \phi_{n_{\nu}} \right) \right] d\Pi \left( P \right) \tag{28}
$$

$$
\leq \int_{P \in B} \left\{ \sum_{\delta_{\nu} \in \Delta_{\nu}} \mathbb{E}_{P_{0}} \left[ \prod_{\ell=1}^{n_{\nu}} \left[ \frac{p}{p_{0}} \left( \mathbf{X}_{\ell} \right) \right]^{\frac{1}{\pi_{\nu \ell}}} \left( 1 - \phi_{n_{\nu}} \right) \middle| \delta_{\nu} \right] P_{P_{\nu}} \left( \delta_{\nu} \right) \right\} d\Pi \left( P \right) \tag{29}
$$

$$
\leq \int_{P \in B} \max_{\delta_{\nu} \in \Delta_{\nu}} \mathbb{E}_{P_0} \left[ \prod_{\ell=1}^{n_{\nu}} \left[ \frac{p}{p_0} \left( \mathbf{X}_{\ell} \right) \right]^{\frac{1}{\pi_{\nu \ell}}} \left( 1 - \phi_{n_{\nu}} \right) \middle| \delta_{\nu} \right] d\Pi \left( P \right) \tag{30}
$$

$$
\leq \int_{P \in B} \mathbb{E}_{P_0} \left[ \prod_{\ell=1}^{n_{\nu}} \left[ \frac{p}{p_0} \left( \mathbf{X}_{\ell} \right) \right]^{\frac{1}{\pi_{\nu \ell}}} \left( 1 - \phi_{n_{\nu}} \right) \middle| \boldsymbol{\delta}_{\nu}^* \right] d\Pi \left( P \right) \tag{31}
$$

<span id="page-25-0"></span>
$$
\leq \int_{P \in B} \mathbb{E}_{P_0} \left[ \prod_{\ell=1}^{n_{\nu}} \left[ \frac{p}{p_0} \left( \mathbf{X}_{\ell} \right) \right] \left( 1 - \phi_{n_{\nu}} \right) \middle| \delta_{\nu}^{*} \right] d\Pi \left( P \right) \tag{32}
$$

$$
\leq \int_{P\in B} P_{\pmb{\delta}_\nu^*} \left(1-\phi_{n_\nu}\right) d\Pi\left(P\right),
$$

where *δ*ν∈Δ<sup>ν</sup> P<sup>P</sup><sup>ν</sup> (*δ*ν) = 1 (S¨arndal et al., [2003](#page-30-15)) and

*δ*∗ <sup>ν</sup> ∈ Δ<sup>ν</sup> = ' {δ<sup>∗</sup> νi}<sup>i</sup>=1,...,N<sup>ν</sup> , δ<sup>∗</sup> νi ∈ {0, 1} ( denotes that sample, drawn from the space of all possible samples, Δν, which maximizes the probability under the

population generating distribution for the event of interest. The inequality in Equation [32](#page-25-0) results from <sup>p</sup> <sup>p</sup><sup>0</sup> <sup>≤</sup> 1 and <sup>1</sup> <sup>π</sup>ν ≥ 1. The conditional expectation of (1 − φn<sup>ν</sup> ) given *δ*<sup>∗</sup> <sup>ν</sup> is denoted by, P*δ*<sup>∗</sup> <sup>ν</sup> (1 − φn<sup>ν</sup> ). If P ∈ P\PN<sup>ν</sup> ,

$$
\mathbb{E}_{P_0, P_{\nu}}\left[\int_{P \in \mathcal{P} \setminus \mathcal{P}_{N_{\nu}}} \prod_{i=1}^{N_{\nu}} \frac{p^{\pi}}{p_0^{\pi}} \left(\mathbf{X}_i \delta_{\nu i}\right) (1 - \phi_{n_{\nu}})\right] d\Pi(P)
$$
\n
$$
\leq \int_{P \in \mathcal{P} \setminus \mathcal{P}_{N_{\nu}}} P_{\delta_{\nu}^{*}} \left(1 - \phi_{n_{\nu}}\right) d\Pi(P) \leq \int_{P \in \mathcal{P} \setminus \mathcal{P}_{N_{\nu}}} d\Pi(P) = \Pi\left(\mathcal{P} \setminus \mathcal{P}_{N_{\nu}}\right),
$$

since (1 − φn<sup>ν</sup> ) ≤ 1.

We next establish a bound for P*<sup>δ</sup>*<sup>∗</sup> <sup>ν</sup> (1 − φ<sup>n</sup><sup>ν</sup> ) on a sieve or slice. Let A<sup>π</sup> <sup>r</sup> = {P ∈ P<sup>N</sup><sup>ν</sup> : r<sup>N</sup><sup>ν</sup> ≤ d<sup>π</sup> <sup>N</sup><sup>ν</sup> (P, P0) ≤ 2r<sup>N</sup><sup>ν</sup> } for integers, r. Under observed **X**1δ<sup>∗</sup> <sup>ν</sup>1,..., **X**<sup>N</sup><sup>ν</sup> δ<sup>∗</sup> νN<sup>ν</sup> ∈ X , by conditions [\(A1\)](#page-8-0) and [\(A4\)](#page-8-0) we have,

$$
\sup_{P \in \mathcal{A}_r^{\pi}} P_{\delta_{\nu}^*} (1 - \phi_{n_{\nu}}) \tag{33}
$$

$$
= \sup_{\{P \in \mathcal{P}_{N_{\nu}} : r\xi \le d_{N_{\nu}}^{\pi}(P, P_0) \le 2r\xi\}} P_{\delta_{\nu}^*} (1 - \phi_{n_{\nu}})
$$
(34)

$$
\leq \sup_{\left\{P \in \mathcal{P}_{N_{\nu}}: \frac{r\xi}{\sqrt{\gamma}} \leq d_{N_{\nu}}(P, P_0) \leq \frac{2r\xi}{\sqrt{\gamma}}\right\}} P_{\delta_{\nu}^*}\left(1 - \phi_{n_{\nu}}\right) \tag{35}
$$

$$
\stackrel{(ii)}{\leq} \exp\left(-\frac{Kn_{\nu}r^2\xi^2}{\gamma}\right),\tag{36}
$$

where the smaller range in (i), P ∈ P<sup>N</sup><sup>ν</sup> : rξ <sup>√</sup><sup>γ</sup> <sup>≤</sup> <sup>d</sup><sup>N</sup><sup>ν</sup> (P, P0) <sup>≤</sup> <sup>2</sup>rξ <sup>√</sup><sup>γ</sup> , increases P*<sup>δ</sup>*<sup>∗</sup> <sup>ν</sup> (1 − φ<sup>n</sup><sup>ν</sup> ). The result in (ii) uses condition [\(A2\)](#page-8-0) to obtain the result of Lemmas 2 and 9 in Ghosal & van der Vaart [\(2007\)](#page-30-10) where we set <sup>ξ</sup> <sup>→</sup> ξ/√γ.

Finally, fixing some value for δ > 0, set r = 2δ for a given, for integers, ≥ 0. Following the approach for bounding the sum over the slices in Wong & Shen [\(1995\)](#page-31-2), let L be the smallest integer such that 2<sup>2</sup><sup>L</sup>δ2ξ<sup>2</sup> > 2γ, since d<sup>π</sup> <sup>N</sup><sup>ν</sup> <sup>&</sup>lt; <sup>√</sup>2<sup>γ</sup> (by our definition of the pseudo Hellinger metric in Section [3.2\)](#page-8-1). Then,

$$
\mathbb{E}_{P_{\theta_0}, P_{\nu}} \left[ \int_{\{P \in \mathcal{P}_{N_{\nu}} : d_{N_{\nu}}^{\pi}(P, P_0) \ge \delta \xi\}} \prod_{i=1}^{N_{\nu}} \frac{p^{\pi}}{p_0^{\pi}} \left( \mathbf{X}_i \delta_{\nu i} \right) d\Pi \left( P \right) \left( 1 - \phi_{n_{\nu}} \right) \right]
$$
(37)

$$
= \sum_{\ell=0}^{L} \mathbb{E}_{P_{\theta_0}, P_{\nu}} \int_{\{P \in \mathcal{P}_{N_{\nu}} : 2^{\ell} \delta \xi \le d_{N_{\nu}}^{\pi}(P, P_0) \le 2^{\ell+1} \delta \xi\}} \prod_{i=1}^{N_{\nu}} \frac{p^{\pi}}{p_0^{\pi}} \left( \mathbf{X}_i \delta_{\nu i} \right) d\Pi \left( P \right) (1 - \phi_{N_{\nu}}) \tag{38}
$$

1704 T. D. Savitsky and D. Toth

$$
\leq \gamma \sum_{\ell=0}^{L} \exp\left(-\frac{2^{2\ell} K n_{\nu} \delta^2 \xi^2}{\gamma}\right) \tag{39}
$$

$$
\leq 2\gamma \exp\left(-\frac{Kn_{\nu}\delta^2\xi^2}{\gamma}\right),\tag{40}
$$

for <sup>n</sup><sup>ν</sup> sufficiently large such that Kn<sup>ν</sup> <sup>δ</sup>2ξ<sup>2</sup> <sup>γ</sup> ≥ 1.

This concludes the proof.

$$
\square
$$

<span id="page-27-0"></span>**Lemma B.2.** For every ξ > 0 and measure Π on the set,

$$
B = \left\{ P : -P_0 \log \left( \frac{p}{p_0} \right) \le \xi^2, P_0 \left( \log \frac{p}{p_0} \right)^2 \le \xi^2 \right\}
$$

under the conditions [\(A2\),](#page-8-0) [\(A3\),](#page-8-0) [\(A4\),](#page-8-0) and [\(A5\),](#page-8-0) we have for every C > 0 and N<sup>ν</sup> sufficiently large,

<span id="page-27-1"></span>
$$
Pr\left\{\int\limits_{P\in\mathcal{P}}\prod_{i=1}^{N_{\nu}}\frac{p^{\pi}}{p_0^{\pi}}\left(\mathbf{X}_i\delta_{\nu i}\right)d\Pi\left(P\right)\leq \exp\left[-(1+C)N_{\nu}\xi^2\right]\right\}\leq \frac{\gamma+C_3}{C^2N_{\nu}\xi^2},\qquad(41)
$$

where the above probability is taken with the respect to P<sup>0</sup> and the sampling generating distribution, Pν, jointly.

The bound of "1" in the numerator of the result for Lemma 8.1 of Ghosal et al. [\(2000](#page-30-11)), is replaced with γ + C<sup>3</sup> for our generalization of this result in Equation [41.](#page-27-1) The sum of positive constants, γ + C3, is greater than 1 and will be larger for sampling designs where the inclusion probabilities, {πνi}, express relatively higher gradients. Observing each finite population in a skewed fashion through the taking of an informative sample may only slow the rate of posterior contraction (as compared to contraction of the posterior distribution defined on the fully observed finite population).

Proof. By Jensen's inequality,

$$
\log \int_{P \in \mathcal{P}} \prod_{i=1}^{N_{\nu}} \frac{p^{\pi}}{p_0^{\pi}} \left( \mathbf{X}_i \delta_{\nu i} \right) d\Pi \left( P \right) \ge \sum_{i=1}^{N_{\nu}} \int_{P \in \mathcal{P}} \log \frac{p^{\pi}}{p_0^{\pi}} \left( \mathbf{X}_i \delta_{\nu i} \right) d\Pi \left( P \right)
$$

$$
= N_{\nu} \cdot \mathbb{P}_{N_{\nu}} \int_{P \in \mathcal{P}} \log \frac{p^{\pi}}{p_0^{\pi}} d\Pi \left( P \right),
$$

where we recall that the last equation denotes the empirical expectation functional taken with respect to the joint distribution over population generating and informative sampling. By Fubini,

$$
\mathbb{P}_{N_{\nu}} \int_{P \in \mathcal{P}} \log \frac{p^{\pi}}{p_0^{\pi}} d\Pi(P) = \int_{P \in \mathcal{P}} \left[ \mathbb{P}_{N_{\nu}} \log \frac{p^{\pi}}{p_0^{\pi}} \right] d\Pi(P)
$$

$$
= \int_{P \in \mathcal{P}} \left[ \mathbb{P}_{N_{\nu}} \frac{\delta_{\nu}}{\pi_{\nu}} \log \frac{p}{p_0} \right] d\Pi(P)
$$

$$
= \int_{P \in \mathcal{P}} \left[ \mathbb{P}_{N_{\nu}}^{\pi} \log \frac{p}{p_0} \right] d\Pi(P)
$$

$$
= \mathbb{P}_{N_{\nu}}^{\pi} \int_{P \in \mathcal{P}} \log \frac{p}{p_0} d\Pi(P),
$$

where we, again, apply Fubini.

Then, the probability statement in the result of Equation [41](#page-27-1) is bounded (from above) by,

$$
\Pr\left\{\mathbb{G}_{N_{\nu}}^{\pi}\int\limits_{P\in\mathcal{P}}\log\frac{p}{p_{0}}d\Pi\left(P\right)\leq-\sqrt{N_{\nu}}\xi^{2}\left(1+C\right)-\sqrt{N_{\nu}}P_{0}\int\limits_{P\in\mathcal{P}}\log\frac{p}{p_{0}}d\Pi\left(P\right)\right\}
$$
\n
$$
=\Pr\left\{\mathbb{G}_{N_{\nu}}^{\pi}\int\limits_{P\in\mathcal{P}}\log\frac{p}{p_{0}}d\Pi\left(P\right)\leq-\sqrt{N_{\nu}}\xi^{2}\left(1+C\right)-\sqrt{N_{\nu}}\int\limits_{P\in\mathcal{P}}P_{0}\log\frac{p}{p_{0}}d\Pi\left(P\right)\right\}
$$
\n
$$
=\Pr\left\{\mathbb{G}_{N_{\nu}}^{\pi}\int\limits_{P\in\mathcal{P}}\log\frac{p}{p_{0}}d\Pi\left(P\right)\leq-\sqrt{N_{\nu}}\xi^{2}\left(1+C\right)+\sqrt{N_{\nu}}\xi^{2}=-\sqrt{N_{\nu}}\xi^{2}C\right\},\right\}
$$

where we have again applied Fubini in the second inequality and also the bound for P<sup>0</sup> log <sup>p</sup> <sup>p</sup><sup>0</sup> <sup>≤</sup> <sup>ξ</sup><sup>2</sup> for <sup>P</sup> on the set <sup>B</sup>.

We now apply Chebyshev and Jensen's inequality to bound the probability,

$$
\Pr\left\{\mathbb{G}_{N_{\nu}}^{\pi} \int_{P \in \mathcal{P}} \log \frac{p}{p_0} d\Pi\left(P\right) \leq -\sqrt{N_{\nu}} \xi^2 C\right\}
$$
\n
$$
\leq \frac{\text{Var}\left[\int_{P \in \mathcal{P}} \mathbb{G}_{N_{\nu}}^{\pi} \log \frac{p}{p_0} d\Pi\left(P\right)\right]}{N_{\nu} \xi^4 C^2} \tag{42a}
$$

<span id="page-28-0"></span>
$$
\leq \frac{\int_{P \in \mathcal{P}} \left[ \text{Var}\left( \mathbb{G}_{N_{\nu}}^{\pi} \log \frac{p}{p_0} \right) \right] d\Pi \left( P \right)}{N_{\nu} \xi^4 C^2} \tag{42b}
$$

$$
\leq \frac{\int\limits_{P \in \mathcal{P}} \left[ \mathbb{E}_{P_0, P_\nu} \left( \mathbb{G}_{N_\nu}^{\pi} \log \frac{p}{p_0} \right)^2 \right] d\Pi \left( P \right)}{N_\nu \xi^4 C^2} \tag{42c}
$$

<span id="page-28-2"></span><span id="page-28-1"></span>
$$
\leq \frac{\int_{P \in \mathcal{P}} \left[ \mathbb{E}_{P_0, P_{\nu}} \left( \sqrt{N_{\nu}} \mathbb{P}_{N_{\nu}}^{\pi} \log \frac{p}{p_0} \right)^2 \right] d\Pi \left( P \right)}{N_{\nu} \xi^4 C^2}, \tag{42d}
$$

1706 T. D. Savitsky and D. Toth

where E<sup>P</sup>0,P<sup>ν</sup> (·) denotes the expectation with respect to the joint distribution over population generation and sampling (from that population) without replacement. We apply Jensen's inequality in Equation [42b](#page-28-0) and use E X<sup>2</sup> > Var (X) in the third inequality, stated in Equation [42c,](#page-28-1) and drop the centering term in Equation [42d.](#page-28-2) We now bound the expectation inside the square brackets on the right-hand side of Equation [42d,](#page-28-2) which is taken with respect to this joint distribution. In the sequel, define A<sup>ν</sup> = σ (**X**1,..., **X**N<sup>ν</sup> ) as the sigma field of information potentially available for the N<sup>ν</sup> units in population, Uν.

$$
\mathbb{E}_{P_0, P_{\nu}} \left( \sqrt{N_{\nu}} \mathbb{P}_{N_{\nu}}^{\pi} \log \frac{p}{p_0} \right)^2
$$
\n
$$
= \frac{1}{N_{\nu}} \sum_{i,j \in U_{\nu}} \mathbb{E}_{P_0, P_{\nu}} \left( \frac{\delta_{\nu i} \delta_{\nu j}}{\pi_{\nu i} \pi_{\nu j}} \log \frac{p}{p_0} (\mathbf{X}_i) \log \frac{p}{p_0} (\mathbf{X}_j) \right)
$$
\n
$$
= \frac{1}{N_{\nu}} \sum_{i=j \in U_{\nu}} \mathbb{E}_{P_0} \left[ \mathbb{E}_{P_{\nu}} \left\{ \left( \frac{\delta_{\nu i}}{\pi_{\nu i}^2} \left( \log \frac{p}{p_0} (\mathbf{X}_i) \right)^2 \right) \middle| \mathcal{A}_{\nu} \right\} \right]
$$
\n
$$
+ \frac{1}{N_{\nu}^2} \sum_{i \neq j \in U_{\nu}} \mathbb{E}_{P_0} \left[ \frac{\mathbb{E}_{P_{\nu}} \left[ \delta_{\nu i} \delta_{\nu j} \middle| \mathcal{A}_{\nu} \right]}{\pi_{\nu i} \pi_{\nu j}} \log \frac{p}{p_0} (\mathbf{X}_i) \log \frac{p}{p_0} (\mathbf{X}_j) \right]
$$
\n
$$
= \frac{1}{N_{\nu}} \sum_{i=j \in U_{\nu}} \mathbb{E}_{P_0} \left[ \left( \frac{1}{\pi_{\nu i}} \right) \left( \log \frac{p}{p_0} (\mathbf{X}_i) \right)^2 \right]
$$
\n
$$
+ \frac{1}{N_{\nu}} \sum_{i \neq j \in U_{\nu}} \mathbb{E}_{P_0} \left[ \frac{\pi_{\nu ij}}{\pi_{\nu i} \pi_{\nu j}} \log \frac{p}{p_0} (\mathbf{X}_i) \log \frac{p}{p_0} (\mathbf{X}_j) \right]
$$
\n
$$
\leq \xi^2 \sup_{\nu} \left[ \frac{1}{\min_{i \in U_{\nu}} \pi_{\nu i}} \right] + \xi^2 (N_{\nu} - 1) \sup_{\nu} \max_{i \neq j \in U_{\nu}} \left[ \left| \frac{\pi_{
$$

for sufficiently large Nν, where we have applied the condition for P ∈ B for the first term of the last two inequalities and conditions and [\(A4\)](#page-8-0) and [\(A5\)](#page-8-0) for the last inequality. We additionally note that πνij = πνj when i = j, i, j ∈ Uν. This concludes the proof.

# **References**

- <span id="page-29-1"></span>Barnard, J., McCulloch, R. & Meng, X.-L. (2000), 'Modeling covariance matrices in terms of standard deviations and correlations, with application to shrinkage', Statistica Sinica **10**(4), 1281–1311.
- <span id="page-29-0"></span>Breslow, N. E. & Wellner, J. A. (2007), 'Weighted likelihood for semiparametric models and two-phase stratified samples, with application to cox regression',

Scandinavian Journal of Statistics **34**(1), 86–102. URL: [http://EconPapers.](http://EconPapers.repec.org/RePEc:bla:scjsta:v:34:y:2007:i:1:p:86-102) [repec.org/RePEc:bla:scjsta:v:34:y:2007:i:1:p:86-102](http://EconPapers.repec.org/RePEc:bla:scjsta:v:34:y:2007:i:1:p:86-102)

- <span id="page-30-2"></span>Chambers, R. & Skinner, C. (2003), Analysis of Survey Data, Wiley Series in Survey Methodology, Wiley. URL: [http://books.google.com/books?](http://books.google.com/books?id=4pYGz69d-LkC) [id=4pYGz69d-LkC](http://books.google.com/books?id=4pYGz69d-LkC)
- <span id="page-30-12"></span>Dawid, A. (1981), 'Some matrix-variate distribution theory: Notational considerations and a Bayesian application', Biometrika **68**(1), 265–274.
- <span id="page-30-4"></span>Dong, Q., Elliott, M. R. & Raghunathan, T. E. (2014), 'A nonparametric method to generate synthetic populations to adjust for complex sampling design features', Survey Methodology **40**(1), 29–46.
- <span id="page-30-0"></span>Dunson, D. B. (2010), 'Nonparametric bayes applications to biostatistics', Bayesian nonparametrics **28**, 223–273.
- <span id="page-30-11"></span>Ghosal, S., Ghosh, J. K. & Vaart, A. W. V. D. (2000), 'Convergence rates of posterior distributions', Ann. Statist pp. 500–531.
- <span id="page-30-10"></span>Ghosal, S. & van der Vaart, A. (2007), 'Convergence rates of posterior distributions for noniid observations', Ann. Statist. **35**(1), 192–223. URL: [http://](http://dx.doi.org/10.1214/009053606000001172) [dx.doi.org/10.1214/009053606000001172](http://dx.doi.org/10.1214/009053606000001172)
- <span id="page-30-13"></span>Hoff, P. D. (2011), 'Separable covariance arrays via the tucker product, with applications to multivariate relational data', Bayesian Anal. **6**(2), 179–196. URL: <http://dx.doi.org/10.1214/11-BA606>
- <span id="page-30-1"></span>Holt, D., Smith, T. & Winter, P. (1980), 'A nonparametric method to generate synthetic populations to adjust for complex sampling design features', Journal of the Royal Statistical Society. Series A (General) **143**, 474–487.
- <span id="page-30-5"></span>Kunihama, T., Herring, A. H., Halpern, C. T. & Dunson, D. B. (2014), Nonparametric bayes modeling with sample survey weights, Technical report, Submitted to Biometrika.
- <span id="page-30-3"></span>Little, R. J. (2004), 'To model or not to model? Competing modes of inference for finite population sampling', Journal of the American Statistical Association **99**(466), 546–556.
- <span id="page-30-7"></span>Malec, D., Davis, W. W. & Cao, X. (1999), 'Model-based small area estimates of overweight prevalence using sample selection adjustment', Statistics in Medicine **18**, 3189–3200.
- <span id="page-30-14"></span>Murray, I., Adams, R. P. & MacKay, D. J. (2010), 'Elliptical slice sampling', JMLR: W&CP **9**, 541–548.
- <span id="page-30-8"></span>Pfeffermann, D., Da Silva Moura, F. A. & Do Nascimento Silva, P. L. (2006), 'Multi-level modelling under informative sampling', Biometrika **93**(4), 943– 959.
- <span id="page-30-9"></span>Pfeffermann, D. & Sverchkov, M. (2009), Inference under informative sampling, in D. Pfeffermann & C. Rao, eds, 'Handbook of statistics 29B: sample surveys: inference and analysis', Elsevier Science Ltd., pp. 455–487.
- <span id="page-30-6"></span>Rao, J. N. K. & Wu, C. (2010), 'Bayesian pseudo-empirical-likelihood intervals for complex surveys', Journal of the Royal Statistical Society Series B **72**(4), 533–544. URL: [http://EconPapers.repec.org/RePEc:bla:jorssb:](http://EconPapers.repec.org/RePEc:bla:jorssb:v:72:y:2010:i:4:p:533-544) [v:72:y:2010:i:4:p:533-544](http://EconPapers.repec.org/RePEc:bla:jorssb:v:72:y:2010:i:4:p:533-544)
- <span id="page-30-15"></span>S¨arndal, C.-E., Swensson, B. & Wretman, J. (2003), 'Model assisted survey sampling (springer series in statistics)'.
- <span id="page-31-3"></span>Savitsky, T. (2015), growfunctions: Bayesian Non-Parametric Dependent Models for Time-Indexed Functional Data. R package version 0.12. URL: [https://](https://CRAN.R-project.org/package=growfunctions) [CRAN.R-project.org/package=growfunctions](https://CRAN.R-project.org/package=growfunctions)
- <span id="page-31-4"></span>Savitsky, T. D. (2014), 'Bayesian Non-parametric Mixture Estimation for Timeindexed Functional Data in R', To appear in Journal of Statistical Software.
- <span id="page-31-0"></span>Savitsky, T. D. & Dalal, S. R. (2013), 'Bayesian non-parametric analysis of multirater ordinal data, with application to prioritizing research goals for prevention of suicide', Journal of the Royal Statistical Society: Series C (Applied Statistics). URL: <http://dx.doi.org/10.1111/rssc.12049>
- <span id="page-31-1"></span>Si, Y., Pillai, N. S. & Gelman, A. (2015), 'Bayesian nonparametric weighted sampling inference', Bayesian Anal. **10**(3), 605–625. URL: [http://dx.doi.](http://dx.doi.org/10.1214/14-BA924) [org/10.1214/14-BA924](http://dx.doi.org/10.1214/14-BA924)
- <span id="page-31-2"></span>Wong, W. H. & Shen, X. (1995), 'Probability inequalities for likelihood ratios and convergence rates of sieve mles', Ann. Statist. **23**(2), 339–362. URL: <http://dx.doi.org/10.1214/aos/1176324524>