*International Statistical Review* (2021), 89, 1, 72[–107](#page--1-0) doi:10.1111/insr.12376

# **Uncertainty Estimation for Pseudo-Bayesian Inference Under Complex Sampling**

## **Matthew R. Williams1 [a](http://orcid.org/0000-0001-8894-1240)nd Terrance D. Savitsky2**

<sup>1</sup>*National Center for Science and Engineering Statistics, National Science Foundation, Alexandria, Virginia 22314, USA* <sup>2</sup>*Office of Survey Methods Research, U.S. Bureau of Labor Statistics, Washington, DC 20212-0001, USA E-mail: matthew.dunn.williams@gmail.com*

#### **Summary**

**Social and economic studies are often implemented as complex survey designs. For example, multistage, unequal probability sampling designs utilised by federal statistical agencies are typically constructed to maximise the efficiency of the target domain level estimator (e.g. indexed by geographic area) within cost constraints for survey administration. Such designs may induce dependence between the sampled units; for example, with employment of a sampling step that selects geographically indexed clusters of units. A sampling-weighted pseudo-posterior distribution may be used to estimate the population model on the observed sample. The dependence induced between coclustered units inflates the scale of the resulting pseudo-posterior covariance matrix that has been shown to induce under coverage of the credibility sets. By bridging results across Bayesian model misspecification and survey sampling, we demonstrate that the scale and shape of the asymptotic distributions are different between each of the pseudo-maximum likelihood estimate (MLE), the pseudo-posterior and the MLE under simple random sampling. Through insights from survey-sampling variance estimation and recent advances in computational methods, we devise a correction applied as a simple and fast postprocessing step to Markov chain Monte Carlo draws of the pseudo-posterior distribution. This adjustment projects the pseudo-posterior covariance matrix such that the nominal coverage is approximately achieved. We make an application to the National Survey on Drug Use and Health as a motivating example and we demonstrate the efficacy of our scale and shape projection procedure on synthetic data on several common archetypes of survey designs.**

*Key words*: cluster sampling; credible set; Markov chain Monte Carlo; multistage sampling; pseudo-posterior distribution; sampling weights; survey sampling.

#### **1 INTRODUCTION**

We focus on the task of the data analyst to estimate a Bayesian model, P-<sup>0</sup> , that they suppose generates values for a random variable, Y , for units of a population, U D .1;:::;N/, from an observed sample, S D .1;:::;n - N /, drawn from that population under a complex sampling design governed by distribution, P. Common designs include one-stage samples of businesses and multistage samples of households (HHs). These two major classes of sampling designs cover the majority of survey designs and represent a standard data collection method for social and economic data. While our methods and simulations will apply to both of these and to other types, we will focus on multistage designs for our motivating example.

In this section, we introduce the concepts and current state of the literature for implementing Bayesian models for survey data. We also provide an introduction to variance estimation for the survey-weighted maximum likelihood estimator and to model misspecification in Bayesian models. As we will demonstrate, combining results across these different areas provides an opportunity to innovate and improve on current methods. In particular, we can formally demonstrate a deficiency in the coverage of Bayesian credible intervals using the popular survey-weighted pseudo-posterior as well as provide an efficient and effective correction, achieving nominal coverage asymptotically and in practice for moderate sample sizes. With this adjustment, the pseudo-posterior now possesses many desirable theoretical properties. Unlike other methods in the literature, the pseudo-posterior can be implemented across a very wide variety of survey designs and analyst-specified population models. Even though the adjustment is asymptotic, the resulting posterior still retains its small sample characteristics (i.e. is not forced to be asymptotically normal).

#### *1.1 Informative Sampling Designs*

Multistage sampling designs are created to achieve efficient (low variance) estimation of a desired simple quantile, mean, or total estimator for a collection of domains within constraints on cost to administer the survey. A first stage of the sampling design often collects contiguous geographic areas from the population into clusters, where a subset of the clusters are randomly selected into the sample for this first stage. The contiguity of areas within each cluster is defined for convenience and cost to collect the sample, but it induces a dependence among units nested in areas within each cluster. Dependencies among sampled units may be additionally promulgated through the drawing of a fixed-sized sample, without replacement, in any stage of the sampling design, for example, by constructing a systematic sampling step with a fixed interval using a random starting point, as is common when sampling HHs on the same street.

The sampling design distribution is induced by specifying marginal inclusion probabilities at each stage. Survey agencies, such as federal statistical agencies, publish marginal inclusion probabilities for last-stage sampled units, <sup>i</sup> D P r fı<sup>i</sup> D 1 j Ag 2 .0; 1, for (observed) units, i D 1;:::;n, sampled in the last stage of the sampling design, where n denotes the number of units in the observed sample, ı<sup>i</sup> 2 f0; 1g specifies a unit inclusion indicator and A represents 'prior' informatio[n \(Isaki & Fuller, 198](#page-25-0)2) about the population that is then used to create the sample design (e.g. firm size or geographic information). In the traditional perspective where the population is assumed to be fixed, the conditional inclusion probabilities are simply treated as the unconditional or marginal inclusion probabilities P r fı<sup>i</sup> D 1 j Ag D P r fı<sup>i</sup> D 1g. To maintain generality, we retain the A notation throughout.

Efficiency of the population estimator, g.Y /, is enhanced through designing the inclusion probabilities, .<sup>i</sup>/iD1;:::;N , to be correlated with .yi/iD1;:::;N , where N denotes the size of population, U; an example is the use of a proportion-to-size sampling design in the Current Employment Statistics survey of business establishments, administered by the U.S. Bureau of Labor Statistics, for the purpose of measuring total employment by geographic area and industry. Higher unit inclusion probabilities are assigned to larger employers because they drive the variance of the resulting total employment estimator. Sampling designs that induce this correlation are termed 'informative' and the balance of information in the sample is different from that in the population.

Analyzing the data 'as is' assuming a simple random sampling (SRS) mechanism will then lead to bias. The most common approach in the survey literature is to use a survey-weighted likelihood. This is a plug-in estimator that formulates a sampling-weighted pseudo-likelihood density by exponentiating each (last stage) unit-indexed likelihood contribution by a sampling weight constructed to be inversely proportional to the unit marginal inclusion probability, w<sup>i</sup> / 1=<sup>i</sup> , where <sup>i</sup> D P .ı<sup>i</sup> D 1 j A/; for units, i D 1;:::;n, where n denotes the number of units in the observed sample. An approximate, weight-exponentiated pseudo-likelihood for the population, <sup>Q</sup><sup>n</sup> <sup>i</sup>D<sup>1</sup> p .yij/ <sup>w</sup><sup>i</sup> , is constructed from the n units observed in the sample. The survey weights mitigate the estimation bias, but the uncertainty distribution (covariance structure) also needs to be estimated and adjusted.

#### *1.2 Variance Estimation for Survey-Weighted Maximum Likelihood*

When establishing the consistency of the survey-weighted pseudo-maximum likelihood estimate (MLE) O ML, it is almost ubiquitous in the literature (I[saki & Fuller, 1982; P](#page-25-0)feffermann *et al.*, 1998; [C](#page-25-1)[hambers & Skinner, 2003\) to](#page-24-0) [assume](#page-25-1) [that](#page-25-1) [the](#page-25-1) [pairwise](#page-25-1) [inclusion](#page-25-1) [probab](#page-25-1)ilities ii<sup>0</sup> D P r fı<sup>i</sup> D 1; ıi<sup>0</sup> D 1 j Ag 2 .0; 1 factor or are independent asymptotically (ii<sup>0</sup> D i<sup>i</sup><sup>0</sup>). However, the typical assumption for variance estimation (of the same models) is to assume an arbitrary amount of within cluster dependence both in the sampling design and the populationgenerating model (Heeringa *et al.*[, 2010; R](#page-25-2)ao *et al.*[, 1992\). R](#page-25-3)ecent work (Williams & Savitsky, 2018a) [provides](#page-25-4) [the](#page-25-4) [theoretical](#page-25-4) [conditions](#page-25-4) [and](#page-25-4) [motivating](#page-25-4) [examples](#page-25-4) [to](#page-25-4) [exte](#page-25-4)nd our understanding of consistency to these real-world situations for which there are already well-established procedures for estimating variance.

Wang *et al.* [\(2017\) s](#page-25-5)uggest estimating variance based on the Horvitz–Thompson estimator for weighted sums. Theoretically this is appealing, but in practice, it requires knowledge of second order inclusion probabilities ii<sup>0</sup>. These joint probabilities are rarely calculated and are disseminated even less frequently. Instead, the de facto approach for variance estimation is based on the approximate sampling independence of the primary sampling units (PSUs) (Heeringa *et al.*, 2010). [Variance](#page-25-2) [estimation](#page-25-2) [can](#page-25-2) [be](#page-25-2) [in](#page-25-2) [the](#page-25-2) [form](#page-25-2) [of](#page-25-2) [Taylor](#page-25-2) [linearisation](#page-25-2) [or](#page-25-2) [replicat](#page-25-2)ion-based methods. We provide a high-level overview of the two below, but a variety of implementations are available ([Binder, 1996; R](#page-24-1)ao *et al.*[, 1992\). L](#page-25-3)et yij , Xij and wij be the observed data for individual i in cluster j of the sample. Assume the parameter is a vector of dimension d with population model value 0.

- 1. Taylor Linearisation
	- (a) Approximate an estimate O , or a 'residual' .O 0/, as a weighted sum: <sup>O</sup> <sup>P</sup> i;j wij ´ij . / where ´ij is a function evaluated at the current values of yij , Xij and O . (c) Compute the variance between clusters: Var. <sup>2</sup><sup>O</sup>
	- (b) Compute the weighted components for each primary cluster (e.g. PSUs): <sup>O</sup> <sup>P</sup> <sup>j</sup> D <sup>i</sup> wij ´ij . /. (d) For stratified designs, compute <sup>O</sup><sup>s</sup> and Var. <sup>2</sup><sup>O</sup><sup>s</sup> / within strata and sum Var. <sup>2</sup><sup>O</sup>
	- / D <sup>1</sup> Jd PJ <sup>j</sup>D1.<sup>O</sup> O <sup>j</sup> /.O O <sup>j</sup> /<sup>T</sup>
	- / <sup>D</sup> <sup>P</sup> s Var. 2Os /.
- 2. Replication
	- (a) Through randomisation (bootstrap), leave-one-out (jackknife) or orthogonal contrasts (balanced repeated replicates), create a set of K replicate weights .wi/<sup>k</sup> for all i 2 S and for every k D 1;:::;K.
- (b) Each set of weights has a modified value (usually 0) for a subset of clusters and typically has a weight adjustment to the other clusters to compensate: <sup>P</sup> i2S P .wi/<sup>k</sup> D <sup>i</sup>2<sup>S</sup> w<sup>i</sup> for every k. (d) Compute the variance between replicates: Var. <sup>2</sup><sup>O</sup>
- (c) Estimate O <sup>k</sup> for each replicate k 2 1;:::;K.
- / D <sup>1</sup> P<sup>K</sup> <sup>k</sup>D1.<sup>O</sup> O k/.O O k/<sup>T</sup> .
- Kd (e) For stratified designs, generate replicates such that each strata is represented in every replicate.

There are two notable challenges associated with these methods:

- For Taylor linearisation, the value of O is typically only computed *once* and then used in a plug-in value for ´i. /. Whereas for the replication methods, the estimate O <sup>k</sup> must be computed K times. This may lead to a sizable differences in *computational effort* for models of moderate complexity and a moderate number of replicates K.
- For the replication methods, no additional derivatives are needed. In contrast, the Taylor linearisation method often requires the calculation of a gradient based on the estimating equation for O to derive the analytical form of the first order approximation ´i. /. This poses significant *analytical challenges* for all but the simplest models.

In this work, we present a third option which can be implemented as a hybrid of the two approaches. This approach can then be applied to any sampling design, as long as the analyst has either the replication weights or the cluster and stratification information. The implementation (Secti[on](#page-13-0) 4) is made possible by recent advances in algorithmic differentiation (Margossian, 2018), which allows us to specify the model as a log density but only treat the gradient in the abstract *without* specifying it analytically.

#### *1.3 Bayesian Models for Survey Data*

#### *1.3.1 Generalised estimating equations and method of moments*

[Yin \(200](#page-25-6)9) proposes using the normal distribution with parameters from the solution of generalised method of moments to generate Markov chain Monte Carlo (MCMC) posterior draws. An analogous approach with survey-weighted estimating equations was proposed and implemented [by Shah](#page-25-7) *et al.* (2000) to estimate a logistic mixed-effects model for survey data, where the fixed effects were sampled from a normal distribution with mean estimates from a survey-weighted likelihood maximisation and variance estimates from a Taylor series linearisation approac[h. Wang](#page-25-5) *et al.* (2017) formalises the use of the generalised method of moments approach for survey data and provides a Bernstein–von Mises result, demonstrating that the resulting posterior intervals achieve correct frequentist coverage asymptotically. These methods assume asymptotic normality when making posterior draws. While large sample performance is generally good, small sample properties and performance are uncertain. By contrast, our approach addresses Bayesian estimation performed by the data analyst under a model of their choosing. The data analyst will extract numerical draws from the joint distribution for their model parameters, from which a variety of statistics of interest may be computed (e.g. quantiles or probabilities for events). Our method performs a rescaling of draws from the marginal or joint posterior distributions obtained from the model and MCMC algorithm specified by the data analyst to produce credibility sets with asymptotically correct frequentist coverage. The approaches [of Shah](#page-25-7) *et al.* (2000) a[nd Wang](#page-25-5) *et al.* (2017), by contrast, require the availability of an argmax of the posterior (maximum a posterior probability) point estimate obtained from

an estimating equation, which may not be reliably computable for a complicated hierarchical Bayesian model specification. These approaches also circumvent the estimation of the Bayesian model specified by the data analyst, so they may not address the small sample borrowing of strength provided by a Bayesian specification as does our method.

#### *1.3.2 Likelihood approaches*

Savitsky and Toth (2[016\) p](#page-25-8)roposed a plug-in estimator that formulates a sampling-weighted pseudo-posterior density that is analogous to the weighted pseudo-likelihood. However, the sampling weights, .wi/, are normalised, to .wQ <sup>i</sup>/, to control the amount of estimated posterior uncertainty. Savitsky and Toth (2[016\) d](#page-25-8)efault to normalizing, <sup>P</sup><sup>n</sup> <sup>i</sup>D<sup>1</sup> wQ <sup>i</sup> D n. León-Novelo & Savitsky (2019) demonstrate that the pseudo-posterior estimator constructed from weights normalised to n generally produce credibility intervals that fail to contract on frequentist confidence sets by under covering because they do not account for dependencies among units induced by the joint distribution .P-<sup>0</sup> ; P/. L[eón-Novelo & Savitsky \(2019\) d](#page-25-9)evelop an alternative approach to the pseudo-posterior distribution that multiplicatively adjusts the likelihood to accomplish asymptotically unbiased estimation of the population model on the observed informative sample. This extends the formulation of the observed likelihood by Pfeffermann *et al.* (1998) to a fully Bayesian implementation by specifying a conditional population model, p.<sup>i</sup> j yi;/, for the inclusion probabilities, .<sup>i</sup>/i2<sup>U</sup> . P[feffermann](#page-25-10) *et al.* (2006) also extend P[feffermann](#page-25-1) *et al.* (1998) to a partially Bayesian estimation, but they treat .<sup>i</sup>/i2<sup>U</sup> as *fixed*. So the approach of P[feffermann](#page-25-10) *et al.* (2006) may be viewed to be *not* fully Bayesian, because it does not specify a joint or conditional model for p.<sup>i</sup>jyi;/.

L[eón-Novelo & Savitsky \(2019\) s](#page-25-9)how that credible intervals estimated from their adjusted, fully Bayes posterior achieves correct coverage in the case of a simple, single stage proportion to size sampling design. Their likelihood adjustment, however, requires a different MCMC sampler than that developed for the population model and the adjusted likelihood includes an integration that must be numerically computed in each MCMC draw. So the fully Bayesian estimator lacks the broad ease-of-implementation of the pseudo-posterior approach.

R[ao & Wu \(2010\) a](#page-25-11)lso address the under coverage of the pseudo-posterior in the specific case of formulating Pg.Y / as an empirical likelihood for the purpose of estimating a total or mean, g.Y / <sup>1</sup>. They replace <sup>n</sup> as the normaliser for .w<sup>Q</sup> <sup>i</sup>/ with <sup>n</sup> <sup>D</sup> n=DEFFg.Y / <sup>b</sup>, where DEFFg.Y / <sup>b</sup> <sup>D</sup> VarP- Œg.Y /= <sup>1</sup> VarSRSŒg.Y / <sup>1</sup> denotes the design effect (DEFF), defined as the variance induced under sampling design distribution, P, divided by that under SRS. Their approach improves the coverage properties for estimation of simple statistics, rather than some of interest to the data analyst for a general P- . In addition, the simultaneous modelling of multiple outcomes or parameters would require multiple DEFFs to be used, which is not possible if DEFF is only incorporated via scaling the sample size n.

In the sequel, we demonstrate that a postprocessing adjustment step applied to the pseudoposterior MCMC samples corrects the under coverage demonstrated by León-Novelo & Savitsky (2019). The fully Bayes approach will tend to produce more efficient credible sets, however, under the requirement to specify a conditional population model for the inclusion probabilities that is *assumed* to be correctly specified. In practice, sample designs are often algorithmically defined, becoming quite complex. The fully Bayes approach has not been applied to multistage cluster designs. The impact of clustering on the effective sample size may still be a challenge. In this work, we demonstrate that the survey-weighted pseudo-posterior can be adjusted to give correct inference even under complex survey designs that include within-cluster dependence.

#### *1.4 Misspecification for Bayesian models*

When a Bayesian model has a misspecified likelihood, the resulting posterior distribution contracts on an alternative distribution that is the minimum Kullback–Leibler distance from the true generating distributio[n \(Kleijn & van der Vaart, 201](#page-25-12)2). However Kleijn & van der Vaart (2012) demonstrate that Bayesian credible sets from the posterior are not valid confidence sets. In other words, under misspecification, the MLE and the posterior distribution for the misspecified model have different limiting distributions. This suggests that an adjustment is needed to achieve valid confidence sets asymptoticall[y. Ribatet](#page-25-13) *et al.* (2012) motivate a similar sandwich form of an adjustment of the asymptotic covariance of the pseudo-posterior distribution under specification of a composite weight-exponentiated pseudo-likelihood, where their pseudo-likelihood is employed to approximate a likelihood that is not able to be specified. Ribatet *et al.* (2012) redesign the MCMC sampler to accomplish the adjustment, such that their approach requires the development of a specialised MCMC sampler, distinct from the sampler developed for the population, P-0 .

Our survey-sampling formulation assumes existence of a population model, P-<sup>0</sup> , which, though unknown, has a tractable form that allows consistency of our estimator, P- . Even though consistency is achieved, the survey-weighted pseudo-posterior is still misspecified because the exponentially weighted likelihood is a noisy approximation to the true likelihood of the joint distribution .P-<sup>0</sup> ; P/. The sampling-weighted pseudo-posterior arises out of a random sampling mechanism to approximate the information in the population using a partially observed sample taken from that population. In contra[st, Ribatet](#page-25-13) *et al.* (2012) do not compute expectations with respect to the joint distribution, .P-<sup>0</sup> ; P/, to develop their adjustment because they do not contemplate a random sampling process governed by P. We provide theoretical results for the form of the asymptotic sampling-weighted pseudo-MLE covariance matrix under the joint distribution for population generation and the taking of a sample. We also provide a clean post-estimation adjustment that allows for minimal changing of a user's intended MCMC implementation.

#### *1.5 Adjusting the Distribution of the Pseudo-Posterior*

The current work constructs a simple postprocessing step that adjusts the scale and shape of sampling-weighted, pseudo-posterior parameter credibility sets that we show in the sequel achieves approximately correct coverage under a broad class of generally-used sampling designs. Our procedure applies an adjustment step to the posterior draws to achieve an asymptotic sandwich form for the pseudo-posterior covariance that is the same as that for the sampling-weighted pseudo-MLE. We accomplish the adjustment by computing the variance of the score function and the expectation of the square of its gradient under the joint distribution, .P-<sup>0</sup> ; P/. The variance of the gradient is estimated via a hybrid approach combining the principles of Taylor linearisation and replication. The implementation (Secti[on](#page-13-0) 4) is through the use of algorithmic differentiatio[n \(Margossian, 201](#page-25-14)8) and randomised replication sampling [\(Preston, 200](#page-25-15)9).

The adjustment step is applied, numerically, by resampling the observed data, y1;:::;yn, under an empirical distribution approximation for .P-<sup>0</sup> ; P/. The resampling step is implemented by simply drawing blocks of units from the existing sample at those stages where dependence is induced within the blocks. All units nested within each re-sampled block are included in each re-sample; for example, if the multistage design includes a clustering step, we use the known cluster memberships of the last stage units and just re-sample the clusters. The population-generating distribution, P-<sup>0</sup> , is estimated, once, on our original sample, and the adjustment is evaluated using the best available estimate for , the posterior mean. Our adjustment is, therefore, computationally fast and achieves nearly correct coverage for . The pseudo-posterior MCMC sampler, used for estimation of , requires only a simple edit to the population posterior sampler (to insert sampling weights) because the same posterior geometry is employed. Our adjustment procedure requires *no* change to the MCMC sampler for the pseudo-posterior, which preserves its ease of use.

#### <span id="page-6-0"></span>**2 Motivating Multistage Cluster Design: The National Survey on Drug Use and Health**

Our motivating survey design is the National Survey on Drug Use and Health (NSDUH), sponsored by the Substance Abuse and Mental Health Services Administration. NSDUH is the primary source for statistical information on illicit drug use, alcohol use, substance use disorders, mental health issues and their co-occurrence for the civilian, noninstitutionalised population of the United States. The NSDUH employs a multistage state-based design (Morton *et al.*, 2016), with the earlier stages defined by geography within each state in order to select HHs (and group quarters) nested within these geographically-defined PSUs. Williams & Savitsky (2018a) provides conditions for asymptotic consistency for the pseudo-posterior for designs like the NSDUH, which are characterised by:

- Cluster sampling, such as selecting only one unit per cluster, or selecting multiple individuals from a dwelling unit.
- Population information such as socio-economic indicators used to sort sampling units along gradients.

Both features are common, in practice, and create sampling dependencies that do not attenuate even if the population grows. For simplicity of exposition, we examine the relationship between two measures, current (past month) smoking of cigarettes and past year major depressive episode for adults, through a two-parameter logistic regression model. Both cigarette smoking and depression may be clustered geographically and within HHs. For example, rates for each tend to vary by age, urban/rural status, education and other demographics that typically cluster geographically and within HH (Center for Behavioral Health Statistics and Quality, 2015b; [2015a\).](#page-24-2)

#### <span id="page-6-1"></span>**3 Asymptotic Covariance of the Pseudo-Posterior Distribution**

#### *3.1 The Pseudo-Posterior Framework*

We suppose random variables of the population are generated, **X** D .**X**1;:::; **X**N- / ind P-0 where <sup>0</sup> 2 R<sup>d</sup> and we perform inference on 2 ' of the population model from the sample of size, n. A sampling design imposes a *known* distribution (given information A) on a vector of random inclusion indicators, ı D .ı1;:::;ıN- /, and on units composing a population, U. The sampling distribution takes an *observed* random sample, S U, of size n - N from U. Our conditions for the main results are based on marginal unit inclusion probabilities, i D Prfıi D 1 j Ag for all i 2 U and the second order pairwise probabilities, ij D Prfıi D 1 \ ıj D 1 j Ag for i; j 2 U, which are obtained from the joint distribution over .ı1;:::;ıN-/. We denote the sampling distribution by P, which governs the taking of samples from the population. P is implicitly conditionally defined given realisations from P-0 In other words, the joint distribution for ı can depend on some population information from **X**, which we note by the use of A to indicate information about the population that is available to the survey designer.

We denote the observed sample of size n as f**X**; ıg D .f**X**1; ı1g;:::; f**X**N- ; ıN g/ , following Savitsky and Tot[h \(201](#page-25-8)[6\), Savitsky & Srivastava \(201](#page-25-16)8) and Williams & Savitsky (2018a), where ıi D 0 indicates unit i is *not* included in the sample removes the associated **X**i . It is a notational convention that emphasises the dependence of generated samples on both P (which governs ıi) and P-<sup>0</sup> (which governs the generation of population values, **X**i ). Because the ı are random with respect to P, f**X**; ıg is jointly random with respect to .P; P-0 /.

The inclusion probabilities are formulated to depend on the finite population data values, **X**, so that we employ the pseudo-posterior estimator to approximate the population likelihood from the observed sample with,

<span id="page-7-0"></span>
$$
p^{\pi}\left(\mathbf{X}_{\nu i},\delta_{\nu i}\right) := p\left(\mathbf{X}_{\nu i}\right)^{\delta_{\nu i}/\pi_{\nu i}}, \quad i \in U_{\nu},\tag{1}
$$

which weights each density contribution, p.**X**i/, by the inverse of its marginal inclusion probability (Savitsky and Tot[h, 201](#page-25-8)6). When ıi D 0, the pseudo-posterior likelihood contribution for unit i under conditional independence (given ) is removed. This approximation for the population likelihood produces the associated pseudo-posterior density,

$$
p^{\pi}(\theta \mid \mathbf{X}_{\nu}, \boldsymbol{\delta}_{\nu}) = \frac{\prod_{i \in U_{\nu}} p_{\theta}^{\pi}(\mathbf{X}_{\nu i}, \delta_{\nu i}) \pi(\theta)}{\int_{\Theta} \prod_{i \in U_{\nu}} p_{\theta}^{\pi}(\mathbf{X}_{\nu i}, \delta_{\nu i}) \pi(\theta) d\theta},
$$
\n(2)

where **X**; ı D .**X**1; ı1;:::; **X**N- ; ıN- / denotes the observed sample of size, n. The pseudo-posterior mass placed on subset B ' becomes

$$
\Pi^{\pi} (B | \mathbf{X}_{\nu}, \delta_{\nu}) = \int_{\theta \in B} p^{\pi} (\mathbf{X}_{\nu}, \delta_{\nu} | \theta) \pi(\theta) d\theta.
$$
 (3)

In typical application[s \(Savitsky & Srivastava, 201](#page-25-16)8), sampling weights are normalised to satisfy <sup>P</sup> i2S- --1 i D n, which regulates the scale of uncertainty in the estimated pseudoposterior distribution. In practice, dependencies induced by informative, multistage sampling designs produce a smaller effective sample size than n, such that the typical procedure underestimates posterior uncertainty. In addition, the shape (geometry) of the pseudo-posterior distribution is impacted by the dependence induced in each stage of the sampling design such that the asymptotic covariance matrix will not be the same as that for the MLE obtained under SRS. We proceed to derive the form of the limiting covariance matrix for the pseudo-MLE under informative sampling, which we define as the MLE of Equatio[n \(](#page-7-0)1). We demonstrate that the covariance matrix of the pseudo-MLE is different from that for the MLE under SRS but that the latter is a special case of the former. We next demonstrate that the limiting covariance matrix of the pseudo-posterior distribution differs from the pseudo-MLE under informative sampling (because of the failure of Bartlett's second identity) such that resulting credibility intervals would not be expected to contract on valid frequentist confidence intervals, absent adjustment.

The difference between the limiting covariance matrix for the pseudo-posterior distribution, on the one hand, from that for the MLE under SRS, on the other hand, may only be partly driven by informativeness of the sampling design. The dependencies induced under employment of a multistage sampling design, such as the within cluster dependence of units, will also impact the scale of the limiting covariance matrix of the pseudo-posterior distribution, even absent

sampling informativeness. In other words, even where sampling inclusion probabilities, .i/, are not required to provide unbiased estimation of 2 ' (i.e. the design is 'ignorable'), the resulting limiting covariance matrix of the posterior distribution under multistage sampling would be different from that for the MLE under SRS.

Our main result is achieved in the limit as " 1, under the countable set of successively larger-sized populations, fUg2Z<sup>C</sup>. The asymptotics under our construction is controlled by 2 N to map to the process where we fix a , construct an associated finite population of size, N, generate random variables **X**1;:::; **X**Nind P-<sup>0</sup> , construct unit marginal sample inclusion probabilities, .-1;:::;-N / (and other design features such as cluster and strata identifiers) under P, and then draw a sample, f1;:::;ng from that population. The process is repeated for each increment of . We define the associated stochastic rates of convergences notations, A D o<sup>P</sup> .B/ to denote that A D YB where Y !<sup>P</sup> 0 and A D O<sup>P</sup> .B/ denotes A D YB where Y D O<sup>P</sup> .1/. For deterministic sequences, A and B, the notations reduce to the usual o and O.

#### <span id="page-8-1"></span>*3.2 Review of Survey-Weighted Empirical Functionals*

We use the empirical distribution approximation for the joint distribution over population generation and the draw of an informative sample that produces our observed data. Our empirical distribution construction follows B[reslow & Wellner \(2007\) a](#page-24-3)nd incorporates inverse inclusion probability weights, f1=i giD1;:::;N-, to account for the informative sampling design,

$$
\mathbb{P}_{N_{\nu}}^{\pi} = \frac{1}{N_{\nu}} \sum_{i=1}^{N_{\nu}} \frac{\delta_{\nu i}}{\pi_{\nu i}} \delta\left(\mathbf{X}_{\nu i}\right),\tag{4}
$$

where ı .**X**i/ denotes the Dirac delta function, with probability mass 1 on **X**i and we recall that N D jUj denotes the size of the finite population. This construction contrasts with the usual empirical distribution, PN- D <sup>1</sup> Nv PN<sup>i</sup>D<sup>1</sup> ı .**X**i/.

We will construct asymptotic distributions for the sequence of centered and scaled random quantities,

$$
h_{N_v} = \sqrt{N_v} \left( \theta - \theta_0 \right), \tag{5}
$$

for specific estimators. Let O ;N- D arg max-P i2Uıi <sup>i</sup> logfp- .Xi/g denote the MLE of the pseudo-likelihood in Equation (1[\) \(](#page-7-0)that we denote as the pseudo-MLE). It is the MLE of the logarithm of the sample-weighted likelihood for an observed sample (where ıi D 1 for those units included in the observed sample). The pseudo-MLE defines the sequence,

<span id="page-8-0"></span>
$$
\hat{h}_{N_{\nu}}^{\pi} = \sqrt{N_{\nu}} \left( \hat{\theta}_{\pi, N_{\nu}} - \theta_0 \right), \tag{6}
$$

as contrasted with a centered and scaled sequence for the MLE, O N- , for the population (as if fully-observed),

$$
\hat{h}_{N_{\nu}} = \sqrt{N_{\nu}} \left( \hat{\theta}_{N_{\nu}} - \theta_0 \right).
$$
\n(7)

Define the log-likelihood, `- D log p- D log p-0ChN- = <sup>p</sup>N and the associated score function, P `- D r- `- . Equation (6[\) i](#page-8-0)s scaled by N because we later sum empirical expectations and variances over N sampling-weighted units in the population, where ı indexes all possible realisable samples (B[reslow & Wellner, 2007\).](#page-24-3)

We follow the notational convention [of Ghosal](#page-25-17) *et al.* (2000) and define the associated expectation functionals with respect to these empirical distributions by P N f D 1 N-P<sup>N</sup>iD1 ıi <sup>i</sup> f .**X**i/. Similarly, P<sup>N</sup> f D <sup>1</sup> N-P<sup>N</sup><sup>i</sup>D<sup>1</sup> f .**X**i/. Lastly, we use the associated centered empirical processes, G N- <sup>D</sup> <sup>p</sup>N -P N- P<sup>0</sup> and G<sup>N</sup>- <sup>D</sup> <sup>p</sup>N .P<sup>N</sup>-P0/.

We construct two variance expressions, starting with Fisher's information:

$$
H_{\theta_0} = -\frac{1}{N_{\nu}} \sum_{i \in U_{\nu}} \mathbb{E}_{P_{\theta_0}} \ddot{\ell}_{\theta_0}(\mathbf{X}_{\nu i}), \tag{8}
$$

whose inverse provides the asymptotic covariance of the pseudo-posterior under our Bernstein Von–Mises result that follows. Next, we define

$$
J_{\theta_0} = \frac{1}{N_{\nu}} \sum_{i \in U_{\nu}} \mathbb{E}_{P_{\theta_0}} \dot{\ell}_{\theta_0}(\mathbf{X}_{\nu i}) \dot{\ell}_{\theta_0}(\mathbf{X}_{\nu i})^T, \tag{9}
$$

which is the middle term in the asymptotic variance of the MLE. Under the population model (and an SRS subsample), the likelihood is properly specified, so J-<sup>0</sup> D H-0 .

Because our pseudo-posterior framework arises from a random sampling process governed by P,

$$
H_{\theta_0}^{\pi} = -\mathbb{E}_{P_{\theta_0}, P_{\nu}} \left[ \mathbb{P}_{N_{\nu}}^{\pi} \ddot{\ell}_{\theta_0} \right]
$$
  
\n
$$
= -\frac{1}{N_{\nu}} \sum_{i \in U_{\nu}} \mathbb{E}_{P_{\theta_0}} \left[ \mathbb{E}_{P_{\nu}} \frac{[\delta_{\nu i}|A_{\nu}]}{\pi_{\nu i}} \ddot{\ell}_{\theta_0}(\mathbf{X}_{\nu i}) \right]
$$
  
\n
$$
= -\frac{1}{N_{\nu}} \sum_{i \in U_{\nu}} \mathbb{E}_{P_{\theta_0}} \ddot{\ell}_{\theta_0}(\mathbf{X}_{\nu i})
$$
  
\n
$$
= H_{\theta_0},
$$

where A denotes the sigma field of information in U (i.e. the information available to the survey designer). We note that this equivalence between H -<sup>0</sup> and H-<sup>0</sup> does not hold for the weighted composite likelihood [of Ribatet](#page-25-13) *et al.* (2012), where the weights are arbitrary and arise from a deterministic process to approximate an intractable likelihood for the population, U.

Our main results in the following section are anchored in the observation that the surveyweighted J -<sup>0</sup> D EP<sup>0</sup> ;Ph P N-R `-0 R `T -0 i ¤ J-<sup>0</sup> because of the misspecification from using a noisy approximation to the likelihood for .P-<sup>0</sup> ; P/.

#### <span id="page-9-0"></span>*3.3 Main Results*

The following conditions guarantee three results on the forms for asymptotic covariance matrices of the distributions for the pseudo-MLE estimator and the pseudo-posterior. The first theorem extends Theorem 5.23 [of van der Vaart \(199](#page-25-18)8) to derive the asymptotic expansion of the centered and scaled pseudo-MLE. The second theorem specifies the form of the associated sandwich covariance matrix for the (asymptotic expansion of the) pseudo-MLE. The third theorem extends similar theorems [in Kleijn & van der Vaart \(201](#page-25-12)2) a[nd van der Vaart \(199](#page-25-18)8) that specify the covariance matrix of the asymptotic Gaussian form for the pseudo-posterior distribution. We observe that the asymptotic covariance matrices are different for each of the MLE, the pseudo-MLE and the pseudo-posterior, which sets up our proposed scale and shape adjustment, introduced in the sequel.

*(A1)* (Continuity) For each 2 ' 2 R<sup>d</sup> (an open subset of Euclidean space), `-<sup>0</sup> .**x**/ be a measurable function (of **x**) and differentiable at <sup>0</sup> for P-<sup>0</sup> almost every **<sup>x</sup>** [with derivative, <sup>P</sup> `-<sup>0</sup> .**x**/], such that for every <sup>1</sup> and <sup>2</sup> in a neighborhood of <sup>0</sup> with E- P `-<sup>0</sup> .**x**/ <sup>P</sup> `-<sup>0</sup> .**x**/ <sup>T</sup> < 1, we have a Lipschitz condition:

$$
\left|\ell_{\theta_1}(\mathbf{x})-\ell_{\theta_2}(\mathbf{x})\right|\leq \dot{\ell}_{\theta_0}(\mathbf{x})\left\|\theta_1-\theta_2\right\| \text{a.s. } P_{\theta_0}.
$$

*(A2)* (Local quadratic expansion) The Kullback–Liebler divergence with respect to P-<sup>0</sup> has a second order Taylor expansion about 0,

$$
\mathbb{E}_{P_{\theta_0}} \log \frac{p_{\theta}}{p_{\theta_0}} = \frac{1}{2} \left( \theta - \theta_0 \right)^T H_{\theta_0} \left( \theta - \theta_0 \right) + o \left( \|\theta - \theta_0\|^2 \right),
$$

where H-<sup>0</sup> is a d d positive definite matrix. *(A3)* (Bartlett's first identity)

$$
\mathbb{E}_{P_{\theta_0}}\dot{\ell}_{\theta_0}=0
$$

*(A4)* (Consistency of the MLE for the population)

$$
\mathbb{P}_{N_{\nu}} \ell_{\hat{\theta}_{N_{\nu}}} \geq \sup_{\theta} \mathbb{P}_{N_{\nu}} \ell_{\theta} - o_{P_{\theta_{0}}}\left(N_{\nu}^{-1}\right)
$$

and O N-P<sup>0</sup> ! 0:

*(A5)* (Nonzero inclusion probabilities)

$$
\sup_{\nu} \left[ \frac{1}{\min_{i \in U_{\nu}} |\pi_{\nu i}|} \right] \leq \gamma, \text{with } P_{\theta_0} - \text{probability 1.}
$$

*(A6)* (Growth of dependence is restricted)

For every U there exists a binary partition fS1; S2g of the set of all pairs S D ffi; j g W i ¤ j 2 Ug such that

$$
\limsup_{\nu\uparrow\infty} |S_{\nu 1}| \leq \mathcal{O}(N_{\nu}),
$$

and

$$
\limsup_{\nu \uparrow \infty} \max_{i,j \in S_{\nu^2}} \left| \frac{\pi_{\nu ij}}{\pi_{\nu i} \pi_{\nu j}} - 1 \right| \le \mathcal{O}\left(N_{\nu}^{-1}\right), \text{with } P_{\theta_0} - \text{probability 1}
$$

*(A7)* (Constant sampling fraction) For some constant, f 2 .0; 1/, that we term the 'sampling fraction',

$$
\limsup_{\nu} \left| \frac{n_{\nu}}{N_{\nu}} - f \right| = \mathcal{O}(1), \text{with } P_0 - \text{probability 1.}
$$

We note that Conditions ([A4\)–\(A7\) a](#page-9-0)re necessary to produce consistency of the sampleweighted pseudo-posterior estimator and, by extension, the MLE of the sample-weighted likelihood,

$$
\mathbb{P}_{N_{\nu}}^{\pi} \ell_{\hat{\theta}_{\pi, N_{\nu}}} \geq \sup_{\theta, \nu} \mathbb{P}_{N_{\nu}}^{\pi} \ell_{\theta} - o_{P_{\theta_{0}}, P_{\nu}} (N_{\nu}^{-1})
$$

and O ;N-P<sup>0</sup> ;P-! 0, which was previously shown by Savitsky and Toth (2[016\) a](#page-25-8)nd Williams & Savitsky (2018a). This paper focuses on achieving correct uncertainty quantification from

the pseudo-posterior by adjusting credibility sets achieved from Bayesian hierarchical model specifications that are guaranteed to asymptotically contract on correct frequentist confidence intervals, which will rely or build upon this consistency result. Bounding the supremum of over all of the inverse of inclusion probabilities from above in Conditi[on \(A5](#page-9-0)) is equivalent to bounding *all* of the inclusion probabilities away from 0. This assumption is used by Savitsky and Tot[h \(201](#page-25-8)6) a[nd Williams & Savitsky \(2018](#page-25-4)a) to achieve their consistency result by ensuring that no portion of the population is systemically excluded from the sample such that the portion may never be sampled, which could otherwise lead to systematic bias. This requirement is used in essentially all consistency results in the literature; see, for example, Pfeffermann *et al.* (1998). We note that Conditi[on \(A6](#page-9-0)) defines two sets containing pairs of population units. The set S<sup>1</sup> contains units, such as those within clusters, where sampling dependence among units is asymptotically *unattenuated*, so long as the *number* of dependent pairs of units is O .N/. The set, S2, by contrast, contains pairs of units (e.g. pairs of units where each unit of the pair resides in a *different* cluster or PSU from the other) where dependence *is* required to asymptotically attenuate to 0. This condition relaxes the usual assumption of asymptotic independence among all units (e.g. requiring that all units are in S2) that has been typically used to guarantee the consistency result for O ;N-P<sup>0</sup> ;P-! 0. While the more restrictive condition is met by nearly all single-stage designs used, in practice (including SRS), it does not apply to multistage clustered sampling designs [which our Conditi[on \(A6](#page-9-0)) *does* cover]. For example, in a two-stage sampling design where the first stage conducts a sampling of clusters of units, the set S<sup>1</sup> captures pairs of units *within* PSU, whose dependence will not asymptotically attenuate to 0, while S<sup>2</sup> contains pairs of units *between* PSUs, where asymptotic independence is required; that is, the PSUs are required to be asymptotically independent. See Savitsky and Tot[h \(201](#page-25-8)6) and Williams & Savitsky (2018a) for extensive details and proofs. Lastly, Conditi[on \(A7](#page-9-0)) is utilised by Savitsky and Tot[h \(201](#page-25-8)6) a[nd Williams & Savitsky \(2018](#page-25-4)a) and relaxes the assumption of an asymptotically 0 sampling fraction used [in Pfeffermann \(1993\). Pfeffermann \(199](#page-25-19)3) assumes an asymptotically 0 sampling fraction in order to approximate the variance of their statistic of interest with respect to the joint distribution, .P; P-<sup>0</sup> /, with just the marginal population-generating distribution, P-<sup>0</sup> ; in other words, they *ignore* the sampling design distribution, P. Savitsky and Toth [\(201](#page-25-8)6) show on page 27 that Conditi[on \(A7](#page-9-0)) allows the replacement of f N with n for sufficiently large to specify the rate of contraction for their consistency result. Their results also go through in the special case that the sampling fraction asymptotes to exactly 0.

<span id="page-11-1"></span>**Theorem 1** (Asymptotic normality of the pseudo-MLE)**.** *Suppose conditio[ns \(A1\)–\(A7](#page-9-0)) hold. Then*

$$
\sqrt{N_{\nu}}\left(\hat{\theta}_{\pi,N_{\nu}}-\theta_{0}\right)=-H_{\theta_{0}}^{-1}\frac{1}{\sqrt{N_{\nu}}}\sum_{i=1}^{N_{\nu}}\frac{\delta_{\nu i}}{\pi_{\nu i}}\dot{\ell}_{\theta_{0}}(\mathbf{X}_{\nu i})+o_{P_{\theta_{0}},P_{\nu}}(1)
$$
(10a)

$$
= -H_{\theta_0}^{-1} \sqrt{N_{\nu}} \mathbb{P}_{N_{\nu}}^{\pi} \dot{\ell}_{\theta_0} + o_{P_{\theta_0}, P_{\nu}}(1)
$$
 (10b)

$$
= -H_{\theta_0}^{-1} \mathbb{G}_{N_v}^{\pi} \dot{\ell}_{\theta_0} + o_{P_{\theta_0}, P_v}(1).
$$
 (10c)

<span id="page-11-0"></span>**Theorem 2** (Asymptotic variance of the pseudo-MLE)**.** *Suppose condition[s \(A1\)–\(A5](#page-9-0)) hold. Then*

$$
Var_{P_{\theta_0}, P_{\nu}} \{-H_{\theta_0}^{-1} \sqrt{N_{\nu}} \mathbb{P}_{N_{\nu}}^{\pi} \dot{\ell}_{\theta_0}\} = H_{\theta_0}^{-1} J_{\theta_0}^{\pi} H_{\theta_0}^{-1}
$$
(11a)

<span id="page-12-1"></span>
$$
=H_{\theta_0}^{-1}\left[J_{\theta_0}+\frac{1}{N_{\nu}}\sum_{i=1}^{N_{\nu}}\mathbb{E}_{P_{\theta_0}}\left\{\left[\frac{1}{\pi_{\nu i}}-1\right]\dot{\ell}_{\theta_0}(\mathbf{X}_{\nu i})\dot{\ell}_{\theta_0}(\mathbf{X}_{\nu i})^T\right\}\right]H_{\theta_0}^{-1}
$$
(11b)

<span id="page-12-0"></span>
$$
\leq \gamma H_{\theta_0}^{-1} J_{\theta_0} H_{\theta_0}^{-1} = \gamma H_{\theta_0}^{-1}.
$$
\n(11c)

The upper bound in Equation (1[1c\) d](#page-12-0)emonstrates a multiplicative injury to <sup>p</sup>N convergence rate achieved for the MLE (under simple random sample of size, N, the population size) in the case of the pseudo-MLE. The larger is , the more varied will be information in the samples around that for the population, which indicates a decreasing efficiency of the sampling design. The amount of injury would be higher for less efficient sampling designs. The maximum penalty paid is a uniformly inflated scale, which will produce wider confidence regions. Theorem 2 [do](#page-11-0)es not restrict the possibility that some designs may be *more* efficient than an SRS or that the efficiency varies by model parameter.

Equation (1[1b\) d](#page-12-1)emonstrates that the shape or geometry of the limiting distribution will be impacted in the case of unequal sampling inclusion probabilities. This 'warping' effect would be expected to be more pronounced in a highly-skewed proportion-to-size sampling design than in an unequally-weighted stratified sampling design with relatively few strata. We demonstrate both of these (warping and scaling) phenomena via simulations in Section 5[.](#page-15-0)

<span id="page-12-2"></span>**Theorem 3** (Asymptotic distribution of the pseudo-posterior)**.** *Suppose conditions ([A1\)–\(A7\)](#page-9-0) hold. Then*

$$
\sup_{B \in \Theta} \left| \Pi_{N_{\nu}}^{\pi} \left( \theta \in B \mid \mathbf{X}_{\nu} \delta_{\nu} \right) - \mathcal{N}_{\hat{\theta}_{\pi, N_{\nu}}, N_{\nu}^{-1} H_{\theta_{0}}^{-1}}(B) \right| \stackrel{P_{\theta_{0}}, P_{\nu}}{\rightarrow} 0, \tag{12}
$$

*where* O ;N*may be the pseudo-MLE or the pseudo-posterior mean.*

The different forms of the asymptotic covariance matrices for the pseudo-MLE (Theorem 2[\),](#page-11-0) on the one hand, and the pseudo-posterior (Theorem 3[\),](#page-12-2) on the other hand, are driven by the failure of Bartlett's second identity under informative sampling J -<sup>0</sup> ¤ J-<sup>0</sup> D H-<sup>0</sup> . This difference motivates our postprocessing step, which we next introduce, that performs multiplicative adjustments to draw from the pseudo-posterior distribution such that their covariance is approximately equal to that of the pseudo-MLE. Please see Section [B of](#page-6-0) the Appendix for detailed proofs of Theorems 1–3. Our theoretical results assume that sources of unit dependence in the population are parameterised. It is the typical or common case under Bayesian modelling to specify sources of dependence in the population through parameterisation of the hierarchical model such that the observed data are conditionally independent given the parameters; an example is the employment of spatially indexed random effects. Our theoretical results (that require consistency of the pseudo-posterior distribution) reflect this case. It may happen, however, that some source of population level dependence is not fully parameterised in the model. For example, the population model response values would be conditionally dependent within sampling clusters. In this case, the theoretical exposition becomes highly complicated and technical. However, this case is common for a survey statistician estimating means and simple models with complex survey data. We will demonstrate in our simulation study (Section 5[\),](#page-15-0) that our proposed adjustment is robust and captures this residual dependence of **X** within sampling units.

#### <span id="page-13-0"></span>**4 Postprocessing the Pseudo-Posterior**

From Sectio[n](#page-6-1) 3, we see that the asymptotic covariance of the pseudo-MLE or pseudoposterior mean is H-1 -<sup>0</sup> <sup>J</sup> -0 H-1 -<sup>0</sup> , yet the asymptotic covariance of our samples drawn from the pseudo-posterior is H-1 -<sup>0</sup> . This is analogous to the differences observed [in Ribatet](#page-25-13) *et al.* (2012), though our formulation for J -<sup>0</sup> (and also <sup>H</sup> -0 ) arises from a random-sampling mechanism, which we leverage in the sequel to perform a post hoc adjustment to draws from the pseudoposterior. Let O <sup>m</sup> represent the sample from the pseudo-posterior for m D 1;:::;M draws with sample mean N . Define the adjusted sample:

<span id="page-13-1"></span>
$$
\hat{\theta}_m^a = \left(\hat{\theta}_m - \bar{\theta}\right) R_2^{-1} R_1 + \bar{\theta},\tag{13}
$$

where R<sup>0</sup> <sup>1</sup>R<sup>1</sup> D H-1 -<sup>0</sup> <sup>J</sup> -0 H-1 -<sup>0</sup> and <sup>R</sup><sup>0</sup> <sup>2</sup>R<sup>2</sup> D H-1 -<sup>0</sup> . We may loosely think of <sup>R</sup>-1 <sup>2</sup> R<sup>1</sup> as a multivariate 'design effect' adjustment (For the SRS sample, we expect Barlett's second identity to hold, and thus, H-1 -<sup>0</sup> J-0H-1 -<sup>0</sup> <sup>D</sup> <sup>H</sup>-1 -<sup>0</sup> , which is the same asymptotic variance as the unadjusted pseudo-posterior). Because O <sup>m</sup> a N.0; N -<sup>1</sup> H-1 -<sup>0</sup> /, we now have <sup>O</sup> a m a N.0; N -<sup>1</sup> H-1 -<sup>0</sup> <sup>J</sup> -0 H-1 -<sup>0</sup> /, which is the asymptotic distribution of the MLE under the pseudolikelihood. Unli[ke Ribatet](#page-25-13) *et al.* (2012), who pre-compute the MLE and change the geometry of their posterior sampler, our implementation is applied as a post hoc projection of the pseudoposterior sample, leaving the initial Monte Carlo sampler intact. So the data analyst may use the Monte Carlo sampler that they designed for population model estimation (under SRS).

For composite likelihoo[ds, Ribatet](#page-25-13) *et al.* (2012) calculate VarP<sup>0</sup> P `-<sup>0</sup> D J-<sup>0</sup> analytically. However, we have an additional distribution P for the sampling design that is unlikely to be in analytic form. In rare cases, the full nn matrix of pairwise inclusion probabilities ij may be available or approximated. In that case, an expression for variance of a weighted sum such as the score function P `-<sup>0</sup> is readily available and could be used for estimating J -<sup>0</sup> conditioned on the data **X**. For example, s[ee Wang](#page-25-5) *et al.* (2017). In practice, the design is often algorithmically defined; for example, designs may use the sorting and clustering of population units in addition to unequal probabilities of selection. Rather than assuming a simplifying model for this distribution, we instead approximate the joint distribution .P-<sup>0</sup> ; P/ with the empirical distribution by resampling the units and associated response values.

Under a multistage sampling design with PSUs constructed as blocks (e.g. geographic regions or HHs) of last stage units (e.g. persons), we would re-sample a subset of the PSUs that contain dependent last-stage units, followed by including all last-stage units with each PSU. We use information about the PSU memberships of each last stage unit in the observed sample in order to conduct the resampling. The data analyst is expected to have this information about the structure of the sampling design, in addition to possessing the sampling weights for the last stage units (e.g. persons). It is necessary when conducting the resampling to explicitly re-sample blocks of units, such as PSUs, when member units express dependence. Such a procedure preserves the dependence structure within the replicate re-samples. This resampling procedure ensures our adjustment properly estimates the scale inflation of the pseudo-posterior distribution induced by the dependent step(s). We use an SRS without replacement procedure to re-sample the PSUs because they are nearly independent from one another, in practice. Equivalently, the survey producer may issue sets of replicates weights that are created internally by using strata and cluster information. The data analyst can then skip directly to the estimation with the provided replicates (Step 13). This variance estimation approach is a hybrid because it uses the Taylor linear expansion to create transformed variables

$$
(\hat{\psi} - \psi_0) = H_{\theta_0}(\hat{\theta} - \theta_0) \approx \sum_{i \in S} w_i \dot{\ell}_{\hat{\theta}}(\mathbf{X}_i) = \sum_{i \in S} w_i z_i(\hat{\theta}).
$$

Variance estimations methods (Taylor linearisation or replication methods) are then applied to P <sup>i</sup>2<sup>S</sup> <sup>w</sup><sup>i</sup> <sup>P</sup> ` O - .**X**i/ where the 'total' is the estimate, and O is a plug-in, calculated only once.

Algorithm 1 provides a simple and computationally efficient resampling approach to estimate VarP<sup>0</sup> ;Ph P N-P `-0 i D J -0 . We recall from Section 3[.2 th](#page-8-1)at H-<sup>0</sup> D EP<sup>0</sup> R `-<sup>0</sup> and H -<sup>0</sup> D EP<sup>0</sup> ;Ph P N-R `-0 i D H-<sup>0</sup> . Therefore, consistent estimates of H-<sup>0</sup> are available without Algorithm 1. Both the plug-in estimate <sup>P</sup> <sup>i</sup>2<sup>S</sup> <sup>w</sup><sup>i</sup> <sup>R</sup> ` N - .**X**i/ and the posterior average <sup>1</sup> M PM mD1 P <sup>i</sup>2<sup>S</sup> <sup>w</sup><sup>i</sup> <sup>R</sup> ` O m.**X**i/ using the original sample S will provide consistent estimates of H-<sup>0</sup> . (We drop the '' subscript from **X** for readability). In our <sup>R</sup> implementation (Appendix [A\),](#page--1-1) we use the plug-in estimate. Estimating HO-<sup>0</sup> within each replication in Algorithm 1 is also possible: HO-<sup>0</sup> D <sup>1</sup> R P<sup>R</sup> <sup>r</sup>D<sup>1</sup> <sup>h</sup><sup>r</sup> with <sup>h</sup><sup>r</sup> <sup>D</sup> <sup>P</sup> <sup>l</sup>2S<sup>r</sup> wQ <sup>r</sup> l R ` N - .**X**<sup>r</sup> <sup>l</sup> /. However, the estimation of <sup>J</sup>O -<sup>0</sup> cannot be performed without estimating across-PSU (or across-replicate) variance. For simplicity, we use half the PSUs from the sample in each replicat[e \(Preston, 200](#page-25-15)9). Other resampling without replacement approaches should be effectiv[e \(Rao](#page-25-3) *et al.*, 1992). However, sampling the PSUs *with* replacement underestimates the variance when the number of PSUs nested within strata is very small because with replacement sampling inaccurately reproduces the sampling design of PSUs from the population. For example, the NSDUH sample only has two PSUs available per strata.

#### <span id="page-15-0"></span>**5 Simulation Study**

We construct a population model to address our inferential interest of a binary outcome y with a linear predictor .

$$
y_i \mid \mu_i \stackrel{\text{ind}}{\sim} \text{Bern}[F_l(\mu_i)], \ i = 1, \dots, N \tag{14}
$$

where F<sup>l</sup> is the cumulative distribution function for the logistic distribution. The first set of simulations (Secti[on 5.1](#page-16-0).1) is based on equal probability sampling. We let depend on a single predictor x1. The second set of simulations (Secti[on 5.1](#page-17-0).2) is based on unequal probability sampling. We let depend on two predictors x<sup>1</sup> and x2, where x<sup>2</sup> is a size variable to set the selection probabilities into the sample. The third set of simulations (Secti[on 5.1](#page-18-0).3) is also based on unequal probability sampling, but we let depend on three predictors x1, x<sup>2</sup> and ´2, where the latter is a random cluster effect at the PSU level. The quantity of inferential interest for all of our simulations is the estimation of the population model coefficients (intercept and slope) for x1. The .x2; ´2/ are nuisance.

The variable x<sup>1</sup> represents the observed information available for analysis, whereas x<sup>2</sup> represents auxiliary information available for setting inclusion probabilities used to conduct sampling, which is either ignored or not available for analysis. The x<sup>1</sup> and x<sup>2</sup> distributions are N .0; 1/ and E.r D 1=5/ with rate r, where N . / and E. / represent normal and exponential distributions, respectively. The cluster effect ´<sup>2</sup> is neither a design variable used for sampling nor part of the analytical model but is a nuisance representing unknown and unmodelled dependence between units within the same cluster (PSU). We choose ´<sup>2</sup> E.1=5/ for a skewed distribution.

We formulate the logarithm of the sampling-weighted pseudo-likelihood for estimating .; / from our observed data for the n -N sampled units,

<span id="page-15-1"></span>
$$
\log \left[ \prod_{i=1}^{n} p(y_i | x_{1i}, \beta_0, \beta_1)^{\tilde{w}_i} \right] = \sum_{i=1}^{n} \tilde{w}_i \log p(y_i | x_{1i}, \beta_0, \beta_1)
$$
  
= 
$$
\sum_{i=1}^{n} \tilde{w}_i y_i \log [F_l(\beta_0 + x_{1i}\beta_1)]
$$
  
+ 
$$
\tilde{w}_i (1 - y_i) \log [1 - F_l(\beta_0 + x_{1i}\beta_1)],
$$
 (15)

where D .ˇ0; ˇ1/, <sup>i</sup> D ˇ<sup>0</sup> C x1iˇ<sup>1</sup> and the sampling weights, wQ <sup>i</sup> , are normalised such that the sum of the weights equals the sample size <sup>P</sup><sup>n</sup> <sup>i</sup>D<sup>1</sup> wQ <sup>i</sup> D n.

Finally, we estimate the joint posterior distribution using Equatio[n \(1](#page-15-1)5), coupled with our prior distribution assignments, using the No U-Turn Sampler Hamiltonian Monte Carlo algorithm implemented in STAN ([Carpenter, 2015; S](#page-24-4)[tan Development Team, 2016\).](#page-25-20) All computations were performed in R (R [Core Team, 2017\).](#page-25-21) See Appendix [A fo](#page--1-1)r example code for fitting and adjusting a STAN model. Generating samples from the different sample designs can be implemented in several ways. We chose to use the 'sampling' package (T[illé & Matei, 2016\)](#page-25-22) for the probability proportional to size (PPS) method of [Brewer \(1975\) a](#page-24-5)nd for systematic sampling.

#### *5.1 Simulation Designs*

In the following subsections, we discuss how we construct sampling design distributions, P, that will induce dependence and skewed information about the population in the observed sample as a means of assessing the performance of our postprocessing adjustment procedure specified in Algorithm 1. In Section 5[.2, w](#page-19-0)e will assess whether the adjustments performed to the posterior draws generate credibility sets that achieve nominal frequentist coverage. We recall from Section 1 [th](#page--1-1)at the survey-sampling literature defines the DEFF as the ratio of the variance of a estimate for the population mean YN under a complex survey design compared with the variance under SRS: DEFFY<sup>N</sup> D VarP- .Y /= bN VarSRS.Y / bN . In addition to nominal coverage, we are also interested in comparing our model-based DEFFs to the standard DEFFY<sup>N</sup> output of designbased survey software, such as the R 'survey' package (L[umley, 2016\) th](#page-25-23)at implements surveyweighted maximum likelihood for point estimation and by default uses Taylor linearisation methods for variance estimation. We estimate the marginal DEFF for each parameter: DEFF- D diagfH-1 - J - H-1 g=diagfH-1 g. These parameter-specific DEFF provide an estimate of the marginal rescaling induced by the complex sample design relative to a simple random sample.

#### <span id="page-16-0"></span>*5.1.1 Equal Probability Dependent Designs (DE)*

For these designs, we induce dependence in the observed samples by clustering units; for example, by aggregating individuals in the population by geographically-indexed domains. This type of clustering or grouping of units is performed by the sampling designers, in practice, in order to control the costs (in this case, travel and labor costs) of administering the survey. It is typically the case that the clustering structure will be coincident with a dependence structure in the population variables of interest; for example, geographically-indexed domains capture a spatial dependence among measures for individuals induced by similarities in culture and economic factors. The effect is that individuals are sampled in dependent groups or clusters, which is expected to lower the amount of information about the population in a realised random sample under this type of sampling design as compared with an SRS of individuals taken from the same population. Even if a sampling design distribution, P, is not informative, the design will induce a scale inflation in the asymptotic covariance of the posterior distribution if the design includes a stage that samples dependent clusters. Our theoretical results do not directly address this possibility and assume that the analyst model is correctly specified to include all population level dependence. Instead, our theoretical results focus on warping and scale adjustments because of approximation errors of the pseudo-posterior induced by unequal weighting. However, we demonstrate in the sequel that the postprocessing adjustment procedure of Algorithm 1, nevertheless, adjusts the scale of the posterior distribution under this scenario of model misspecification to achieve nominal coverage.

The population-generating model is

$$
\mu_i = 0.0 + 1.0x_{1i}
$$

where the intercept was chosen such that the median of is 0; therefore, the median of Fl./ is 0.5.

The first design (DE1), is a one-stage cluster design where clusters of size 5 are selected according to SRS. All individuals have responses that are unconditionally independent. In other words, the clustering membership is randomised and uninformative. Under this scenario, the pseudo-likelihood reduces to the true likelihood with correctly specified independence between units. Therefore, both the unadjusted MCMC samples O <sup>m</sup> and the adjusted MCMC samples O a m should ideally have similar coverage.

The second design (DE5) is also a one-stage SRS design, except that all five members of each cluster have complete dependence. Both the y and the x<sup>1</sup> have identical values within each cluster: yij D yi0<sup>j</sup> and x1ij D x1i0<sup>j</sup> for all individuals i ¤ i<sup>0</sup> in cluster j . Under this scenario, the pseudo-likelihood is again reduced to the simple likelihood. While the likelihood is correctly specified for any given individual, joint cluster dependence is misspecified as independence. Effectively, the sum of the (equal) weights should really sum to n=5 rather than n. Under this scenario, the unadjusted MCMC samples O <sup>m</sup> should have intervals that are too narrow by a factor of <sup>p</sup>5 while the adjusted intervals for <sup>O</sup> a <sup>m</sup> should be longer and achieve the nominal coverage. This idealised example, in which within cluster dependence is both unspecified in the analyst's model and complete, demonstrates the *sensitivity* of the posterior (and pseudo-posterior) to the misspecification of the effective sample size n and the *robustness* of Algorithm 1 to correct for this.

#### <span id="page-17-0"></span>*5.1.2 One stage unequal probability designs (PPS1)*

For these next designs, we have no dependence induced by the clustering of units. Instead, we use an informative design P that uses information from the population to sample units with unequal probabilities of selection; for example, selecting larger businesses with higher probability than smaller businesses in the Current Employment Statistics survey, administrated by the U.S. Bureau of Labor Statistics. In practice, these designs control costs because large businesses contribute proportionately more to estimates for industry totals, such as total production or number of employees. Further refinements to the design, such as stratification of units into size classes, also create statistical efficiencies by reducing the possibility of extreme sample outcomes (such as selecting a sample composed entirely of small businesses). Our theoretical results directly address these informative designs that lead to warping and scale effects because of the approximation error of the pseudo-posterior induced by unequal weighting. We demonstrate that our postprocessing adjustment via Algorithm 1 achieves nominal coverage under these informative sampling designs.

The population-generating model is now

$$
\mu_i = -1.88 + 1.0x_{1i} + 0.5x_{2i}
$$

where the intercept was chosen such that the median of is approximately 0; therefore, the median of Fl./ is approximately 0.5. The size measure used for sample selection is xQ2<sup>i</sup> D x2<sup>i</sup> mini.x2i/ C 1.

Even though the population response y was simulated with D f .x1; x2/, we estimate the marginal models at the population level for D f .x1/. This exclusion of x<sup>2</sup> is analogous to the situation in which an analyst does not have access to all the sample design information and ensures that our sampling design instantiates informativeness (where y is correlated with the selection variable, x2, that defines inclusion probabilities). In particular, we estimate the models under informative design scenarios and compare the population fitted models, D f .x1/, to those from the samples. The first unequally weighted design is a one-stage probability proportional to size design (PPS1), where probabilities of selection are proportional to the size measure <sup>i</sup> / Qx2<sup>i</sup> . For the same population model, we also create a stratified design (SPPS1). We add this additional design because stratification is expected to *improve* the efficiency of the sampling design as compared with SRS because it will—on average—produce samples that are more informationally representative of the population, such that DEFF may be less than 1. We demonstrate that our scale adjustment adapts to efficient refers to the sampling designs not the adjustment more efficient, as well as less efficient, sampling designs. The population is sorted by size measure xQ<sup>2</sup> and then placed into 10 strata. We then select n=10 units from each strata k with ik / Qx2ik.

#### <span id="page-18-0"></span>*5.1.3 Three-stage unequal probability designs (PPS3)*

The last set of designs combines feature of the first two sets. In practice, multistage designs such as the NSDUH first select geographic PSUs (such as states, counties, census tracts, etc.) in proportion to a measure of population size. This provides both cost savings (collecting data in geographic clusters) and statistical efficiencies (higher population areas represent more of the population total), especially when combined with geographic-based stratification (e.g. by state). The final stages for multistage surveys are often the HH and individual. The effect is that individuals within each PSU cluster may likely have outcome measures related to others in their HH and geographic cluster. The within PSU dependence and the unequal probabilities of selection will induce both a scale inflation and a warping in the asymptotic covariance of the posterior distribution. Our theoretical results do not directly address the rescaling because of within cluster dependence; however, our postprocessing adjustment procedure of Algorithm 1, nevertheless, adjusts both the scale and shape of the posterior distribution under this scenario to achieve nominal coverage.

The population-generating model is now

$$
\mu_{ij} = -1.88 + 1.0x_{1ij} + 0.25x_{2ij} + 0.25z_{2j}
$$

where ´2<sup>j</sup> E.1=5/ is the random effect for PSU j . The median of is still close to 0, and the median of Fl./ is still close to 0.5. The size measure used for sample selection is xQ2<sup>i</sup> D x2<sup>i</sup> mini.x2i/ C 1. Compared with the population model for PPS1 and SPPS1, the relationship between y and the size variable x<sup>2</sup> is weaker (0.25 vs. 0.50). This is often the case for HH surveys compared with establishment surveys, because the amount of information available to the sample designer is much greater for establishments than for HHs.

The next design is a three-stage PPS design (PPS3), analogous to a HH survey in which a geographic area is selected as a PSU, followed by a HH and an individual. We employ a simplified, but broadly representative, version of the design used for NSDUH where we first select the PSU based on the size xQ<sup>2</sup> aggregated up to the PSU level. We next select 5 out of 10 HHs within each PSU, where the HHs are sorted based on an aggregate size measure from xQ<sup>2</sup> and sampled systematically (i.e. every other one along the rank sorted list). Finally, one of three individuals are selected within each HH in proportion to the individual size measure xQ2. The nested sampling within PSU, the systematic sampling of HHs and the mutually exclusive sampling of individuals within HHs creates a sampling dependence that does not attenuate (i.e. factor). See [Williams & Savitsky \(2018a\) f](#page-25-4)or a richer discussion of the sources of sampling dependence.

We include a PSU level random effect ´2<sup>j</sup> to allow for the possibility of unmodelled population level dependence that coincides with the sample design induced dependence and together reduce the effective sample size. For example, geographic covariates such as state or census tract may be related to the outcome of interest, but like x2, they are unavailable to the analyst of a public use file because of confidentiality protections. We expect the unadjusted MCMC sample O <sup>m</sup> to under cover both because of the warping effect from unequal weighting and because

<span id="page-19-1"></span>Table 1. *Summary of coverage, average width and design effect estimates for simulations based on 90% posterior intervals. Based on* M D 100 *realisations with sample size* n D 200*,* R D 100 *replications and population sizes* N D 5000; 5000; 6000 *for the SRS (DE1, DE5), one-stage PPS (PPS1, SPPS1) and three-stage (PPS3, SPPS3) designs, respectively, where S denotes the nesting within a stratified sampling stage*

PPS, probability proportional to size; SRS, simple random sampling..

of the overestimation of the effective sample size from the nuisance PSU dependence. We expect the adjusted MCMC sample O a <sup>m</sup> to capture this dependence, leading to wider uncertainty intervals with closer to nominal coverage.

Lastly, we include a stratified version of the design (SPPS3) in which the aggregate size variable for the PSUs is used to sort the clusters into 10 strata, which are then sampled in a three-stage design. Because the size variable x<sup>2</sup> has a weaker relationship with the outcome, the impact of stratification will be weaker for SPPS3 compared with SPPS1. This example is the closest to our motivating NSDUH design and provides insight into the performance of Algorithm 1 when resampling PSUs nested with strata.

#### <span id="page-19-0"></span>*5.2 Results*

Ta[ble](#page-19-1) 1 provides a summary of results for 100 Monte Carlo realisations for each of the six designs based on a target nominal coverage of 90%. A separate population with a specific formulation for the mean of the linear predictor is generated for each of DE1, DE5, (PPS1, SPPS1) and (PPS3, SPPS3) scenarios in each Monte Carlo realisation. A separate sample is taken from the population for each sampling design (of the six total sampling designs). Estimation of points and intervals were conducted for each sample. Total sample sizes of n D 200 were used to explore performance for moderate sample siz[es. Williams & Savitsky \(2018](#page-25-4)a) demonstrate good convergence for this sample size and similar model settings. In other words, the bias for the posterior mean and the MLE is negligible, so we can focus only on the coverage. While discussion of the mean squared error may be interesting, that is a property of the point estimate, not the uncertainty distribution. Because the unadjusted and the adjusted MCMC samples have identical means and are very close to the MLE, the bias and the MSE of the three are essentially the same.

For each Monte Carlo realisation, we create R D 100 replicates to perform an adjustment via Algorithm 1. We consider coverage estimates from 85% to 95% to be reasonably close to the nominal 90% given the simulation noise from the 100 realisations. Marginal coverage is assessed from the two-sided intervals from sample quantiles .q05; q95/. For simplicity, joint coverage is assessed by comparing the Mahalanobis distance .O /<sup>0</sup> Var.O /.O / to the 90% quantile of a <sup>2</sup> <sup>2</sup> distribution. Figu[re](#page-20-0) 1 displays one realisation from each of the design simulations before and after adjustment to visually demonstrate the rescaling and rotation (to undo warping from unequal probability informative sampling) of the adjustment. It also compares the two pseudo-posterior ellipses with those from the pseudo-MLE. This pseudo-MLE ellipse is equivalent to the joint region from numerical samples in Wang *et al.* (2017). Figu[re](#page-21-0) 2 displays the marginal distributions, medians and 90% quantiles for the unad-

![](_page_20_Figure_1.jpeg)

<span id="page-20-0"></span>**Figure 1.** *Joint pseudo-posterior sample for the intercept (horizontal) and slope (vertical) for one realisation of each of six sample designs. Unadjusted (red circles) and adjusted (blue triangles) with approximate 90% density ellipses. Asymptotic normal 90% ellipse for pseudo-maximum likelihood estimation (MLE) (dashed). Created with 'ggplot2'* ([Wickham, 2009\)](#page-25-24)*. PPS, probability proportional to size; SPPS, stratified PPS [Colour figure can be viewed at* w[ileyonlinelibrary.com](wileyonlinelibrary.com)*]*

justed and adjusted pseudo-posterior distribution. The reference lines display the MLE-based median and 90% intervals that are equivalent to the marginal regions from numerical samples in Wang *et al.* [\(2017\).](#page-25-5)

#### *5.2.1 Coverage*

The one-stage equal probability designs (DE1 and DE5) demonstrate that Algorithm 1 is effective across the entire range of within cluster independence to complete cluster dependence. DE1 serves as a control, in which no adjustment should be needed. The marginal coverage and interval widths for the adjusted sample O a <sup>m</sup> are slightly lower than for the unadjusted <sup>O</sup> m, but the joint elliptical coverage is about as good. DE5 serves as an extreme example under which O <sup>m</sup> is clearly under covering and O a <sup>m</sup> performs much better. Figures 1 [an](#page-20-0)d 2 [sh](#page-21-0)ow one realisation in which the densities mostly overlap for DE1. For DE5, we see the adjusted density is much more diffused and indicates some design-induced dependence between the parameters. This may explain why the joint coverage for O <sup>m</sup> is even worse than the marginal coverage and it suggests that a naive rescaling of the weights by five might not lead to correct joint coverage as postulated in Section 5[.1.1. C](#page-16-0)omparisons to the asymptotic MLE distribution reflect good but not perfect alignment between the adjusted distribution and the asymptotic distribution. We

![](_page_21_Figure_1.jpeg)

<span id="page-21-0"></span>**Figure 2.** *Marginal pseudo-posterior sample for the intercept (*-0*) and slope (*-<sup>1</sup>*) for one realisation of each of six sample designs. Unadjusted (left) and adjusted (right) with median and 90% quantiles (solid bars). Asymptotic normal mean (dotted) and 90% interval (dashed) for pseudo-maximum likelihood estimation (MLE). Created with 'ggplot2'* [\(Wickham, 200](#page-25-24)9)*. PPS, probability proportional to size; SPPS, stratified PPS [Colour figure can be viewed at* [wileyonlinelibrary.co](wileyonlinelibrary.com)m*]*

expect this because the adjusted pseudo-posterior, while having improved asymptotic coverage, still maintains its small sample properties. For example, the marginal distributions are not perfectly symmetric or unimodal.

The one-stage unequal probability designs (PPS1 and SPPS1) demonstrate the warping effect without the presence of within cluster dependence. PPS1 demonstrates improvements for both marginal and joint coverage. The stratified version (SPPS1) shows that the unadjusted sample O <sup>m</sup> is over covering, particularly for the joint region. For the moderate sample size n D 200, the adjusted coverage shows a decrease for the intercept but a much closer to nominal coverage for the joint region (88% vs. 99%). Figur[es](#page-20-0) 1 a[nd](#page-21-0) 2 show a similar pattern. The increase in dispersion for the PPS1 design is reduced and offset by stratification in SPPS1. These designs are similar to establishment surveys such as the CES, which may use frame data to form efficient strata and samples for businesses. Again, the adjusted pseudo-posterior compares well to the asymptotic MLE distribution, while still maintaining small sample properties.

The three-stage designs with PSU level dependence (PPS3 and SPPS3) show similar results. The unequal selection is weaker in the three stage designs than in the one-stage. Therefore, the stratification does not lead to much gain in efficiency. Both designs show an improvement in coverage for both marginal and joint coverage and match well (but not identically) to the asymptotic MLE distribution. This is consistent with results from the one-stage designs, but combines the unequal weighting, within PSU dependence and stratification into a single design, similar to HH surveys such as the NSDUH.

#### *5.2.2 Design effects*

We next compare the parameter-specific DEFF to the DEFFY<sup>N</sup> based on Taylor linearisation (L[umley, 2004\).](#page-25-25) Table 1 [sh](#page-19-1)ows that the DEFF for the intercept <sup>0</sup> is very similar to the overall DEFF for y, where the latter is computed from Taylor linearisation methods. This is not surprising because an intercept is very similar to a mean. Examining the DEFF for the slope 1, we see that the effect of the design is typically less dramatic than for the intercept but still notably different from 1. We remind the reader that these estimates for DEFFs assume that the bias has been removed because of incorporation of the weights and do not suggest that equally weighted likelihoods will lead to estimates for slopes that have correct coverage. For comparisons between consistent weighted estimates and biased unweighted estimates, see Savitsky and Toth (2[016\) a](#page-25-8)nd [Williams & Savitsky \(2018a\).](#page-25-4)

### <span id="page-22-1"></span>**6 Application: NSDUH**

A simple logistic model relating current (past month) smoking status to the presence of a past year major depressive episode was fit via the survey-weighted pseudo-posterior as described in Section 5 [us](#page-15-0)ing probability-based analysis weights for adults from the 2014 NSDUH public use data set. Overall, the DEFF for yN is 1.87, and the parameter specific DEFFs are 1.88 for the

![](_page_22_Figure_7.jpeg)

<span id="page-22-0"></span>**Figure 3.** *Joint pseudo-posterior sample for the intercept (horizontal) and slope (vertical) for a logistic regression modelling current cigarette smoking by past year major depressive episode based the 2014 National Survey on Drug Use and Health. Unadjusted (red circles) and adjusted (blue triangles) draws with approximate 90% density ellipses. Asymptotic normal 90% ellipse for pseudo-maximum likelihood estimation (MLE) (dashed). Created with 'ggplot2'* ([Wickham, 2009\)](#page-25-24) *[Colour figure can be viewed at* w[ileyonlinelibrary.com](wileyonlinelibrary.com)*]*

![](_page_23_Figure_1.jpeg)

<span id="page-23-0"></span>**Figure 4.** *Marginal pseudo-posterior sample for the intercept (top) and slope (bottom) for a logistic regression modelling current cigarette smoking by past year major depressive episode based the 2014 National Survey on Drug Use and Health. Unadjusted (left) and adjusted (right) with median and 90% quantiles (solid bars). Asymptotic normal mean (dotted) and 90% interval (dashed) for pseudo-MLE. Created with 'ggplot2'* [\(Wickham, 20](#page-25-24)09) *[Colour figure can be viewed at* [wileyonlinelibrary.co](wileyonlinelibrary.com)m*]*

intercept and 1.12 for the slope. In addition to the marginal rescaling, Figu[re](#page-22-0) 3 demonstrates the presence of a joint warping effect from the complex sample design of the NSDUH.

The rates for both smoking and depression vary by age, urban/rural status, education and other demographics. Some of these factors are related to the sample inclusion probabilities; thus, weighting is needed to mitigate bias. For example s[ee Williams & Savitsky \(2018](#page-25-4)a) for a comparison of the weighted and unweighted estimates. Some of these features (geographic and HH) also correspond to nested clusters of sampling, creating the potential for intracluster correlations. Together, these design features contributed to an almost doubling of the variance associated with the mean and intercept, but a relatively smaller increase in the variance associated with the slope.

Marginal estimates for agree closely when comparing with the survey-weighted MLE (Figur[e](#page-23-0) 4). The covariance structure also matches when comparing the adjusted MCMC samples O a <sup>m</sup> to the pseudo-MLE estimates (Figu[re](#page-22-0) 3). Given the large sample size of approximately 42 000 adults, close agreement is expected. However, we still note potential deviations from asymptotic normality in the pseudo-posterior, which may serve as a tool for model diagnostics that are not available with asymptotic normality methods such as in Wang *et al.* (2017).

#### **7 Discussion**

This work is motivated by the need to apply Bayesian models to survey data. Previous works [\(Savitsky & Toth, 201](#page-25-8)[6; Williams & Savitsky, 2018](#page-25-4)a) have demonstrated consistency of the survey-weighted pseudo-posterior for a large class of population models and complex survey designs. However, L[eón-Novelo & Savitsky \(2019\) o](#page-25-9)bserve that the resulting posterior intervals can have poor frequentist performance. Insights from the composite likelihood (Ribatet *et al.*, 2012) and model misspecification ([Kleijn & van der Vaart, 2012\) l](#page-25-12)iterature motivated the development of the theory and adjustment of the asymptotic covariance of the survey-weighted pseudo-posterior. This resulting adjusted pseudo-posterior can then be used for inference in the same manner as the posterior distribution from a simple random sample. It also achieves the same asymptotic properties as the pseudo-likelihood under 'design-based' frequentist inference methods. While the results match up well with asymptotic normal methods based on the surveyweighted MLE (Wang *et al.*[, 2017\),](#page-25-5) the adjusted pseudo-posterior provides more information with respect to small sample properties. These results allow for modellers to better incorporate informative sample design features into their own analysis models while allowing survey statisticians to incorporate more complex modelling approaches into their analysis and production of official statistics, for example quantile regression and penalised splines (Williams & Savitsky, 2018b) and multivariate latent variable models for count data (Savitsky and Toth, 2[016\).](#page-25-8)

Adjustment 1[3 im](#page-13-1)plemented via Algorithm 1 provides a simple, computationally inexpensive and effective approach to quantifying and adjusting for the warping of the pseudo-posterior because of unequal weighting and complex sampling dependence between sampling units. Our resampling algorithm eliminates the need to analytically integrate VarP<sup>0</sup> P `-<sup>0</sup> and thus can be applied to the composite pseudo-likelihood as a more flexible alternative to the modified MCMC approaches presented in Ribatet *et al.* [\(2012\). I](#page-25-13)ts implementation (Section [A\) i](#page--1-1)s straight-forward by leveraging existing software for Bayesian estimation, algorithmic differentiation and variance estimation via survey replicates. We note that Adjustment 1[3 is](#page-13-1) a projection, but does not force the pseudo-posterior variance to equal that of the pseudo-MLE exactly for small-to-moderate samples. Instead, it provides an asymptotic adjustment that allows the analyst to base inference on the sample distribution of the posterior (adjusted by the DEFF) rather than using the asymptotic MLE covariance. If the latter is desired, benchmarking to force the posterior samples to exactly match specified mean and covariance can be performed in closed form using a constrained linear projection ([Ghosh, 1992; D](#page-25-26)atta *et al.*[, 2011\) o](#page-24-6)r via an iterative Newton–Raphson approach for other constrained projections ([Williams & Berg, 2013\).](#page-25-27) This benchmarked pseudo-posterior would exactly match the mean and covariance of the samples from Wang *et al.* [\(2017\) b](#page-25-5)ut would still preserve some small sample properties by not forcing the sampling distribution to be normal.

#### **References**

- <span id="page-24-1"></span>Binder, D. A. (1996). Linearization methods for single phase and two-phase samples: A cookbook approach. *Surv. Method.*, **22**, 17–22.
- <span id="page-24-3"></span>Breslow, N. E. & Wellner, J. A. (2007). Weighted likelihood for semiparametric models and two-phase stratified samples, with application to cox regression. *Scand. J. Stat.*, **34**(1), 86–102.
- <span id="page-24-5"></span>Brewer, K. (1975). A simple procedure for pswor. *Aust. J. Stat.*, **17**, 166–172.
- <span id="page-24-4"></span>Carpenter, B. (2015). Stan: A probabilistic programming language. *J. Stat. Softw.*
- <span id="page-24-2"></span>Center for Behavioral Health Statistics and Quality. (2015a). Section 1: Adult mental health tables. In *2014 National Survey on Drug Use and Health: Mental Health Detailed Tables*. Rockville, MD: Substance Abuse and Mental Health Services Administration.
- Center for Behavioral Health Statistics and Quality. (2015b). Section 2: Tobacco product and alcohol use tables. In *2014 National Survey on Drug Use and Health: Detailed Tables*. Rockville, MD: Substance Abuse and Mental Health Services Administration.

<span id="page-24-0"></span>Chambers, R. & Skinner, C. (2003). *Analysis of Survey Data*: Wiley Series in Survey Methodology. Wiley.

<span id="page-24-6"></span>Datta, G. S., Ghosh, M., Steorts, R. & Maples, J. (2011). Bayesian benchmarking with applications to small area estimation. *TEST*, **20**(3), 574–588.

- <span id="page-25-26"></span><span id="page-25-17"></span>Ghosal, S., Ghosh, J. K. & Vaart, A. W. V. D. (2000). Convergence rates of posterior distributions. *Ann. Stat*, 500–531. Ghosh, M. (1992). Constrained bayes estimation with applications. *J. Am. Stat. Assoc.*, **87**(418), 533–540.
- <span id="page-25-2"></span>Heeringa, S. G., West, B. T. & Berglund, P. A. (2010). *Applied Survey Data Analysis*: Chapman and Hall/CRC.
- <span id="page-25-0"></span>Isaki, C. T. & Fuller, W. A. (1982). Survey design under the regression superpopulation model. *J. Am. Stat. Assoc.*, **77**, 89–96.
- <span id="page-25-12"></span>Kleijn, B. & van der Vaart, A. (2012). The Bernstein-von-Mises theorem under misspecification. *Electron. J. Stat.*, **6**, 354–381.
- <span id="page-25-9"></span>León-Novelo, L. G. & Savitsky, T. D. (2019). Fully Bayesian estimation under informative sampling. *Electron. J. Stat.*, **13**(1), 1608–1645.
- <span id="page-25-25"></span>Lumley, T. (2004). Analysis of complex survey samples. *J. Stat. Softw.*, **9**(1), 1–19. R package verson 2.2.
- <span id="page-25-23"></span>Lumley, T. (2016). *survey: analysis of complex survey samples*. R package version 3.32.
- <span id="page-25-14"></span>Margossian, C. C. (2018). *A review of automatic differentiation and its efficient implementation*. CoRR abs/1811.05031.
- Morton, K. B., Aldworth, J., Hirsch, E. L., Martin, P. C. & Shook-Sa, B. E. (2016). Section 2, sample design report. In *2014 National Survey on Drug Use and Health: Methodological Resource Book*. Rockville, MD: Center for Behavioral Health Statistics and Quality, Substance Abuse and Mental Health Services Administration.
- <span id="page-25-19"></span>Pfeffermann, D. (1993). The role of sampling weights when modeling survey data. *Int. Stat. Rev. / Rev. Int. de Statistique*, **61**(2), 317–337.
- <span id="page-25-10"></span>Pfeffermann, D., Da Silva Moura, F. A. & Do Nascimento Silva, P. L. (2006). Multi-level modelling under informative sampling. *Biometrika*, **93**(4), 943–959.
- <span id="page-25-1"></span>Pfeffermann, D., Krieger, A. M. & Rinott, Y. (1998). Parametric distributions of complex survey data under informative probability sampling. *Statistica Sinica*, 1087–1114.
- <span id="page-25-15"></span>Preston, J. (2009). Rescaled bootstrap for stratified multistage sampling. *Surv. Methodol.*, **35**(2), 227–234.
- <span id="page-25-21"></span>R Core Team. (2017). *R: A Language and Environment for Statistical Computing* Vienna, Austria: R Foundation for Statistical Computing.
- <span id="page-25-11"></span>Rao, J. N. K. & Wu, C. F. J. (2010). Bayesian pseudo-empirical-likelihood intervals for complex surveys. *J. R. Stat. Soc. Ser. B*, **72**, 533–544.
- <span id="page-25-3"></span>Rao, J. N. K., Wu, C. F. J. & Yue, K. (1992). Some recent work on resampling methods for complex surveys. *Surv. Methodol.*, **18**, 209–217.
- <span id="page-25-13"></span>Ribatet, M., Cooley, D. & Davison, A. C. (2012). Bayesian inference from composite likelihoods, with an application to spatial extremes. *Stat. Sin.*, **22**(2), 813–845.
- <span id="page-25-16"></span>Savitsky, T. D. & Srivastava, S. (2018). Scalable bayes under informative sampling. *Scand. J. Stat.*, **72**, 533–544. [https://doi.org/10.1111/sjos.1231](https://doi.org/10.1111/sjos.12312)2.
- <span id="page-25-8"></span>Savitsky, T. D. & Toth, D. (2016). Bayesian estimation under informative sampling. *Electron. J. Stat.*, **10**(1), 1677– 1708.
- <span id="page-25-7"></span>Shah, B., Bamwell, B., Folsom, R. & Vaish, A. (2000) Design consistent small area estimates using Gibbs algorithm for logistic models.
- <span id="page-25-20"></span>Stan Development Team. (2016). *RStan: the R interface to Stan*. R package version 2.14.1.
- <span id="page-25-22"></span>Tillé, Y. & Matei, A. (2016). *sampling: Survey Sampling*. R package version 2.8.
- <span id="page-25-18"></span>van der Vaart, A. W. (1998): Cambridge University Press.
- <span id="page-25-5"></span>Wang, Z., Kim, J. K. & Yang, S. (2017). Approximate Bayesian inference under informative sampling. *Biometrika*, **105**(1), 91–102.
- <span id="page-25-24"></span>Wickham, H. (2009). *ggplot2: Elegant Graphics for Data Analysis.* New York: Springer-Verlag.
- <span id="page-25-27"></span>Williams, M. & Berg, E. (2013). Incorporating user input into optimal constraining procedures for survey estimates. *J. Off. Stat.*, **29**(3), 375–396.
- <span id="page-25-4"></span>Williams, M. R. & Savitsky, T. D. (2018a). Bayesian estimation under informative sampling with unattenuated dependence. *Bayesian Analysis*. Advance publication.
- Williams, M. R. & Savitsky, T. D. (2018b). Bayesian pairwise estimation under dependent informative sampling. *Electron. J. Stat.*, **12**(1), 1631–1661.
- <span id="page-25-6"></span>Yin, G. (2009). Bayesian generalized method of moments. *Bayesian Analysis*, **4**, 191–207.

#### **A Example Code**

We present a working R [\(R Core Team, 201](#page-25-21)7) implementation and the code for the NSDUH example in Sectio[n](#page-22-1) 6. The function 'cs\_sampling' is a wrapper that takes a STAN model [\(Carpenter, 201](#page-24-4)5), computes Markov chain Monte Carlo (MCMC) draws from the (pseudo) posterior, extracts the gradient function via RSTAN [\(Stan Development Team, 201](#page-25-20)6), creates a replicate design and estimates the variance of the gradient via the 'survey' package (Lumley, 2016). We then compute and apply the sandwich adjustment in Equation (1[3\).](#page-13-1) The resampling method to estimate J in Algorithm 1 corresponds to a special case of the 'mrbbootstrap' replication option (P[reston, 2009\).](#page-25-15) For the weighted logistic STAN model, see appendix B of supplementary information of [Williams & Savitsky \(2018a\): d](#page-25-4)[oi:10.1214/18-BA1143SUPP](doi:10.1214/18-BA1143SUPP)

#### *A1 cs\_sampling*

```
 17515823, 2021, 1, Downloaded from https://onlinelibrary.wiley.com/doi/10.1111/insr.12376 by University of North Carolina at Chapel Hill, Wiley Online Library on [14/10/2025]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License
```
#### **B Proofs**

#### *B1 Proof of Theorem 1 [\(a](#page-11-1)symptotic normality of the pseudo-maximum likelihood estimator)*

The proof strategy closely follows Theorem 5.23 of v[an der Vaart \(1998\) w](#page-25-18)here we update the centered and scaled empirical process, GN to its sampling-weighted extension, G N- . For every random sequence, hN-, we extend v[an der Vaart \(1998\) L](#page-25-18)emma 19.31 to achieve,

<span id="page-30-0"></span>
$$
\mathbb{G}_{N_{\nu}}^{\pi} \left[ \sqrt{N_{\nu}} \left( \ell_{\theta_{0} + \frac{h_{N_{\nu}}}{\sqrt{N_{\nu}}}} - \ell_{\theta_{0}} \right) - h_{N_{\nu}}^{T} \dot{\ell}_{\theta_{0}} \right] \stackrel{P_{\theta_{0}, P_{\nu}}}{\rightarrow} 0. \tag{B1}
$$

Conditions ([A1\) an](#page-9-0)d ([A3\), a](#page-9-0)long with

$$
\mathbb{E}_{P_{\nu}}\left[\mathbb{P}_{N_{\nu}}^{\pi}\ell_{\theta}\right] = \mathbb{E}_{P_{\nu}}\left[\frac{1}{N_{\nu}}\sum_{i=1}^{N_{\nu}}\frac{\delta_{\nu i}}{\pi_{\nu i}}\ell_{\theta}(\mathbf{X}_{i})\right]
$$
(B2)

$$
= \mathbb{P}_{N_v} \ell_{\theta} \tag{B3}
$$

produces a 0 mean for the random sequence of Equation ([B1\) w](#page-30-0)ith respect to the joint distribution, .P- ; P/. By the boundedness requirement for sequence .--1 i / in Condition ([A5\), th](#page-9-0)e Lipschitz condition in Condition ([A1\) an](#page-9-0)d the dominated convergence theorem, their variance converges to 0, and the result in Equation (B[1\) i](#page-30-0)s achieved.

Conditio[ns \(A1\), \(A4](#page-9-0)) a[nd \(A5](#page-9-0)) and Corollary 5.53 [of van der Vaart \(199](#page-25-18)8), the sequence h<sup>N</sup>-<sup>D</sup> <sup>p</sup>N . 0/ is bounded in probability.

We may rewrite Equatio[n \(B](#page-30-0)1) as,

$$
N_{\nu} \mathbb{P}_{N_{\nu}}^{\pi} \log \frac{p_{\theta_{0} + \frac{h_{N_{\nu}}}{\sqrt{N_{\nu}}}}}{p_{\theta_{0}}} - h_{N_{\nu}}^{T} \mathbb{G}_{N_{\nu}}^{\pi} \dot{\ell}_{\theta_{0}} - N_{\nu} \mathbb{E}_{P_{\theta_{0}}} \log \frac{p_{\theta_{0} + \frac{h_{N_{\nu}}}{\sqrt{N_{\nu}}}}}{p_{\theta_{0}}} = o_{P_{\theta_{0}},P_{\nu}}(1)
$$

From Conditi[on \(A2](#page-9-0)), we have,

$$
\mathbb{E}_{P_{\theta_0}} \log \frac{p_{\theta_0 + \frac{h_{N_v}}{\sqrt{N_v}}} }{p_{\theta_0}} - \frac{1}{2N_v} h_{N_v}^T H_{\theta_0} h_{N_v} = o_{P_{\theta_0}}(1)
$$

Substituting this expression above yields,

<span id="page-31-0"></span>
$$
N_{\nu} \mathbb{P}_{N_{\nu}}^{\pi} \log \frac{p_{\theta_0 + \frac{h_{N_{\nu}}}{\sqrt{N_{\nu}}}}}{p_{\theta_0}} = \frac{1}{2} h_{N_{\nu}}^T H_{\theta_0} h_{N_{\nu}} + h_{N_{\nu}}^T \mathbb{G}_{N_{\nu}}^{\pi} \dot{\ell}_{\theta_0} + o_{P_{\theta_0}, P_{\nu}}(1),
$$
(B4)

that we recognise as the local asymptotic normality condition [of Kleijn & van der Vaart \(201](#page-25-12)2) (which we will later use to derive the form for the asymptotic covariance matrix of the pseudoposterior distribution). Equatio[n \(B](#page-31-0)4) is true for both hO N and <sup>h</sup>Q N- D H-1 -<sup>0</sup> <sup>G</sup> N-P `-<sup>0</sup> by Conditi[on \(A4](#page-9-0)). The remainder of the proof exactly follo[ws van der Vaart \(199](#page-25-18)8) where we separately plug in each of hO N and <sup>h</sup>Q N for hN into Equatio[n \(B](#page-31-0)4) to achieve two equivalent equations [up to oP<sup>0</sup> ;P- .1/]. We take the difference between the two equations and complete the square, which produces the result of the theorem.

#### *B2 Proof of Theore[m](#page-11-0) 2 (asymptotic variance of the pseudo-maximum likelihood estimation)*

We begin by constructively expanding the variance with respect to the joint distribution, .P-<sup>0</sup> ; P/,

$$
\text{Var}_{P_{\theta_0}, P_{\nu}} \{-H_{\theta_0}^{-1} \sqrt{N_{\nu}} \mathbb{P}_{N_{\nu}}^{\pi} \dot{\ell}_{\theta_0}\} = N_{\nu} H_{\theta_0}^{-1} \text{Var}_{P_{\theta_0}, P_{\nu}} \mathbb{P}_{N_{\nu}}^{\pi} \dot{\ell}_{\theta_0} H_{\theta_0}^{-1}.
$$
 (B5)

We proceed to apply the total variance decomposition to the variance of the random sequence in the middle of the above expression,

<span id="page-31-2"></span>
$$
\text{Var}_{P_{\theta_0}, P_{\nu}} \mathbb{P}_{N_{\nu}}^{\pi} \dot{\ell}_{\theta_0} = \text{Var}_{P_{\theta_0}} \mathbb{E}_{P_{\nu}} \left[ \mathbb{P}_{N_{\nu}}^{\pi} \dot{\ell}_{\theta_0} \mid \mathcal{A}_{\nu} \right] + \mathbb{E}_{P_{\theta_0}} \text{Var}_{P_{\nu}} \left[ \mathbb{P}_{N_{\nu}}^{\pi} \dot{\ell}_{\theta_0} \mid \mathcal{A}_{\nu} \right], \tag{B6}
$$

where, fixing , A denotes the sigma field of information in the population, U, which encompasses the information about the realised population that is available to the survey designer. Next, we constructively evaluate each of the two terms.

$$
\text{Var}_{P_{\theta_0}} \mathbb{E}_{P_{\nu}} \left[ \mathbb{P}_{N_{\nu}}^{\pi} \dot{\ell}_{\theta_0} \mid \mathcal{A}_{\nu} \right] = \text{Var}_{P_{\theta_0}} \left\{ \frac{1}{N_{\nu}} \sum_{i=1}^{N_{\nu}} \frac{\mathbb{E}_{\nu} \left[ \delta_{\nu i} \mid \mathcal{A}_{\nu} \right]}{\pi_{\nu i}} \dot{\ell}_{\theta_0}(\mathbf{X}_i) \right\}
$$
(B7a)

<span id="page-31-1"></span>
$$
= \text{Var}_{P_{\theta_0}} \left\{ \frac{1}{N_v} \sum_{i=1}^{N_v} \dot{\ell}_{\theta_0}(\mathbf{X}_i) \right\}
$$
 (B7b)

17515823, 2021, 1, Downloaded from https://onlinelibrary.wiley.com/doi/10.1111/insr.12376 by University of North Carolina at Chapel Hill, Wiley Online Library on [14/10/2025]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License

*International Statistical Review* (2021), **89**, 1, 72–107

Published 2020. This article is a U.S. Government work and is in the public domain in the USA.

$$
= \frac{1}{N_v^2} \mathbb{E}_{P_{\theta_0}} \left\{ \sum_{i=1}^{N_v} \dot{\ell}_{\theta_0}(\mathbf{X}_i) \right\}^2 \tag{B7c}
$$

$$
= \frac{1}{N_{\nu}^{2}} \left[ \sum_{i=1}^{N_{\nu}} \mathbb{E}_{P_{\theta_{0}}} \dot{\ell}_{\theta_{0}}(\mathbf{X}_{i}) \dot{\ell}_{\theta_{0}}(\mathbf{X}_{i})^{T} + \sum_{i \neq j \in U_{\nu}} \mathbb{E}_{P_{\theta_{0}}} \dot{\ell}_{\theta_{0}}(\mathbf{X}_{i}) \dot{\ell}_{\theta_{0}}(\mathbf{X}_{j})^{T} \right]
$$
(B7d)

$$
= \frac{1}{N_{\nu}^2} \sum_{i=1}^{N_{\nu}} \mathbb{E}_{P_{\theta_0}} \dot{\ell}_{\theta_0}(\mathbf{X}_i) \dot{\ell}_{\theta_0}(\mathbf{X}_i)^T, \tag{B7e}
$$

<span id="page-32-3"></span>where the second term in the second equation from the bottom results because **X**<sup>i</sup> ? **X**<sup>j</sup> under P-<sup>0</sup> and by Condition ([A3\).](#page-9-0)

$$
\mathbb{E}_{P_{\theta_0}} \text{Var}_{P_{\nu}} \left[ \mathbb{P}_{N_{\nu}}^{\pi} \dot{\ell}_{\theta_0} \mid \mathcal{A}_{\nu} \right] = \mathbb{E}_{P_{\theta_0}} \left\{ \frac{1}{N_{\nu}^2} \text{Var}_{P_{\nu}} \left[ \sum_{i=1}^{N_{\nu}} \frac{\delta_{\nu i}}{\pi_{\nu i}} \dot{\ell}_{\theta_0} \right] \mid \mathcal{A}_{\nu} \right\}
$$
(B8a)

$$
= \frac{1}{N_{\nu}^{2}} \mathbb{E}_{P_{\theta_{0}}} \left\{ \sum_{i=1}^{N_{\nu}} \text{Var}_{P_{\nu}} \left[ \frac{\delta_{\nu i} \mid \mathcal{A}_{\nu}}{\pi_{\nu i}} \right] \dot{\ell}_{\theta_{0}}(\mathbf{X}_{i}) \dot{\ell}_{\theta_{0}}(\mathbf{X}_{i})^{T} + \sum_{i \neq j \in U_{\nu}} \text{Cov}_{P_{\nu}} \left[ \frac{\delta_{\nu i} \delta_{\nu j} \mid \mathcal{A}_{\nu}}{\pi_{\nu i} \pi_{\nu j}} \right] \dot{\ell}_{\theta_{0}}(\mathbf{X}_{i}) \dot{\ell}_{\theta_{0}}(\mathbf{X}_{j})^{T} \right\}
$$
(B8b)

<span id="page-32-0"></span>
$$
= \frac{1}{N_{\nu}^{2}} \sum_{i=1}^{N_{\nu}} \mathbb{E}_{P_{\theta_{0}}} \left\{ \left[ \frac{1}{\pi_{\nu i}} - 1 \right] \dot{\ell}_{\theta_{0}}(\mathbf{X}_{i}) \dot{\ell}_{\theta_{0}}(\mathbf{X}_{i})^{T} \right\} + \frac{1}{N_{\nu}^{2}} \sum_{i \neq j \in U_{\nu}} \mathbb{E}_{P_{\theta_{0}}} \left\{ \left[ \frac{\pi_{\nu ij}}{\pi_{\nu i} \pi_{\nu j}} - 1 \right] \dot{\ell}_{\theta_{0}}(\mathbf{X}_{i}) \dot{\ell}_{\theta_{0}}(\mathbf{X}_{j})^{T} \right\}
$$
(B8c)

<span id="page-32-1"></span>
$$
\leq \frac{1}{N_v^2} \sum_{i=1}^{N_v} \mathbb{E}_{P_{\theta_0}} \left\{ \left[ \frac{1}{\pi_{vi}} - 1 \right] \dot{\ell}_{\theta_0}(\mathbf{X}_i) \dot{\ell}_{\theta_0}(\mathbf{X}_i)^T \right\} + \max\{1, \gamma - 1\} \frac{1}{N_v^2} \sum_{i \neq j \in U_v} \left| \mathbb{E}_{P_{\theta_0}} \left\{ \dot{\ell}_{\theta_0}(\mathbf{X}_i) \dot{\ell}_{\theta_0}(\mathbf{X}_j)^T \right\} \right|
$$
(B8d)

<span id="page-32-2"></span>
$$
= \frac{1}{N_v^2} \sum_{i=1}^{N_v} \mathbb{E}_{P_{\theta_0}} \left\{ \left[ \frac{1}{\pi_{vi}} - 1 \right] \dot{\ell}_{\theta_0}(\mathbf{X}_i) \dot{\ell}_{\theta_0}(\mathbf{X}_i)^T \right\}.
$$
 (B8e)

The sequence, 1 h ij i<sup>j</sup> 1 i , in Equation (B[8c\) is](#page-32-0) bounded from above by <sup>h</sup> <sup>1</sup> <sup>i</sup> 1 i -. 1/ by Condition ([A5\). S](#page-9-0)ee [Williams & Savitsky \(2018a\) f](#page-25-4)or more details. The second expression in Equation ([B8d\) e](#page-32-1)xactly equals 0 by the independence of **X**<sup>i</sup> and **X**<sup>j</sup> (8 i ¤ j 2 U) under P-<sup>0</sup> and by Condition ([A3\). B](#page-9-0)ecause the second expression in Equation ([B8d\) i](#page-32-1)s bounded from above by 0, it exactly equals 0 (for all 2 ZC; i ¤ j 2 U), producing the equality in Equatio[n \(B8](#page-32-2)e). Equatio[n \(B8](#page-32-0)c) results from the following computations:

$$
\operatorname{Var}_{P_{\nu}}\left[\frac{\delta_{\nu i} \mid A_{\nu}}{\pi_{\nu i}}\right] = \mathbb{E}_{P_{\nu}}\left[\frac{\delta_{\nu i} \mid A_{\nu}}{\pi_{\nu i}}\right]^{2} - \left[\mathbb{E}_{P_{\nu}}\frac{\delta_{\nu i} \mid A_{\nu}}{\pi_{\nu i}}\right]^{2}
$$

$$
= \frac{1}{\pi_{\nu i}} - 1
$$

$$
\operatorname{Cov}_{P_{\nu}}\left[\frac{\delta_{\nu i}\delta_{\nu j} \mid A_{\nu}}{\pi_{\nu i}\pi_{\nu j}}\right] = \mathbb{E}_{P_{\nu}}\left[\frac{\delta_{\nu i}\delta_{\nu j} \mid A_{\nu}}{\pi_{\nu i}\pi_{\nu j}}\right] - \mathbb{E}_{P_{\nu}}\left[\frac{\delta_{\nu i} \mid A_{\nu}}{\pi_{\nu i}}\right]\mathbb{E}_{P_{\nu}}\left[\frac{\delta_{\nu j} \mid A_{\nu}}{\pi_{\nu j}}\right]
$$

$$
= \frac{\pi_{\nu ij}}{\pi_{\nu i}\pi_{\nu j}} - 1
$$

We plug in the results for Equation[s \(B](#page-31-1)7) an[d \(B](#page-32-3)8) back into Equatio[n \(B](#page-31-2)6),

$$
N_{\nu} \text{Var}_{P_{\theta_{0}}, P_{\nu}} \mathbb{P}_{N_{\nu}}^{\pi} \dot{\ell}_{\theta_{0}} = \frac{1}{N_{\nu}} \sum_{i=1}^{N_{\nu}} \mathbb{E}_{P_{\theta_{0}}} \dot{\ell}_{\theta_{0}}(\mathbf{X}_{i}) \dot{\ell}_{\theta_{0}}(\mathbf{X}_{i})^{T} + \frac{1}{N_{\nu}} \sum_{i=1}^{N_{\nu}} \mathbb{E}_{P_{\theta_{0}}} \left\{ \left[ \frac{1}{\pi_{\nu i}} - 1 \right] \dot{\ell}_{\theta_{0}}(\mathbf{X}_{i}) \dot{\ell}_{\theta_{0}}(\mathbf{X}_{i})^{T} \right\} \leq \gamma J_{\theta_{0}},
$$
\n(B9b)

and the result is achieved.

#### *B3 Proof of Theore[m](#page-12-2) 3 (asymptotic distribution of the pseudo-posterior)*

The proof strategy is the same [as Kleijn & van der Vaart \(201](#page-25-12)2) where they first prove the assertion on any two sets of random sequences, .hN- ; gN- / 2 K, where K 2 R<sup>d</sup> is an arbitrary compact set. They then extend the result to a sequence of balls, KN- , centered on 0 with increasing radii, MN- " 1. We extend their strategy by updating notation to incorporate the .ıi ; i/, where i D Prfıi D 1g, governed by the sampling design distribution, P, such that our result applies for .P-<sup>0</sup> ; P/, jointly. Recall that we have the local asymptotic normality result,

<span id="page-33-0"></span>
$$
N_{\nu} \mathbb{P}_{N_{\nu}}^{\pi} \log \frac{p_{\theta_0 + \frac{h_{N_{\nu}}}{\sqrt{N_{\nu}}}}}{p_{\theta_0}} = \frac{1}{2} h_{N_{\nu}}^T H_{\theta_0} h_{N_{\nu}} + h_{N_{\nu}}^T \mathbb{G}_{N_{\nu}}^{\pi} \dot{\ell}_{\theta_0} + \mathcal{O}_{P_{\theta_0}, P_{\nu}}(1),
$$
(B10)

from the proof of Theore[m](#page-11-1) 1 under Conditio[ns \(A1\), \(A2\), \(A3](#page-9-0)) an[d \(A5](#page-9-0)). Define the sampleweighted empirical log-likelihood ratio,

$$
s_{N_v}^{\pi}(h) = N_v \mathbb{P}_{N_v}^{\pi} \log \frac{P_{\theta_0} + \frac{h}{\sqrt{N_v}}}{P_{\theta_0}},
$$
\n(B11)

and let N-;-<sup>0</sup> <sup>D</sup> <sup>H</sup>-1 -<sup>0</sup> <sup>G</sup> N-P `-<sup>0</sup> . Plugging into Equatio[n \(B1](#page-33-0)0), we achieve,

$$
s_{N_{\nu}}^{\pi}(h_{N_{\nu}})=h_{N_{\nu}}^T H_{\theta_0} \Delta_{N_{\nu},\theta_0}^{\pi}-\frac{1}{2}h_{N_{\nu}}^T H_{\theta_0} h_{N_{\nu}}+o_{P_{\theta_0},P_{\nu}}(1).
$$

Published 2020. This article is a U.S. Government work and is in the public domain in the USA. *International Statistical Review* (2021), **89**, 1, 72–107

Let <sup>N</sup> denote the normal distribution, N - N-;-0 ; H-1 -0 and define the sequence of random functions,

<span id="page-34-0"></span>
$$
f_{N_{\nu}}^{\pi}(g_{N_{\nu}},h_{N_{\nu}})=\left(1-\frac{\phi_{N_{\nu}}(h_{N_{\nu}})s_{N_{\nu}}^{\pi}(g_{N_{\nu}})\pi_{N_{\nu}}(g_{N_{\nu}})}{\phi_{N_{\nu}}(g_{N_{\nu}})s_{N_{\nu}}^{\pi}(h_{N_{\nu}})\pi_{N_{\nu}}(h_{N_{\nu}})}\right)_{+}.
$$
\n(B12)

Plugging into the logarithm of Equation ([B12\) f](#page-34-0)or s N- . / and <sup>N</sup>- . /, where for any .h<sup>N</sup>- ; g<sup>N</sup>- / 2 K, the prior ratio, -N- .g<sup>N</sup>- /=-N- .h<sup>N</sup>-/ ! 1 as " 1, we achieve:

$$
\log \left( \frac{\phi_{N_{\nu}}(h_{N_{\nu}}) s_{N_{\nu}}^{\pi}(g_{N_{\nu}}) \pi_{N_{\nu}}(g_{N_{\nu}})}{\phi_{N_{\nu}}(g_{N_{\nu}}) s_{N_{\nu}}^{\pi}(h_{N_{\nu}}) \pi_{N_{\nu}}(h_{N_{\nu}})} \right) = \tag{B13a}
$$

$$
= (g_{N_{\nu}} - h_{N_{\nu}})^T H_{\theta_0} \Delta_{N_{\nu}, \theta_0}^{\pi} + \frac{1}{2} h_{N_{\nu}}^T H_{\theta_0} h_{N_{\nu}} - \frac{1}{2} g_{N_{\nu}}^T H_{\theta_0} g_{N_{\nu}} + o_{P_{\theta_0}, P_{\nu}}(1)   
- \frac{1}{2} \left( h_{N_{\nu}} - \Delta_{N_{\nu}, \theta_0}^{\pi} \right)^T H_{\theta_0} \left( h_{N_{\nu}} - \Delta_{N_{\nu}, \theta_0}^{\pi} \right) + \frac{1}{2} \left( g_{N_{\nu}} - \Delta_{N_{\nu}, \theta_0}^{\pi} \right)^T H_{\theta_0} \left( g_{N_{\nu}} - \Delta_{N_{\nu}, \theta_0}^{\pi} \right)
$$
\n(B13b)

$$
=o_{P_{\theta_0},P_v}(1),\tag{B13c}
$$

as " 1. Conclude that

<span id="page-34-1"></span>
$$
\sup_{g,h\in K} f_{N_{\nu}}^{\pi}(g,h) \stackrel{P_{\theta_0},P_{\nu}}{\to} 0,
$$
\n(B14)

as " 1. Define <sup>N</sup> as the event that ˘ N-.K/ > 0. Define

> ˘;K N- .B j **X**; ı/ D ˘ N- .h 2 B j **X**; ı/ =˘ N-.K j **X**; ı/

to the posterior mass truncated to the compact space, K, and similarly for ˚<sup>K</sup> N- . Fix (any) > 0 and define the sequence of events, ˝N- D n supg;h2<sup>K</sup> f N- .g; h/ - o . Construct the inequality,

$$
\mathbb{E}_{P_{\theta_0}, P_{\nu}} \left\| \Pi_{N_{\nu}}^{\pi, K} - \Phi_{N_{\nu}}^{K} \right\| \mathbf{1}_{\mathcal{B}_{N_{\nu}}} \leq \mathbb{E}_{P_{\theta_0}, P_{\nu}} \left\| \Pi_{N_{\nu}}^{\pi, K} - \Phi_{N_{\nu}}^{K} \right\| \mathbf{1}_{\Omega_{N_{\nu}} \cap \mathcal{B}_{N_{\nu}}} + 2 \mathbb{E}_{P_{\theta_0}, P_{\nu}} \left( \mathcal{B}_{N_{\nu}} \backslash \Omega_{N_{\nu}} \right), \tag{B15}
$$

where the total variation normal, k k, is bounded above by two and the second on the right-hand side is o .1/ from Equation ([B14\). B](#page-34-1)ecause kP Qk D 2 R .1 p=q/CdQ, we may expand the first term on the right-hand side,

$$
\frac{1}{2} \mathbb{E}_{P_{\theta_0}, P_{\nu}} \| \Pi_{N_{\nu}}^{\pi, K} - \Phi_{N_{\nu}}^K \mathbf{1}_{\Omega_{N_{\nu}}} \cap \mathcal{Z}_{N_{\nu}}
$$
(B16)

$$
\leq \mathbb{E}_{P_{\theta_0}, P_{\nu}} \int \left(1 - \frac{\phi_{N_{\nu}}(h) s_{N_{\nu}}^{\pi}(g) \pi_{N_{\nu}}(g)}{\phi_{N_{\nu}}(g) s_{N_{\nu}}^{\pi}(h) \pi_{N_{\nu}}(h)}\right)_{+} d \phi_{N_{\nu}}^{K}(g) d \Pi_{N_{\nu}}^{\pi, K}(h) \mathbf{1}_{\Omega_{N_{\nu}} \cap \Xi_{N_{\nu}}} \tag{B17}
$$

$$
\leq \mathbb{E}_{P_{\theta_0}, P_{\nu}} \int \sup_{g, h \in K} f_{N_{\nu}}^{\pi}(g, h) \mathbf{1}_{\Omega_{N_{\nu}} \cap \mathcal{B}_{N_{\nu}}} d\Phi_{N_{\nu}}^K(g) d\Pi_{N_{\nu}}^{\pi, K}(h) \leq \eta.
$$
 (B18)

The proof next follo[ws Kleijn & van der Vaart \(201](#page-25-12)2) to expand the result on a compact K to compact sets, .K<sup>N</sup>- / of balls with radii M<sup>N</sup>- " 1, which provides the result for R<sup>d</sup> in the limit of . From Theore[m](#page-11-1) 1, we have:

$$
\hat{h}_{\pi, N_{\nu}} = \sqrt{N}_{\nu} \left( \hat{\theta}_{\pi, N_{\nu}} - \theta_0 \right) = -\mathbb{A}_{N_{\nu}, \theta_0}^{\pi} + o_{P_{\theta_0}, P_{\nu}}(1),
$$
\n(B19)

and the stated result is achieved with a rescaling and shift because the total variation norm is invariant to rescalings and shifts.