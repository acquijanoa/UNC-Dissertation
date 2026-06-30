![](_page_0_Picture_0.jpeg)

**Journal of the American Statistical Association**

**ISSN: 0162-1459 (Print) 1537-274X (Online) Journal homepage: [www.tandfonline.com/journals/uasa20](https://www.tandfonline.com/journals/uasa20?src=pdf)**

# **Hypotheses Testing from Complex Survey Data Using Bootstrap Weights: A Unified Approach**

**Jae Kwang Kim, J. N. K. Rao & Zhonglei Wang**

**To cite this article:** Jae Kwang Kim, J. N. K. Rao & Zhonglei Wang (2024) Hypotheses Testing from Complex Survey Data Using Bootstrap Weights: A Unified Approach, Journal of the American Statistical Association, 119:546, 1229-1239, DOI: [10.1080/01621459.2023.2183130](https://www.tandfonline.com/action/showCitFormats?doi=10.1080/01621459.2023.2183130)

**To link to this article:** <https://doi.org/10.1080/01621459.2023.2183130>

© 2023 The Author(s). Published with license by Taylor & Francis Group, LLC.

![](_page_0_Picture_9.jpeg)

Published online: 03 Apr 2023.

[Submit your article to this journal](https://www.tandfonline.com/action/authorSubmission?journalCode=uasa20&show=instructions&src=pdf) 

![](_page_0_Figure_13.jpeg)

![](_page_0_Picture_14.jpeg)

View related [articles](https://www.tandfonline.com/doi/mlt/10.1080/01621459.2023.2183130?src=pdf)

View [Crossmark](http://crossmark.crossref.org/dialog/?doi=10.1080/01621459.2023.2183130&domain=pdf&date_stamp=03%20Apr%202023) data

![](_page_0_Picture_17.jpeg)

Citing [articles:](https://www.tandfonline.com/doi/citedby/10.1080/01621459.2023.2183130?src=pdf) 4 View citing articles

# <span id="page-1-4"></span>**Hypotheses Testing from Complex Survey Data Using Bootstrap Weights: A Unified Approach**

Jae Kwang Ki[m](#page-1-0) [a](#page-1-0) , J. N. K. Ra[ob](#page-1-1), and Zhonglei Wan[g](#page-1-0) [c](#page-1-2)

<span id="page-1-0"></span>a Department of Statistics, Iowa State University, Ames, IA; bSchool of Mathematics and Statistics, Carleton University, Ottawa, ON, Canada; <sup>c</sup> Wang Yanan Institute for Studies in Economics (WISE) and School of Economics, Xiamen University, Xiamen, Fujian, P.R. China

#### <span id="page-1-1"></span>**ABSTRACT**

Standard statistical methods without taking proper account of the complexity of a survey design can lead to erroneous inferences when applied to survey data due to unequal selection probabilities, clustering, and other design features. In particular, the Type I error rates of hypotheses tests using standard methods can be much larger than the nominal significance level. Methods incorporating design features in testing hypotheses have been proposed, including Wald tests and quasi-score tests that involve estimated covariance matrices of parameter estimates. In this article, we present a unified approach to hypothesis testing without requiring estimated covariance matrices or design effects, by constructing bootstrap approximations to quasi-likelihood ratio statistics and quasi-score statistics and establishing its asymptotic validity. The proposed method can be easily implemented without a specific software designed for complex survey sampling. We also consider hypothesis testing for categorical data and present a bootstrap procedure for testing simple goodness of fit and independence in a two-way table. In simulation studies, the Type I error rates of the proposed approach are much closer to their nominal significance level compared with the naive likelihood ratio test and quasi-score test. An application to an educational survey under a logistic regression model is also presented. Supplementary materials for this article are available online.

#### <span id="page-1-2"></span>**ARTICLE HISTORY**

Received March 2019 Accepted January 2023

#### **KEYWORDS**

Quasi-likelihood-ratio test; Quasi-score test; Wald test; Wilks' theorem

# <span id="page-1-3"></span>**1. Introduction**

Testing statistical hypotheses is one of the fundamental problems in statistics. In the parametric model approach, hypothesis testing can be implemented using Wald test, likelihood ratio test, or score test. In each case, a test statistic is calculated and then compared with the 100*(*1 − *α)*%-quantile of the reference distribution, which is the limiting distribution of the test statistic under the null hypothesis, where *α* is the nominal significance level. The limiting distribution is often a Chi-squared distribution (Shao [2003\)](#page-11-0).

However, statistical inference with survey data involves additional steps incorporating the sampling design features. Korn and Graubard [\(1999\)](#page-10-0) and Chambers and Skinner [\(2003\)](#page-10-1) provided comprehensive overviews of the methods for analyzing survey data from complex sampling designs. In hypothesis testing with survey data, the limiting distribution of the test statistic is generally not a Chi-squared distribution. Rather, it can be expressed as a weighted sum of several independent random variables from *χ*2*(*1*)*, which is a Chi-square distribution with one degree of freedom, and the weights depend on unknown parameters. To handle such problems, one may consider some corrections to the test statistics to obtain a Chi-square limiting distribution approximately. Such an approach usually involves computing "design effects" associated with test statistics (Lumley and Scott [2014\)](#page-10-2).

In this article, we use a different approach to computing the limiting distribution using a bootstrap approximation. Beaumont and Bocci [\(2009\)](#page-10-3) investigated a general weighted bootstrap method for hypothesis testing under complex sampling designs in the context of Wald tests for linear regression analysis. In this article, we generalize the bootstrap testing idea of Beaumont and Bocci [\(2009\)](#page-10-3) and present a unified bootstrap weight approach to obtain the limiting distribution of test statistics under complex sampling designs, including stratified multistage sampling with clusters selected with replacement.

The proposed method is a bootstrap weight method (Mashreghi, Haziza, and Léger [2016\)](#page-10-4) under complex survey sampling. Although Praestgaard and Wellner [\(1993\)](#page-10-5) and Chatterjee and Bose [\(2005\)](#page-10-6) rigorously investigate bootstrap weight methods under regularity conditions, both may provide erroneous inferences under complex survey sampling. Praestgaard and Wellner [\(1993\)](#page-10-5) established central limit theorems for the weighted bootstrap empirical process, but they assumed that the sample is independently and identically distributed (iid). Although it is common to assume that the finite population is iid, elements in a sample may not be due to complex sampling designs; see Pfeffermann [\(1993\)](#page-10-7) for details. In addition, Praestgaard and Wellner [\(1993\)](#page-10-5) only focused on mean estimation, but we are interested in hypothesis testing for parameters that solve a certain estimating equation. Chatterjee and Bose [\(2005\)](#page-10-6) proposed a bootstrap weight method for estimating functions,

**CONTACT** Zhonglei Wang [wangzl@xmu.edu.cn](mailto:wangzl@xmu.edu.cn) Wang Yanan Institute for Studies in Economics (WISE) and School of Economics, Xiamen University, Xiamen, Fujian 361005, P.R. China.

Supplementary materials for this article are available online. Please go to [www.tandfonline.com/r/JASA](http://www.tandfonline.com/r/JASA).

<sup>© 2023</sup> The Author(s). Published with license by Taylor & Francis Group, LLC.

This is an Open Access article distributed under the terms of the Creative Commons Attribution-NonCommercial-NoDerivatives License (<http://creativecommons.org/licenses/by-nc-nd/4.0/>), which permits non-commercial re-use, distribution, and reproduction in any medium, provided the original work is properly cited, and is not altered, transformed, or built upon in any way. The terms on which this article has been published allow the posting of the Accepted Manuscript in a repository by the author(s) or with their consent.

<span id="page-2-2"></span>and their theoretical results are valid even for correlated random variables. However, the bootstrap weight method of Chatterjee and Bose [\(2005\)](#page-10-6) is not applicable to informative sampling designs, since survey weights are not taken into consideration; see Pfeffermann [\(1993\)](#page-10-7) for details.

There exist bootstrap weight methods under complex survey sampling, but the proposed bootstrap method differs from them in the following aspects. Rao, Wu, and Yue [\(1992\)](#page-10-8) proposed a bootstrap weight method to estimate the variance of a function of population total estimators under stratified random sampling, but did not investigate the theoretical properties of their bootstrap method. Chipperfield and Preston [\(2007\)](#page-10-9) proposed a without replacement scaled bootstrap to achieve the same goal as Rao, Wu, and Yue [\(1992\)](#page-10-8) under stratified random sampling, but their method is only applicable when the parameter of interest is a smooth function of population totals. Moreover, neither method is applicable to other complex survey sampling. Bertail and Combris [\(1997\)](#page-10-10) proposed a bootstrap weight method to make inferences for population means under general sampling designs, and the corresponding theoretical properties were investigated. However, their results are not applicable to parameters obtained by solving an estimating equation. Although Beaumont and Patak [\(2012\)](#page-10-11) developed a general bootstrap weight method for parameter solving an estimating equation in general sampling designs, they focused only on variance estimation without a rigorous discussion about hypothesis testing. Antal and Tillé [\(2011\)](#page-10-12) proposed a one-one bootstrap method for variance estimation under general sampling designs. Similar to the one-one bootstrap method, Antal and Tillé [\(2014\)](#page-10-13) developed a doubled half-bootstrap to estimate the variance of the estimator for population means. However, neither Antal and Tillé [\(2011\)](#page-10-12) nor Antal and Tillé [\(2014\)](#page-10-13) considered parameters solving an estimating equation. Moreover, existing work did not establish limiting distributions for the corresponding bootstrap estimator, especially for the one solving an estimating equation. Thus, up to our knowledge, there is no theoretical guarantee for existing methods when used for hypothesis testing.

Since the use of bootstrap for hypothesis testing is not fully investigated in the literature, our goal is to fill this important gap under general sampling designs. To do this, we first establish a bootstrap central limit theorem under regularity conditions. The sampling design is allowed to be informative in the sense that the sampling mechanism depends on the response variables being sampled (Pfeffermann [1993\)](#page-10-7). Once the bootstrap central limit theorem is established, the proposed bootstrap method can be applied to the quasi-likelihood-ratio test and the quasi-score test. As long as the bootstrap central limit theorem holds for the sampling design, the bootstrap replicates of the test statistics can be used to approximate their distributions under the null hypothesis. Besides, this article also has a practical contribution to survey sampling, since bootstrap weights are commonly provided in agencies, including Statistics Canada.

The article is organized as follows. In [Section 2,](#page-2-0) the basic setup is introduced. The proposed bootstrap method is presented in [Section 3.](#page-4-0) In [Section 4,](#page-5-0) the quasi-likelihood ratio test and the quasi-score test are introduced. Test for categorical survey data is briefly discussed in [Section 5.](#page-6-0) The results of two simulation studies are presented in [Section 5.](#page-6-0) An application to data from an educational survey under a logistic regression model is presented in [Section 6.](#page-7-0) Concluding remarks are made in [Section 7.](#page-8-0)

# <span id="page-2-0"></span>**2. Basic Setup**

Suppose that a finite population F*<sup>N</sup>* = {*(zi*, *xi*, *yi)* : *i* = 1, *...* , *N*} is a random sample of size *N* from a super-population model (Isaki and Fuller [1982;](#page-10-14) Pfeffermann [1993;](#page-10-7) Sverchkov and Pfeffermann [2004;](#page-11-1) Rubin-Bleuer and Schiopu-Kratina [2005\)](#page-10-15), where *yi* is the response of interest, and *x<sup>i</sup>* is the corresponding auxiliary vector of length *d*, and *z<sup>i</sup>* is the vector of design variables such that the first-order inclusion probability for unit *i* can be obtained by *π<sup>i</sup>* = *π(zi*; *φN)*, where *φ<sup>N</sup>* is a design parameter determined by {*z*1, *...* , *zN*}. In other literature, the inclusion probability is expressed as *π<sup>i</sup>* = *π(*F*N*, *i)*; see Sverchkov and Pfeffermann [\(2004\)](#page-11-1) for details.

We are interested in hypothesis testing for a parameter *θ* <sup>0</sup> = arg max*θ*∈ *E*{*l(θ*; *X*, *Y)*}, where ⊂ R*<sup>p</sup>* is the parameter space, *(X*, *Y)* is the random variable associated with elements of F*N*, the expectation is taken with respect to the super-population model, and *l(θ*; *x*, *y)* is a pre-defined objective function. For example, if the finite population is a random sample generated from a super-population model with a density function *f(x*, *y*; *θ* <sup>0</sup>*)* and if *θ* <sup>0</sup> is the parameter of interest, then *l(θ*; *x*, *y)* = log *f(x*, *y*; *θ)*. For regression problems, *θ* <sup>0</sup> can be the parameter in the conditional expectation of the response of interest given the auxiliary vectors.

From the finite population F*N*, a probability sample {*(xi*, *yi)* : *i* ∈ *A*} of size *n* is selected by a probability sampling design with the first-order inclusion probability *πi*, where *A* is the set of sample indexes. The design variables *z<sup>i</sup>* and the design parameter *φ<sup>N</sup>* are not available at the sample level. From the sample, we can compute a pseudo *M*-estimator of *θ* <sup>0</sup> by

$$
\widehat{\boldsymbol{\theta}} = \arg \max_{\boldsymbol{\theta} \in \Theta} l_{w}(\boldsymbol{\theta}),
$$

where *lw(θ)* = *N*−<sup>1</sup> *<sup>i</sup>*∈*<sup>A</sup> wil(θ*; *xi*, *yi)* is the survey-weighted objective function and *wi* = *π*−<sup>1</sup> *<sup>i</sup>* is the sampling weight for *i* ∈ *A*. Note that *lw(θ)* is design-unbiased for *lN(θ)* = *N*−<sup>1</sup> *<sup>N</sup> <sup>i</sup>*=<sup>1</sup> *l(θ*; *xi*, *yi)* in the sense that *E*{*lw(θ)* | F*N*} = *lN(θ)*, where the conditional expectation conditioning on F*<sup>N</sup>* is with respect to the sampling mechanism. To be rigorous, a subscript *N* should also be used to index the sample, the sampling weights, the estimator and other quantities, but it is suppressed in the sequel for simplicity.

*Remark 1.* When discussing a finite population in the preceding paragraph, we implicitly assume that elements of F*<sup>N</sup>* are iid. Such an assumption is widely adopted in survey sampling; see Corollary 1.3.2.1 and Section 2.2.1 of Fuller [\(2009\)](#page-10-16) for details. Besides, under a linear regression model,

<span id="page-2-1"></span>
$$
Y = \mathbf{X}^{\mathrm{T}} \boldsymbol{\beta}_0 + \epsilon \tag{1}
$$

with *E(* | *x)* = 0, the finite population {*(xi*, *yi)* : *i* = 1, *...* , *N*} satisfies the iid assumption; see Särndal, Swensson, and Wretman [\(1989,](#page-11-2) sec. 5.1), Breidt and Opsomer [\(2000\)](#page-10-17), Särndal, Swensson, and Wretman [\(2003,](#page-11-3) sec. 6.4), Kim et al. [\(2006,](#page-10-18) sec. 5), Toth and Eltinge [\(2011\)](#page-11-4), Wu and Thompson [\(2020,](#page-11-5) sec. 5.1) and Han and Wellner [\(2021\)](#page-10-19) for more general <span id="page-3-11"></span>setups about the super-population model. If we are interested in the regression parameter *β*<sup>0</sup> in [\(1\)](#page-2-1), an objective function is *l(θ*; *x*, *y)* = −*(y* − *x*T*β)*2.

Although we implicitly make an iid assumption for the finite population, elements in the sample *A* may not be iid due to complex survey designs. That is, the informative sampling no longer preserves the iid structure in the sample; see Pfeffermann [\(1993\)](#page-10-7) for details. Under certain sampling designs, we can relax the iid assumption, and the super-population model can involve cluster random effects (Rubin-Bleuer and Schiopu-Kratina [2005\)](#page-10-15); refer to Section S4 of supplementary material for an extension to stratified multi-stage sampling with cluster random effects. We focus on element sampling mainly in this article.

Since the parameter of interest is the superpopulation parameter, the final sample can be viewed as a realization of the twophase sampling design, where the finite population is the firstphase sample and the final sample is the second-phase sample. Chen and Rao [\(2007\)](#page-10-20) and Rubin-Bleuer and Schiopu-Kratina [\(2005\)](#page-10-15) present rigorous asymptotic setups for two-phase sampling.

Often, the pseudo *M*-estimator*θ* can be obtained by solving

<span id="page-3-0"></span>
$$
\widehat{\mathbf{S}}_{w}(\boldsymbol{\theta}) = \frac{\partial}{\partial \boldsymbol{\theta}} l_{w}(\boldsymbol{\theta}) = \frac{1}{N} \sum_{i \in A} w_{i} \mathbf{S}(\boldsymbol{\theta}; \mathbf{x}_{i}, y_{i}) = \mathbf{0}, \qquad (2)
$$

where *S(θ*; *x*, *y)* = *∂l(θ*; *x*, *y)/∂θ* is the tangent function of *θ*. To discuss the asymptotic properties of the pseudo M-estimator, we assume the following regularity conditions for a sequence of finite populations and samples (Isaki and Fuller [1982\)](#page-10-14).

- <span id="page-3-1"></span>C1. The parameter space is an open and convex set containing *θ* <sup>0</sup> as an interior point. For *θ* ∈ , *l(θ*; *x*, *y)* is concave and twice-continuously differentiable with respect to *θ*, *E*{|*l(θ*; *X*, *Y)*|} *<* ∞, and *E*{*l(θ*; *X*, *Y)*} is uniquely maximized at *θ* 0.
- <span id="page-3-2"></span>C2. For *θ* ∈ , *V*{*lw(θ)* | F*N*} → 0 in probability with respect to the super-population model as *n* → ∞, where *V*{*lw(θ)* | F*N*} is the design variance of *lw(θ)* conditional on the finite population F*N*. In addition, there exists a nonnegative design variance estimator *V* -{*lw(θ)* | F*N*}, such that *E*[*V* -{*lw(θ)* | F*N*} | F*N*] = *V*{*lw(θ)* | F*N*} almost surely.
- <span id="page-3-3"></span>C3. There exists a compact set B ⊂ containing *θ* <sup>0</sup> as an interior point, such that sup*θ*∈B*nV*{-*Sw(θ)* | F*N*} − *<sup>S</sup>(θ)* → 0 in probability with respect to the super-population model as *n* → ∞, where *V*{-*Sw(θ)* | F*N*} is the design variance of -*Sw(θ)* conditional on the finite population, *M* = sup*x*∈R*m*,*x*=1*Mx* is the operator norm of an *m* × *m* matrix *M* associated with the Euclidean vector norm ·, and  *<sup>S</sup>(θ)* is a positive-definitive nonstochastic matrix.
- <span id="page-3-5"></span>C4. A central limit theorem holds for the survey-weighted score function-*Sw(θ* <sup>0</sup>*)* in [\(2\)](#page-3-0). Specifically,

$$
n^{1/2}\widehat{\mathbf{S}}_{w}(\boldsymbol{\theta}_{0}) \rightarrow N(\mathbf{0}, \Sigma_{S}(\boldsymbol{\theta}_{0}))
$$

in distribution as *n* → ∞ with respect to the superpopulation model and the sampling mechanism.

<span id="page-3-4"></span>C5. The weak convergence of-*Iw(θ)* = *N*−<sup>1</sup> *<sup>i</sup>*∈*<sup>A</sup> wiI(θ*; *xi*, *yi)* is uniform in *<sup>θ</sup>* <sup>∈</sup> <sup>B</sup>, where *<sup>I</sup>(θ*; *<sup>x</sup>*, *<sup>y</sup>)* = −*∂***S***(θ*; *<sup>x</sup>*, *<sup>y</sup>)/∂θ*T. Specifically, sup*θ*∈B-*Iw(θ)*−I*(θ)* → 0 in probability with respect to the super-population model and the sampling mechanism as *n* → ∞, where I*(θ)* = *E*{*I(θ*; *X*, *Y)*}, and I*(θ* <sup>0</sup>*)* is invertible.

Conditions [C1–](#page-3-1)[C2](#page-3-2) are sufficient conditions for the weak consistency of *θ*. Conditions [C3–](#page-3-3)[C5](#page-3-4) are used to establish the asymptotic normality of *θ*. In Conditions [C3–](#page-3-3)[C4,](#page-3-5) we consider the case where *V*[*E*{-*Sw(θ)* | F*N*}] is asymptotically negligible compared with *E*[*V*{-*Sw(θ)* | F*N*}], and it is often the case when *n/N* → 0 as *n* → ∞. That is, we can safely ignore the variability for generating the finite population and focus on the design variance of -*Sw(θ)* given the finite population. This assumption is a building block of the proposed bootstrap method in the next section, and a similar assumption is also proposed by Beaumont and Bocci [\(2009\)](#page-10-3); see Assumption 5 and its discussion of Beaumont and Bocci [\(2009\)](#page-10-3) for details.

As mentioned by Han and Wellner [\(2021\)](#page-10-19), "establishing a finite-dimensional CLT ... can be a nontrivial problem for general sampling designs," and it is beyond our scope to prove Condition [C4](#page-3-5) under general sampling designs. In fact, there is no guarantee that the central limit theorem holds for general sampling designs. Instead, we simply consider the set of sampling designs such that the central limit theorem holds; see Hájek [\(1964\)](#page-10-21), Rosén [\(1972a\)](#page-10-22), Rosén [\(1972b\)](#page-10-23), Berger [\(1998\)](#page-10-24), Rubin-Bleuer and Schiopu-Kratina [\(2005\)](#page-10-15), sections 1.3.3–1.3.4 of Fuller [\(2009\)](#page-10-16), Chauvet [\(2015\)](#page-10-25), and Han and Wellner [\(2021\)](#page-10-19) for details. For simplicity, we use the same notation, "·", to denote the matrix norm and the vector norm. Detailed comments on these regularity conditions are presented in Section S1 of supplementary material. For the proposed conditions, we implicitly assume that *N* → ∞ as *n* → ∞. The following theorem establishes the asymptotic normality of*θ*. Its proof can be found in Section S2 of supplementary material.

<span id="page-3-7"></span>*Theorem 1.* Under Conditions [C1–](#page-3-1)[C5,](#page-3-4) we have

<span id="page-3-9"></span>
$$
\sqrt{n}(\widehat{\boldsymbol{\theta}} - \boldsymbol{\theta}_0) \longrightarrow N(\mathbf{0}, \Sigma_{\boldsymbol{\theta}})
$$
 (3)

in distribution as *n* → ∞ with respect to the joint distribution of the super-population model and the sampling mechanism, where  *<sup>θ</sup>* <sup>=</sup> *<sup>I</sup>(<sup>θ</sup>* <sup>0</sup>*)*−<sup>1</sup> *<sup>S</sup>(<sup>θ</sup>* <sup>0</sup>*)*{*I(<sup>θ</sup>* <sup>0</sup>*)*−1}T.

Now, by the second-order Taylor expansion, we obtain

$$
l_{w}(\theta_{0}) = l_{w}(\widehat{\theta}) + \widehat{S}_{w}^{T}(\widehat{\theta})(\theta_{0} - \widehat{\theta}) - \frac{1}{2}(\widehat{\theta} - \theta_{0})^{T}\widehat{I}_{w}(\widehat{\theta})(\widehat{\theta} - \theta_{0})
$$
  
+  $o_{p}(n^{-1})$   
=  $l_{w}(\widehat{\theta}) - \frac{1}{2}(\widehat{\theta} - \theta_{0})^{T}\widehat{I}_{w}(\widehat{\theta})(\widehat{\theta} - \theta_{0}) + o_{p}(n^{-1}),$  (4)

where-*Iw(θ)* = *N*−<sup>1</sup> *<sup>i</sup>*∈*<sup>A</sup> wiI(θ*; *xi*, *yi)*. Define

<span id="page-3-10"></span><span id="page-3-6"></span>
$$
W(\boldsymbol{\theta}_0) = -2n\{l_w(\boldsymbol{\theta}_0) - l_w(\boldsymbol{\theta})\}.
$$
 (5)

By [\(4\)](#page-3-6), we obtain

$$
W(\boldsymbol{\theta}_0) = n(\widehat{\boldsymbol{\theta}} - \boldsymbol{\theta}_0)^{\mathrm{T}} \widehat{\boldsymbol{I}}_w(\widehat{\boldsymbol{\theta}})(\widehat{\boldsymbol{\theta}} - \boldsymbol{\theta}_0) + o_p(1).
$$

Thus, by [Theorem 1,](#page-3-7) we have

<span id="page-3-8"></span>
$$
W(\boldsymbol{\theta}_0) \longrightarrow \mathcal{G} = \sum_{i=1}^p c_i Z_i^2 \tag{6}
$$

<span id="page-4-6"></span>in distribution as *n* → ∞, and the reference distribution is the joint distribution of the super-population model and the sampling mechanism, where *c*1, *...* ,*cp* are the eigenvalues of *D* =  *<sup>θ</sup>I(θ* <sup>0</sup>*)*, and *Z*1, *...* , *Zp* are *p* independent random variables from the standard normal distribution. Result [\(6\)](#page-3-8) was established by Rao and Scott [\(1984\)](#page-10-26) and Lumley and Scott [\(2014\)](#page-10-2), and it can be regarded as a version of the Wilks' theorem for survey sampling. For *p* = 1, we can use *c* −1 <sup>1</sup> *W(θ* <sup>0</sup>*)* as the test statistic with *χ*2*(*1*)* distribution as the limiting distribution under the null hypothesis, where *c*<sup>1</sup> = *σ*<sup>2</sup> *<sup>S</sup> (θ*0*)/I(θ*0*)*, and *<sup>σ</sup>*<sup>2</sup> *<sup>S</sup> (θ)* and *I(θ)* are the counterparts of  *<sup>s</sup>(θ)* and *I(θ)*, respectively. Unless the sampling design is simple random sampling and the sampling fraction is negligible, the limiting distribution does not reduce to the standard Chi-squared distribution with *p* degrees of freedom.

# <span id="page-4-0"></span>**3. Bootstrap Calibration**

We propose using a bootstrap weight method to approximate the limiting distribution in [\(6\)](#page-3-8) and other test statistics under complex survey sampling. Such a bootstrap calibration is very attractive in practice since there is no need to derive the analytic form of the limiting distribution of the test statistic. That is, standard software that takes account of the weights can be used to implement bootstrap hypothesis testing from data files reporting weights and associated sets of bootstrap weights, as done in Statistics Canada and some other agencies. Given the sample *A*, the proposed bootstrap method is as follows:

- Step 1 Generate a rescaling vector *r*<sup>∗</sup> of size *n* = |*A*| from a bootstrap distribution with *E*∗*(r*∗*)* = **1***n*, and the covariance of *r*<sup>∗</sup> is determined by the sampling design, where |*A*| is the cardinality of set *A*, *E*∗*(*·*)* is the expectation with respect to the bootstrap distribution conditional on the sample *A*, and **1***<sup>n</sup>* is a vector of 1 with length *n*. The bootstrap weights are chosen in such a way that *E*∗{*S*<sup>∗</sup> *<sup>w</sup>(θ)*} = -*Sw(θ)* and *V*∗{*S*<sup>∗</sup> *<sup>w</sup>(θ)*} = *V* -{-*Sw(θ)* | F*N*}, where -*S* ∗ *<sup>w</sup>(θ)* <sup>=</sup> *<sup>i</sup>*∈*<sup>A</sup> <sup>w</sup>*<sup>∗</sup> *<sup>i</sup> <sup>S</sup>(θ*; *<sup>x</sup>i*, *yi)*, *<sup>w</sup>*<sup>∗</sup> *<sup>i</sup>* <sup>=</sup> *<sup>r</sup>*<sup>∗</sup> *<sup>i</sup> wi*, *V*∗*(*·*)* is the bootstrap variance conditional on the sample *A*, and *V* -{-*Sw(θ)* | F*N*} is a nonnegative design variance estimator of -*Sw(θ)* conditional on the finite population F*N*.
- Step 2 Obtain a bootstrap replicate of *θ*, denoted as *θ* ∗ , by solving-*S* ∗ *<sup>w</sup>(θ)* = **0**.

Different distributions may be used to generate *r*<sup>∗</sup> under different sampling designs, but the generation of*r*<sup>∗</sup> only depends on the inclusion probabilities and not on the specific form of the estimating function for parameter estimation. Thus, once the bootstrap weights are generated, they can be used for various inference problems based on the sample *A*. However, the commonly used linearization technique is based on the specific form of the estimating function, so different linearization is required for different hypothesis testing problems. Thus, compared with the linearization technique, the proposed bootstrap method is more appealing to the practitioners.

We now present additional regularity conditions to validate the proposed bootstrap method.

<span id="page-4-1"></span>C6. For *θ* ∈ , the bootstrap weights are nonnegative and satisfy *E*∗{-*S* ∗ *<sup>w</sup>(θ)*} = -*Sw(θ)* and *V*∗{-*S* ∗ *<sup>w</sup>(θ)*} = *V* -{-*Sw(θ)*} conditional on the sample *A*, where *V* -{-*Sw(θ)*} satisfies

$$
n \sup_{\theta \in \mathcal{B}} \|\widehat{V} \{\widehat{S}_w(\theta)\} - V \{\widehat{S}_w(\theta)\}\| \to 0
$$

in probability with respect to the super-population model and the sampling mechanism as *n* → ∞.

<span id="page-4-2"></span>C7. A central limit theorem holds for the bootstrap weighted score function-*S* ∗ *<sup>w</sup>(θ)* uniformly over B. Specifically,

$$
n^{1/2}\{\widehat{\mathbf{S}}_{w}^{*}(\boldsymbol{\theta})-\widehat{\mathbf{S}}_{w}(\boldsymbol{\theta})\} | A \longrightarrow N(\mathbf{0}, \Sigma_{S}(\boldsymbol{\theta})) \qquad (7)
$$

in distribution uniformly over B, and the reference distribution is the bootstrap sampling distribution conditional on the sample *A*.

<span id="page-4-3"></span>C8. sup*θ*∈B-*I* ∗ *<sup>w</sup>(θ)* −-*Iw(θ)* → 0 in probability with respect to the bootstrap sampling distribution conditional on the sample *A* as *n* → ∞, where-*I* ∗ *<sup>w</sup>(θ)* <sup>=</sup> *<sup>N</sup>*−<sup>1</sup> *<sup>i</sup>*∈*<sup>A</sup> <sup>w</sup>*<sup>∗</sup> *<sup>i</sup> I(θ*; *xi*, *yi)*.

Condition [C6](#page-4-1) guarantees that the scaled bootstrap variance of -*S* ∗ *<sup>w</sup>(θ)* converges to  *<sup>S</sup>(θ)* in Condition [C3,](#page-3-3) and they are also used to show the bootstrap central limit theorem in Condition [C7.](#page-4-2) Condition [C8](#page-4-3) is needed to establish the asymptotic normality of the bootstrap estimator*θ* ∗ . To check Condition [C7,](#page-4-2) we may use a bootstrap version of the Lyapounov condition. As discussed in Condition [C4,](#page-3-5) we do not claim that Condition [C7](#page-4-2) holds under general sampling designs, but we can validate both under certain sampling designs. For example, in Section S4 of the supplementary material, we use a Lyapounov-type condition to prove the two central limit theorems under stratified multistage sampling. Detailed comments on these regularity conditions are presented in Section S1 of the supplementary material. The high-level conditions are also verifiable under other commonly used sampling designs, including Poisson sampling.

<span id="page-4-5"></span>*Theorem 2.* Let the assumptions for [Theorem 1](#page-3-7) hold. Then, under Conditions [C6–](#page-4-1)[C8,](#page-4-3) we have

<span id="page-4-4"></span>
$$
\sqrt{n}(\widehat{\boldsymbol{\theta}}^* - \widehat{\boldsymbol{\theta}}) \mid A \longrightarrow N(\mathbf{0}, \Sigma_{\boldsymbol{\theta}})
$$
 (8)

in distribution as *n* → ∞, and the reference distribution on the left side of [\(8\)](#page-4-4) is the bootstrap sampling distribution conditional on the sample *A*, where  *<sup>θ</sup>* is defined in [Theorem 1.](#page-3-7)

[Theorem 2](#page-4-5) is proved under general sampling designs; see Section S3 of the supplementary material for details. [Theorem 2](#page-4-5) establishes the bootstrap central limit theorem for *θ* ∗ , and the limiting distribution coincides with the one in [\(3\)](#page-3-9) for *θ*. Thus, the proposed bootstrap method can be used to approximate the sampling distribution of*θ*.

For estimating the mean of a finite population, Chao and Lo [\(1985\)](#page-10-27) first proved a bootstrap central limit theorem under simple random sampling without replacement, and Bickel and Freedman [\(1984\)](#page-10-28) proved the bootstrap central limit theorem under stratified random sampling without replacement. We now present stratified multi-stage sampling with clusters selected with replacement as an example to validate the bootstrap central limit theorem in [Theorem 2.](#page-4-5) Stratified multi-stage sampling with clusters selected with replacement is extensively used in <span id="page-5-2"></span>socio-economic surveys when the first-stage sampling fractions within strata are negligible, as is often the case with large-scale sample surveys.

*Example 1.* Under stratified multi-stage sampling, a popular approach is to assume that *nh(*≥ 2*)* clusters are selected with replacement from the *h*th stratum, with the selection probabilities *phi* satisfying *Nh <sup>i</sup>*=<sup>1</sup> *phi* = 1 for *h* = 1, *...* , *H* (Rao and Wu [1988\)](#page-10-29), where *Nh* is the number of clusters in the *h*th stratum. As mentioned in Section 2, the selection probabilities can be random with respect to the super-population model. Within each sampled cluster, denote-*Sw*,*hi(θ)* as a design-unbiased estimator of the cluster mean under sampling at the second and subsequent stages. For simplicity, we still use-*Sw*,*hi(θ)*for the case when the *(hi)*th cluster is fully observed. Then, the weighted score function is -*Sw(θ)* = *<sup>H</sup> <sup>h</sup>*=<sup>1</sup> *Wh* -*Sw*,*h(θ)*, where *Wh* = *Mh/M*, *M* = *<sup>H</sup> <sup>h</sup>*=<sup>1</sup> *Mh*, *Mh* <sup>=</sup> *Nh <sup>i</sup>*=<sup>1</sup> *Mhi*, *Mhi* is the size of the *(hi)*th cluster, -*<sup>S</sup>w*,*h(θ)* <sup>=</sup> *<sup>n</sup>*−<sup>1</sup> *h nh i*=1 *Sw*,*hi(θ)*, *Sw*,*hi(θ)* = *da(hi)p*−<sup>1</sup> *a(hi)* -*Sw*,*a(hi)(θ)*, *da(hi)* = *Ma(hi)/Mh*, and *a(hi)* is the index of the *i*th selected cluster in the *h*th stratum.

The bootstrap replicate of -*Sw(θ)* can be obtained in the following way. For the *h*th stratum, generate *n*<sup>∗</sup> *<sup>h</sup>* <sup>=</sup> *(n*<sup>∗</sup> *<sup>h</sup>*1, *...* , *<sup>n</sup>*<sup>∗</sup> *hnh )*T by a multinomial distribution with *nh* − 1 trials and a success probability vector *n*−<sup>1</sup> *<sup>h</sup> (*1, *...* , 1*)*<sup>T</sup> of length *nh*. Then, the bootstrap weighted score function is -*S* ∗ *<sup>w</sup>(θ)* <sup>=</sup> *<sup>H</sup> <sup>h</sup>*=<sup>1</sup> *Wh* -*S* ∗ *<sup>w</sup>*,*h(θ)*, where -*S* ∗ *<sup>w</sup>*,*h(θ)* <sup>=</sup> *<sup>n</sup>*−<sup>1</sup> *h nh <sup>i</sup>*=<sup>1</sup> *<sup>r</sup>*<sup>∗</sup> *hi Sw*,*hi(θ)*, *<sup>r</sup>*<sup>∗</sup> *hi* <sup>=</sup> *khn*<sup>∗</sup> *hi*, and *kh* = *nh/(nh* −1*)*; Rao, Wu, and Yue [\(1992\)](#page-10-8) proposed a similar rescaling bootstrap method under stratified multi-stage sampling, but they did not investigate the corresponding theoretical properties with respect to hypothesis testing of super-population model parameters. By Lemma S7 in Section S4 of supplementary material, we conclude that *E*∗{-*S* ∗ *<sup>w</sup>(θ)*} = -*Sw(θ)* and *V*∗{-*S* ∗ *<sup>w</sup>(θ)*} = *V* -{-*Sw(θ)* | F*N*}, where *V* -{-*<sup>S</sup>w(θ)* <sup>|</sup> <sup>F</sup>*N*} = *<sup>H</sup> <sup>h</sup>*=<sup>1</sup> *W*<sup>2</sup> *<sup>h</sup>*{*nh(nh* − <sup>1</sup>*)*}−1*Sh*,vec*(θ)(Inh* <sup>−</sup> <sup>P</sup>**1**,*nh )S*<sup>T</sup> *<sup>h</sup>*,vec*(θ)*, <sup>P</sup>**1**,*nh* <sup>=</sup> **<sup>1</sup>***nh (***1**<sup>T</sup> *nh* **<sup>1</sup>***nh )*−1**1**<sup>T</sup> *nh* is the projection matrix to the linear space spanned by **1***nh* , **1***nh* = *(*1, *...* , 1*)*<sup>T</sup> is a vector consisting of *nh* 1's,*Inh* is an *nh*×*nh* identity matrix, and

$$
\widetilde{S}_{h, \text{vec}}(\theta) = (d_{a(h1)} p_{a(h1)}^{-1} \widehat{S}_{w,a(h1)}(\theta), \dots, d_{a(hn_h)} p_{a(hn_h)}^{-1} \widehat{S}_{w,a(hn_h)}(\theta)).
$$

In Section S4 of supplementary material, we show that the bootstrap central limit theorem in [Theorem 2](#page-4-5) holds with *<sup>n</sup>* <sup>=</sup> *<sup>H</sup> <sup>h</sup>*=<sup>1</sup> *nh*.

Under stratified simple random sampling without replacement (Bickel and Freedman [1984\)](#page-10-28), a bootstrap method is proposed in Section S7 of the supplementary material, and it is essentially the same as the one considered by Rao and Wu [\(1988\)](#page-10-29) when estimating the population mean; see Rao and Wu [\(1988\)](#page-10-29) and Beaumont and Patak [\(2012\)](#page-10-11) for details.

To establish a bootstrap version of Wilks' Theorem in [\(6\)](#page-3-8), we use

$$
W^*(\widehat{\pmb{\theta}}) = -2n\{l_w^*(\widehat{\pmb{\theta}}) - l_w^*(\widehat{\pmb{\theta}}^*)\}
$$

as the bootstrap version of *W(θ* <sup>0</sup>*)* in [\(5\)](#page-3-10). Under the assumptions for [Theorem 2,](#page-4-5) we can show that

<span id="page-5-1"></span>
$$
W^*(\widehat{\boldsymbol{\theta}}) \mid A \longrightarrow \mathcal{G} \tag{9}
$$

in distribution as *n* → ∞, and the reference distribution is the bootstrap sampling distribution conditional on the sample *A*, where G is the limiting distribution [\(6\)](#page-3-8) of the quasi-likelihood ratio test statistic *W(θ* <sup>0</sup>*)* in [\(5\)](#page-3-10) under the null hypothesis.

By [\(9\)](#page-5-1), we conclude that the corresponding test statistic generated by the proposed bootstrap method has the same asymptotic distribution as in [\(5\)](#page-3-10). Thus, we can use the bootstrap distribution of *<sup>W</sup>*∗*(θ)* to approximate the sampling distribution of *W(θ* <sup>0</sup>*)*.

*Remark 2.* As mentioned in [Section 1,](#page-1-3) most existing bootstrap methods were proposed for variance estimation under survey sampling. There indeed exist some papers discussing bootstrap confidence intervals under survey sampling, but these methods do not apply to general hypothesis testing problems. For example, Rao, Wu, and Yue [\(1992\)](#page-10-8) proposed a bootstrap-*t* confidence interval for the parameter of interest, but their bootstrap confidence interval does not apply when the parameter of interest is multi-dimensional. In addition, their bootstrap method is only valid under the stratified simple random sampling. Beaumont and Patak [\(2012\)](#page-10-11) and Bertail and Combris [\(1997\)](#page-10-10) discussed bootstrap confidence intervals in their simulation studies, but it is not clear how the corresponding intervals are constructed. Also, neither Bertail and Combris [\(1997\)](#page-10-10) nor Beaumont and Patak [\(2012\)](#page-10-11) investigated their bootstrap confidence intervals theoretically. [Theorem 2,](#page-4-5) on the other hand, is a building block to show that the proposed bootstrap method can be used for hypothesis testing under general sampling designs. Although we adopt almost the same technique to prove [Theorems 1–](#page-3-7)[2,](#page-4-5) we are the first to establish that the bootstrap estimator*θ* <sup>∗</sup> has the same asymptotic distribution as the original *θ*. Then, we can show that the corresponding bootstrap statistic*W*∗*(θ)*shares the same asymptotic distribution as the original one *W(θ* <sup>0</sup>*)* under the null hypothesis, so the proposed bootstrap method can be used for hypothesis testing without deriving the corresponding variance estimator as in Lumley and Scott [\(2014\)](#page-10-2). In addition, the proposed bootstrap method can be widely applied in commonly used sampling designs; see Sections S4–S5 of the supplementary material for details.

*Remark 3.* In many practical situations, the design weights *wi* = *π*−<sup>1</sup> *<sup>i</sup>* are further modified to incorporate population-level constraints (as in calibration weighting) or to handle unit nonresponse (as in nonresponse weighting adjustment). In this case, the bootstrap design weights can be modified in the same way to obtain the bootstrap final weights. See Bessonneau et al. [\(2021\)](#page-10-30) for a detailed illustration of constructing bootstrap final weights in household surveys. When the test statistic is computed using the final weights rather than the design weights, the bootstrap final weights can be applied to the same test statistics to obtain the bootstrap distribution of the test statistic.

# <span id="page-5-0"></span>**4. Quasi-Likelihood-Ratio Test**

Without loss of generality, we denote *θ* as the true parameter *θ* <sup>0</sup> in this section. Let <sup>0</sup> ⊂ be the parameter space under the null hypothesis, so the null hypothesis can be written as *H*<sup>0</sup> : *<sup>θ</sup>* <sup>∈</sup> 0. In this section, we consider *<sup>H</sup>*<sup>0</sup> : *<sup>θ</sup>* <sup>2</sup> <sup>=</sup> *<sup>θ</sup>(*0*)* <sup>2</sup> for some <span id="page-6-7"></span>known vector *θ(*0*)* <sup>2</sup> , where *θ* = *(θ*<sup>T</sup> <sup>1</sup> , *θ*<sup>T</sup> <sup>2</sup> *)*T. Thus, we have <sup>0</sup> = {*<sup>θ</sup>* <sup>∈</sup> ; *<sup>θ</sup>* <sup>2</sup> <sup>=</sup> *<sup>θ</sup>(*0*)* <sup>2</sup> }. Let *θ (*0*)* <sup>1</sup> be the profile pseudo *M*-estimator of *<sup>θ</sup>* <sup>1</sup> under *<sup>H</sup>*<sup>0</sup> : *<sup>θ</sup>* <sup>2</sup> <sup>=</sup> *<sup>θ</sup>(*0*)* <sup>2</sup> , which is obtained by maximizing *lw(<sup>θ</sup>* 1, *<sup>θ</sup>(*0*)* <sup>2</sup> *)* with respect to *θ* 1.

The quasi-likelihood ratio test statistic for testing *H*<sup>0</sup> : *θ* <sup>2</sup> = *θ(*0*)* <sup>2</sup> is defined as

<span id="page-6-1"></span>
$$
W(\boldsymbol{\theta}_2^{(0)}) = -2n\{l_w(\widehat{\boldsymbol{\theta}}^{(0)}) - l_w(\widehat{\boldsymbol{\theta}})\},\tag{10}
$$

where *θ (*0*)* <sup>=</sup> *(θ (*0*)* 1 T, *θ(*0*)* 2 <sup>T</sup>*)*T. Under simple random sampling with replacement, *W(θ(*0*)* <sup>2</sup> *)* in [\(10\)](#page-6-1) is asymptotically distributed as *χ*2*(q)* with *q* = *p* − *p*<sup>0</sup> and *p*<sup>0</sup> = dim*(*0*)*, and the reference distribution is the joint distribution of the superpopulation model and the sampling mechanism, where dim*()* is the dimension of . The following lemma, proved by Lumley and Scott [\(2014\)](#page-10-2), presents the limiting distribution of the quasilikelihood ratio test statistic in [\(10\)](#page-6-1) for general sampling designs.

*Lemma 1.* Under *<sup>H</sup>*<sup>0</sup> : *<sup>θ</sup>* <sup>2</sup> <sup>=</sup> *<sup>θ</sup>(*0*)* <sup>2</sup> and the assumptions for [Theorem 2,](#page-4-5)

<span id="page-6-2"></span>
$$
W(\boldsymbol{\theta}_2^{(0)}) \longrightarrow \mathcal{G}_1 = \sum_{i=1}^q c_i Z_i^2 \tag{11}
$$

in distribution as *n* → ∞, and the reference distribution is the joint distribution of the super-population model and the sampling mechanism, where *c*<sup>1</sup> ≥ *c*<sup>2</sup> ≥ ··· ≥ *cq <sup>&</sup>gt;* 0 are the eigenvalues of *<sup>P</sup>* <sup>=</sup>  *<sup>θ</sup>*,2*I*22·1*(<sup>θ</sup>* 1, *<sup>θ</sup>(*0*)* <sup>2</sup> *)*, *Z*1, *...* , *Zq* are *q* independent random variables from the standard normal distribution,  *<sup>θ</sup>*,2 is the asymptotic variance of *<sup>n</sup>*1*/*2*θ* <sup>2</sup> in [\(3\)](#page-3-9), *I*22·1*(θ* 1, *θ* <sup>2</sup>*)* = *I*22*(θ* 1, *θ* <sup>2</sup>*)* − *<sup>I</sup>*21*(<sup>θ</sup>* 1, *<sup>θ</sup>* <sup>2</sup>*)*{*I*11*(<sup>θ</sup>* 1, *<sup>θ</sup>* <sup>2</sup>*)*}−1*I*12*(<sup>θ</sup>* 1, *<sup>θ</sup>* <sup>2</sup>*)*, and *<sup>I</sup>ij(<sup>θ</sup>* 1, *<sup>θ</sup>* <sup>2</sup>*)* <sup>=</sup> *E*[−*∂*2*l*{*(θ* 1, *θ* <sup>2</sup>*)*; *X*, *Y*}*/(∂θi∂θ*<sup>T</sup> *<sup>j</sup> )*] for *i*, *j* = 1, 2.

Lumley and Scott [\(2014\)](#page-10-2) proposed to estimate the limiting distribution in [\(11\)](#page-6-2) using a design-based estimator of *P* =  *<sup>θ</sup>*,2*I*22·1*(<sup>θ</sup>* 1, *<sup>θ</sup>(*0*)* <sup>2</sup> *)*.

We consider an alternative test using a novel application of the bootstrap method discussed in [Section 3.](#page-4-0) To do this, a bootstrap version of the quasi-likelihood ratio test statistic in [\(10\)](#page-6-1) is obtained as

<span id="page-6-3"></span>
$$
W^*(\widehat{\boldsymbol{\theta}}_2) = -2n\{l_w^*(\widehat{\boldsymbol{\theta}}_1^{*(0)}, \widehat{\boldsymbol{\theta}}_2) - l_w^*(\widehat{\boldsymbol{\theta}}^*)\},\tag{12}
$$

where *θ* ∗*(*0*)* <sup>1</sup> = arg max*<sup>θ</sup>* <sup>1</sup> *l* ∗ *<sup>w</sup>(θ* 1,*θ* <sup>2</sup>*)*, and *θ* <sup>2</sup> is the second component of*θ* = *(θ*T <sup>1</sup> ,*θ*T <sup>2</sup> *)*T. Under the conditions of [Theorem 2,](#page-4-5) we can show that the limiting distribution for *<sup>W</sup>*∗*(θ* <sup>2</sup>*)* is the same as that in [\(11\)](#page-6-2). Thus, the bootstrap distribution using *<sup>W</sup>*∗*(θ* <sup>2</sup>*)* in [\(12\)](#page-6-3) can be used to obtain the reference distribution of the test statistic *W(θ(*0*)* <sup>2</sup> *)* in [\(10\)](#page-6-1) under *<sup>H</sup>*<sup>0</sup> : *<sup>θ</sup>* <sup>2</sup> <sup>=</sup> *<sup>θ</sup>(*0*)* 2 .

*Remark 4.* In addition to the quasi-likelihood ratio test, the proposed bootstrap method is also applicable to quasi-score tests (Rao, Scott, and Skinner [1998\)](#page-10-31). To develop a quasi-score test for *<sup>H</sup>*<sup>0</sup> : *<sup>θ</sup>* <sup>2</sup> <sup>=</sup> *<sup>θ</sup>(*0*)* <sup>2</sup> , we first decompose the survey-weighted score function as

<span id="page-6-4"></span>
$$
\widehat{\mathbf{S}}_{w}(\boldsymbol{\theta}) = \left(\begin{matrix} \widehat{\mathbf{S}}_{w1}(\boldsymbol{\theta}) \\ \widehat{\mathbf{S}}_{w2}(\boldsymbol{\theta}) \end{matrix}\right) = \left(\begin{array}{c} N^{-1} \sum_{i \in A} w_i \boldsymbol{u}_1(\boldsymbol{\theta}; \boldsymbol{x}_i, y_i) \\ N^{-1} \sum_{i \in A} w_i \boldsymbol{u}_2(\boldsymbol{\theta}; \boldsymbol{x}_i, y_i) \end{array}\right) \quad (13)
$$

with *uj(θ*; *xi*, *yi)* = *(yi* − *μi)*{*V*0*(μi)*}−<sup>1</sup> *∂μi/∂θ<sup>j</sup>* for *j* = 1, 2, where *μ<sup>i</sup>* = *μ(xi*; *θ)*, *μ(x*; *θ)* = *E(Y* | *x)* is a known function with an unknown parameter *θ*, and *V*0*(μi)* = *V(Y* | *x)* is a known function of *μi*. Note that quasi-score tests clearly specify the model without requiring a correctly specified density function, so they are more important in practice. Let*θ (*0*)* <sup>1</sup> be the solution to-*Sw*<sup>1</sup> *<sup>θ</sup>* 1, *<sup>θ</sup>(*0*)* 2 = **0**, where-*Sw*1*(θ)* is defined in [\(13\)](#page-6-4). The quasi-score test statistic for *<sup>H</sup>*<sup>0</sup> : *<sup>θ</sup>* <sup>2</sup> <sup>=</sup> *<sup>θ</sup>(*0*)* <sup>2</sup> is

<span id="page-6-6"></span>
$$
X_{\text{QS}}^2\left(\boldsymbol{\theta}_2^{(0)}\right) = \widehat{\mathbf{S}}_{w2}^{\text{T}}\left(\widehat{\boldsymbol{\theta}}^{(0)}\right)\left\{\widehat{\boldsymbol{I}}_{w22\cdot 1}\left(\widehat{\boldsymbol{\theta}}^{(0)}\right)\right\}^{-1}\widehat{\mathbf{S}}_{w2}\left(\widehat{\boldsymbol{\theta}}^{(0)}\right), (14)
$$

where

<span id="page-6-5"></span>
$$
\widehat{I}_{w22\cdot 1}(\boldsymbol{\theta}) = \widehat{I}_{w22}(\boldsymbol{\theta}) - \widehat{I}_{w21}(\boldsymbol{\theta}) \{\widehat{I}_{w11}(\boldsymbol{\theta})\}^{-1} \widehat{I}_{w12}(\boldsymbol{\theta}) \qquad (15)
$$

and

$$
\begin{pmatrix}\n\widehat{I}_{w11}(\boldsymbol{\theta}) & \widehat{I}_{w12}(\boldsymbol{\theta}) \\
\widehat{I}_{w21}(\boldsymbol{\theta}) & \widehat{I}_{w22}(\boldsymbol{\theta})\n\end{pmatrix} = - \begin{pmatrix}\n\widehat{S}_{w1}(\boldsymbol{\theta})/\partial \theta_1' & \widehat{S}_{w1}(\boldsymbol{\theta})/\partial \theta_2' \\
\widehat{S}_{w2}(\boldsymbol{\theta})/\partial \theta_1' & \widehat{S}_{w2}(\boldsymbol{\theta})/\partial \theta_2'\n\end{pmatrix}.
$$

The bootstrap replicates of *X*<sup>2</sup> *QS(θ(*0*)* <sup>2</sup> *)* can be constructed by replacing the sampling weights *wi* and *θ (*0*)* , respectively, by the bootstrap weights *w*<sup>∗</sup> *<sup>i</sup>* and the bootstrap estimator that solves -*S* ∗ *<sup>w</sup>(θ)* <sup>=</sup> *<sup>N</sup>*−<sup>1</sup> *<sup>i</sup>*∈*<sup>A</sup> <sup>w</sup>*<sup>∗</sup> *<sup>i</sup>* {*yi* <sup>−</sup> *pi(θ)*}*(*1, *<sup>x</sup>*<sup>T</sup> *<sup>i</sup> )*<sup>T</sup> <sup>=</sup> **<sup>0</sup>** subject to *θ* <sup>2</sup> = *θ* 2. Then, the bootstrap version of the quasi-score test statistic *X*<sup>2</sup> *QS(θ(*0*)* <sup>2</sup> *)* is

$$
X_{\text{QS}}^{2*}(\widehat{\boldsymbol{\theta}}_2) = \widehat{\mathbf{S}}_{w2}^{*T}(\widehat{\boldsymbol{\theta}}^{*(0)}) \left\{ \widehat{\boldsymbol{I}}_{w22\cdot 1}^{*}(\widehat{\boldsymbol{\theta}}^{*(0)}) \right\}^{-1} \widehat{\mathbf{S}}_{w2}^{*}(\widehat{\boldsymbol{\theta}}^{*(0)}), \quad (16)
$$

where *θ* <sup>∗</sup>*(*0*)* <sup>=</sup> *(θ* ∗*(*0*)* 1 T,*θ* 2 <sup>T</sup>*)*T, and -*I* ∗ *w*22·1*(θ* ∗*(*0*) )* is computed similarly to [\(15\)](#page-6-5) using bootstrap weights.

# <span id="page-6-0"></span>**5. Test for Categorical Survey Data**

#### *5.1. Goodness-of-Fit Test*

Suppose that a finite population *U* is partitioned into *K* categories with *U* = *U(*1*)* ∪···∪ *U(K)* . Define *pk* = *Nk/N* to be the population proportion of the *k*th category, where *Nk* is the size of *U(k)* . Thus, we implicitly assume that the finite population is an iid realization of the super-population model following a multinomial distribution with *K* categories.

Suppose that we are interested in the simple goodness-of-fit test *H*<sup>0</sup> : *pk* = *p (*0*) <sup>k</sup>* for *k* = 1, *...* , *K*, where *(p (*0*)* <sup>1</sup> , *...* , *p (*0*) <sup>K</sup> )*<sup>T</sup> is a pre-specified vector satisfying *<sup>K</sup> <sup>k</sup>*=<sup>1</sup> *p (*0*) <sup>k</sup>* = 1. From the sample *A*, we compute *pk* = *N <sup>k</sup>/N* as an estimator of *pk*, where *N* - *<sup>k</sup>* = *<sup>i</sup>*∈*Ak wi* is a design-unbiased estimator of *Nk*, *Ak* <sup>=</sup> *<sup>A</sup>* <sup>∩</sup> *<sup>U</sup>(k)* , and *N* - = *<sup>K</sup> <sup>k</sup>*=<sup>1</sup> *N <sup>k</sup>* = *<sup>i</sup>*∈*<sup>A</sup> wi*. The estimated proportions can be obtained by *p* = arg max*p*∈ *<sup>i</sup>*∈*<sup>A</sup> wil(p*; *xi)*, where *l(p*; *xi)* = *<sup>K</sup> <sup>k</sup>*=<sup>1</sup> *xik* log*(pk)*, with *xik* = 1 if the *i*th element belongs to the *k*th category and 0 otherwise, *p* = *(p*1, *...* , *pK)*T.

Then, the Pearson Chi-squared goodness-of-fit test statistic for *H*<sup>0</sup> is

$$
X^{2}(\boldsymbol{p}^{(0)}) = n \sum_{k=1}^{K} (\widehat{p}_{k} - p_{k}^{(0)})^{2} / p_{k}^{(0)},
$$

<span id="page-7-4"></span>where *p(*0*)* = *(p (*0*)* <sup>1</sup> , *...* , *p (*0*) <sup>K</sup>*−1*)*T. Under regularity conditions, according to Rao and Scott [\(1981\)](#page-10-32), we have

<span id="page-7-1"></span>
$$
X^{2}(\boldsymbol{p}^{(0)}) \longrightarrow \mathcal{G}_{2} = \sum_{k=1}^{K-1} \lambda_{k} Z_{k}^{2} \qquad (17)
$$

in distribution as *n* → ∞ under *H*0, and the reference distribution is the joint distribution of the super-population model and the sampling mechanism, where *λ*<sup>1</sup> ≥ *λ*<sup>2</sup> ≥ ··· ≥ *λK*−<sup>1</sup> *>* 0 are the eigenvalues of the design effect matrix *D* = *P*−<sup>1</sup> <sup>0</sup>  *<sup>p</sup>*, *P*<sup>0</sup> = diag*(p(*0*) )* − *p(*0*) (p(*0*) )*T,  *<sup>p</sup>* is the asymptotic variance of √*np*, and *Z*1, *...* , *ZK*<sup>−</sup><sup>1</sup> are *K*−1 independent random variables from the standard normal distribution. Under simple random sampling with replacement, the limiting distribution in [\(17\)](#page-7-1) is a Chi-squared distribution with *K* − 1 degrees of freedom since *P*<sup>0</sup> =  *<sup>p</sup>*.

We now apply the proposed bootstrap method to approximate the limiting distribution in [\(17\)](#page-7-1), without computing the estimated covariance matrix  *<sup>p</sup>*. To describe the proposed method, let *p*<sup>∗</sup> be the estimator of the population proportion *p* based on the bootstrap weights *w*<sup>∗</sup> *<sup>i</sup>* . The proposed bootstrap statistic of *X*2*(p(*0*) )* is

$$
X^{2*}(\widehat{\boldsymbol{p}}) = n \sum_{i=1}^{K} (\widehat{p}_i^* - \widehat{p}_i)^2 / \widehat{p}_i.
$$

We use*pi* in place of *p (*0*) <sup>i</sup>* in the bootstrap test statistics. Under the conditions of [Theorem 2,](#page-4-5) we can show that the limiting distribution for *<sup>X</sup>*2∗*(p)* is the same as that in [\(17\)](#page-7-1).

#### *5.2. Test of Independence in a Two-Way Table*

The proposed bootstrap can also be used to test independence in a two-way count table. Let *pij* = *Nij/N* be the population proportion for cell *(i*, *j)* with margins *pi*<sup>+</sup> and *p*+*<sup>j</sup>* for *i* = 1, *...* , *R* and *j* = 1, *...* , *C*, where *R* and *C* are the numbers of rows and columns, and {*Nij*; *i* = 1, *...* , *R*, *j* = 1, *...* , *C*} is the set of population counts with margins *Ni*<sup>+</sup> and *N*+*j*. Let *N ij* be a design-unbiased estimator of *Nij* and*pij* = *N ij/N* -. By assuming a multinomial distribution for the super-population model, the Chi-squared statistic for testing independence *H*<sup>0</sup> : *pij* = *pi*+*p*+*<sup>j</sup>* for all *i* and *j* is

<span id="page-7-2"></span>
$$
X_I^2 = n \sum_{i=1}^R \sum_{j=1}^C \frac{(\widehat{p}_{ij} - \widehat{p}_{i+} \widehat{p}_{+j})^2}{\widehat{p}_{i+} \widehat{p}_{+j}}.
$$
 (18)

The limiting distribution of *X*<sup>2</sup> *<sup>I</sup>* is a mixture of *<sup>χ</sup>*2*(*1*)* distributions (Rao and Scott [1981\)](#page-10-32).

To develop bootstrap replicates of the test statistics under *H*0, let*p*∗ *ij* be the bootstrap cell proportion computed using the bootstrap weights, *p*∗ *<sup>i</sup>*<sup>+</sup> <sup>=</sup> *<sup>C</sup> j*=1*p*∗ *ij*, and *p*∗ <sup>+</sup>*<sup>j</sup>* <sup>=</sup> *<sup>R</sup> i*=1*p*∗ *ij*. The proposed bootstrap version of *X*<sup>2</sup> *<sup>I</sup>* is

<span id="page-7-3"></span>
$$
X_I^{2*} = n \sum_{i=1}^R \sum_{j=1}^C \frac{\{(\widehat{p}_{ij}^* - \widehat{p}_{i+}^* \widehat{p}_{+j}^*) - (\widehat{p}_{ij} - \widehat{p}_{i+} \widehat{p}_{+j})\}^2}{(\widehat{p}_{i+} \widehat{p}_{+j})}.
$$
 (19)

Under *H*0, the terms in the numerator of *X*<sup>2</sup> *<sup>I</sup>* are identical to {*(pij* <sup>−</sup> *pi*<sup>+</sup>*p*+*j)* − *(pij* − *pi*+*p*+*j)*}2. That is, the bootstrap test statistic is calculated by simply replacing {*pij*,*pi*+,*p*+*j*} and {*pij*, *pi*+, *<sup>p</sup>*+*j*} with {*p*∗ *ij*,*p*∗ *i*+,*p*∗ +*j* } and {*pij*,*pi*+,*p*+*j*}, respectively. It can be shown that the limiting distribution of *X*2<sup>∗</sup> *<sup>I</sup>* conditional on the sample is the same as that of *X*<sup>2</sup> *I* .

# <span id="page-7-0"></span>**6. Simulation Study**

#### *6.1. Single-Stage Sampling*

In this section, we test the performance of the proposed bootstrap method under the PPS (probability proportional to size) sampling with replacement. A finite population of size *N* is generated as follows. For *i* = 1, *...* , *N*, *xi* ∼ U*(*−5, 5*)* and *yi* = *θ*<sup>1</sup> + *θ*2*xi* + *i*, where U*(a*, *b)* is a uniform distribution on the interval *(a*, *b)*, *(θ*1, *θ*2*)* = *(*1, 1*)*, and *<sup>i</sup>* ∼ *N(*0, 22*)*. To make the sampling design informative, the PPS sampling with replacement is carried out to generate a sample of size *n* with selection probabilities *pi* ∝ *(*6*Ii* + 1*)* and *<sup>N</sup> <sup>i</sup>*=<sup>1</sup> *pi* = 1, where *Ii* = 1 if *xi<sup>i</sup> >* 0 and 0 otherwise; see Section S6 of the supplementary material for details about the probability proportional to size sampling and the bootstrap weights. The sampling design is informative (Pfeffermann [1993\)](#page-10-7), since the selection probabilities are partially determined by the residuals {*<sup>i</sup>* : *i* = 1, *...* , *N*}. A similar model setup was also used by Verret, Rao, and Hidiroglou [\(2015\)](#page-11-6) in the context of small area estimation. We consider two scenarios: *(N*, *n)* = *(*2,000,200*)* and *(N*, *n)* = *(*15,000,500*)*.

We are interested in testing *<sup>H</sup>*<sup>0</sup> : *<sup>θ</sup>*<sup>2</sup> <sup>=</sup> *<sup>θ</sup>(*0*)* <sup>2</sup> with the nominal significance level of *α* = 0.05 and consider three different values for *θ(*0*)* <sup>2</sup> ∈ {1, 1.1, 1.2}. The following testing methods are compared:

- 1. Naive likelihood ratio method with *W(θ(*0*)* <sup>2</sup> *)*in [\(10\)](#page-6-1) and *χ*2*(*1*)* as test statistic and reference distribution, respectively.
- 2. Naive quasi-score method with test statistic [\(14\)](#page-6-6) and *χ*2*(*1*)* as the reference distribution.
- 3. Lumley and Scott [\(2014\)](#page-10-2) method. The test statistic is *W(θ(*0*)* <sup>2</sup> *)/δ*, where*δ* = *nV* -*(θ*2*)* -*Iw*,22·1,-*Iw*,22·<sup>1</sup> is defined in [\(15\)](#page-6-5), and *V* -*(θ*2*)* is obtained from -*I*−<sup>1</sup> *<sup>w</sup> (θ)V* -{-*Sw(θ)*}{-*I*−<sup>1</sup> *<sup>w</sup> (θ)*}T. The reference distribution is *F*1,*k*, where *Fν*1,*ν*<sup>2</sup> is an *F* distribution with parameters *ν*<sup>1</sup> and *ν*2, *k* is the degrees of freedom of the variance estimator based on the sampling design, and *k* is obtained by subtracting the number of parameters from the effective sample size associated with the sampling design.
- 4. Bootstrap quasi-likelihood-ratio method with *W(θ(*0*)* <sup>2</sup> *)* in [\(10\)](#page-6-1) being the test statistic, and the reference distribution is approximated by the empirical distribution of *<sup>W</sup>*∗*(θ*2*)*in [\(12\)](#page-6-3).
- 5. Bootstrap quasi-score method with *X*<sup>2</sup> *QS(θ(*0*)* <sup>2</sup> *)* in [\(14\)](#page-6-6) being the test statistic, and the reference distribution is approximated by the empirical distribution of the bootstrap test statistics,

When calculating the likelihood-ratio test statistics, we assume a normal distribution with a constant but unknown variance for the conditional distribution *yi* given *xi*. For each scenario, we generate 1000 Monte Carlo samples, and we consider *M* = 500 and *M* = 1000 iterations for both bootstrap methods. [Table 1](#page-8-1) summarizes the simulation results. For both scenarios, the naive likelihood-ratio method and the naive quasi-score

<span id="page-8-3"></span><span id="page-8-1"></span>**Table 1.** Test power for the hypothesis test <sup>H</sup><sup>0</sup> : *<sup>θ</sup>*<sup>2</sup> <sup>=</sup> *<sup>θ</sup>(*0*)* <sup>2</sup> based on 1 000 Monte Carlo simulations, the true parameter is *θ*<sup>2</sup> = 1.00, and the nominal significance level is 0.05.

NOTE: NLR, naive likelihood-ratio method; NQS, naive quasi-score method; LS, Lumley and Scott [\(2014\)](#page-10-2) method; BLR, bootstrap quasi-likelihood-ratio method; BQS: bootstrap quasi-score method. The numbers of bootstrap repetitions are set to be M = 500 and M = 1000.

method have a significantly inflated Type I error rate when the null hypothesis is true. The corresponding design effects (Wu and Thompson [2020\)](#page-11-5) are more than 2 for estimating *θ*<sup>2</sup> under both scenarios; see Section S8.1 of the supplementary material for details. Thus, the distributions of *W(θ(*0*)* <sup>2</sup> *)* and the quasiscore test statistic [\(14\)](#page-6-6) are more "dispersed" than *χ*2*(*1*)* under the null hypothesis, and due to the lack of design-effect corrections, the two naive methods are not theoretically valid. The Type I error rates of the Lumley and Scott [\(2014\)](#page-10-2) method and the two proposed bootstrap testing methods are similar under *H*0, regardless of the population and sample sizes. The Type I error rates of the three methods come closer to the nominal truth as the population and sample sizes increase. Furthermore, the powers of the two bootstrap testing methods are reasonable compared to the Lumley-Scott method. The naive methods have a slightly larger power at the expense of elevated Type I error rates. In addition, we get similar test results regardless of the number of bootstrap repetitions.

Besides, Section S9 of the supplementary material reports results of an additional simulation study for a stratified two-stage sampling design.

#### *6.2. Test of Independence*

In this section, we consider test of independence in a 3× 3 table of counts to check the performance of the proposed bootstrap testing method. For a finite population *U* = {*y<sup>i</sup>* : *i* = 1, *...* , *N*}, *y<sup>i</sup>* is generated by a multinomial distribution with one trial and a success probability vector *p*, where *p* = *(p*11, *...* , *pij*, *...* , *p*33*)*T, and *pij* is the success probability for the cell in the *i*th row and *j*th column for *i*, *j* = 1, *...* , 3. That is, *y<sup>i</sup>* is a dummy variable consisting of eight 0's and a 1. For the success probability vector *p*, we consider three cases:

Case I : *p*<sup>11</sup> = 1*/*4, *p*<sup>12</sup> = *p*<sup>13</sup> = *p*<sup>21</sup> = *p*<sup>31</sup> = 1*/*8, *p*<sup>22</sup> = *p*<sup>23</sup> = *p*<sup>32</sup> = *p*<sup>33</sup> = 1*/*16.

Case II : *p*<sup>11</sup> = 1*/*4, *p*<sup>12</sup> = *p*<sup>13</sup> = *(*1.4*)/*8, *p*<sup>21</sup> = *p*<sup>31</sup> = *(*0.6*)/*8, *p*<sup>22</sup> = *p*<sup>33</sup> = 1*/*16, *p*<sup>23</sup> = *(*1.4*)/*16, *p*<sup>32</sup> = *(*0.6*)/*16.

Case III : *p*<sup>11</sup> = *p*2,3 = *p*3,2 = 1*/*6, *p*<sup>12</sup> = *p*<sup>13</sup> = *p*<sup>21</sup> = *p*<sup>31</sup> = *p*<sup>22</sup> = *p*<sup>33</sup> = 1*/*12.

Case I satisfies independence for the two-way table of counts, but Cases II–III do not. The level of nonindependence can be expressed by *γ* = <sup>3</sup> *i*=1 3 *<sup>j</sup>*=<sup>1</sup> *(pij* − *pi*+*p*+*j)*2*/(pi*+*p*+*j)*, and

<span id="page-8-2"></span>**Table 2.** Power of the test procedures for independence based on 1000 Monte Carlo simulation samples, Case I corresponds to the independence scenario, and the nominal significance level is 0.05.

NOTE: NP, naive Pearson method; NLR: naive likelihood-ratio method; RS, Rao-Scott method using second-order correction; BP\_I, bootstrap Pearson method with M = 500; BP\_II, bootstrap Pearson method with M = 1000. The numbers of bootstrap repetitions are set to be M = 500 and M = 1000.

*γ* = 0 corresponds to independence. The values of *γ* are 0, 0.017 and 0.125 for Cases I–III, respectively.

For each *y<sup>i</sup>* , we generate an auxiliary variable *xi* = *β*T*y<sup>i</sup>* , where *β* = *(β*1, *...* , *β*9*)*T, *β<sup>j</sup>* = 0.5 + *ej* for *j* = 1, *...* , 9, *ej* ∼ Ex*(*1*)*. Probability proportional to size sampling with replacement is used to generate a sample of size *n* with selection probability proportional to *xi*. We consider two scenarios for the population and sample sizes: *(N*, *n)* = *(*2000, 75*)* and *(N*, *n)* = *(*10,000, 150*)*.

We are interested in testing independence in the two-way table with *α* = 0.05 nominal significance level. For each sample, we compare the following test methods:

- 1. Naive Pearson method based on *X*<sup>2</sup> *<sup>I</sup>* in [\(18\)](#page-7-2) with *<sup>χ</sup>*2*(*4*)* being the reference distribution.
- 2. The second-order Rao-Scott correction method using the Satterthwaite approximation (Thomas and Rao [1987\)](#page-11-7). This method is widely used in analyzing two-way contingency table and is implemented by the command svychisq in the R package survey (Lumley [2020\)](#page-10-33).
- 3. Bootstrap Pearson method using *X*<sup>2</sup> *<sup>I</sup>* , and its distribution is approximated by that of *X*2<sup>∗</sup> *<sup>I</sup>* in [\(19\)](#page-7-3).

For each scenario, we generate 1000 Monte Carlo samples, and we consider *M* = 500 and *M* = 1000 iterations for both bootstrap methods. [Table 2](#page-8-2) summarizes the simulation results. For Case I when the null hypothesis is true, the Type I error rate of the naive Pearson method is much larger than the nominal significance level 0.05 for different sample sizes, since the asymptotic distribution of the test statistic *X*<sup>2</sup> *<sup>I</sup>* in [\(18\)](#page-7-2) is more "dispersed" than a *χ*2*(*4*)* distribution under the null hypothesis; see Section S8.2 of the supplementary material for details. The performance of the Rao-Scott method with secondorder correction is similar to the bootstrap testing method in terms of Type I error and power. The power of the bootstrap testing method increases with the value of *γ* . Similar to the previous simulation, the power of the naive method is larger at the expense of elevated Type I error rates. In addition, we get similar test results regardless of the number of bootstrap repetitions.

# <span id="page-8-0"></span>**7. Application**

We present an analysis of the 2011 Private Education Expenditure Survey (PEES) in South Korea using the proposed bootstrap <span id="page-9-2"></span>testing methods. The purpose of this survey is to study the relationship between private education expenditure and academic performance of students before entering college.

A stratified two-stage sampling design was used for the 2011 PEES, and the strata consist of 16 first-tier administrative divisions, including most provinces and metropolitan cities in South Korea. For each stratum, PPS sampling with replacement was conducted in the first stage, and the primary sampling unit was the school. The students were randomly selected from the sample school in the second stage. There are about 1000 sample schools and 45,000 students involved in this survey.

For the student *i* in the sample *A*, let *yi* be the academic performance assessed by the teacher, and it takes a value from 1 through 3 corresponding to low, middle and high academic performance, respectively. Associated with *yi*, let *x<sup>i</sup>* be the vector of covariates of interest. As discussed by Kim, Park, and Lee [\(2017\)](#page-10-34), we consider the following covariates: after-school education, hours taking lessons provided by the school after regular classes in a month; father's education, 1 for college or higher and 0 otherwise; gender, 1 for female and 0 for male; household income per month; mother's education, 1 for college or higher and 0 otherwise; private education, hours taking private lessons in a month.

In this section, we study the academic performance of students in middle school and high school separately and are interested in studying the conditional probability of achieving high academic performance. Specifically, consider the following logistic model

<span id="page-9-1"></span>
$$
logit{pr(Y = 1 | x)} = (1, x<sup>T</sup>)\theta,
$$
\n(20)

where logit*(x)* = log*(x)*−log*(*1−*x)*for *x* ∈ *(*0, 1*)*, *Y* = 1 if high academic performance is achieved and 0 otherwise, *x* is a vector of six covariates, and *θ* = *(θ*0, *θ*1, *...* , *θ*6*)*T. We are interested in testing the null hypotheses *H*0,*<sup>i</sup>* : *θ<sup>i</sup>* = 0 for *i* = 1, *...* , 6 with *α* = 0.05 nominal significance level. Since the Wald method is widely used in practice, the naive likelihood ratio method, naive quasi-score method, bootstrap quasi-likelihood ratio method, bootstrap quasi-score method with 1000 iterations and a twosided Wald test are compared. The *p*-values for the two naive methods and the Wald test are obtained using reference distributions for simple random sampling. Specifically, the reference distribution for the two naive methods is *χ*2*(*1*)*, and that of the Wald test is a normal distribution with estimated variance (Fuller [2009,](#page-10-16) sec. 1.2.8) for*θ* using the sandwich formula. The *p*values for the proposed bootstrap testing methods are obtained by bootstrap empirical distributions of the corresponding test statistics.

The results of our analysis are summarized in [Table 3.](#page-9-0) The two bootstrap testing methods and the Wald test perform approximately the same in terms of *p*-values. The *p*-values of the two naive methods are approximately the same, but they differ from those of the bootstrap testing methods, especially for "after school education" at the middle school level and most covariates at the high school level, as the naive methods may not properly reflect the intra-cluster correlation due to cluster sampling.

Based on the two bootstrap testing methods in [Table 3,](#page-9-0) we have the following conclusions under *α* = 0.05 nominal significance level. Controlling other covariates, the probability

<span id="page-9-0"></span>**Table 3.** Estimates (Est) and the <sup>p</sup>-values (Unit: 10<sup>−</sup>2) for testing <sup>H</sup>0,<sup>i</sup> : *<sup>θ</sup>*<sup>i</sup> <sup>=</sup> 0, where i = 1, *...* , 6, in [\(20\)](#page-9-1) for the middle school and high school.

NOTE: Cov, covariate; Est, estimate; NLR, naive likelihood-ratio method; NQS, naive quasi-score method; BLR, bootstrap quasi-likelihood-ratio method; BQS: bootstrap quasi-score method; Wald, two-sided Wald test.

that female students achieve high academic performance is significantly higher than that of male students in middle school, but the gender effect is not significant in high school. Hours spent on private education and after-school education can significantly increase the probability of achieving high academic performance in both middle school and high school. The household income has a significantly positive influence on their child's academic performance. However, the level of education of the father and mother has a significant influence during the middle school period only.

# **8. Concluding Remarks**

Many statistical agencies provide microdata files containing survey weights and several sets of associated replication weights, in particular bootstrap weights. Standard statistical packages often permit the use of survey-weighted test statistics, and we have shown how to approximate their distributions under the null hypothesis by their bootstrap analogues computed from the bootstrap weights supplied in the data file. Using bootstrap weights, we can easily apply hypothesis testing without using specialized software designed for complex survey data, as noted by Beaumont and Bocci [\(2009\)](#page-10-3). The proposed bootstrap method is applied to quasi-likelihood ratio tests and quasi-score tests.

Our theories depend on establishing bootstrap central limit theorems under specified sampling designs, including stratified multi-stage sampling with clusters selected with replacement. Under some sampling designs, the central limiting theorem does not necessarily hold and the proposed method is not applicable. In the case of sampling designs that allow for the central limit theorem, the bootstrap central limit theorems can be established under some additional assumptions. Under the assumption that the bootstrap central limit theorem holds in the sampling design, the proposed hypothesis testing methods are applicable. The hypotheses we have discussed are special cases of the one considered in Section 4.2 of Jiang [\(2019\)](#page-10-35), so the theoretical properties of the proposed bootstrap method can also be investigated using the general method developed in Section 4.2 of Jiang [\(2019\)](#page-10-35). The proposed method can be applied to the case of categorical data by developing bootstrap procedures for testing simple goodness

<span id="page-10-36"></span>of fit and independence in a two-way table. We plan to extend our bootstrap method for categorical data to test hypotheses from multi-way tables of weighted counts or proportions, using a log-linear model approach proposed by Rao and Scott [\(1984\)](#page-10-26). Moreover, the bootstrap method can also be applied to nonresponse and weighting adjustments for variance estimation under regularity conditions, and extending the proposed bootstrap method for weight adjustments for hypothesis testing is also a promising project.

# **Supplementary Materials**

Supplementary Material contains the following: comments on the regularity conditions (S1), proof of Theorem 1 (S2), proof of Theorem 2 (S3), verification of the proposed bootstrap method under stratified multistage sampling (S4), bootstrap method for Poisson sampling (S5) bootstrap method for probability proportional to size sampling with replacement (S6), stratified simple random sampling (S7), design effect (S8) and an additional simulation study: stratified two-stage sampling (S9).

# **Acknowledgments**

We are grateful to the associate editor and three anonymous referees for their constructive comments, which have greatly improved the article.

# **Funding**

Kim is partially supported by NSF (OAC-1931380). Rao is supported by a grant from the Natural Sciences and Engineering Research Council of Canada. Wang is partially supported by NSFC (11901487, 71988101, 72033002, 12231011) and National Key R&D Program of China (2022YFA10038002).

# **References**

- <span id="page-10-12"></span>Antal, E., and Tillé, Y. (2011), "A Direct Bootstrap Method for Complex Sampling Designs from a Finite Population," *Journal of the American Statistical Association*, 106, 534–543. [\[1230\]](#page-2-2)
- <span id="page-10-13"></span>(2014), "A New Resampling Method for Sampling Designs Without Replacement: The Doubled Half Bootstrap,"*Computational Statistics*, 29, 1345–1363. [\[1230\]](#page-2-2)
- <span id="page-10-3"></span>Beaumont, J.-F., and Bocci, C. (2009), "A Practical Bootstrap Method for Testing Hypotheses from Survey Data," *Survey Methodology*, 35, 25–35. [\[1229,](#page-1-4)[1231](#page-3-11)[,1237\]](#page-9-2)
- <span id="page-10-11"></span>Beaumont, J.-F., and Patak, Z. (2012), "On the Generalized Bootstrap for Sample Surveys with Special Attention to Poisson Sampling," *International Statistical Review*, 80, 127–148. [\[1230](#page-2-2)[,1233\]](#page-5-2)
- <span id="page-10-24"></span>Berger, Y. G. (1998), "Rate of convergence to normal distribution for the Horvitz-Thompson estimator," *Journal of Statistical Planning and Inference*, 67, 209–226. [\[1231\]](#page-3-11)
- <span id="page-10-10"></span>Bertail, P., and Combris, P. (1997), "Bootstrap Généralisé d'un Sondage," *Annales d'Economie et de Statistique*, 46, 49–83. [\[1230,](#page-2-2)[1233\]](#page-5-2)
- <span id="page-10-30"></span>Bessonneau, P., Brilhaut, G., Chauvet, G., and Garcia, C. (2021), "With-Replacement Bootstrap Variance Estimation for Household Surveys: Principles, Examples and Implementation," *Survey Methodology*, 47, 313–347. [\[1233\]](#page-5-2)
- <span id="page-10-28"></span>Bickel, P. J., and Freedman, D. A. (1984), "Asymptotic Normality and the Bootstrap in Stratified Sampling," *Annals of Statistics*, 12, 470–482. [\[1232,](#page-4-6)[1233\]](#page-5-2)
- <span id="page-10-17"></span>Breidt, F. J., and Opsomer, J. D. (2000), "Local Polynomial Regression Estimators in Survey Sampling," *Annals of Statistics*, 28, 1026–1053. [\[1230\]](#page-2-2)
- <span id="page-10-1"></span>Chambers, R., and Skinner, C. J. (2003),*Analysis of Survey Data*, Chichester: Wiley. [\[1229\]](#page-1-4)
- <span id="page-10-27"></span>Chao, M.-T., and Lo, S.-H. (1985), "A Bootstrap Method for Finite Population," *Sankhya: Series A*, 47, 399–405. [\[1232\]](#page-4-6)
- <span id="page-10-6"></span>Chatterjee, S., and Bose, A. (2005), "Generalized Bootstrap for Estimating Equations," *Annals of Statistics*, 33, 414–436. [\[1229](#page-1-4)[,1230\]](#page-2-2)
- <span id="page-10-25"></span>Chauvet, G. (2015), "Coupling Methods for Multistage Sampling," *Annals of Statistics*, 43, 2484–2506. [\[1231\]](#page-3-11)
- <span id="page-10-20"></span>Chen, J., and Rao, J. (2007), "Asymptotic Normality Under Two-Phase Sampling Designs," *Statistica Sinica*, 17, 1047–1064. [\[1231\]](#page-3-11)
- <span id="page-10-9"></span>Chipperfield, J., and Preston, J. (2007), "Efficient Bootstrap for Business Surveys," *Survey Methodology*, 33, 167–172. [\[1230\]](#page-2-2)
- <span id="page-10-16"></span>Fuller, W. A. (2009), *Sampling Statistics*, New Jersey: Wiley. [\[1230](#page-2-2)[,1231,](#page-3-11)[1237\]](#page-9-2)
- <span id="page-10-21"></span>Hájek, J. (1964), "Asymptotic Theory of Rejective Sampling with Varying Probabilities from a Finite Population,"*Annals of Mathematical Statistics*, 35, 1491–1523. [\[1231\]](#page-3-11)
- <span id="page-10-19"></span>Han, Q., and Wellner, J. A. (2021), "Complex Sampling Designs: Uniform Limit Theorems and Applications," *Annals of Statistics*, 49, 459–485. [\[1230](#page-2-2)[,1231\]](#page-3-11)
- <span id="page-10-14"></span>Isaki, C. T., and Fuller, W. A. (1982), "Survey Design under the Regression Superpopulation Model," *Journal of the American Statistical Association*, 77, 89–96. [\[1230,](#page-2-2)[1231\]](#page-3-11)
- <span id="page-10-35"></span>Jiang, J. (2019), *Robust Mixed Model Analysis*, Singapore: World Scientific. [\[1237\]](#page-9-2)
- <span id="page-10-18"></span>Kim, J. K., Michael Brick, J., Fuller, W. A., and Kalton, G. (2006), "On the Bias of the Multiple-Imputation Variance Estimator in Survey Sampling," *Journal of the Royal Statistical Society*, Series B, 68, 509–521. [\[1230\]](#page-2-2)
- <span id="page-10-34"></span>Kim, J. K., Park, S., and Lee, Y. (2017), "Statistical Inference using Generalized Linear Mixed Models under Informative Cluster Sampling," *Canadian Journal of Statistics*, 45, 479–497. [\[1237\]](#page-9-2)
- <span id="page-10-0"></span>Korn, E. L., and Graubard, B. I. (1999),*Analysis of Health Surveys*, New York: Wiley. [\[1229\]](#page-1-4)
- <span id="page-10-33"></span>Lumley, T. (2020), *Analysis of Complex Survey Samples*, R package version 3.37. [\[1236\]](#page-8-3)
- <span id="page-10-2"></span>Lumley, T., and Scott, A. J. (2014), "Tests for Regression Models Fitted to Survey Data," *Australian & New Zealand Journal of Statistics*, 56, 1–14. [\[1229](#page-1-4)[,1232,](#page-4-6)[1233,](#page-5-2)[1234](#page-6-7)[,1235](#page-7-4)[,1236\]](#page-8-3)
- <span id="page-10-4"></span>Mashreghi, Z., Haziza, D., and Léger, C. (2016), "A Survey of Bootstrap Methods in Finite Population Sampling," *Statistics Surveys*, 10, 1–52. [\[1229\]](#page-1-4)
- <span id="page-10-7"></span>Pfeffermann, D. (1993), "The Role of Sampling Weights when Modeling Survey Data," *International Statistical Review*, 61, 317–337. [\[1229](#page-1-4)[,1230,](#page-2-2)[1231,](#page-3-11)[1235\]](#page-7-4)
- <span id="page-10-5"></span>Praestgaard, J., and Wellner, J. A. (1993), "Exchangeably Weighted Bootstraps of the General Empirical Process," *Annals of Probability*, 21, 2053– 2086. [\[1229\]](#page-1-4)
- <span id="page-10-32"></span>Rao, J. N. K., and Scott, A. J. (1981), "The Analysis of Categorical Data from Complex Sample Surveys: Chi-Squared Tests for Goodness of Fit and Independence in Two-Way Tables," *Journal of the American Statistical Association*, 76, 221–230. [\[1235\]](#page-7-4)
- <span id="page-10-26"></span>(1984), "On Chi-Squared Tests for Multiway Contingency Tables with Cell Proportions Estimated from Survey Data," *Annals of Statistics*, 12, 46–60. [\[1232,](#page-4-6)[1238\]](#page-10-36)
- <span id="page-10-31"></span>Rao, J. N. K., Scott, A. J., and Skinner, C. J. (1998), "Quasi-Score Tests with Survey Data," *Statistica Sinica*, 8, 1059–1070. [\[1234\]](#page-6-7)
- <span id="page-10-29"></span>Rao, J. N. K., and Wu, C. F. J. (1988), "Resampling Inference with Complex Survey Data," *Journal of the American Statistical Association*, 83, 231– 241. [\[1233\]](#page-5-2)
- <span id="page-10-8"></span>Rao, J. N. K., Wu, C. F. J., and Yue, K. (1992), "Some Recent Work on Resampling Methods for Complex Surveys," *Survey Methodology*, 18, 209–217. [\[1230](#page-2-2)[,1233\]](#page-5-2)
- <span id="page-10-22"></span>Rosén, B. (1972a), "Asymptotic Theory for Successive Sampling with Varying Probabilities without Replacement, I," *Annals of Mathematical Statistics*, 43, 373–397. [\[1231\]](#page-3-11)
- <span id="page-10-23"></span>(1972b), "Asymptotic Theory for Successive Sampling with Varying Probabilities without Replacement, II," *The Annals of Mathematical Statistics*, 43, 748–776. [\[1231\]](#page-3-11)
- <span id="page-10-15"></span>Rubin-Bleuer, S., and Schiopu-Kratina, I. (2005), "On the Two-Phase Framework for Joint Model and Design-based Inference," *Annals of Statistics*, 33, 2789–2810. [\[1230,](#page-2-2)[1231\]](#page-3-11)
- <span id="page-11-3"></span>Särndal, C.-E., Swensson, B., and Wretman, J. (2003), *Model Assisted Survey Sampling*, New York: Springer. [\[1230\]](#page-2-2)
- <span id="page-11-2"></span>Särndal, C.-E., Swensson, B., and Wretman, J. H. (1989), "The Weighted Residual Technique for Estimating the Variance of the General Regression Estimator of the Finite Population Total," *Biometrika*, 76, 527–537. [\[1230\]](#page-2-2)
- <span id="page-11-0"></span>Shao, J. (2003), *Mathematical Statistics*, New York: Springer. [\[1229\]](#page-1-4)
- <span id="page-11-1"></span>Sverchkov, M., and Pfeffermann, D. (2004), "Prediction of Finite Population Totals based on the Sample Distribution," *Survey Methodology*, 30, 79– 92. [\[1230\]](#page-2-2)
- <span id="page-11-7"></span>Thomas, D. R., and Rao, J. N. K. (1987), "Small-Sample Comparisons of Level and Power for Simple Goodness-of-Fit Statistics under Cluster

Sampling," *Journal of the American Statistical Association*, 82, 630–636. [\[1236\]](#page-8-3)

- <span id="page-11-4"></span>Toth, D., and Eltinge, J. L. (2011), "Building Consistent Regression Trees From Complex Sample Data," *Journal of the American Statistical Association*, 106, 1626–1636. [\[1230\]](#page-2-2)
- <span id="page-11-6"></span>Verret, F., Rao, J. N. K., and Hidiroglou, M. A. (2015), "Model-based Small Area Estimation Under Informative Sampling," *Survey Methodology*, 41, 333–347. [\[1235\]](#page-7-4)
- <span id="page-11-5"></span>Wu, C., and Thompson, M. E. (2020), *Sampling Theory and Practice*, Gewerbestrasse: Springer. [\[1230](#page-2-2)[,1236\]](#page-8-3)