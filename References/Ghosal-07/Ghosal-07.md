## **CONVERGENCE RATES OF POSTERIOR DISTRIBUTIONS FOR NONIID OBSERVATIONS**

## BY SUBHASHIS GHOSAL<sup>1</sup> AND AAD VAN DER VAART

## *North Carolina State University and Vrije Universiteit Amsterdam*

We consider the asymptotic behavior of posterior distributions and Bayes estimators based on observations which are required to be neither independent nor identically distributed. We give general results on the rate of convergence of the posterior measure relative to distances derived from a testing criterion. We then specialize our results to independent, nonidentically distributed observations, Markov processes, stationary Gaussian time series and the white noise model. We apply our general results to several examples of infinite-dimensional statistical models including nonparametric regression with normal errors, binary regression, Poisson regression, an interval censoring model, Whittle estimation of the spectral density of a time series and a nonlinear autoregressive model.

**1. Introduction.** Let *(*X*(n),* A*(n), P(n) <sup>θ</sup>* : *θ* ∈ *-)* be a sequence of statistical experiments with observations *X(n)*, where the parameter set  is arbitrary and *n* is an indexing parameter, usually the sample size. We put a prior distribution *n* on *θ* ∈  and study the rate of convergence of the posterior distribution *n(*·|*X(n))* under *P(n) <sup>θ</sup>*<sup>0</sup> , where *θ*<sup>0</sup> is the "true value" of the parameter. The rate of this convergence can be measured by the size of the smallest shrinking balls around *θ*<sup>0</sup> that contain most of the posterior probability. For parametric models with independent and identically distributed (i.i.d.) observations, it is well known that the posterior distribution converges at the rate *n*−1*/*2. When  is infinite-dimensional, but the observations are i.i.d., Ghosal, Ghosh and van der Vaart [\[14\]](#page-30-0) obtained rates of convergence in terms of the size of the model (measured by the metric entropy or existence of certain tests) and the concentration rate of the prior around *θ*<sup>0</sup> and computed the rate of convergence for a variety of examples. A similar result was obtained by Shen and Wasserman [\[27\]](#page-31-0) under stronger conditions.

Little is known about the asymptotic behavior of the posterior distribution in infinite-dimensional models when the observations are not i.i.d. For independent, nonidentically distributed (i.n.i.d.) observations, consistency has recently been ad-

Received January 2004; revised March 2006.

<sup>1</sup>Supported in part by NSF Grant number DMS-03-49111.

*[AMS 2000 subject classifications.](http://www.ams.org/msc/)* Primary 62G20; secondary 62G08.

*Key words and phrases.* Covering numbers, Hellinger distance, independent nonidentically distributed observations, infinite dimensional model, Markov chains, posterior distribution, rate of convergence, tests.

dressed by Amewou-Atisso, Ghosal, Ghosh and Ramamoorthi [\[1\]](#page-29-0) and Choudhuri, Ghosal and Roy [\[7\]](#page-30-0). The main purpose of the present paper is to obtain a theorem on rates of convergence of posterior distributions in a general framework not restricted to the setup of i.i.d. observations. We specialize this theorem to several classes of non-i.i.d. models including i.n.i.d. observations, Gaussian time series, Markov processes and the white noise model. The theorem applies in every situation where it is possible to test the true parameter versus balls of alternatives with exponential error probabilities and it is not restricted to any particular structure on the joint distribution. The existence of such tests has been proven in many special cases by Le Cam [\[20–22\]](#page-30-0) and Birgé [\[3–5\]](#page-30-0), who used them to construct estimators with optimal rates of convergence, determined by the (local) metric entropy or "Le Cam dimension" of the model. Our main theorem uses the same metric entropy measure of the complexity of the model and combines this with a measure of prior concentration around the true parameter to obtain a bound on the posterior rate of convergence, generalizing the corresponding result of Ghosal, Ghosh and van der Vaart [\[14\]](#page-30-0). We apply these results to obtain posterior convergence rates for linear regression, nonparametric regression, binary regression, Poisson regression, interval censoring, spectral density estimation and nonlinear autoregression. van der Meulen, van der Vaart and van Zanten [\[30\]](#page-31-0) have extended the approach of this paper to several types of diffusion models.

The organization of the paper is as follows. In the next section, we describe our main theorem in an abstract framework. In Sections [3,](#page-4-0) [4,](#page-5-0) [5](#page-8-0) and [6,](#page-9-0) we specialize to i.n.i.d. observations, Markov chains, the white noise model and Gaussian time series, respectively. In Section [7,](#page-11-0) we discuss a large number of more concrete applications, combining models of various types with many types of different priors, including priors based on the Dirichlet process, mixture representations or sequence expansions on spline bases, priors supported on finite sieves and conjugate Gaussian priors. Technical proofs, including the proofs of the main results, are collected in Section [8.](#page-27-0)

The notation will be used to denote inequality up to a constant that is fixed throughout. The notation *Pf* will abbreviate *f dP* . The symbol *x* will stand for the greatest integer less than or equal to *x*. Let *h(f, g)* = *( (f* <sup>1</sup>*/*<sup>2</sup> − *g*1*/*2*)*<sup>2</sup> *dµ)*1*/*<sup>2</sup> and *K(f, g)* = *f* log*(f/g) dµ* stand for the Hellinger distance and Kullback–Leibler divergence, respectively, between two nonnegative densities *f* and *g* relative to a measure *µ*. Furthermore, we define additional discrepancy measures by *Vk(f, g)* = *f* |log*(f/g)*| *<sup>k</sup> dµ* and *Vk,*0*(f, g)* = *f* |log*(f/g)* − *K(f, g)*| *<sup>k</sup> dµ*, *k >* 1. The index *k* = 2 of *V*<sup>2</sup> and *V*2*,*<sup>0</sup> may be omitted and these simply written as *V* and *V*0, respectively. The symbols N and R will denote the sets of natural and real numbers, respectively. The *ε*-covering number of a set for a semimetric *d*, denoted by *N (ε, -, d)*, is the minimal number of *d*-balls of radius *ε* needed to cover the set *-*; see, for example, [\[31\]](#page-31-0).

<span id="page-2-0"></span>**2. General theorem.** For each *n* ∈ N and *θ* ∈ *-*, let *P(n) <sup>θ</sup>* admit densities *<sup>p</sup>(n) θ* relative to a *σ* -finite measure *µ(n)*. Assume that *(x, θ )* → *p(n) <sup>θ</sup> (x)* is jointly measurable relative to A ⊗ B, where B is a *σ* -field on *-*. By Bayes' theorem, the posterior distribution is given by

(2.1) 
$$
\Pi_n(B|X^{(n)}) = \frac{\int_B p_\theta^{(n)}(X^{(n)}) d\Pi_n(\theta)}{\int_{\Theta} p_\theta^{(n)}(X^{(n)}) d\Pi_n(\theta)}, \qquad B \in \mathcal{B}.
$$

Here, *X(n)* is an "observation," which, in our setup, will be understood to be generated according to *P(n) <sup>θ</sup>*<sup>0</sup> for some given *θ*<sup>0</sup> ∈ *-*.

For each *n*, let *dn* and *en* be semimetrics on  with the property that there exist universal constants *ξ >* 0 and *K >* 0 such that for every *ε >* 0 and for each *θ*<sup>1</sup> ∈ with *dn(θ*1*, θ*0*)>ε*, there exists a test *φn* such that

$$
(2.2) \tP_{\theta_0}^{(n)}\phi_n \leq e^{-Kn\varepsilon^2}, \t\sup_{\theta \in \Theta: e_n(\theta,\theta_1) < \varepsilon_\xi} P_{\theta}^{(n)}(1-\phi_n) \leq e^{-Kn\varepsilon^2}.
$$

Typically, we have *dn* ≤ *en* and in many cases we choose *dn* = *en*, but using two semimetrics provides some added flexibility. Le Cam [\[20–22\]](#page-30-0) and Birgé [\[3–5\]](#page-30-0) showed that the rate of convergence, in a minimax sense, of the best estimators of *θ* relative to the distance *dn* can be understood in terms of the *Le Cam dimension* or *local entropy* function of the set  relative to *dn*. For our purposes, this dimension is a function whose value at *ε >* 0 is defined to be log*N (εξ ,*{*θ* : *dn(θ , θ*0*)* ≤ *ε*}*, en)*, that is, the logarithm of the minimum number of *dn*-balls of radius *εξ* needed to cover an *en*-ball of radius *ε* around the true parameter *θ*0. Birgé [\[3, 4\]](#page-30-0) and Le Cam [\[20–22\]](#page-30-0) showed that there exist estimators *θ*ˆ *<sup>n</sup>* = *θ*ˆ *n(X(n))* such that *dn(θ*ˆ *n, θ*0*)* = *OP (εn)* under *P(n) <sup>θ</sup>*<sup>0</sup> , where

(2.3) 
$$
\sup_{\varepsilon > \varepsilon_n} \log N(\varepsilon \xi, \{\theta : d_n(\theta, \theta_0) \le \varepsilon\}, e_n) \le n\varepsilon_n^2.
$$

Further, under certain conditions *εn* is the best rate obtainable, given the model, and hence gives a minimax rate.

As in the i.i.d. case, the behavior of posterior distributions depends on the size of the model measured by (2.3) and the concentration rate of the prior *n* at *θ*0. For a given *k >* 1, let

$$
B_n(\theta_0, \varepsilon; k) = \{ \theta \in \Theta : K(p_{\theta_0}^{(n)}, p_{\theta}^{(n)}) \le n\varepsilon^2, V_{k,0}(p_{\theta_0}^{(n)}, p_{\theta}^{(n)}) \le n^{k/2} \varepsilon^k \}.
$$

An appropriate condition will appear as a lower bound on *n(Bn(θ*0; *ε, k))* with *k* = 2 being good enough to establish convergence in mean. For almost sure convergence, or convergence of the posterior mean, better control may be needed (through a larger value of *k*), depending on the rate of convergence.

The following result, generalizing Theorem 2.4 of Ghosal, Ghosh and van der Vaart [\[14\]](#page-30-0) for the i.i.d. case, bounds the rate of posterior convergence.

<span id="page-3-0"></span>THEOREM 1. *Let dn and en be semimetrics on for which tests satisfying* [\(2.2\)](#page-2-0) *exist*. *Let εn >* 0, *εn* → 0, *(nε*<sup>2</sup> *n)*<sup>−</sup><sup>1</sup> <sup>=</sup> *O(*1*)*, *k >* 1, *and <sup>n</sup>* ⊂  *be such that for every sufficiently large j* ∈ N,

(2.4) 
$$
\sup_{\varepsilon > \varepsilon_n} \log N\left(\frac{1}{2}\varepsilon\xi, \{\theta \in \Theta_n : d_n(\theta, \theta_0) < \varepsilon\}, e_n\right) \le n\varepsilon_n^2,
$$

$$
(2.5) \qquad \frac{\Pi_n(\theta \in \Theta_n : j\varepsilon_n < d_n(\theta, \theta_0) \le 2j\varepsilon_n)}{\Pi_n(B_n(\theta_0, \varepsilon_n; k))} \le e^{Kn\varepsilon_n^2 j^2/2}.
$$

*Then for every Mn* → ∞, *we have that*

(2.6) 
$$
P_{\theta_0}^{(n)} \Pi_n (\theta \in \Theta_n : d_n(\theta, \theta_0) \geq M_n \varepsilon_n |X^{(n)}| \to 0.
$$

The theorem uses the fact that *<sup>n</sup>* ⊂  to alleviate the entropy condition (2.4), but returns an assertion about the posterior distribution on *<sup>n</sup>* only. The complementary assertion *P(n) θ*0 *n(-* \ *<sup>n</sup>*|*X(n))* → 0 may be handled either by a direct argument or by the following analog of Lemma 5 of [\[2\]](#page-30-0).

LEMMA 1. If 
$$
\frac{\Pi_n(\Theta \setminus \Theta_n)}{\Pi_n(B_n(\theta_0, \varepsilon_n; k))} = o(e^{-2n\varepsilon_n^2})
$$
 for some  $k > 1$ , then  $P_{\theta_0}^{(n)} \Pi_n(\Theta \setminus \Theta_n | X^{(n)}) \to 0$ .

The choice *<sup>n</sup>* = *-*, which makes the condition of Lemma 1 trivial, imposes a much stronger restriction on (2.4) and is generally unattainable when  is not compact.

The following theorem extends the convergence in Theorem 1 to almost sure convergence and yields a rate for the convergence under slightly stronger conditions.

THEOREM 2. *In the situation of Theorem* 1,

(i) *if all X(n) are defined on a fixed sample space and εn n*−*<sup>α</sup> for some α* ∈ *(*0*,* 1*/*2*) such that k(*1 − 2*α) >* 2, *then the convergence* (2.6) *also holds in the almost sure sense*;

(ii) *if εn n*−*<sup>α</sup> for some α* ∈ *(*0*,* 1*/*2*) such that k(*1 − 2*α) >* 4*α*, *then the left side of* (2.6) *is O(ε*<sup>2</sup> *n)*.

If  is a convex set and *d*<sup>2</sup> *<sup>n</sup>* is a convex function in one argument keeping the other fixed and is bounded above by *B*, then for *θ*ˆ *<sup>n</sup>* = *θ dn(θ*|*X(n))*, we have, by Jensen's inequality, that

$$
d_n^2(\hat{\theta}_n, \theta_0) \le \int d_n^2(\theta, \theta_0) d \Pi_n(\theta|X^{(n)}) \le \varepsilon_n^2 + B^2 \Pi_n(d_n(\theta, \theta_0) \ge \varepsilon_n|X^{(n)}).
$$

This yields the rate *εn* for the point estimator *θ*ˆ *<sup>n</sup>* under the conditions of Theorem 1.

<span id="page-4-0"></span>The complicated-looking condition [\(2.5\)](#page-3-0) can often be simplified in infinitedimensional cases, where, typically, *nε*<sup>2</sup> *<sup>n</sup>* → ∞. Because the numerator in [\(2.5\)](#page-3-0) is trivially bounded by one, a sufficient condition for [\(2.5\)](#page-3-0) is that *n(Bn(θ*0*, εn, k)) e*−*cnε*<sup>2</sup> *<sup>n</sup>* . The local entropy in condition [\(2.4\)](#page-3-0) can also often be replaced by the global entropy log*N (εξ/*2*, n, en)* without affecting rates. Also, if the prior is such that the minimax rate given by [\(2.3\)](#page-2-0) satisfies [\(2.5\)](#page-3-0) and the condition of Lemma [1,](#page-3-0) then the posterior convergence rate attains the minimax rate.

Entropy conditions, however, may not always be appropriate to ensure the existence of tests. Ad hoc tests may sometimes be more conveniently constructed. A more general theorem on convergence rates, which is formulated directly in terms of tests and stated below, may be proven in a similar manner.

THEOREM 3. *Let dn be a semimetric on -*, *εn* → 0, *(nε*<sup>2</sup> *n)*<sup>−</sup><sup>1</sup> <sup>=</sup> *O(*1*)*, *k >* 1, *K >* 0, *<sup>n</sup>* ⊂ *and φn be a sequence of test functions such that*

$$
(2.7) \t P_{\theta_0}^{(n)} \phi_n \to 0, \t \sup_{\theta \in \Theta_n : j \varepsilon_n < d_n(\theta, \theta_0) \leq 2j \varepsilon_n} P_{\theta}^{(n)} (1 - \phi_n) \lesssim e^{-Kj^2 n \varepsilon_n^2}
$$

*and* [\(2.5\)](#page-3-0) *holds*. *Then for every Mn* → ∞, *we have that P(n) θ*0 *n(θ* ∈ *n* : *dn(θ , θ*0*)* ≥ *Mnεn*|*X(n))* → 0.

**3. Independent observations.** In this section, we consider the case where the observation *X(n)* is a vector *X(n)* = *(X*1*, X*2*,...,Xn)* of independent observations *Xi*. Thus, we take the measures *P(n) θ* of Section [2](#page-2-0) equal to product measures *<sup>n</sup> <sup>i</sup>*=<sup>1</sup> *Pθ ,i* on a product measurable space *<sup>n</sup> <sup>i</sup>*=1*(*X*i,* A*i)*. We assume that the distribution *Pθ ,i* of the *i*th component *Xi* possesses a density *pθ ,i* relative to a *σ* -finite measure *µi* on *(*X*i,* A*i)*, *i* = 1*,...,n*. In this case, tests can be constructed relative to the semimetric *dn*, whose square is given by

(3.1) 
$$
d_n^2(\theta, \theta') = \frac{1}{n} \sum_{i=1}^n \int (\sqrt{p_{\theta, i}} - \sqrt{p_{\theta', i}})^2 d\mu_i.
$$

Thus, *d*<sup>2</sup> *<sup>n</sup>* is the average of the squares of the Hellinger distances for the distributions of the individual observations.

The following lemma, due to Birgé (cf. [\[22\]](#page-30-0), page 491, or [\[4\]](#page-30-0), Corollary 2 on page 149), guarantees the existence of tests satisfying the conditions of [\(2.2\)](#page-2-0).

LEMMA 2. *If P(n) <sup>θ</sup> are product measures and dn is defined by* (3.1), *then there exist tests φn such that P(n) <sup>θ</sup>*<sup>0</sup> *φn* <sup>≤</sup> *<sup>e</sup>*−<sup>1</sup> <sup>2</sup> *nd*<sup>2</sup> *n(θ*0*,θ*1*) and P(n) <sup>θ</sup> (*1−*φn)* <sup>≤</sup> *<sup>e</sup>*−<sup>1</sup> <sup>2</sup> *nd*<sup>2</sup> *n(θ*0*,θ*1*) for all θ* ∈  *such that dn(θ , θ*1*)* ≤ <sup>1</sup> <sup>18</sup> *dn(θ*0*, θ*1*)*.

The Kullback–Leibler divergence between product measures is equal to the sum of the Kullback–Leibler divergences between the individual components. <span id="page-5-0"></span>Furthermore, as a consequence of the Marcinkiewiz–Zygmund inequality (e.g., [\[9\]](#page-30-0), page 356), the mean *Y*¯ *<sup>n</sup>* of *n* independent random variables satisfies E|*Y*¯ *<sup>n</sup>* − E*Y*¯ *n*| *<sup>k</sup>* ≤ *Ckn*<sup>−</sup>*k/*<sup>2</sup> <sup>1</sup> *n <sup>n</sup> <sup>i</sup>*=<sup>1</sup> E|*Yi*| *<sup>k</sup>* for *k* ≥ 2, where *Ck* is a constant depending only on *k*. Therefore, the set *Bn(θ*0*, ε*; *k)* contains the set

$$
B_n^*(\theta_0, \varepsilon; k) = \left\{\theta \in \Theta : \frac{1}{n} \sum_{i=1}^n K_i(\theta_0, \theta) \le \varepsilon^2, \ \frac{1}{n} \sum_{i=1}^n V_{k,0;i}(\theta_0, \theta) \le C_k \varepsilon^k \right\},\
$$

where *Ki(θ*0*,θ)* = *K(Pθ*0*,i, Pθ ,i)* and *Vk,*0;*i(θ*0*,θ)* = *Vk,*0*(Pθ*0*,i, Pθ ,i)*. Thus, we can work with a "ball" around *θ*<sup>0</sup> relative to the average Kullback–Leibler divergence and the average *k*th order moments, as in the preceding display, and simplify Theorem [1](#page-3-0) to the following result:

THEOREM 4. *Let P(n) <sup>θ</sup> be product measures and dn be defined by* [\(3.1\)](#page-4-0). *Suppose that for a sequence εn* → 0 *such that nε*<sup>2</sup> *<sup>n</sup> is bounded away from zero*, *some k >* 1, *all sufficiently large j and sets <sup>n</sup>* ⊂ *-*, *the following conditions hold*:

(3.2) 
$$
\sup_{\varepsilon > \varepsilon_n} \log N(\varepsilon/36, \{\theta \in \Theta_n : d_n(\theta, \theta_0) < \varepsilon\}, d_n) \le n\varepsilon_n^2;
$$

(3.3) 
$$
\frac{\Pi_n(\Theta \setminus \Theta_n)}{\Pi_n(B_n^*(\theta_0, \varepsilon_n; k))} = o(e^{-2n\varepsilon_n^2});
$$

(3.4) 
$$
\frac{\Pi_n(\theta \in \Theta_n : j\varepsilon_n < d_n(\theta, \theta_0) \leq 2j\varepsilon_n)}{\Pi_n(B_n^*(\theta_0, \varepsilon_n; k))} \leq e^{n\varepsilon_n^2 j^2/4}.
$$

*Then P(n) θ*0 *n(θ* : *dn(θ , θ*0*)* ≥ *Mnεn*|*X(n))* → 0 *for every Mn* → ∞.

The average Hellinger distance is not always the most natural choice. It can be replaced by any other distance *dn* that satisfies (3.2)–(3.3) and for which the conclusion of Lemma [2](#page-4-0) holds. Often, we set *k* = 2 and work with the smaller neighborhood

$$
(3.5) \qquad \bar{B}_n(\theta_0, \varepsilon) = \left\{ \theta : \frac{1}{n} \sum_{i=1}^n K_i(\theta_0, \theta) \le \varepsilon^2, \frac{1}{n} \sum_{i=1}^n V_{2;i}(\theta_0, \theta) \le \varepsilon^2 \right\}.
$$

**4. Markov chains.** For *θ* ranging over a set *-*, let *(x, y)* → *pθ (y*|*x)* be a collection of transition densities from a measurable space *(*X*,* A*)* into itself, relative to some reference measure *ν*. Thus, for each *θ* ∈ *-*, the map *(x, y)* → *pθ (y*|*x)* is measurable and for each *x*, the map *y* → *pθ (y*|*x)* is a probability density relative to *µ*. Let *X*0*, X*1*,...* be a stationary Markov chain generated according to the transition density *pθ* , where it is assumed that there exists a stationary distribution *Qθ* with *µ*-density *qθ* . Let *P(n) <sup>θ</sup>* be the law of *(X*0*, X*1*,...,Xn)*.

<span id="page-6-0"></span>Tests satisfying the conditions of [\(2.2\)](#page-2-0) can be obtained from results of Birgé [\[4\]](#page-30-0), which are more refined versions of his own results in [\[3\]](#page-30-0). A special case is presented as Lemma 3 below. Actually, Birgé's result ([\[4\]](#page-30-0), Theorem 3, page 155) is much more general in that it also applies to nonstationary chains and allows different upper and lower bounds, as seen in the following display.

Assume that there exists a finite measure *ν* on *(*X*,* A*)* such that, for some *k,l* ∈ N, every *θ* ∈ and every *x* ∈ X and *A* ∈ A,

(4.1) 
$$
P_{\theta}(X_l \in A | X_0 = x) \lesssim \nu(A) \lesssim \frac{1}{k} \sum_{j=1}^k P_{\theta}(X_j \in A | X_0 = x),
$$

where P*<sup>θ</sup>* is the generic notation for any probability law governed by *θ*. For instance, if there exists a *µ*-integrable function *r* such that *r(y) pθ (y*|*x) r(y)* for every *(x, y)*, then (4.1) holds with the measure *ν* given by *dν(y)* = *r(y) dµ(y)*. Define the square of a semidistance *d* by

(4.2) 
$$
d^2(\theta, \theta') = \iint \left[ \sqrt{p_\theta(y|x)} - \sqrt{p_{\theta'}(y|x)} \right]^2 d\mu(y) d\nu(x).
$$

LEMMA 3. *If there exist k*, *l and a measure ν such that* (4.1) *holds*, *then there exist a constant K depending only on (k, l) and tests φn such that*

$$
P_{\theta_0}^{(n)} \phi_n \le e^{-Knd^2(\theta_0, \theta_1)}, \qquad \sup_{\theta \in \Theta: d(\theta, \theta_1) \le d(\theta_0, \theta_1)/8} P_{\theta}^{(n)}(1 - \phi_n) \le e^{-Knd^2(\theta_0, \theta_1)}.
$$

The preceding lemma is also true if the chain is not started at stationarity. If, as we assume, *X*<sup>0</sup> is generated from a stationary distribution under *θ*0, then the Kullback–Leibler divergence of *P(n) θ*0 and *P(n) <sup>θ</sup>* satisfies

(4.3) 
$$
K(P_{\theta_0}^{(n)}, P_{\theta}^{(n)}) = n \int K(p_{\theta_0}(\cdot|x), p_{\theta}(\cdot|x)) dQ_{\theta_0}(x) + K(q_{\theta_0}, q_{\theta}).
$$

To handle the neighborhoods *Bn(θ*0*, ε*; 2*)*, we need a bound on *V (P(n) <sup>θ</sup>*<sup>0</sup> *, P(n) <sup>θ</sup> )*, which will also be of the order of *n* times an expression depending only on individual observations, under a variety of conditions. In the following lemma, we use an *α*-mixing assumption. For a sequence {*Xn*}, let the *α*-mixing coefficient be given by *αh* = sup{| Pr*(X*<sup>0</sup> ∈ *A,Xh* ∈ *B)* − Pr*(X*<sup>0</sup> ∈ *A)* Pr*(Xh* ∈ *B)*| : *A,B* ∈ B*(*R*)*}.

LEMMA 4. *Suppose that the Markov chain X*0*, X*1*,... is α-mixing under θ*0, *with mixing coefficients αh*. *Then for every s >* 2, *V (p(n) <sup>θ</sup>*<sup>0</sup> *, p(n) <sup>θ</sup> ) is bounded by*

$$
\frac{8sn}{s-2}\sum_{h=0}^{\infty}\alpha_h^{1-2/s}\left(\int\!\!\int\!\!\left|\log\frac{p_{\theta_0}(y|x)}{p_{\theta}(y|x)}\right|^s p_{\theta_0}(y|x)\,d\mu(y)\,dQ_{\theta_0}(x)\right)^{2/s}+2V(q_{\theta_0},q_{\theta}).
$$

<span id="page-7-0"></span>PROOF. We can write

(4.4) 
$$
\log \frac{p_{\theta_0}^{(n)}}{p_{\theta}^{(n)}} = \sum_{i=1}^n \log \frac{p_{\theta_0}(X_i|X_{i-1})}{p_{\theta}(X_i|X_{i-1})} + \log \frac{q_{\theta_0}(X_0)}{q_{\theta}(X_0)} =: n\bar{Y}_n + Z_0,
$$

where *Yi* = log*(pθ*<sup>0</sup> *(Xi*|*Xi*<sup>−</sup>1*)/pθ (Xi*|*Xi*<sup>−</sup>1*))* and *Z*<sup>0</sup> = log*(qθ*<sup>0</sup> *(X*0*)/qθ (X*0*))*. Then *Y*1*, Y*2*,...* are *α*-mixing with mixing coefficients *αh*−1. Therefore, the variance of the left-hand side of (4.4) is bounded above by *n(*E|*Yi*| *s)*2*/s* × 4*s(s* − 2*)*−<sup>1</sup> <sup>∞</sup> *<sup>h</sup>*=<sup>1</sup> *α* 1−2*/s <sup>h</sup>*−<sup>1</sup> , by the bound of Ibragimov [\[18\]](#page-30-0).

Let *-*<sup>1</sup> ⊂  be the set of parameter values such that *K(qθ*<sup>0</sup> *, qθ )* and *V (qθ*<sup>0</sup> *, qθ )* are bounded by 1. Then from [\(4.3\)](#page-6-0) and Lemma [4,](#page-6-0) it follows that for large *n* and *ε*<sup>2</sup> ≥ 2*/n*, the set *Bn(θ*0*, ε*; 2*)* contains the set *B*∗*(θ*0*, ε*;*s)* defined by

$$
\left\{\theta\in\Theta_1:\mathrm{P}_{\theta_0}\log\bigg(\frac{p_{\theta_0}}{p_{\theta}}(X_1|X_0)\bigg)\leq\frac{1}{2}\varepsilon^2,\mathrm{P}_{\theta_0}\Big|\mathrm{log}\,\frac{p_{\theta_0}}{p_{\theta}}(X_1|X_0)\Big|^s\leq C_s\varepsilon^s\right\},\right\}
$$

where the power *s* must be chosen sufficiently large to ensure that the mixing coefficients satisfy <sup>∞</sup> *<sup>h</sup>*=<sup>0</sup> *α* 1−2*/s <sup>h</sup> <sup>&</sup>lt;* <sup>∞</sup> and where *<sup>C</sup>*−2*/s <sup>s</sup>* <sup>=</sup> <sup>16</sup>*s(*<sup>2</sup> <sup>−</sup> *s)*−<sup>1</sup> <sup>∞</sup> *<sup>h</sup>*=<sup>0</sup> *α* 1−2*/s <sup>h</sup>* . The contributions of *Qθ*<sup>0</sup> *(*log*(qθ*<sup>0</sup> */qθ ))* and *Qθ*<sup>0</sup> *(*log*(qθ*<sup>0</sup> */qθ ))*<sup>2</sup> may also be incorporated into the bound.

The above facts may be combined to obtain the following result.

THEOREM 5. *Let P(n) <sup>θ</sup> be the distribution of (X*0*, X*1*,...,Xn) for a stationary Markov chain X*0*, X*1*,... with transition densities pθ (y*|*x) and stationary density qθ satisfying* [\(4.1\)](#page-6-0) *and let d be defined by* [\(4.2\)](#page-6-0). *Assume*, *further*, *that the chain is αmixing with coefficients αh satisfying* <sup>∞</sup> *<sup>h</sup>*=<sup>0</sup> *α* 1−1*/s <sup>h</sup> <* ∞ *for some s >* 2. *Suppose that for a sequence εn* → 0 *such that nε*<sup>2</sup> *<sup>n</sup>* ≥ 2, *some s >* 2, *every sufficiently large j and sets <sup>n</sup>* ⊂ *-*, *the following conditions are satisfied*:

(4.5) 
$$
\sup_{\varepsilon > \varepsilon_n} \log N(\varepsilon/16, \{\theta \in \Theta_n : d(\theta, \theta_0) < \varepsilon\}, d) \le n\varepsilon_n^2;
$$

(4.6) 
$$
\frac{\Pi_n(\Theta \setminus \Theta_n)}{\Pi_n(B^*(\theta_0, \varepsilon_n; s))} = o(e^{-2n\varepsilon_n^2});
$$

(4.7) 
$$
\frac{\Pi_n(\theta \in \Theta_n : (j-1)\varepsilon_n < d(\theta, \theta_0) \leq j\varepsilon_n)}{\Pi_n(B^*(\theta_0, \varepsilon_n; s))} \leq e^{Kn\varepsilon_n^2 j^2/8},
$$

*for the constant K of Lemma* [3.](#page-6-0) *Then P(n) θ*0 *n(θ* : *d*∗*(θ , θ*0*)* ≥ *Mnεn*|*X(n))* → 0 *for every Mn* → ∞.

A Markov chain with *n*-step transition probability *Pn(x,*·*)* = Pr*(Xn* ∈ *A*|*X*<sup>0</sup> = *x)* and stationary measure *Q* is called *uniformly ergodic* if *Pn(x,*·*)* − *Q* → 0 as *n* → ∞, uniformly in *x*, where · is the total variation norm. It can be <span id="page-8-0"></span>shown that the convergence is then automatically exponentially fast (cf. [\[23\]](#page-30-0), Theorem 16.0.2). Thus, the *α*-mixing coefficients are exponentially decreasing and hence satisfy <sup>∞</sup> *<sup>h</sup>*=<sup>0</sup> *α* 1−2*/s <sup>h</sup> <* ∞ for every *s >* 2. Hence, it suffices to verify [\(4.7\)](#page-7-0) with some arbitrary fixed *s >* 2. If sup{ |*pθ*<sup>0</sup> *(y*|*x*1*)* − *pθ*<sup>0</sup> *(y*|*x*2*)*| *dµ(y)* : *x*1*, x*<sup>2</sup> ∈ R} *<* 2, then integrating out *x*<sup>2</sup> relative to the stationary measure *qθ*<sup>0</sup> , we see that Condition (16.8) of Theorem 16.0.2 of [\[23\]](#page-30-0) holds and hence the chain is uniformly ergodic.

**5. White noise model.** Let *-* ⊂ *L*2[0*,* 1] and for *θ* ∈ *-*, let *P(n) <sup>θ</sup>* be the distribution on *C*[0*,* 1] of the stochastic process *X(n)* = *(X(n) <sup>t</sup>* : 0 ≤ *t* ≤ 1*)* defined structurally as *X(n) <sup>t</sup>* <sup>=</sup> *<sup>t</sup>* <sup>0</sup> *θ (s) ds* + <sup>√</sup> 1 *<sup>n</sup>Wt* for a standard Brownian motion *W*. This is the standard white noise model, which is known to arise as an approximation of many particular sequences of experiments. An equivalent experiment is obtained by the one-to-one correspondence of *X(n)* with the sequence defined by *Xn,i* = *X(n), ei*, where ·*,*· is the inner product of *L*2[0*,* 1] and {*e*1*, e*2*,...*} is a given orthonormal basis of *L*2[0*,* 1]. The variables *Xn,*1*, Xn,*2*,...* are independent and normally distributed, with means *θ,ei* and variance *n*<sup>−</sup>1. In the following, we use this concrete representation and abuse notation by identifying *X(n)* with the sequence *(Xn,*1*, Xn,*2*, . . .)* and *θ* ∈  with the sequence *(θ*1*, θ*2*, . . .)* defined by *θi* = *θ,ei*. In the latter representation, we have that *-* ⊂ 2, the space of square summable sequences. Let *θ*<sup>2</sup> = <sup>1</sup> <sup>0</sup> *<sup>θ</sup>* <sup>2</sup>*(s) ds* <sup>=</sup> <sup>∞</sup> *<sup>i</sup>*=<sup>1</sup> *θ* <sup>2</sup> *<sup>i</sup>* denote the squared *L*2 norm.

Tests satisfying the conditions of [\(2.2\)](#page-2-0) can easily be found explicitly, namely, as the likelihood ratio test for *θ*<sup>0</sup> versus *θ*1, where we can use the *L*2-norm for both *dn* and *en*. Furthermore, the Kullback–Leibler divergence and discrepancy *V*2*,*<sup>0</sup> also turn out to be multiples of the *L*2-norm.

LEMMA 5. *The test φn* = 1{2*θ*<sup>1</sup> − *θ*0*, X(n) > θ*1<sup>2</sup> − *θ*02} *satisfies P(n) <sup>θ</sup>*<sup>0</sup> *φn* <sup>≤</sup> <sup>1</sup> <sup>−</sup> *(*√*nθ*<sup>1</sup> <sup>−</sup> *<sup>θ</sup>*0*/*2*) and <sup>P</sup>(n) <sup>θ</sup> (*<sup>1</sup> <sup>−</sup> *φn)* <sup>≤</sup> <sup>1</sup> <sup>−</sup> *(*√*nθ*<sup>1</sup> <sup>−</sup> *<sup>θ</sup>*0*/*4*) for any θ* ∈ *such that θ* − *θ*1≤*θ*<sup>1</sup> − *θ*0*/*4.

LEMMA 6. *For every θ,θ*<sup>0</sup> ∈ *-* ⊂ *L*2[0*,* 1], *we have K(P(n) <sup>θ</sup>*<sup>0</sup> *, P(n) <sup>θ</sup> )* <sup>=</sup> <sup>1</sup> <sup>2</sup>*nθ* − *θ*0<sup>2</sup> *and V*2*,*0*(P(n) <sup>θ</sup>*<sup>0</sup> *, P(n) <sup>θ</sup> )* = *nθ* − *θ*02. *Consequently*, *we have Bn(θ*0*, ε*; 2*)* = {*θ* ∈ *-*: *θ* − *θ*0 ≤ *ε*}.

PROOF OF LEMMA 5. The test rejects the null hypothesis for positive values of the statistic *Tn* = *θ*<sup>1</sup> − *θ*0*, X(n)* − <sup>1</sup> <sup>2</sup> *θ*1<sup>2</sup> <sup>+</sup> <sup>1</sup> <sup>2</sup> *θ*02, which, under *<sup>θ</sup>*, is distributed as *θ*<sup>1</sup> − *θ*0*, θ* − *θ*1 + <sup>1</sup> <sup>2</sup> *θ*<sup>1</sup> <sup>−</sup> *<sup>θ</sup>*0<sup>2</sup> <sup>+</sup> <sup>√</sup> 1 *<sup>n</sup> θ*<sup>1</sup> <sup>−</sup> *<sup>θ</sup>*0*,W*. The variable *θ*<sup>1</sup> − *θ*0*,W* is normally distributed with mean zero and variance *θ*<sup>1</sup> − *θ*02. Under *θ* = *θ*0, the mean of the test statistic is equal to −<sup>1</sup> <sup>2</sup> *θ*<sup>0</sup> <sup>−</sup> *<sup>θ</sup>*12, whereas for *θ* − *θ*1 ≤ *ξθ*<sup>1</sup> − *θ*0 and *ξ* ∈ *(*0*,* <sup>1</sup> <sup>2</sup> *)*, the mean of the statistic under *θ* is bounded <span id="page-9-0"></span>below by *(* <sup>1</sup> <sup>2</sup> <sup>−</sup>*ξ )θ*0−*θ*12, in view of the Cauchy–Schwarz inequality. The lemma follows upon choosing *ξ* = 1*/*4.

PROOF OF LEMMA [6.](#page-8-0) We write log*(p(n) <sup>θ</sup>*<sup>0</sup> */p(n) <sup>θ</sup> )* <sup>=</sup> *<sup>n</sup>θ*<sup>0</sup> <sup>−</sup>*θ,X(n)* − *<sup>n</sup>* <sup>2</sup> *θ*0<sup>2</sup> <sup>+</sup> *<sup>n</sup>* <sup>2</sup> *θ*2, whence the mean and variance under *<sup>θ</sup>*<sup>0</sup> are easily obtained.

In the preceding lemmas, no restriction on the parameter set *-* ⊂ *L*2[0*,* 1] was imposed. The lemmas lead to the following theorem, which gives bounds on the rate of convergence in terms of quantities involving the *L*2-norm only.

THEOREM 6. *Let P(n) <sup>θ</sup> be the distribution on C*[0*,* 1] *of the solution of the diffusion equation dXt* = *θ (t) dt* +*n*−1*/*<sup>2</sup> *dWt with X*<sup>0</sup> = 0. *Suppose that for εn* → 0, *(nε*<sup>2</sup> *n)*<sup>−</sup><sup>1</sup> <sup>=</sup> *O(*1*) and -*⊂ *L*2[0*,* 1], *the following conditions are satisfied*:

(5.1) 
$$
\sup_{\varepsilon > \varepsilon_n} \log N(\varepsilon/8, \{\theta \in \Theta : \|\theta - \theta_0\| < \varepsilon\}, \|\cdot\|) \le n\varepsilon_n^2;
$$

*for every j* ∈ N

(5.2) 
$$
\frac{\Pi_n(\theta \in \Theta : \|\theta - \theta_0\| \leq j\varepsilon_n)}{\Pi_n(\theta \in \Theta : \|\theta - \theta_0\| \leq \varepsilon_n)} \leq e^{n\varepsilon_n^2 j^2/64}.
$$

*Then P(n) θ*0 *n(θ* ∈ *-*: *θ* − *θ*0 ≥ *Mnεn*|*X(n))* → 0 *for every Mn* → ∞.

In Section [7.6,](#page-20-0) we shall calculate the rate of convergence for a conjugate prior.

**6. Gaussian time series.** Suppose that *X*1*, X*2*,...* is a stationary Gaussian process with mean zero and spectral density *f* , which is known to belong to a model F . Let *γh(f )* = *<sup>π</sup>* <sup>−</sup>*<sup>π</sup> eihλf (λ) dλ* be the corresponding autocovariance function. Let *P(n) <sup>f</sup>* be the distribution of *(X*1*,...,Xn)*.

For this situation, we can derive the following lemma from [\[3\]](#page-30-0). Let *f* <sup>2</sup> and *f* ∞ be the *L*2-norm relative to Lebesgue measure and the uniform norm of a function *f* : *(*−*π,π*] → R, respectively.

LEMMA 7. *Suppose that there exist constants and M such that* log *f* ∞ ≤ *and* <sup>∞</sup> *<sup>h</sup>*=−∞ |*h*|*γ* <sup>2</sup> *<sup>h</sup> (f )* ≤ *M for every f* ∈ F . *Then there exist constants ξ and K depending only on and M such that for every ε* - 1*/* <sup>√</sup>*<sup>n</sup> and every <sup>f</sup>*0*, f*<sup>1</sup> <sup>∈</sup> <sup>F</sup> *with f*<sup>1</sup> − *f*0<sup>2</sup> ≥ *ξε*, *we have*

(6.1) 
$$
P_{f_0}^{(n)}\phi_n \vee \sup_{f \in \mathcal{F}: \|f-f_1\|_{\infty} \leq \xi \varepsilon} P_f^{(n)}(1-\phi_n) \leq e^{-Kn\varepsilon^2}.
$$

PROOF. It follows from the assumptions that <sup>|</sup>*h*|*>n/*<sup>2</sup> *γ* <sup>2</sup> *<sup>h</sup> (f )* ≤ 2*M/n*. This is bounded by *<sup>ε</sup>*<sup>2</sup> for *<sup>ε</sup>* <sup>≥</sup> <sup>√</sup>2*M/n*. The assertion follows from Proposition 5.5, page 222 of [\[3\]](#page-30-0), with *φn* = 1{log*(p(n) <sup>f</sup>*<sup>1</sup> */p(n) <sup>f</sup>*<sup>0</sup> *)* ≥ 0}.

The preceding lemma shows that tests satisfying the conditions of [\(2.2\)](#page-2-0) exist when *dn* is the *L*2-distance and when *en* is the uniform distance, leading to conditions in terms of *N (εξ ,*{*f* ∈ F : *f* − *f*0<sup>2</sup> *< ε*}*,* ·∞*)*. We do not know if the *L*∞-distance can be replaced by the *L*2-distance. The uniform bound on log *f* ∞ is not unreasonable as it is known that the structure of the time series changes dramatically if the spectral density approaches zero. The following lemma allows the neighborhoods *Bn(f*0*, ε*; 2*)* to be dealt with entirely in terms of balls for the *L*2 norm.

LEMMA 8. *Suppose that there exists constant such that* log *f* ∞ ≤ *for every f* ∈ F . *Then there exists a constant C depending only on such that for every f, g* ∈ F , *we have P(n) <sup>f</sup> (*log*(p(n) <sup>f</sup> /p(n) <sup>g</sup> )) Cnf* − *g*<sup>2</sup> <sup>2</sup> *and* var*P(n) f (*log*(p(n) <sup>f</sup> /p(n) <sup>g</sup> )) Cnf* − *g*<sup>2</sup> 2.

PROOF. The *(k, l)*th element of the covariance matrix *Tn(f )* of *X(n)* = *(X*1*,...,Xn)*, given the spectral density *f* , is given by *<sup>π</sup>* <sup>−</sup>*<sup>π</sup> <sup>e</sup>iλ(k*<sup>−</sup>*l)f (λ) dλ* for 1 <sup>≤</sup> *k*, *l* ≤ *n*. Using the matrix identities det*(AB*−1*)* = det*(I* + *B*−1*/*2*(A* − *B)B*−1*/*2*)* and *A*−<sup>1</sup> − *B*−<sup>1</sup> = *A*−1*(A* − *B)B*<sup>−</sup>1, we can write

$$
\log \frac{p_f^{(n)}}{p_g^{(n)}} = -\frac{1}{2} \log \det(I + T_n(g)^{-1/2} T_n(f - g) T_n(g)^{-1/2})
$$

$$
- \frac{1}{2} (X^{(n)})^T T_n(f)^{-1} T_n(g - f) T_n(g)^{-1} X^{(n)}.
$$

For a random vector *X* with mean zero and covariance matrix , we have E*(X<sup>T</sup> AX)* = tr*(A)* and var*(X<sup>T</sup> AX)* = tr*(AA)* + tr*(AAT )*. Hence,

$$
P_f^{(n)}\left(\log \frac{p_f^{(n)}}{p_g^{(n)}}\right) = -\frac{1}{2} \log \det(I + T_n(g)^{-1/2} T_n(f - g) T_n(g)^{-1/2}) - \frac{1}{2} \text{tr}(T_n(g - f) T_n(g)^{-1}),
$$

$$
4 \operatorname{var}_{P_f^{(n)}} \left( \log \frac{p_f^{(n)}}{p_g^{(n)}} \right) = \operatorname{tr}(T_n(g - f) T_n(g)^{-1} T_n(g - f) T_n(g)^{-1}) + \operatorname{tr}(T_n(g - f) T_n(g)^{-1} T_n(f) T_n(g)^{-1} T_n(g - f) T_n(f)^{-1}).
$$

Define matrix norms by *A*<sup>2</sup> = *k l a*2 *k,l* = tr*(AAT )* and |*A*| = sup{*Ax* : *x* = 1}, where *x* is the Euclidean norm. Then tr*(A*2*)* ≤ *A*<sup>2</sup> and *AB* ≤ |*A*|*B*. Furthermore, as a result of the inequalities −<sup>1</sup> <sup>2</sup>*µ*<sup>2</sup> <sup>≤</sup> log*(*<sup>1</sup> <sup>+</sup> *µ)* <sup>−</sup> *µ* ≤ 0, for all *µ* ≥ 0, we have for any nonnegative definite matrix *A* that

<span id="page-11-0"></span>−<sup>1</sup> <sup>2</sup> tr*(A*2*)* <sup>≤</sup> log det*(I* <sup>+</sup> *A)* <sup>−</sup> tr*(A)* <sup>≤</sup> 0. In view of the identities *<sup>x</sup><sup>T</sup> Tn(f )x* <sup>=</sup> *<sup>k</sup> xkeikλ* 2 *f (λ) dλ* and *x<sup>T</sup> Tn(*1*)x* = 2*πx*2, we also have that |*Tn(f )*| ≤ 2*πf* ∞ and |*Tn(f )*<sup>−</sup>1| ≤ *(*2*π )*−11*/f* ∞. To see the validity of the second inequality, we use the fact that *A*−1 ≤ *c*−<sup>1</sup> if *Ax* ≥ *cx* for all *x*. For *f* ∈ F , *f* ∞ *<* ∞ and 1*/f* ∞ *<* ∞. Furthermore,

(6.2) 
$$
||T_n(f)||^2 = \sum_{|h| < n} (n - |h|) \gamma_h^2(f) \le 2\pi n \int_{-\pi}^{\pi} f^2(\lambda) d\lambda.
$$

Using the preceding inequalities and the identity tr*(AB)* = tr*(BA)*, it is straightforward to obtain the desired bounds on the mean and variance of log*(p(n) <sup>f</sup> /p(n) <sup>g</sup> )*. 

The preceding lemmas can be combined to obtain the following theorem, where the constants *ξ* and *K* are those introduced in Lemma [7.](#page-9-0)

THEOREM 7. *Let P(n) <sup>f</sup> be the distribution of (X*1*,...,Xn) for a stationary Gaussian time series* {*Xt* : *t* = 0*,*±1*,...*} *with spectral density f* ∈ F . *Assume that there exist constants and M such that* log *f* ∞ ≤ *and <sup>h</sup>* |*h*|*γ* <sup>2</sup> *<sup>h</sup> (f )* ≤ *M for every f* ∈ F . *Let εn* ≥ 1*/* <sup>√</sup>*<sup>n</sup> satisfy*, *for every <sup>j</sup>* <sup>∈</sup> <sup>N</sup>,

$$
\sup_{\varepsilon > \varepsilon_n} \log N(\xi \varepsilon/2, \{ f \in \mathcal{F} : \| f - f_0 \|_2 \le \varepsilon \}, \| \cdot \|_\infty) \le n\varepsilon_n^2,
$$

$$
\frac{\Pi(f: \|f - f_0\|_2 \leq j\varepsilon)}{\Pi(f: \|f - f_0\|_2 \leq \varepsilon)} \lesssim e^{Kn\varepsilon_n^2 j^2/8}.
$$

*Then P(n) <sup>f</sup>*<sup>0</sup> *(f* : *f* − *f*0<sup>2</sup> ≥ *Mnεn*|*X*1*,...,Xn)* → 0 *for every Mn* → ∞.

**7. Applications.** In this section, we present a number of examples of application of the general results obtained in the preceding sections. The examples concern combinations of a variety of models with various prior distributions.

7.1. *Finite sieves.* Consider the setting of independent, nonidentically distributed observations of Section [3.](#page-4-0) We construct sequences of priors, each supported on finitely many points such that the posterior distribution converges at a rate equal to the solution of an equation involving bracketing entropy numbers. Because bracketing entropy numbers are often close to metric entropy numbers, this construction exhibits priors for which the prior mass condition [\(2.5\)](#page-3-0) is automatically satisfied. The construction is similar to that for the i.i.d. case given by Ghosal, Ghosh and van der Vaart [\[14\]](#page-30-0), Section 3. However, in this case, some extra care is needed to appropriately define the bracketing numbers in the product space of densities. In the following, we consider a componentwise bracketing.

Consider a sequence of models P *(n)* = {*P(n) <sup>θ</sup>* : *θ* ∈ *-*} of *n*-fold product measures *P(n) <sup>θ</sup>* , where each measure is given by a density *(x*1*,...,xn)* → *<sup>n</sup> <sup>i</sup>*=<sup>1</sup> *pθ ,i(xi)* relative to a product-dominating measure *<sup>n</sup> <sup>i</sup>*=<sup>1</sup> *µi*. For a given *n* and *ε >* 0, we define the *componentwise Hellinger upper bracketing number* for  to be the smallest number *N* such that there are integrable nonnegative functions *uj,i* for *j* = 1*,* 2*,...,N* and *i* = 1*,* 2*,...,n*, with the property that for any *θ* ∈ *-*, there exists some *j* such that *pθ ,i* ≤ *uj,i* for all *i* = 1*,* 2*,...,n* and *<sup>n</sup> <sup>i</sup>*=<sup>1</sup> *h*<sup>2</sup>*(pθ ,i, uj,i)*<sup>2</sup> ≤ *nε*2. We shall denote this by *Nn*<sup>⊗</sup> ] *(ε, -, dn)*.

Given a sequence of sets *<sup>n</sup>* ↑  and *εn* → 0 such that log*Nn*<sup>⊗</sup> ] *(εn, n, dn)* ≤ *nε*<sup>2</sup> *<sup>n</sup>*, let *(uj,i* : *j* = 1*,* 2*,...,N, i* = 1*,* 2*, . . . , n)* be a componentwise Hellinger upper bracketing for *<sup>n</sup>* [where *N* = *Nn*<sup>⊗</sup> ] *(εn, n, dn)*]. From this bracketing, we construct a prior distribution *n* on the collection of densities of product measures, by defining *n* to be the measure that assigns mass *N*−<sup>1</sup> to each of the joint densities *p(n) <sup>j</sup>* = ⊗*<sup>n</sup> <sup>i</sup>*=1*(uj,i/ uj,i dµi)*, *j* = 1*,* 2*,...,N*. The collection P*<sup>n</sup>* = {*p(n) <sup>j</sup>* : *j* = 1*,* 2*,...,N*} forms a sieve for the models P *(n)* and can be considered as the parameter space for a given *n*. Although it is possible for the spaces P*<sup>n</sup>* to not be embedded in a fixed space, Theorem [4](#page-5-0) still applies and implies the following result.

THEOREM 8. *Let <sup>n</sup>* ↑  *and θ*<sup>0</sup> ∈ *-*. *Assume that* log*Nn*<sup>⊗</sup> ] *(εn, n, dn)* ≤ *nε*<sup>2</sup> *<sup>n</sup> for some sequence εn* → 0 *with nε*<sup>2</sup> *<sup>n</sup>* → ∞. *Let n be the uniform measure on the renormalized collection of upper product brackets*, *as indicated previously*. *Then for all sufficiently large M*,

$$
(7.1) \tP_{\theta_0}^{(n)} \Pi_n(p^{(n)}: d_n^2(p_{\theta_0}^{(n)}, p^{(n)}) \ge M \varepsilon_n^2 | X_1, X_2, \dots, X_n) \to 0.
$$

PROOF. As P*<sup>n</sup>* consists of finitely many points, its covering number with respect to any metric is bounded by its cardinality. Thus, [\(3.2\)](#page-5-0) holds and [\(3.3\)](#page-5-0) holds trivially.

Let *θ*<sup>0</sup> ∈ *<sup>n</sup>* for all *n>n*0. For a given *n>n*0, let *j*<sup>0</sup> be the index for which *pθ*0*,i* ≤ *uj*0*,i* and *<sup>n</sup> <sup>i</sup>*=<sup>1</sup> *h*<sup>2</sup>*(pθ*0*,i, uj*0*,i)* ≤ *nε*<sup>2</sup> *<sup>n</sup>*. If *p* is a probability density, *u* is an integrable function such that *u* ≥ *p* and *v* = *u/ u*, then because 2*ab* ≤ *(a*<sup>2</sup> + *b*2*)*, it easily follows that *h*2*(p, v)* ≤ *( u dµ)*−1*/*2*h*2*(p, u)*.

For any two probability densities *p* and *q*, we have (see, e.g., Lemma 8 of [\[17\]](#page-30-0))

$$
K(p,q) \lesssim h^2(p,q) \bigg(1 + \log \bigg\| \frac{p}{q} \bigg\|_{\infty}\bigg), \qquad V(p,q) \lesssim h^2(p,q) \bigg(1 + \log \bigg\| \frac{p}{q} \bigg\|_{\infty}\bigg)^2.
$$

Together with the elementary inequalities 1 + log *x* ≤ 2 <sup>√</sup>*<sup>x</sup>* and *(*<sup>1</sup> <sup>+</sup> log *x)*<sup>2</sup> <sup>≤</sup> *(*4*x*1*/*4*)*<sup>2</sup> = 16*x*1*/*<sup>2</sup> for all *x* ≥ 1, the bounds imply that

$$
K(p,q) \lesssim h^2(p,q) \left\| \frac{p}{q} \right\|_{\infty}^{1/2}, \qquad V(p,q) \lesssim h^2(p,q) \left\| \frac{p}{q} \right\|_{\infty}^{1/2}
$$

*.*

Because *(pθ*0*,i/vj*0*,i)* ≤ *uj*0*,i dµ*, it follows that *n*−<sup>1</sup> *<sup>n</sup> <sup>i</sup>*=<sup>1</sup> *K(pθ*0*,i, vj*0*,i) ε*<sup>2</sup> *n* and *n*−<sup>1</sup> *<sup>n</sup> <sup>i</sup>*=<sup>1</sup> *V (pθ*0*,i, vj*0*,i) ε*<sup>2</sup> *<sup>n</sup>*. Thus, *<sup>n</sup> <sup>i</sup>*=<sup>1</sup> *vj*0*,i* gets prior probability equal to *N*−<sup>1</sup> ≥ *e*−*nε*<sup>2</sup> *<sup>n</sup>* and hence relation [\(3.4\)](#page-5-0) also holds for a multiple of the present *εn*. Thus, the posterior converges at the rate *εn* with respect to the metric *dn*.

7.1.1. *Nonparametric Poisson regression.* Let *X*1*, X*2*,...* be independent Poisson-distributed random variables with parameters *ψ(z*1*), ψ(z*2*), . . . ,* where *ψ* : R → *(*0*,*∞*)* is an unknown increasing link function and *z*1*, z*2*,...* are onedimensional covariates. We assume that *L* ≤ *ψ* ≤ *U* for some constants 0 *<L< U <* ∞.

If *l* ≤ *ψ* ≤ *u*, then for any *z* and *x*, we have *e*<sup>−</sup>*ψ(z)(ψ(z))x/x*! ≤ *e*<sup>−</sup>*l(z)(u(z))x/ x*!. For a pair of link functions *l* ≤ *u*, let *ql,u(x, z)* = *e*<sup>−</sup>*l(z)(u(z))x/x*! and put *f (n) l,u (x*1*, x*2*,...,xn)* <sup>=</sup> *<sup>n</sup> <sup>i</sup>*=<sup>1</sup> *ql,u(xi, zi)*. For any constants *L<λ*1*, λ*2*, µ*1*, µ*<sup>2</sup> *< U*, we have

$$
\sum_{x=0}^{\infty} \left( \left( e^{-\lambda_1} \frac{\mu_1^x}{x!} \right)^{1/2} - \left( e^{-\lambda_2} \frac{\mu_2^x}{x!} \right)^{1/2} \right)^2
$$
  
=  $(e^{-(\lambda_1 + \mu_1)/2} - e^{-(\lambda_2 + \mu_2)/2})^2 + 2e^{-(\lambda_1 + \lambda_2)/2} (e^{(\mu_1 + \mu_2)/2} - e^{\sqrt{\mu_1 \mu_2}})$   
 $\le \left( \frac{1}{2} + \frac{1}{4} L^{-1} \right) e^{U - L} (|\lambda_1 - \lambda_2|^2 + |\mu_1 - \mu_2|^2).$ 

Let *l*<sup>1</sup> ≤ *u*<sup>1</sup> and *l*<sup>2</sup> ≤ *u*<sup>2</sup> be two pairs of link functions taking their values in the interval [*L,U*]. Therefore, with P*<sup>z</sup> <sup>n</sup>* <sup>=</sup> *<sup>n</sup>*−<sup>1</sup> *<sup>n</sup> <sup>i</sup>*=<sup>1</sup> *δzi* being the empirical distributions of *z*1*, z*2*,...,zn*, we have that *d*<sup>2</sup> *n(f (n) l*1*,u*1 *, f (n) l*2*,u*2 *) (*|*l*<sup>1</sup> − *l*2| <sup>2</sup> + |*u*<sup>1</sup> − *u*2| <sup>2</sup>*) d*P*<sup>z</sup> n*. Hence, an *ε*-bracketing of the link functions with respect to the *L*2*(*P*<sup>z</sup> n)*-metric yields a componentwise Hellinger upper bracketing whose size is a multiple of *ε*. Now the *ε*-bracketing entropy numbers of the above class are bounded by a multiple of *ε*<sup>−</sup>1, relative to any *L*2-metric (cf. Theorem 2.7.5 of [\[31\]](#page-31-0)). Equating this with *nε*2, we obtain the rate *n*−1*/*<sup>3</sup> for posterior convergence, which is also the minimax rate, relative to *dn*.

In this example, the normalized upper brackets for the densities are also Poisson mass functions corresponding to the link functions equal to the upper brackets. Hence, the prior can be viewed as charging the space of link functions and the distance *dn* can also be induced on this space. This makes interpretations of the prior and the posterior, as well as the posterior convergence rate, more transparent. Further, as the space of link functions is a fixed space, proceeding as in Theorem 3.1 of [\[14\]](#page-30-0), a fixed prior not depending on *n* may be constructed such that the posterior converges at the same *n*−1*/*<sup>3</sup> rate.

7.2. *Linear regression with unknown error distribution.* Let *X*1*,...,Xn* be independent regression response variables satisfying *Xi* = *α* + *βzi* + *εi*, *i* = 1*,* 2*,...,n*, where the *zi*'s are nonrandom one-dimensional covariates lying in [−*L,L*] for some *L* and the errors *εi* are i.i.d. with density *f* following some prior . Amewou-Atisso et al. [\[1\]](#page-29-0) studied posterior consistency under this setup. Here, we refine the result to a posterior convergence rate. Assume that |*f (x)*| ≤ *C* for all *x* and all *f* in the support of . The priors for *α* and *β* are assumed to be compactly supported with positive densities in the interiors of their supports and all the parameters are assumed to be a priori independent. Let the true value of *(f, α, β)* be *(f*0*, α*0*, β*0*)*, an interior point in the support of the prior.

Let *H (ε)* be a bound for the Hellinger *ε*-entropy of the support of and suppose that *f*0*(x)/f (x)* ≤ *M(x)* for all *x*, where *Mδf*<sup>0</sup> *<* ∞, *δ >* 0. Then by Theorem 5 of [\[33\]](#page-31-0), it follows that max{*K(f*0*, f ), V (f*0*,f)*} *h*2*(f*0*,f)* × log2*(*1*/h(f*0*, f ))*. Let *a(ε)* = −log *(h(f*0*,f)* ≤ *ε)*. The posterior convergence rate for density estimation is then *εn*, given by

(7.2) 
$$
\max\{H(\varepsilon_n), a(\varepsilon_n/(\log \varepsilon_n^{-1}))\} \leq n\varepsilon_n^2.
$$

The following theorem shows that Euclidean parameters do not affect the rate.

THEOREM 9. *Under the above setup*, *if f*0*(x* − *α*<sup>0</sup> − *β*0*z)/f (x* − *α* − *βz)* ≤ *M(x) for all x, z, α, β*, *then the joint posterior of (α, β, f ) concentrates around (α*0*, β*0*, f*0*) at the rate εn defined by* (7.2), *with respect to dn*.

PROOF. We have, by *(a* + *b)*<sup>2</sup> ≤ 2*(a*<sup>2</sup> + *b*2*)*, that *h*2*(f*1*(*· − *α*<sup>1</sup> − *β*1*z), f*2*(*· − *α*<sup>2</sup> − *β*2*z))* ≤ 2*h*2*(f*1*, f*2*)* + 4*C*2|*α*<sup>1</sup> − *α*2| <sup>2</sup> + 4*C*2*L*2|*β*<sup>1</sup> − *β*2| 2, which leads to

$$
d_n^2(P_{f_1,\alpha_1,\beta_1}^{(n)}, P_{f_2,\alpha_2,\beta_2}^{(n)}) \lesssim h^2(f_1,f_2) + |\alpha_1 - \alpha_2|^2 + |\beta_1 - \beta_2|^2
$$

and hence the *dn*-entropy of the parameter space is bounded by a multiple of *H (ε)* + log <sup>1</sup> *<sup>ε</sup> H (ε)*.

To lower bound the prior probability of *B*¯*n((f*0*, α*0*, β*0*), ε*; 2*)* defined by [\(3.5\)](#page-5-0), by Theorem 5 of [\[33\]](#page-31-0) with *h* = *h(f*0*(*· − *α*<sup>0</sup> − *β*0*z), f (*· − *α* − *βz))*, we have that *K(f*0*(*· −*α*<sup>0</sup> −*β*0*z), f (*· −*α* −*βz)) h*<sup>2</sup> log <sup>1</sup> *<sup>h</sup>* and *V (f*0*(*· −*α*<sup>0</sup> −*β*0*z), f (*· − *α* −*βz)) h*<sup>2</sup> log<sup>2</sup> <sup>1</sup> *<sup>h</sup>* . Thus, a multiple of *<sup>ε</sup>*−2*e*−*ca(ε/* log *<sup>ε</sup>*−1*)* lower bounds the prior probability of [\(3.5\)](#page-5-0) and the first factor can be absorbed into the second, where *c* is a suitable positive constant. Thus, Theorem 4 implies that the posterior convergence rate with respect to *dn* is *εn*.

More concretely, if the prior is a Dirichlet mixture of normals (or its symmetrization) with the scale parameter lying between two positive numbers and the base measure having compact support, and if the true error density is also a normal mixture of this type, then by Ghosal and van der Vaart [\[16\]](#page-30-0), it follows that the convergence rate is *(*log *n)/*√*n*. The assumption of compact support of the base measure can be relaxed by using sieves. Compactness of the support of the prior for *α* and *β* may be relaxed by using sieves |*α*| ≤ *c* log *n* if these priors have

sub-Gaussian tails. Also, it is straightforward to extend the result to a multidimensional regressor. For more general error densities, one has to allow arbitrarily small scale parameters and apply the results of Ghosal and van der Vaart [\[17\]](#page-30-0) to obtain a slower rate.

Often, only the Euclidean part is of interest and an *n*−1*/*<sup>2</sup> rate of convergence is generally obtained in the classical context. The posterior of the Euclidean part is also expected to converge at an *n*−1*/*<sup>2</sup> rate and the Bernstein–von Mises theorem may hold; see [\[26\]](#page-31-0) for some results. However, as we consider *(f, α, β)* together and obtain global convergence rates, it seems unlikely that our methods will yield these improved convergence rates for the Euclidean portion of the parameter.

7.3. *Whittle estimation of the spectral density.* Let {*Xt* : *t* ∈ Z} be a second order stationary time series with mean zero and autocovariance function *γr* = E*(XtXt*+*r)*. The spectral density of the process is defined (under the assumption that *<sup>r</sup>* |*γr*| *<* ∞) by *f (λ)* = <sup>1</sup> 2*π* <sup>∞</sup> *<sup>r</sup>*=−∞ *γre*<sup>−</sup>*irπλ*, *<sup>λ</sup>* ∈ [0*,* <sup>1</sup>]; here, we have changed the original domain [−*π,π*] of spectral density to [0*,* 1] by using symmetry and then rescaling. Let *In(λ)* = *(*2*πn)*−1| *<sup>n</sup> <sup>t</sup>*=<sup>1</sup> *Xt <sup>e</sup>*−*itπλ*<sup>|</sup> 2, *λ* ∈ [0*,* 1], denote the periodogram. Because the likelihood is complicated, Whittle [\[32\]](#page-31-0) proposed as an approximate likelihood that of a sample *U*1*,...,Uν* of independent exponential variables with means *f (*2*j/n)*, *j* = 1*,...,ν*, evaluated with *Uj* = *In(*2*j/n)*, where *ν* = *n/*2. The Whittle likelihood is motivated by the fact that if *λn,i* → *λi*, *i* = 1*,...,m*, then under reasonable conditions such as mixing conditions, *(In(λn,*1*), . . . , In(λn,m))* converges weakly to a vector of independent exponential variables with mean vector *(f (λ*1*), . . . , f (λm))*; see, for instance, Theorem 10.3.2 of Brockwell and Davis [\[6\]](#page-30-0). Dahlhaus [\[10\]](#page-30-0) applied the technique of Whittle likelihood to estimating the spectral density by the minimum contrast method. A consistent Bayesian nonparametric method has been proposed by Choudhuri, Ghosal and Roy [\[7\]](#page-30-0). Below, we indicate how to obtain a rate of convergence using Theorem [4.](#page-5-0)

As in the proof of consistency, we use the contiguity result of Choudhuri, Ghosal and Roy [\[8\]](#page-30-0), which shows that for a Gaussian time series, the sequence of laws of *(In(*2*/n), . . . , In(*2*ν/n))* and the sequence of approximating exponential distributions of *(U*1*,...,Uν )* are contiguous. Thus, a rate of convergence of the posterior distribution under the actual distribution follows from a rate of convergence under the assumption that *U*1*,...,Uν* are exactly independent and exponentially distributed with means *f (*2*/n), . . . , f (*2*ν/n)*, to which Theorem [4](#page-5-0) can be applied.

Let *d*¯<sup>2</sup> *n(f*1*, f*2*)* <sup>=</sup> *<sup>ν</sup>*−<sup>1</sup> *<sup>ν</sup> <sup>i</sup>*=1*(f*1*(*2*i/n)* − *f*2*(*2*i/n))*2. If *f*<sup>1</sup> and *f*<sup>2</sup> are spectral densities with *m* ≤ *f*1*, f*<sup>2</sup> ≤ *M* pointwise, then it follows that

$$
(7.3) \quad \frac{1}{4M^2}\bar{d}_n^2(f_1, f_2) \le d_n^2(f_1, f_2) \le \frac{1}{4m^2}\bar{d}_n^2(f_1, f_2) \le \frac{1}{4m^2}||f_1 - f_2||_\infty^2,
$$

where *dn* is given by [\(3.1\)](#page-4-0) and ·∞ is the uniform distance. If the spectral densities are Lipschitz continuous, then a rate for the discretized *L*2-distance *dn* will imply a rate for the ordinary *L*2-distance ·<sup>2</sup> by the relation *f*<sup>1</sup> − *f*2<sup>2</sup> *d*¯ *n(f*1*, f*2*)* + *(L* + *M)/n*, where *L* and *M* are the Lipschitz constant and uniform bound, respectively. To see this, note that <sup>√</sup>*n/νd*¯ *n(f,* 0*)* = *fn*2, where *fn* = *<sup>ν</sup> <sup>j</sup>*=<sup>1</sup> *f (*2*j/n)*1*((*2*j*−2*)/n,*2*j/n*] and hence

$$
|\|f\|_2 - \sqrt{n/\nu} \bar{d}_n(f,0)| \lesssim \|f - f_n\|_2 \lesssim \frac{\|f\|_{\mathrm{Lip}}}{n} + \left(1 - \frac{2\nu}{n}\right) \|f\|_{\infty}.
$$

It follows that for the verification of [\(3.2\)](#page-5-0), we may always replace *dn* by *d*¯ *<sup>n</sup>* and if the spectral densities are restricted to Lipschitz functions with Lipschitz constant *Ln* and where *εn Ln/n*, then we may also replace *dn* by the *L*2-norm ·2.

Now, by easy calculations, for all spectral densities *f, f*<sup>0</sup> taking values in [*m,M*], we have that *ν*−<sup>1</sup> *<sup>ν</sup> <sup>i</sup>*=<sup>1</sup> *K(Pf*0*,i, Pf,i) <sup>d</sup>*¯<sup>2</sup> *n(f*0*,f) f* − *f*0<sup>2</sup> <sup>∞</sup> and *ν*−<sup>1</sup> *<sup>ν</sup> <sup>i</sup>*=<sup>1</sup> *<sup>V</sup>*2*,*0*(Pf*0*,i, Pf,i) <sup>d</sup>*¯<sup>2</sup> *n(f*0*,f) f* −*f*0<sup>2</sup> <sup>∞</sup>, hence it suffices to estimate the prior probability of sets of the form {*f* : *f* − *f*0∞ ≤ *ε*}. Alternatively, if the spectral densities under consideration are Lipschitz, then we may estimate the prior mass of an *L*2-ball around *f*0.

As a concrete prior, we consider the prior used by Choudhuri, Ghosal and Roy [\[7\]](#page-30-0), namely *f* = *τ q*, where *τ* = var*(Xt)* has a nonsingular prior density and *q*, a probability density on [0*,* 1], is given the Dirichlet–Bernstein prior of Petrone [\[24\]](#page-30-0). We then restrict the prior to the set K = {*f* : *m<f <M*}. The order of the Bernstein polynomial, *k*, has prior mass function *ρ*, which is assumed to satisfy *e*−*β*1*<sup>k</sup>* log *<sup>k</sup> ρ(k) e*−*β*2*k*. Let denote the resulting prior.

Clearly, as *f*<sup>0</sup> ∈ K, restricting the prior to K can only increase the prior probability of {*f* : *f* −*f*0∞ *< ε*}. Therefore, following Ghosal [\[12\]](#page-30-0), *(f* −*f*0∞ *< ε) e*−*cε*−<sup>1</sup> log *<sup>ε</sup>*−<sup>1</sup> . Hence, *εn* of the order *n*−1*/*3*(*log *n)*1*/*<sup>3</sup> satisfies [\(3.4\)](#page-5-0).

Consider a sieve F*<sup>n</sup>* for the parameter space K, which consists solely of Bernstein polynomials of order *kn* or less. All of these functions have Lipschitz constant at most *k*<sup>2</sup> *<sup>n</sup>* and are uniformly bounded away from zero and infinity by construction. The *ε*-entropy of F*<sup>n</sup>* relative to *d*¯ *<sup>n</sup>* can be bounded above by that of the simplex, which is further bounded above by *k* log *k* + *k* log *ε*<sup>−</sup>1. Hence, by choosing *kn* of the order *n*1*/*3*(*log *n)*2*/*3, the convergence rate at *f*<sup>0</sup> on F*<sup>n</sup>* with respect to *dn* is given by max*(n*−1*/*2*k* <sup>1</sup>*/*<sup>2</sup> *<sup>n</sup> (*log *n)*1*/*2*, n*−1*/*3*(*log *n)*1*/*3*, k*<sup>2</sup> *n/n)* <sup>=</sup> *<sup>n</sup>*−1*/*3*(*log *n)*4*/*3. Now, *(*F *<sup>c</sup> <sup>n</sup> )* <sup>=</sup> *ρ(k > kn) <sup>e</sup>*−*β*2*kn* <sup>=</sup> *<sup>e</sup>*−*βn*1*/*3*(*log *n)*2*/*<sup>3</sup> = *e*−*βn(n*−1*/*3*(*log *n)*1*/*3*)*<sup>2</sup> . Thus, the posterior probability of F *<sup>c</sup> <sup>n</sup>* goes to zero by Lemma [1](#page-3-0) and hence the convergence rate on K is also *n*−1*/*3*(*log *n)*1*/*3. The minimax rate *n*−2*/*<sup>5</sup> may be obtained, for instance, by using splines, which have better approximation properties.

7.4. *Nonlinear autoregression.* Consider the nonlinear autoregressive model in which we observe the elements *X*1*,...,Xn* of a stationary time series {*Xt* : *t* ∈ Z} satisfying

(7.4) 
$$
X_i = f(X_{i-1}) + \varepsilon_i, \qquad i = 1, 2, ..., n,
$$

where *f* is an unknown function and *ε*1*, ε*2*,...,εn* are i.i.d. *N (*0*, σ*2*)*. For simplicity, we assume that *σ* = 1. Then *Xn* is a Markov chain with transition density *pf (y*|*x)* = *φ(y* − *f (x))*, where *φ(x)* = *(*2*π )*−1*/*2*e*−*x*2*/*2. Assume that *f* ∈ F , a class of functions such that |*f (x)*| ≤ *M* and |*f (x)* − *f (y)*| ≤ *L*|*x* − *y*| for all *x,y* and *f* ∈ F .

Set *r(y)* = <sup>1</sup> <sup>2</sup> *(φ(y* − *M)* + *φ(y* + *M))*. Then *r(y) pf (y*|*x) r(y)* for all *x,y* ∈ R and *f* ∈ F . Further, sup{ |*p(y*|*x*1*)* − *p(y*|*x*2*)*| *dy* : *x*1*, x*<sup>2</sup> ∈ R} *<* 2. Hence, the chain is *α*-mixing with exponentially decaying mixing coefficients and has a unique stationary distribution *Qf* whose density *qf* satisfies *r qf r*. Let *f <sup>s</sup>* = *(* |*f* | *<sup>s</sup> dr)*1*/s*.

Because *h*2*(N (µ*1*,* 1*), N (µ*2*,* 1*))* = 2[1 − exp*(*−|*µ*<sup>1</sup> − *µ*2| <sup>2</sup>*/*8*)*], it easily follows for *f*1*, f*<sup>2</sup> ∈ F , *d* defined in [\(4.2\)](#page-6-0) and *dν* = *r dλ* that *f*<sup>1</sup> − *f*2<sup>2</sup> *d(f*1*, f*2*) f*<sup>1</sup> − *f*22. Thus, we may verify [\(4.5\)](#page-7-0) relative to the *L*2*(r)*-metric. It can also be computed that

$$
P_{f_0} \log \frac{p_{f_0}(X_2|X_1)}{p_f(X_2|X_1)} = \frac{1}{2} \int (f_0 - f)^2 q_{f_0} d\lambda \lesssim ||f - f_0||_2^2,
$$
  

$$
P_{f_0} \left| \log \frac{p_{f_0}(X_2|X_1)}{p_f(X_2|X_1)} \right|^s \lesssim \int |f_0 - f|^s q_{f_0} d\lambda \lesssim ||f - f_0||_s^s.
$$

Thus, *B*∗*(f*0*, ε*;*s)* ⊃ {*f* : *f* − *f*0*<sup>s</sup>* ≤ *cε*} for some constant *c >* 0, where *B*∗*(f*0*, ε*;*s)* is as in Theorem [5.](#page-7-0) Thus, it suffices to verify [\(4.7\)](#page-7-0) with *s >* 2.

7.4.1. *Random histograms.* As a prior on the regression functions *f* , consider a random histogram as follows. For a given number *K* ∈ N, partition a given compact interval in R into *K* intervals *I*1*,...,IK* and let *I*<sup>0</sup> = R \ *<sup>k</sup> Ik*. Let the prior *n* on *f* be induced by the map *α* → *fα* given by *fα* = *<sup>K</sup> <sup>k</sup>*=<sup>1</sup> *αk*1*Ik* , where the coordinates *α*1*,...,αK* of *α* ∈ R*<sup>K</sup>* are chosen to be i.i.d. random variables with the uniform distribution on the interval [−*M,M*] and where *K* = *Kn* is to be chosen later. Let *r(Ik)* = *Ik r dλ*.

The support of *n* consists of all functions with values in [−*M,M*] that are piecewise constant on each interval *Ik* for *k* = 1*,...,K* and which vanish on *I*0. For any pair *fα* and *fβ* of such functions, we have, for any *s* ∈ [2*,*∞], *fα* − *fβ<sup>s</sup>* = *α* − *βs*, where *α<sup>s</sup>* is the *r*-weighted  *<sup>s</sup>*-norm of *α* = *(α*1*,...,αK)* ∈ R*<sup>K</sup>* given by *α<sup>s</sup> <sup>s</sup>* = *<sup>k</sup>* |*αk*| *sr(Ik)*. The dual use of ·*<sup>s</sup>* should not lead to any confusion as it will be clear from the context whether ·*<sup>s</sup>* is a norm on functions or on vectors. The *L*2*(r)*-projection of *f*<sup>0</sup> onto this support is the function *fα*<sup>0</sup> for *α*0*,k* = *Ik <sup>f</sup>*0*r dλ/r(Ik)*, whence, by Pythagoras' theorem, *fα* <sup>−</sup> *<sup>f</sup>*0<sup>2</sup> <sup>2</sup> = *fα* − *fα*<sup>0</sup> <sup>2</sup> <sup>2</sup> + *fα*<sup>0</sup> <sup>−</sup> *<sup>f</sup>*0<sup>2</sup> <sup>2</sup> for any *α* ∈ [−*M,M*] *<sup>K</sup>*. In particular, *fα* − *f*0<sup>2</sup> ≥ *cα* − *α*0<sup>2</sup> for some constant *c* and hence, with F*<sup>n</sup>* denoting the support of *n*,

$$
N(\varepsilon, \{f \in \mathcal{F}_n : \|f - f_0\|_2 \le 16\varepsilon\}, \|\cdot\|_2)
$$
  
\$\le N(\varepsilon, \{\alpha \in \mathbb{R}^K : \|\alpha - \alpha\_0\|\_2 \le 16c\varepsilon\}, \|\cdot\|\_2) \le (80c)^K\$,

<span id="page-18-0"></span>as in Lemma 4.1 of [\[25\]](#page-31-0). Thus, [\(4.5\)](#page-7-0) holds if *nε*<sup>2</sup> *<sup>n</sup>* - *K*. To verify [\(4.7\)](#page-7-0), note that for *λ* = *(λ(I*1*), . . . , λ(IK))*,

$$
\|f_{\alpha_0} - f_0\|_s^s = \int_{I_0} |f_0|^s d\lambda + \sum_k \int_{I_k} |\alpha_{0,k} - f_0|^s r d\lambda \leq M^s r(I_0) + L^s \|\lambda\|_s^s.
$$

Hence, as *f*<sup>0</sup> ∈ F , for every *α* ∈ [−*M,M*] *K*,

$$
||f_{\alpha}-f_0||_s \lesssim ||\alpha-\alpha_0||_s + r(I_0)^{1/s} + ||\lambda||_s \leq ||\alpha-\alpha_0||_{\infty} + r(I_0)^{1/s} + ||\lambda||_s,
$$

where ·∞ is the ordinary maximum norm on R*K*. For *r(I*0*)*1*/s* + *λ<sup>s</sup>* ≤ *ε/*2, we have that {*f* : *f* − *f*0*<sup>s</sup>* ≤ *ε*}⊃{*fα* : *α* − *α*0∞ ≤ *ε/*2}. Using *α* − *α*0<sup>2</sup> ≤ *cfα* − *f*02, for any *ε >* 0 such that *r(I*0*)*1*/s* + *λ<sup>s</sup>* ≤ *ε/*2, we have

$$
\frac{\Pi_n(f: \|f-f_0\|_2 \leq j\varepsilon)}{\Pi_n(f: \|f-f_0\|_s \leq \varepsilon)} \leq \frac{\Pi_n(\alpha: \|\alpha-\alpha_0\|_2 \leq j\varepsilon)}{\Pi_n(\alpha: \|\alpha-\alpha_0\|_{\infty} \leq \varepsilon c/2)}.
$$

We show that the right-hand side is bounded by *eCnε*2*/*<sup>8</sup> for some *C*.

For *<sup>k</sup> Ik*, a regular partition of an interval [−*A,A*], we have that *λ<sup>s</sup>* = 2*A/K* and since *r(Ik)* ≥ *λ(Ik)*inf*x*∈*Ik r(x)* for every *k* ≥ 1, the norm ·<sup>2</sup> is bounded below by <sup>√</sup>2*Aφ(A)/K* - <sup>√</sup>*φ(A)/K* times a multiple of the Euclidean norm. In this case, the preceding display is bounded above by

$$
\frac{(Cj\varepsilon\sqrt{K/\phi(A)}/(2M))^K \text{ vol}_K}{(\varepsilon c/(4M))^K} \sim \left(\frac{j\sqrt{2\pi e}}{\sqrt{\phi(A)}}\right)^K \frac{1}{\sqrt{\pi K}},
$$

by Stirling's approximation, where vol*<sup>K</sup>* is the volume of the *K*-dimensional Euclidean unit ball. The probability *r(I*0*)* is bounded above by 1−2*(A) φ(A)*. Hence, [\(4.7\)](#page-7-0) will hold if *K* log*(*1*/φ(A)) nε*<sup>2</sup> *<sup>n</sup>*, *φ(A) ε<sup>s</sup> <sup>n</sup>* and *A/K εn*. All requirements are met for *εn* equal to a multiple of *<sup>n</sup>*−1*/*3*(*log *n)*1*/*<sup>2</sup> [with *<sup>K</sup>* <sup>∼</sup> <sup>√</sup>log*(*1*/εn)ε*<sup>−</sup><sup>1</sup> *<sup>n</sup>* and *<sup>A</sup>* <sup>∼</sup> <sup>√</sup>log*(*1*/εn)*]. This is only marginally weaker than the minimax rate, which is *n*−1*/*<sup>3</sup> for this problem, provided the autoregression functions are assumed to be only Lipschitz continuous.

The logarithmic factor in the convergence rate appears to be a consequence of the fact that the regression functions are defined on the full real line. The present prior is a special case of a spline-based prior (see, e.g., Section [7.7\)](#page-22-0). If *f* has smoothness beyond Lipschitz continuity, then the use of higher order splines should yield a faster convergence rate.

7.5. *Finite-dimensional i.n.i.d. models.* Theorem [4](#page-5-0) is also applicable to finitedimensional models and yields the usual convergence rate as shown below. The result may be compared with Theorem I.10.2 of [\[19\]](#page-30-0) and Proposition 1 of [\[13\]](#page-30-0).

THEOREM 10. *Let X*1*,...,Xn be i*.*n*.*i*.*d*. *observations following densities pθ ,i*, *where -* ⊂ R*<sup>d</sup>* . *Let θ*<sup>0</sup> *be an interior point of -*. *Assume that there exist* *constants α >* 0 *and* 0 ≤ *ci* ≤ *Ci <* ∞ *with*, *for every θ,θ*1*, θ*<sup>2</sup> ∈ *-*,

(7.5) 
$$
c = \liminf_{n \to \infty} \frac{1}{n} \sum_{i=1}^{n} c_i > 0, \qquad C = \limsup_{n \to \infty} \frac{1}{n} \sum_{i=1}^{n} C_i < \infty
$$

*such that Pθ*0*,i(*log *pθ*0*,i pθ ,i )* <sup>≤</sup> *Ci<sup>θ</sup>* <sup>−</sup> *<sup>θ</sup>*02*α*, *Pθ*0*,i(*log *pθ*0*,i pθ ,i )*<sup>2</sup> <sup>≤</sup> *Ci<sup>θ</sup>* <sup>−</sup> *<sup>θ</sup>*02*<sup>α</sup> and*

$$
(7.6) \t c_i \|\theta_1 - \theta_2\|^{2\alpha} \leq h^2(p_{\theta_1, i}, p_{\theta_2, i}) \leq C_i \|\theta_1 - \theta_2\|^{2\alpha}.
$$

*Assume that the prior measure possesses a density π which is bounded away from zero in a neighborhood of θ*<sup>0</sup> *and bounded above on the entire parameter space*. *Then the posterior converges at the rate n*−1*/(*2*α) with respect to the Euclidean metric*.

For regular families, the above displays are satisfied for *α* = 1 and the usual *n*−1*/*<sup>2</sup> rate is obtained; see [\[19\]](#page-30-0), Chapter III. Nonregular cases, for instance, when the densities have discontinuities depending on the parameter [such as the uniform distribution on *(*0*,θ)*], have *α <* 1 and faster rates are obtained; see [\[19\]](#page-30-0), Chapters V and VI and [\[13\]](#page-30-0).

PROOF OF THEOREM [10.](#page-18-0) By the assumptions (7.5) and (7.6), it suffices to show that the posterior convergence rate with respect to *dn* defined by [\(3.1\)](#page-4-0) is *n*−1*/*2. Now, by Pollard ([\[25\]](#page-31-0), Lemma 4.1),

$$
N(\varepsilon/18, \{\theta \in \Theta : d_n(\theta, \theta_0) < \varepsilon\}, d_n)
$$
\n
$$
\leq N((\varepsilon^2/(36C))^{1/(2\alpha)}, \{\theta \in \Theta : \|\theta - \theta_0\| < (2\varepsilon^2/c)^{1/(2\alpha)}\}, \|\cdot\|)
$$
\n
$$
\leq 6^d \left(\frac{72C}{c}\right)^{d/(2\alpha)},
$$

which verifies [\(3.2\)](#page-5-0). For [\(3.4\)](#page-5-0), note that

$$
\Pi(\theta : d_n(p_\theta, p_{\theta_0}) \le j\varepsilon)
$$
\n
$$
\Pi(\theta : n^{-1} \sum_{i=1}^n K_i(\theta_0, \theta) \le \varepsilon^2, n^{-1} \sum_{i=1}^n V_{2;i}(\theta_0, \theta) \le \varepsilon^2)
$$
\n
$$
\le \frac{\Pi(\theta : \|\theta - \theta_0\| \le (2j^2\varepsilon^2/c)^{1/(2\alpha)})}{\Pi(\theta : \|\theta - \theta_0\| \le (\varepsilon^2/(2C))^{1/(2\alpha)})} \le A j^{d/\alpha}
$$

for sufficiently small *ε >* 0, where *A* is a constant depending on *d*, *c*, *C* and the upper and lower bounds on the prior density. The conclusion follows for *εn* = *M/*√*n*, where *<sup>M</sup>* is a large constant.

The condition that the Hellinger distance is bounded below by a power of the Euclidean distance excludes the possibility of unbounded parameter spaces. This defect may be rectified by applying Theorem [3](#page-4-0) to derive the rate. If there is a uniformly exponentially consistent test for *θ* = *θ*<sup>0</sup> against the complement of a bounded set, then the result holds even if  is not bounded. Often, such tests exist by virtue of bounds on log affinity, as in the case of normal distributions, or by large deviation type inequalities; see [\[20\]](#page-30-0) and [\[14\]](#page-30-0), Section 7. Further, if the prior density is not bounded above, but has a polynomial or subexponential majorant, then the rate calculation also remains valid.

7.6. *White noise with conjugate priors.* In this section, we consider the white noise model of Section [5](#page-8-0) with a conjugate Gaussian prior. This allows us to complement and rederive results of Zhao [\[34\]](#page-31-0) and Shen and Wasserman [\[27\]](#page-31-0) in our framework. Thus, we observe an infinite sequence *X*1*, X*2*,...* of independent random variables, where *Xi* is normally distributed with mean *θi* and variance *n*<sup>−</sup>1.

We consider the prior *n* on the parameter *θ* = *(θ*1*, θ*2*, . . .)* that can be structurally described by saying that *θ*1*,...,θk* are independent with *θi* normally distributed with mean zero and variance *σ*<sup>2</sup> *i,k* and that *θk*<sup>+</sup>1*, θk*<sup>+</sup>2*,...* are set equal to zero. Here, we choose the cutoff *k* dependent on *n* and equal to *k* = *n*1*/(*2*α*+1*)* for some *α >* 0. Zhao [\[34\]](#page-31-0) and Shen and Wasserman [\[27\]](#page-31-0) consider the case where *σ*<sup>2</sup> *i,k* <sup>=</sup> *<sup>i</sup>*−*(*2*α*+1*)* for *<sup>i</sup>* <sup>=</sup> <sup>1</sup>*,...,k* and show that the convergence rate is *εn* <sup>=</sup> *n*−*α/(*2*α*+1*)* if the true parameter *θ*<sup>0</sup> is "*α*-regular" in the sense that <sup>∞</sup> *<sup>i</sup>*=<sup>1</sup> *θ* <sup>2</sup> <sup>0</sup>*,ii*2*<sup>α</sup> <sup>&</sup>lt;* ∞. We shall obtain the same result for any triangular array of variances such that

(7.8) 
$$
\min \{ \sigma_{i,k}^2 i^{2\alpha} : 1 \le i \le k \} \sim k^{-1}.
$$

For instance, for each *k*, the coefficients *θ*1*,...,θk* could be chosen i.i.d. normal with mean zero and variance *k*−<sup>1</sup> or could follow the model of the authors mentioned previously.

THEOREM 11. *If k* ∼ *n*1*/(*2*α*+1*) and* (7.8) *holds*, *then the posterior converges at the rate εn* = *n*−*α/(*2*α*+1*) for any θ*<sup>0</sup> *such that* <sup>∞</sup> *<sup>i</sup>*=<sup>1</sup> *θ* <sup>2</sup> <sup>0</sup>*,ii*2*<sup>α</sup> <sup>&</sup>lt;* <sup>∞</sup>.

PROOF. The support *<sup>n</sup>* of the prior is the set of all *θ* ∈ <sup>2</sup> with *θi* = 0 for *i>k* and can be identified with R*k*. Moreover, the 2-norm · on the support can be identified with the Euclidean norm ·*<sup>k</sup>* on R*k*. Let *Bk(x, ε)* denote the *k*dimensional Euclidean ball of radius *ε* and center *x* ∈ R*<sup>d</sup>* . For any true parameter *θ*<sup>0</sup> ∈ 2, we have *θ* − *θ*0≥**P***θ* − **P***θ*0*k*, where **P** is the projection on *<sup>n</sup>*, and hence

$$
N(\varepsilon/8, \{\theta \in \Theta_n : \|\theta - \theta_0\| \le \varepsilon\}, \|\cdot\|) \le N(\varepsilon/8, B_k(\mathbf{P}\theta_0, \varepsilon), \|\cdot\|_k) \le (40)^k.
$$

It follows that [\(5.1\)](#page-9-0) is satisfied for *nε*<sup>2</sup> *<sup>n</sup> k*, that is, in view of our choice of *k*, *εn n*−*α/(*2*α*+1*)* .

By Pythagoras' theorem, we have that *θ* − *θ*0<sup>2</sup> = **P***θ* − **P***θ*0<sup>2</sup> + *i>k θ* <sup>2</sup> 0*,i* for any *θ* in the support of *n*. Hence, for *i>k θ* <sup>2</sup> <sup>0</sup>*,i* <sup>≤</sup> *<sup>ε</sup>*<sup>2</sup> *n/*2, we have that

$$
\Pi_n(\theta \in \Theta_n : \|\theta - \theta_0\| \le \varepsilon_n) \ge \Pi_n(\theta \in \mathbb{R}^k : \|\theta - \mathbf{P}\theta_0\|_k \le \varepsilon_n/2).
$$

<span id="page-20-0"></span>

By the definition of the prior, the right-hand side involves a quadratic form in Gaussian variables. For the *k* × *k* diagonal matrix with elements *σ*<sup>2</sup> *i,k*, the quotient on the left-hand side of [\(5.2\)](#page-9-0) can be bounded as

$$
\frac{\Pi_n(\theta \in \Theta_n : \|\theta - \theta_0\| \leq j\varepsilon_n)}{\Pi_n(\theta \in \Theta_n : \|\theta - \theta_0\| \leq \varepsilon_n)} \leq \frac{N_k(-\mathbf{P}\theta_0, \Sigma)(B(0, j\varepsilon_n))}{N_k(-\mathbf{P}\theta_0, \Sigma)(B(0, \varepsilon_n/2))}.
$$

The probability in the numerator increases if we center the normal distribution at 0 rather than at −**P***θ*0, by Anderson's lemma. Furthermore, for any *µ* ∈ R*k*,

$$
\frac{dN_k(\mu,\Sigma)}{dN_k(0,\Sigma/2)}(\theta)=\frac{e^{-\sum_{i=1}^k(\theta_i-\mu_i)^2/(2\sigma_{i,k}^2)}}{\sqrt{2}^k e^{-\sum_{i=1}^k\theta_i^2/\sigma_{i,k}^2}}\geq 2^{-k/2}e^{-\sum_{i=1}^k\mu_i^2/\sigma_{i,k}^2}.
$$

Therefore, we may recenter the denominator at 0 at the cost of adding the factor on the right (with *µ* = *θ*0) and dividing the covariance matrix by 2. We obtain that the left-hand side of [\(5.2\)](#page-9-0) is bounded above by

$$
2^{k/2} e^{\sum_{i=1}^k \theta_{0,i}^2/\sigma_{i,k}^2} \frac{N_k(0, \Sigma)(B(0, j\varepsilon_n))}{N_k(0, \Sigma/2)(B(0, \varepsilon_n/2))}
$$
  

$$
\leq 2^{k/2} e^{\sum_{i=1}^k \theta_{0,i}^2/\sigma_{i,k}^2} \left(\frac{\bar{\sigma}_k}{\underline{\sigma}_k}\right)^k \frac{N_k(0, \bar{\sigma}_k^2 I)(B(0, j\varepsilon_n))}{N_k(0, \underline{\sigma}_k^2 I/2)(B(0, \varepsilon_n/2))},
$$

where *σ*¯*<sup>k</sup>* and *σ <sup>k</sup>* denote the maximum and the minimum of *σi,k* for *i* = 1*,* 2*,...,k*. The probabilities on the right-hand side are left tail probabilities of chi-square distributions with *k* degrees of freedom, and can be expressed as integrals. The preceding display is bounded above by

$$
2^{k/2} e^{\sum_{i=1}^k \theta_{0,i}^2/\sigma_{i,k}^2} \left(\frac{\bar{\sigma}_k}{\underline{\sigma}_k}\right)^k \frac{\int_0^{j^2 \epsilon_n^2/\bar{\sigma}_k^2} x^{k/2-1} e^{-x/2} dx}{\int_0^{\epsilon_n^2/(2\underline{\sigma}_k^2)} x^{k/2-1} e^{-x/2} dx}.
$$

The exponential in the integral in the numerator is bounded above by 1 and hence this integral is bounded above by *j kεk n/(kσ*¯ *<sup>k</sup> <sup>k</sup> )*. We now consider two separate cases. If *ε*<sup>2</sup> *n/σ*<sup>2</sup> *<sup>k</sup>* remains bounded, then we can also bound the exponential in the integral in the denominator below by a constant and have that the preceding display is bounded above by a multiple of 4*kj <sup>k</sup>* exp*( <sup>k</sup> <sup>i</sup>*=<sup>1</sup> *θ* <sup>2</sup> 0*,i/σ*<sup>2</sup> *i,k)*. If *ε*<sup>2</sup> *n/σ*<sup>2</sup> *<sup>k</sup>* → ∞, then we bound the integral in the denominator below by *(η/*2*)k/*2−<sup>1</sup> *<sup>η</sup> η/*<sup>2</sup> *<sup>e</sup>*−*x/*<sup>2</sup> *dx* for *<sup>η</sup>* <sup>=</sup> *<sup>ε</sup>*<sup>2</sup> *n/(*2*σ*<sup>2</sup> *k)*. This leads to the upper bound being a multiple of 8*kj <sup>k</sup>* exp*( <sup>k</sup> <sup>i</sup>*=<sup>1</sup> *θ* <sup>2</sup> <sup>0</sup>*,i<sup>σ</sup>* <sup>−</sup><sup>2</sup> *i,k )ε*<sup>2</sup> *nσ* <sup>−</sup><sup>2</sup> *<sup>k</sup>* exp*(ε*<sup>2</sup> *nσ* <sup>−</sup><sup>2</sup> *<sup>k</sup> /*8*)*. By the assumption [\(7.8\)](#page-20-0), we have that *σ*<sup>2</sup> *<sup>k</sup> k*−*(*2*α*+1*)* ∼ *n*<sup>−</sup>1. We also have that *k* ∼ *nε*<sup>2</sup> *<sup>n</sup>*. It follows that *ε*2 *n/σ*<sup>2</sup> *<sup>k</sup> nε*<sup>2</sup> *<sup>n</sup>* and that *<sup>σ</sup>* <sup>−</sup><sup>2</sup> *<sup>k</sup>* is bounded by a polynomial in *k*. We conclude that with our choice of *k* ∼ *n*1*/(*2*α*+1*)* , [\(5.2\)](#page-9-0) is satisfied if *εn* satisfies *<sup>k</sup> <sup>i</sup>*=<sup>1</sup> *θ* <sup>2</sup> 0*,i/σ*<sup>2</sup> *i,k nε*<sup>2</sup> *n* and *i>k θ* <sup>2</sup> <sup>0</sup>*,i* <sup>≤</sup> *<sup>ε</sup>*<sup>2</sup> *n/*2.

It follows that the posterior concentrates at *θ*<sup>0</sup> at the rate *εn* that satisfies these requirements as well as the condition *εn n*−*α/(*2*α*+1*)* . If the true parameter *θ*<sup>0</sup> <span id="page-22-0"></span>satisfies <sup>∞</sup> *<sup>i</sup>*=<sup>1</sup> *θ* <sup>2</sup> <sup>0</sup>*,ii*2*<sup>α</sup> <sup>&</sup>lt;* <sup>∞</sup>, then all three inequalities are satisfied for *εn* a multiple of *n*−*α/(*2*α*+1*)* . The rate *n*−*α/(*2*α*+1*)* is the minimax rate for this problem.

Our prior is dependent on *n*, but with some more effort, it can be seen that the same conclusion can be obtained with a mixture prior of the form *<sup>n</sup> λnn* for suitable *λn*.

7.7. *Nonparametric regression with Gaussian errors.* Consider the nonparametric regression model, where we observe independent random variables *X*1*,...,Xn* distributed as *Xi* = *f (zi)* + *εi* for an unknown regression function *f* , deterministic real-valued covariates *z*1*,...,zn* and normally distributed error variables *ε*1*,...,εn* with zero means and variances *σ*2. For simplicity, we assume that the error variance *σ*<sup>2</sup> is known. We also suppose that the covariates take values in a fixed compact set, which we will take as the unit interval, without loss of generality.

Let *f*<sup>0</sup> denote the true value of the regression function, let *Pf,i* be the distribution of *Xi* and let *P(n) <sup>f</sup>* be the distribution of *(X*1*,...,Xn)*. Thus, *Pf,i* is the normal measure with mean *f (zi)* and variance *σ*2. Let P*<sup>z</sup> <sup>n</sup>* <sup>=</sup> *<sup>n</sup>*−<sup>1</sup> *<sup>n</sup> <sup>i</sup>*=<sup>1</sup> *δzi* be the empirical measure of the covariates and let ·*<sup>n</sup>* denote the norm on *L*2*(*P*<sup>z</sup> n)*.

By easy calculations, *K(Pf*0*,i, Pf,i)* = |*f*0*(zi)* − *f (zi)*| <sup>2</sup>*/(*2*σ*2*)* and *V*2*,*0*(Pf*0*,i, Pf,i)* = |*f*0*(zi)*−*f (zi)*| <sup>2</sup>*/σ*<sup>2</sup> for all *i* = 1*,* 2*,...,n*, whence the average Kullback– Leibler divergence and variance are bounded by a multiple of *f*<sup>0</sup> − *f* <sup>2</sup> *n/σ*<sup>2</sup> and hence it is enough to quantify prior concentration in ·*n*-balls. The average Hellinger distance, as used in Theorem [4,](#page-5-0) is bounded above by ·*n*, but is equivalent to this norm only if the class of regression functions is uniformly bounded, which makes it less attractive. However, it can be verified (cf. [\[5\]](#page-30-0)) that the likelihood ratio test for *f*<sup>0</sup> versus *f*<sup>1</sup> satisfies the conclusion of Lemma [2](#page-4-0) relative to ·*<sup>n</sup>* (instead of *dn* and *θi* = *fi*). Therefore, we may use the norm ·*<sup>n</sup>* instead of the average Hellinger distance throughout.

We shall construct priors based on series representations that are appropriate if *f*<sup>0</sup> ∈ *Cα*[0*,* 1], where *α >* 0 could be fractional. This means that *f*<sup>0</sup> is *α*<sup>0</sup> times continuously differentiable with *f*0*<sup>α</sup> <* ∞, *α*<sup>0</sup> being the greatest integer less than *α* and the seminorm being defined by

(7.9) 
$$
||f||_{\alpha} = \sup_{x \neq x'} \frac{|f^{(\alpha_0)}(x) - f^{(\alpha_0)}(x')|}{|x - x'|^{\alpha - \alpha_0}}.
$$

7.7.1. *Splines.* Fix an integer *q* with *q* ≥ *α*. For a given natural number *K*, which will increase with *n*, partition the interval *(*0*,* 1] into *K* subintervals *((k* − 1*)/K, k/K*] for *k* = 1*,* 2*,...,K*. The space of splines of order *q* relative to this partition is the collection of all functions *f* : *(*0*,* 1] → R that are *q* − 2 times continuously differentiable throughout *(*0*,* 1] and, if restricted to a subinterval *((k* − 1*)/K, k/K*], are polynomials of degree strictly less than *q*. These

<span id="page-23-0"></span>splines form a *J* = *(q* + *K* − 1*)*-dimensional linear space, with a convenient basis *B*1*, B*2*,...,BJ* being the B-splines, as defined in, for example, [\[11\]](#page-30-0). The B-splines satisfy (i) *Bj* ≥ 0, *j* = 1*,* 2*,...,J* , (ii) *<sup>J</sup> <sup>j</sup>*=<sup>1</sup> *Bj* = 1, (iii) *Bj* is supported inside an interval of length *q/K* and (iv) at most *q* of *B*1*(x), . . . , BJ (x)* are nonzero at any given *x*. Let *B(z)* = *(B*1*(z), . . . , BJ (z))<sup>T</sup>* and write *β<sup>T</sup> B* for the function *z* → *<sup>j</sup> βjBj (z)*.

The basic approximation property of splines proved in [\[11\]](#page-30-0), page 170, shows that for some *β*<sup>∞</sup> ∈ R*<sup>J</sup>* (dependent on *J* ),

$$
(7.10) \t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t
$$

Thus, by increasing *J* appropriately with the sample size, we may view the space of splines as a sieve for the construction of the maximum likelihood estimator, as in Stone [\[28, 29\]](#page-31-0), and for Bayes estimates as in [\[14, 15\]](#page-30-0) for the problem of density estimation.

To put a prior on *f* , we represent it as *fβ(z)* = *β<sup>T</sup> B(z)* and induce a prior on *f* from a prior on *β*. Ghosal, Ghosh and van der Vaart [\[14\]](#page-30-0), in the context of density estimation, choose *β*1*,...,βJ* i.i.d. uniform on an interval [−*M,M*], the restriction to a finite interval being necessary to avoid densities with arbitrarily small values. In the present regression situation, a restriction to a compact interval is unnecessary and we shall choose *β*1*,...,βJ* to be a sample from the standard normal distribution.

We need the regressors *z*1*, z*2*,...,zn* to be sufficiently regularly distributed in the interval [0*,* 1]. In view of the spatial separation property of the B-spline functions, the precise condition can be expressed in terms of the covariance matrix *n* = *( BiBj d*P*<sup>z</sup> n)*, namely

$$
(7.11) \t\t\t J^{-1} \|\beta\|^2 \lesssim \beta^T \Sigma_n \beta \lesssim J^{-1} \|\beta\|^2,
$$

where · is the Euclidean norm on R*<sup>J</sup>* .

Under condition (7.11), we have that for all *β*1*, β*<sup>2</sup> ∈ R*<sup>J</sup>* ,

(7.12) 
$$
C \|\beta_1 - \beta_2\| \le \sqrt{J} \|f_{\beta_1} - f_{\beta_2}\|_n \le C' \|\beta_1 - \beta_2\|
$$

for some constants *C* and *C* . This enables us to perform all calculations in terms of the Euclidean norms on the spline coefficients.

THEOREM 12. *Assume that the true density f*<sup>0</sup> *satisfies* (7.10) *for some α* ≥ <sup>1</sup> 2 , *let* (7.11) *hold and let n be priors induced by a NJ (*0*,I) distribution on the spline coefficients*. *If J* = *Jn* ∼ *n*1*/(*1+2*α)*, *then the posterior converges at the minimax rate n*−*α/(*1+2*α) relative to* ·*n*.

PROOF. We verify the conditions of Theorem [4.](#page-5-0) Let *fβn* be the *L*2*(*P*<sup>z</sup> n)* projection of *f*<sup>0</sup> onto the *J* -dimensional space of splines *fβ* = *β<sup>T</sup> B*. Then *fβn* −*fβ<sup>n</sup>* ≤ *f*<sup>0</sup> −*fβ<sup>n</sup>* for every *β* ∈ R*<sup>J</sup>* and hence, by (7.12), for every *ε >* 0, we have {*β* : *fβ* − *f*0*<sup>n</sup>* ≤ *ε*}⊂{*β* : *β* − *βn* ≤ *C* <sup>√</sup>*J ε*}. It follows that the *<sup>ε</sup>*covering numbers of the set {*fβ* : *fβ* − *f*0*<sup>n</sup>* ≤ *ε*} for ·*<sup>n</sup>* are bounded by the *C* <sup>√</sup>*J ε*-covering numbers of a Euclidean ball of radius *<sup>C</sup>* <sup>√</sup>*J ε*, which are of the order *D<sup>J</sup>* for some constant *D*. Thus, the entropy condition [\(3.2\)](#page-5-0) is satisfied, provided that *J nε*<sup>2</sup> *n*.

By the projection property, with *β*∞ as in [\(7.10\)](#page-23-0),

$$
(7.13) \t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t
$$

Combining this with [\(7.12\)](#page-23-0) shows that there exists a constant *C* such that for every *ε* - <sup>2</sup>*<sup>J</sup>* <sup>−</sup>*α*, {*<sup>β</sup>* : *fβ* <sup>−</sup> *<sup>f</sup>*0*<sup>n</sup>* <sup>≤</sup> *<sup>ε</sup>*}⊃{*<sup>β</sup>* : *<sup>β</sup>* <sup>−</sup> *βn* ≤ *<sup>C</sup>*<sup>√</sup>*J ε*}. Together with the inclusion in the preceding paragraph and the definition of the prior, this implies that

$$
\frac{\Pi_n(f: ||f - f_0||_n \le j\varepsilon)}{\Pi_n(f: ||f - f_0||_n \le \varepsilon)} \le \frac{N_J(0, I)(\beta : ||\beta - \beta_n|| \le C'j\sqrt{J}\varepsilon)}{N_J(0, I)(\beta : ||\beta - \beta_n|| \le C''j\sqrt{J}\varepsilon)}\n\le \frac{N_J(0, I)(\beta : ||\beta|| \le C'j\sqrt{J}\varepsilon)}{2^{-J/2}e^{-||\beta_n||^2}N_J(0, I)(\beta : ||\beta|| \le C''j\sqrt{J}\varepsilon/\sqrt{2})}.
$$

In the last step, we use Anderson's lemma to see that the numerator increases if we replace the centering *βn* by the origin, whereas to bound the denominator below, we use the fact that

$$
\frac{dN_J(\beta_n, I)}{dN_J(0, I/2)}(\beta) = \frac{e^{-\|\beta - \beta_n\|^2/2}}{(\sqrt{2})^J e^{-\|\beta\|^2}} \ge 2^{-J/2} e^{-\|\beta_n\|^2}
$$

*.*

√ Here, by the triangle inequality, [\(7.12\)](#page-23-0) and (7.13), we have that *βn <sup>J</sup> fβn <sup>n</sup>* <sup>√</sup>*J (J* <sup>−</sup>*<sup>α</sup>* + *f*0∞*)* <sup>√</sup>*<sup>J</sup>* . Furthermore, the two Gaussian probabilities are left tail probabilities of the chi-square distribution with *J* degrees of freedom. The quotient can be evaluated as

$$
2^{J/2}e^{\|\beta_n\|^2}\frac{\int_0^{(C')^2j^2J\varepsilon^2}x^{J/2-1}e^{-x/2}dx}{\int_0^{(C'')^2J\varepsilon^2/2}x^{J/2-1}e^{-x/2}dx}.
$$

This is bounded above by *(Cj )J* for some constant *<sup>C</sup>* if <sup>√</sup>*J ε* remains bounded. Hence, to satisfy [\(3.4\)](#page-5-0), it again suffices that *nε*<sup>2</sup> *<sup>n</sup>* -*J* .

We conclude the proof by choosing *J* = *Jn* ∼ *n*1*/(*1+2*α)*.

7.7.2. *Orthonormal series priors.* The arguments in the preceding subsection use the special nature of the B-spline basis only through the approximation inequality [\(7.10\)](#page-23-0) and the comparison of norms [\(7.12\)](#page-23-0). Theorem [12](#page-23-0) thus extends to many other possible bases. One possibility is to use a sequence of orthonormal bases with good approximation properties for a given class of regression functions *f*0. Then [\(7.11\)](#page-23-0) should be replaced by

$$
(7.14) \t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t
$$

This is trivially true if the bases are orthonormal in *L*2*(*P*<sup>z</sup> n)*, but this requires that the basis functions change with the design points *z*1*,...,zn*. One possible example is the discrete wavelet bases relative to the design points. All arguments remain valid in this setting.

7.8. *Binary regression.* Let *X*1*,...,Xn* be independent observations with P*(Xi* = 1*)* = 1 − P*(Xi* = 0*)* = *F (α* + *βzi)*, where *zi* is a one-dimensional covariate, *α* and *β* are parameters and *F* is a cumulative distribution. Within the parametric framework, logit regression, where *F (z)* = *(*1 + *e*<sup>−</sup>*z)*<sup>−</sup>1, or probit regression, where *F* is the cumulative distribution function of the standard normal distribution, are usually considered. Recently, there has been interest in link functions of unknown functional form. The parameters *(F , α, β)* are separately not identifiable, unless some suitable restrictions on *F* (such as given values of two quantiles of *F*) are imposed. For Bayesian estimation of *(F , α, β)*, one therefore needs to put a prior on *F* that conforms with the given restriction. However, in practice, one usually puts a Dirichlet process or a similar prior on *F* and, independently of this, a prior on *(α, β)*, and makes inference about, say, *z*0, where *F (α* + *βz*0*)* = 1*/*2. Recently, Amewou-Atisso et al. [\[1\]](#page-29-0) showed that the resulting posterior is consistent. In this section, we obtain the rate of convergence by an application of Theorem [4.](#page-5-0)

Because we directly measure distances between the distributions generating the data, identifiability issues need not concern us. The model and the prior can thus be described in a simpler form. We assume that *X*1*, X*2*,...* are independent Bernoulli variables, *Xi* having success parameter *H (zi)* for an unknown, monotone link function *H*. As a prior on *H*, we use the Dirichlet process prior with base measure *γ ((*· − *α)/β)*, for "hyperparameters" *(α, β)* distributed according to some given prior. This results in a mixture of Dirichlet process priors for *H*. Let the true value of *H* be *H*0, which is assumed to be continuous and nondecreasing.

In practice, *γ* is often chosen to have support equal to the whole of R and *(α, β)* chosen to have support equal to R × *(*0*,*∞*)* so that the conditions on *γ* and *(α, β)* described in the following theorem are satisfied.

THEOREM 13. *Assume that z*1*, z*2*,...,zn lie in an interval* [*a,b*] *strictly within the support of the true link function H*<sup>0</sup> *so that H*0*(a*−*) >* 0 *and H*0*(b) <* 1. *Let H be the given mixture of Dirichlet process priors described previously with γ and (α, β) having densities that are positive and continuous inside their supports*. *Assume that there exists a compact set* K *inside the support of the prior for (α, β) such that whenever (α, β)* ∈ K, *the support of the base measure γ ((*· − *α)/β) strictly contains the interval* [*a,b*]. *Then the posterior distribution of H converges at the rate n*−1*/*3*(*log *n)*1*/*<sup>3</sup> *with respect to the distance dn given by* [\(3.1\)](#page-4-0).

PROOF. Because the Hellinger distance between two Bernoulli distributions with success parameters *p* and *q* is equal to *(p*1*/*<sup>2</sup> − *q*1*/*2*)*<sup>2</sup> + *((*1 − *p)*1*/*<sup>2</sup> − *(*1 − *q)*1*/*2*)*2, we have

$$
d_n^2(H_1, H_2) \le \int |H_1^{1/2} - H_2^{1/2}|^2 d\mathbb{P}_n + \int |(1 - H_1)^{1/2} - (1 - H_2)^{1/2}|^2 d\mathbb{P}_n,
$$

where P*<sup>n</sup>* is the empirical distribution of *z*1*, z*2*,...,zn*. Both the classes {*H*1*/*<sup>2</sup> : *H* is a c.d.f.} and {*(*1−*H )*1*/*<sup>2</sup> : *H* is a c.d.f.} have *ε*-entropy bounded by a multiple of *ε*<sup>−</sup>1, by Theorem 2.7.5 of [\[31\]](#page-31-0). Thus, any *εn n*−1*/*<sup>3</sup> satisfies [\(3.2\)](#page-5-0).

By easy calculations, we have

$$
K_i(H_0, H) = H_0(z_i) \log \frac{H_0(z_i)}{H(z_i)} + (1 - H_0(z_i)) \log \frac{1 - H_0(z_i)}{1 - H(z_i)},
$$

$$
V_{2;i}(H_0, H) \le 2H_0(z_i) \left( \log \frac{H_0(z_i)}{H(z_i)} \right)^2 + 2\left(1 - H_0(z_i)\right) \left( \log \frac{1 - H_0(z_i)}{1 - H(z_i)} \right)^2.
$$

Under the conditions of the theorem, the numbers *H*0*(zi)* are bounded away from 0 and 1. By Taylor's expansion, for any *δ >* 0, there exists a constant *C* (depending on *δ*) such that

$$
\sup_{\delta < p < 1 - \delta} \sup_{q : |q - p| < \varepsilon} \left( p \left( \log \frac{p}{q} \right)^r + (1 - p) \left( \log \frac{1 - p}{1 - q} \right)^r \right) \le C\varepsilon^2, \qquad r = 1, 2.
$$

Therefore, with *H* − *H*0∞ = sup{|*H (z)* − *H*0*(z)*| : *z* ∈ [*a,b*]}, we have max*(n*−<sup>1</sup> *<sup>n</sup> <sup>i</sup>*=<sup>1</sup> *Ki(H*0*, H ), n*−<sup>1</sup> *<sup>n</sup> <sup>i</sup>*=<sup>1</sup> *V*2;*i(H*0*, H )) H* − *H*0<sup>2</sup> ∞. Hence, in order to satisfy [\(3.4\)](#page-5-0), it suffices to lower bound the prior probability of the set {*H* : *H* − *H*0∞ ≤ *ε*}.

For given *α* and *β*, the base measure is *γ ((*·−*α)/β)*. For a given *ε >* 0, partition the line into *N ε*−<sup>1</sup> intervals *E*1*, E*2*,...,EN* such that *H*0*(Ej )* ≤ *ε* and such that the *γ ((*· − *α)/β)*-probability of every set *Ej* (for *j* = 1*,* 2*,...,N*) is between *Aε* and 1 for a given positive constant *A*. Existence of such a partition follows from the continuity of *H*0. It easily follows that for every *H* such that *<sup>N</sup> <sup>j</sup>*=<sup>1</sup> |*H (Ej )* − *H*0*(Ej )*| ≤ *ε*, we have *H* − *H*0∞ *ε*. Furthermore, the conclusion is true even if *(α, β)* varies over K. By Lemma 6.1 of [\[14\]](#page-30-0), the prior probability of the set of all *H* satisfying *<sup>N</sup> <sup>j</sup>*=<sup>1</sup> <sup>|</sup>*H (Ej )* <sup>−</sup> *<sup>H</sup>*0*(Ej )*| ≤ *<sup>ε</sup>* is at least exp*(*−*cε*−<sup>1</sup> log *<sup>ε</sup>*−1*)* for some constant *c*. Furthermore, a uniform estimate works for all *(α, β)* ∈ K. Hence, [\(3.4\)](#page-5-0) holds for *εn*, the solution of *nε*<sup>2</sup> = *ε*−<sup>1</sup> log *ε*<sup>−</sup>1, or for *εn* = *n*−1*/*3*(*log *n)*1*/*3, which is only slightly weaker than the minimax rate *n*−1*/*3.

7.9. *Interval censoring.* Let *T*1*, T*2*,...,Tn* constitute an i.i.d. sample from a life distribution *F* on *(*0*,*∞*)*, which is subject to interval censoring by intervals *(l*1*, u*1*), . . . , (ln, un)*. We assume that the intervals are either nonstochastic or else we work conditionally on the realized values. Putting *(δ*1*, η*1*), . . . , (δn, ηn)*, where *δi* = 1{*Ti* ≤ *li*} and *ηi* = 1{*li < Ti < ui*}, *i* = 1*,* 2*,...,n*, the likelihood is given by *<sup>n</sup> <sup>i</sup>*=1*(F (li))δi(F (ui)* <sup>−</sup> *F (li))ηi(*<sup>1</sup> <sup>−</sup> *F (ui))*<sup>1</sup>−*δi*−*ηi* . We may put the Dirichlet

<span id="page-27-0"></span>process prior on *F*. Under mild assumptions on the true *F*<sup>0</sup> and the base measure, the convergence rate under *dn* turns out to be *n*−1*/*3*(*log *n)*1*/*3, which is the minimax rate, except for the logarithmic factor. Here, we use monotonicity of *F* to bound the *ε*-entropy by a multiple of *ε*−<sup>1</sup> and we estimate prior probability concentration as exp*(*−*cε*−<sup>1</sup> log *ε*−1*)* using methods similar to those used in the previous subsection. The details are omitted.

**8. Proofs.** In this section, we collect a number of technical proofs. For the proofs of the main results, we first present two lemmas.

LEMMA 9. *Let dn and en be semimetrics on for which tests satisfying the conditions of* [\(2.2\)](#page-2-0) *exist*. *Suppose that for some nonincreasing function ε* → *N (ε) and some εn* ≥ 0,

$$
(8.1) \qquad N\bigg(\frac{\varepsilon\xi}{2}, \{\theta \in \Theta : d_n(\theta, \theta_0) < \varepsilon\}, e_n\bigg) \le N(\varepsilon) \qquad \text{for all } \varepsilon > \varepsilon_n.
$$

*Then for every ε>εn*, *there exist tests φn*, *n* ≥ 1, (*depending on ε*) *such that P(n) <sup>θ</sup>*<sup>0</sup> *φn* <sup>≤</sup> *N (ε) <sup>e</sup>*−*Knε*<sup>2</sup> <sup>1</sup>−*e*−*Knε*<sup>2</sup> *and <sup>P</sup>(n) <sup>θ</sup> (*<sup>1</sup> <sup>−</sup> *φn)* <sup>≤</sup> *<sup>e</sup>*−*Knε*2*<sup>j</sup>* <sup>2</sup> *for all θ* ∈  *such that dn(θ , θ*0*)>jε and for every j* ∈ N.

PROOF. For a given *j* ∈ N, choose a maximal set of points in *<sup>j</sup>* = {*θ* ∈ *-* : *jε<dn(θ , θ*0*)* ≤ *(j* + 1*)ε*} with the property that *en(θ , θ ) > jεξ* for every pair of points in the set. Because this set of points is a *jεξ* -net over *<sup>j</sup>* for *en* and because *(j* + 1*)ε* ≤ 2*jε*, this yields a set *- <sup>j</sup>* of at most *N (*2*jε)* points, each at *dn*-distance at least *jε* from *θ*0, and every *θ* ∈ *<sup>j</sup>* is within *en*-distance *jεξ* of at least one of these points. (If *<sup>j</sup>* is empty, we take *- <sup>j</sup>* to be empty also.) By assumption, for every point *θ*<sup>1</sup> ∈ *- <sup>j</sup>* , there exists a test with the properties as in [\(2.2\)](#page-2-0), but with *ε* replaced by *jε*. Let *φn* be the maximum of all tests attached in this way to some point *θ*<sup>1</sup> ∈ *- <sup>j</sup>* for some *<sup>j</sup>* <sup>∈</sup> <sup>N</sup>. Then

$$
P_{\theta_0}^{(n)} \phi_n \le \sum_{j=1}^{\infty} \sum_{\theta_1 \in \Theta_j'} e^{-Knj^2 \varepsilon^2} \le \sum_{j=1}^{\infty} N(2j\varepsilon) e^{-Knj^2 \varepsilon^2} \le N(\varepsilon) \frac{e^{-Kn\varepsilon^2}}{1 - e^{-Kn\varepsilon^2}}
$$

and for every *j* ∈ N,

$$
\sup_{\theta \in \bigcup_{i>j} \Theta_i} P_{\theta}^{(n)}(1 - \phi_n) \le \sup_{i > j} e^{-Kn i^2 \varepsilon^2} \le e^{-Kn j^2 \varepsilon^2},
$$

where we have used the fact that for every *θ* ∈ *<sup>i</sup>*, there exists a test *φ* with *φn* ≥ *φ* and *P(n) <sup>θ</sup> (*<sup>1</sup> <sup>−</sup> *φ)* <sup>≤</sup> *<sup>e</sup>*−*Kni*2*ε*<sup>2</sup> . This concludes the proof.

<span id="page-28-0"></span>LEMMA 10. *For k* ≥ 2, *every ε >* 0 *and every probability measure* ¯ *<sup>n</sup> supported on the set Bn(θ*0*, ε*; *k)*, *we have*, *for every C >* 0,

(8.2) 
$$
P_{\theta_0}^{(n)} \bigg( \int \frac{p_{\theta}^{(n)}}{p_{\theta_0}^{(n)}} d\bar{\Pi}_n(\theta) \leq e^{-(1+C)n\epsilon^2} \bigg) \leq \frac{1}{C^k (n\epsilon^2)^{k/2}}.
$$

PROOF. By Jensen's inequality applied to the logarithm, with *ln,θ* = log*(p(n) <sup>θ</sup> / p(n) <sup>θ</sup>*<sup>0</sup> *)*, we have log *(p(n) <sup>θ</sup> /p(n) <sup>θ</sup>*<sup>0</sup> *) d*¯ *n(θ )* <sup>≥</sup> *ln,θ d*¯ *n(θ )*. Thus, the probability in (8.2) is bounded above by

$$
(8.3) \t P_{\theta_0}^{(n)} \Big( \int (l_{n,\theta} - P_{\theta_0}^{(n)} l_{n,\theta}) d\bar{\Pi}_n(\theta) \le -n(1+C)\varepsilon^2 - \int P_{\theta_0}^{(n)} l_{n,\theta} d\bar{\Pi}_n(\theta) \Big).
$$

For every *θ* ∈ *Bn(θ*0*, ε*; *k)*, we have *P(n) θ*0 *ln,θ* = −*K(p(n) <sup>θ</sup>*<sup>0</sup> *, p(n) <sup>θ</sup> )* ≥ −*nε*2. Consequently, by Fubini's theorem and the assumption that ¯ *<sup>n</sup>* is supported on this set, the expression on the right-hand side of (8.3) is bounded above by −*Cnε*2. An application of Markov's inequality yields the upper bound

$$
\frac{P_{\theta_0}^{(n)}|\int (l_{n,\theta} - P_{\theta_0}^{(n)}l_{n,\theta}) d\bar{\Pi}_n(\theta) \wedge 0|^k}{(Cn\varepsilon^2)^k} \leq \frac{P_{\theta_0}^{(n)}\int |l_{n,\theta} - P_{\theta_0}^{(n)}l_{n,\theta}|^k d\bar{\Pi}_n(\theta)}{(Cn\varepsilon^2)^k},
$$

by another application of Jensen's inequality. The right-hand side is bounded by *C*<sup>−</sup>*k(nε*<sup>2</sup>*)*−*k/*2, by the assumption on ¯ *<sup>n</sup>*. This concludes the proof.

PROOF OF THEOREM [1.](#page-3-0) By Lemma [9,](#page-27-0) applied with *N (ε)* = exp*(nε*<sup>2</sup> *n)* (constant in *ε*) and *ε* = *Mεn* in its assertion, where *M* ≥ 2 is a large constant to be chosen later, there exist tests *φn* that satisfy *P(n) <sup>θ</sup>*<sup>0</sup> *φn* <sup>≤</sup> *<sup>e</sup>nε*<sup>2</sup> *<sup>n</sup> (*1−*e*−*KnM*2*ε*<sup>2</sup> *<sup>n</sup> )*−1*e*−*KnM*2*ε*<sup>2</sup> *n* and *P(n) <sup>θ</sup> (*1−*φn)* <sup>≤</sup> *<sup>e</sup>*−*KnM*2*ε*<sup>2</sup> *nj* <sup>2</sup> for all *θ* ∈ *<sup>n</sup>* such that *dn(θ , θ*0*)>Mεnj* and for every *j* ∈ N. The first assertion implies that if *M* is sufficiently large to ensure that *KM*<sup>2</sup> − 1 *>KM*2*/*2, then as *n* → ∞, for any *J* ≥ 1, we have

$$
(8.4) \tP_{\theta_0}^{(n)}\big[\Pi_n\big(d_n(\theta,\theta_0)\geq J M\varepsilon_n|X^{(n)}\big)\phi_n\big]\leq P_{\theta_0}^{(n)}\phi_n\lesssim e^{-KM^2n\varepsilon_n^2/2}.
$$

Setting *n,j* = {*θ* ∈ *<sup>n</sup>* : *Mεnj<dn(θ , θ*0*)* ≤ *Mεn(j* + 1*)*} and using [\(2.2\)](#page-2-0), we obtain, by Fubini's theorem,

$$
(8.5) \t P_{\theta_0}^{(n)} \bigg[ \int_{\Theta_{n,j}} \frac{p_{\theta}^{(n)}}{p_{\theta_0}^{(n)}} d\Pi_n(\theta) (1 - \phi_n) \bigg] \leq e^{-KnM^2 \varepsilon_n^2 j^2} \Pi_n(\Theta_{n,j}).
$$

Fix some *C >* 0. By Lemma 10, we have, on an event *An* with probability at least 1 − *C*<sup>−</sup>*k(nε*<sup>2</sup> *n)*<sup>−</sup>*k/*2,

$$
\int \frac{p_{\theta}^{(n)}}{p_{\theta_0}^{(n)}} d\Pi_n(\theta) \ge \int_{B_n(\theta_0,\varepsilon_n;k)} \frac{p_{\theta}^{(n)}}{p_{\theta_0}^{(n)}} d\Pi_n(\theta) \ge e^{-(1+C)n\varepsilon_n^2} \Pi_n(B_n(\theta_0,\varepsilon_n;k)).
$$

<span id="page-29-0"></span>Hence, decomposing {*θ* ∈ *-* : *dn(θ , θ*0*) > JMεn*}=∪*j*≥*<sup>J</sup> n,j* and using [\(8.5\)](#page-28-0), the last display and [\(2.5\)](#page-3-0), we have, for every sufficiently large *J* ,

$$
P_{\theta_0}^{(n)} \big[ \Pi_n (\theta \in \Theta_n : d_n(\theta, \theta_0) > J \varepsilon_n M | X^{(n)} \big) (1 - \phi_n) 1_{A_n} \big] \\
\leq \sum_{j \geq J} e^{-n \varepsilon_n^2 (K M^2 j^2 - 1 - C - \frac{1}{2} K M^2 j^2)},
$$

by assumption [\(2.5\)](#page-3-0). This converges to zero as *n* → ∞ for fixed *C* and fixed, sufficiently large *M* and *J* if *nε*<sup>2</sup> *<sup>n</sup>* → ∞; it converges to zero for fixed *M* and *C* as *J* = *Jn* → ∞ if *nε*<sup>2</sup> *<sup>n</sup>* is bounded away from zero.

Combining the preceding results, we have, for sufficiently large *M* and *J* ,

(8.6) 
$$
P_{\theta_0}^{(n)} \Pi_n(\theta \in \Theta : d_n(\theta, \theta_0) > M \varepsilon_n J | X^{(n)})
$$
  
 
$$
\leq \frac{1}{C^k (n \varepsilon_n^2)^{k/2}} + 2e^{-KM^2 n \varepsilon_n^2/2} + \sum_{j \geq J} e^{-n \varepsilon_n^2 (\frac{1}{2}KM^2 j^2 - 1 - C)}.
$$

The rest of the conclusion follows easily; see the proof of Theorem 2.4 of [\[14\]](#page-30-0).

PROOF OF THEOREM [2.](#page-3-0) If *εn n*−*<sup>α</sup>* and *k(*1−2*α) >* 2 for *α* ∈ *(*0*,* 1*/*2*)*, then *nε*<sup>2</sup> *<sup>n</sup>* → ∞ and <sup>∞</sup> *n*=1*(nε*<sup>2</sup> *n)*<sup>−</sup>*k/*<sup>2</sup> *<sup>&</sup>lt;* <sup>∞</sup>. For *<sup>C</sup>* <sup>=</sup> <sup>1</sup>*/*2, the first term on the right-hand side of (8.6) dominates and the sum over *n* of the terms in (8.6) converges. The result (i) follows by the Borel–Cantelli lemma.

For assertion (ii), we note that *εn n*−*<sup>α</sup>* and *k(*1 − 2*α)* ≥ 4*α* together imply that *(nε*<sup>2</sup> *n)*<sup>−</sup>*k/*<sup>2</sup> *<sup>ε</sup>*<sup>2</sup> *<sup>n</sup>*. The other terms are exponentially small.

PROOF OF LEMMA [1.](#page-3-0) Because *P(n) <sup>θ</sup>*<sup>0</sup> *(p(n) <sup>θ</sup> /p(n) <sup>θ</sup>*<sup>0</sup> *)* ≤ 1, Fubini's theorem implies that *P(n) <sup>θ</sup>*<sup>0</sup> [ *-*\*<sup>n</sup> (p(n) <sup>θ</sup> /p(n) <sup>θ</sup>*<sup>0</sup> *)dn(θ )*] ≤ *n(-* \ *n)*. Let the events *An* be as in the proof of Theorem [1,](#page-3-0) so that the denominator of the posterior is bounded below by *e*−*(*1+*C)nε*<sup>2</sup> *nn(Bn(θ*0*, εn*; *k))* on *An*. Combining this with the preceding display gives

$$
P_{\theta_0}^{(n)}[\Pi_n(\theta \notin \Theta_n | X^{(n)})1_{A_n}] \leq \frac{\Pi_n(\Theta \setminus \Theta_n)e^{(1+C)n\epsilon_n^2}}{\Pi_n(B_n(\theta_0, \epsilon_n; k))} \leq o(1)e^{-n\epsilon_n^2(1-C)},
$$

by the assumption on *n(-* \ *n)*. The rest of the proof can be completed along the lines of that of Theorem 2.4 of [\[14\]](#page-30-0).

## REFERENCES

[1] AMEWOU-ATISSO, M., GHOSAL, S., GHOSH, J. K. and RAMAMOORTHI, R. V. (2003). Posterior consistency for semiparametric regression problems. *Bernoulli* **9** 291–312. [MR1997031](http://www.ams.org/mathscinet-getitem?mr=1997031)

- <span id="page-30-0"></span>[2] BARRON, A., SCHERVISH, M. and WASSERMAN, L. (1999). The consistency of posterior distributions in nonparametric problems. *Ann. Statist.* **27** 536–561. [MR1714718](http://www.ams.org/mathscinet-getitem?mr=1714718)
- [3] BIRGÉ, L. (1983). Approximation dans les espaces métriques et théorie de l'estimation. *Z. Wahrsch. Verw. Gebiete* **65** 181–237. [MR0722129](http://www.ams.org/mathscinet-getitem?mr=0722129)
- [4] BIRGÉ, L. (1983). Robust testing for independent non-identically distributed variables and Markov chains. In *Specifying Statistical Models. From Parametric to Non-Parametric. Using Bayesian or Non-Bayesian Approaches*. *Lecture Notes in Statist.* **16** 134–162. Springer, New York. [MR0692785](http://www.ams.org/mathscinet-getitem?mr=0692785)
- [5] BIRGÉ, L. (2006). Model selection via testing: An alternative to (penalized) maximum likelihood estimators. *Ann. Inst. H. Poincaré Probab. Statist.* **42** 273–325 [MR2219712](http://www.ams.org/mathscinet-getitem?mr=2219712)
- [6] BROCKWELL, P. J. and DAVIS, R. A. (1991). *Time Series*: *Theory and Methods*, 2nd ed. Springer, New York. [MR1093459](http://www.ams.org/mathscinet-getitem?mr=1093459)
- [7] CHOUDHURI, N., GHOSAL, S. and ROY, A. (2004). Bayesian estimation of the spectral density of a time series. *J. Amer. Statist. Assoc.* **99** 1050–1059. [MR2109494](http://www.ams.org/mathscinet-getitem?mr=2109494)
- [8] CHOUDHURI, N., GHOSAL, S. and ROY, A. (2004). Contiguity of the Whittle measure for a Gaussian time series. *Biometrika* **91** 211–218. [MR2050470](http://www.ams.org/mathscinet-getitem?mr=2050470)
- [9] CHOW, Y. S. and TEICHER, H. (1978). *Probability Theory. Independence, Interchangeability, Martingales*. Springer, New York. [MR0513230](http://www.ams.org/mathscinet-getitem?mr=0513230)
- [10] DAHLHAUS, R. (1988). Empirical spectral processes and their applications to time series analysis. *Stochastic Process. Appl.* **30** 69–83. [MR0968166](http://www.ams.org/mathscinet-getitem?mr=0968166)
- [11] DE BOOR, C. (1978). *A Practical Guide to Splines*. Springer, New York. [MR0507062](http://www.ams.org/mathscinet-getitem?mr=0507062)
- [12] GHOSAL, S. (2001). Convergence rates for density estimation with Bernstein polynomials. *Ann. Statist.* **29** 1264–1280. [MR1873330](http://www.ams.org/mathscinet-getitem?mr=1873330)
- [13] GHOSAL, S., GHOSH, J. K. and SAMANTA, T. (1995). On convergence of posterior distributions. *Ann. Statist.* **23** 2145–2152. [MR1389869](http://www.ams.org/mathscinet-getitem?mr=1389869)
- [14] GHOSAL, S., GHOSH, J. K. and VAN DER VAART, A. W. (2000). Convergence rates of posterior distributions. *Ann. Statist.* **28** 500–531. [MR1790007](http://www.ams.org/mathscinet-getitem?mr=1790007)
- [15] GHOSAL, S., LEMBER, J. and VAN DER VAART, A. W. (2003). On Bayesian adaptation. *Acta Appl. Math.* **79** 165–175. [MR2021886](http://www.ams.org/mathscinet-getitem?mr=2021886)
- [16] GHOSAL, S. and VAN DER VAART, A. W. (2001). Entropies and rates of convergence for maximum likelihood and Bayes estimation for mixtures of normal densities. *Ann. Statist.* **29** 1233–1263. [MR1873329](http://www.ams.org/mathscinet-getitem?mr=1873329)
- [17] GHOSAL, S. and VAN DER VAART, A. W. (2007). Posterior convergence rates of Dirichlet mixtures at smooth densities. *Ann. Statist.* **35**. To appear. Available at [www4.stat.ncsu.](www4.stat.ncsu.edu/~sghosal/papers.html) [edu/~sghosal/papers.html.](www4.stat.ncsu.edu/~sghosal/papers.html)
- [18] IBRAGIMOV, I. A. (1962). Some limit theorems for stationary processes. *Theory Probab. Appl.* **7** 349–382. [MR0148125](http://www.ams.org/mathscinet-getitem?mr=0148125)
- [19] IBRAGIMOV, I. A. and HAS'MINSKII, R. Z. (1981). *Statistical Estimation: Asymptotic Theory*. Springer, New York. [MR0620321](http://www.ams.org/mathscinet-getitem?mr=0620321)
- [20] LE CAM, L. M. (1973). Convergence of estimates under dimensionality restrictions. *Ann. Statist.* **1** 38–53. [MR0334381](http://www.ams.org/mathscinet-getitem?mr=0334381)
- [21] LE CAM, L. M. (1975). On local and global properties in the theory of asymptotic normality of experiments. In *Stochastic Processes and Related Topics* (M. L. Puri, ed.) 13–54. Academic Press, New York. [MR0395005](http://www.ams.org/mathscinet-getitem?mr=0395005)
- [22] LE CAM, L. M. (1986). *Asymptotic Methods in Statistical Decision Theory.* Springer, New York. [MR0856411](http://www.ams.org/mathscinet-getitem?mr=0856411)
- [23] MEYN, S. P. and TWEEDIE, R. L. (1993). *Markov Chains and Stochastic Stability.* Springer, New York. [MR1287609](http://www.ams.org/mathscinet-getitem?mr=1287609)
- [24] PETRONE, S. (1999). Bayesian density estimation using Bernstein polynomials. *Canad. J. Statist.* **27** 105–126. [MR1703623](http://www.ams.org/mathscinet-getitem?mr=1703623)
- <span id="page-31-0"></span>[25] POLLARD, D. (1990). *Empirical Processes: Theory and Applications.* IMS, Hayward, CA. [MR1089429](http://www.ams.org/mathscinet-getitem?mr=1089429)
- [26] SHEN, X. (2002). Asymptotic normality of semiparametric and nonparametric posterior distributions. *J. Amer. Statist. Assoc.* **97** 222–235. [MR1947282](http://www.ams.org/mathscinet-getitem?mr=1947282)
- [27] SHEN, X. and WASSERMAN, L. (2001). Rates of convergence of posterior distributions. *Ann. Statist.* **29** 687–714. [MR1865337](http://www.ams.org/mathscinet-getitem?mr=1865337)
- [28] STONE, C. J. (1990). Large-sample inference for log-spline models. *Ann. Statist.* **18** 717–741. [MR1056333](http://www.ams.org/mathscinet-getitem?mr=1056333)
- [29] STONE, C. J. (1994). The use of polynomial splines and their tensor products in multivariate function estimation (with discussion). *Ann. Statist.* **22** 118–184. [MR1272079](http://www.ams.org/mathscinet-getitem?mr=1272079)
- [30] VAN DER MEULEN, F., VAN DER VAART, A. W. and VAN ZANTEN, J. H. (2006). Convergence rates of posterior distributions for Brownian semimartingale models. *Bernoulli* **12** 863– 888. [MR2265666](http://www.ams.org/mathscinet-getitem?mr=2265666)
- [31] VAN DER VAART, A. W. and WELLNER, J. A. (1996). *Weak Convergence and Empirical Processes*. *With Applications to Statistics*. Springer, New York. [MR1385671](http://www.ams.org/mathscinet-getitem?mr=1385671)
- [32] WHITTLE, P. (1957). Curve and periodogram smoothing (with discussion). *J. Roy. Statist. Soc. Ser. B* **19** 38–63. [MR0092331](http://www.ams.org/mathscinet-getitem?mr=0092331)
- [33] WONG, W. H. and SHEN, X. (1995). Probability inequalities for likelihood ratios and convergence rates of sieve MLEs. *Ann. Statist.* **23** 339–362. [MR1332570](http://www.ams.org/mathscinet-getitem?mr=1332570)
- [34] ZHAO, L. H. (2000). Bayesian aspects of some nonparametric problems. *Ann. Statist.* **28** 532– 552. [MR1790008](http://www.ams.org/mathscinet-getitem?mr=1790008)

DEPARTMENT OF STATISTICS NORTH CAROLINA STATE UNIVERSITY 2501 FOUNDERS DRIVE RALEIGH, NORTH CAROLINA 27695-8203 USA E-MAIL: [ghosal@stat.ncsu.edu](mailto:ghosal@stat.ncsu.edu)

DIVISION OF MATHEMATICS AND COMPUTER SCIENCE VRIJE UNIVERSITEIT DE BOELELAAN 1081A 1081 HV AMSTERDAM THE NETHERLANDS E-MAIL: [aad@cs.vu.nl](mailto:aad@cs.vu.nl)