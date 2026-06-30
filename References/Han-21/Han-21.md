# **COMPLEX SAMPLING DESIGNS: UNIFORM LIMIT THEOREMS AND APPLICATIONS**

BY QIYANG HAN<sup>1</sup> AND JON A. WELLNER<sup>2</sup>

<sup>1</sup>*Department of Statistics, Rutgers University, [qh85@stat.rutgers.edu](mailto:qh85@stat.rutgers.edu)* <sup>2</sup>*Department of Statistics, University of Washington, [jaw@stat.washington.edu](mailto:jaw@stat.washington.edu)*

In this paper, we develop a general approach to proving global and local uniform limit theorems for the Horvitz–Thompson empirical process arising from complex sampling designs. Global theorems such as Glivenko– Cantelli and Donsker theorems, and local theorems such as local asymptotic modulus and related ratio-type limit theorems are proved for both the Horvitz–Thompson empirical process, and its calibrated version. Limit theorems of other variants and their conditional versions are also established. Our approach reveals an interesting feature: the problem of deriving uniform limit theorems for the Horvitz–Thompson empirical process is essentially no harder than the problem of establishing the corresponding finite-dimensional limit theorems, once the usual complexity conditions on the function class are satisfied. These global and local uniform limit theorems are then applied to important statistical problems including (i) *M*-estimation, (ii) *Z*-estimation and (iii) frequentist theory of pseudo-Bayes procedures, all with weighted likelihood, to illustrate their wide applicability.

## **1. Introduction.**

1.1. *Overview*. Over the past thirty years, uniform limit theorems for the empirical process have proved to be a universal tool in various statistical problems based on independent observations; we only refer readers to the textbooks [\[35,](#page-25-0) [48,](#page-25-0) [77,](#page-26-0) [81\]](#page-26-0) for relevant theoretical developments and various statistical applications.

Our focus here will be uniform limit theorems for the Horvitz–Thompson empirical process arising from complex sampling designs (cf. [\[70\]](#page-26-0)). Such limit theorems provide fundamental probabilistic tools in statistical applications with survey data, for instance, in combination with the functional delta method (see, e.g., [\[4,](#page-23-0) [8,](#page-23-0) [9,](#page-23-0) [27\]](#page-24-0) for applications in econometrics), or in semiparametric modeling (see, e.g., [\[13–15,](#page-24-0) [50,](#page-25-0) [56,](#page-25-0) [57\]](#page-25-0) for applications in biostatistics), just to name a few. Recent years have seen the emergence of interest in further limit theory in this direction (e.g., [\[7,](#page-23-0) [11,](#page-24-0) [16,](#page-24-0) [17,](#page-24-0) [26,](#page-24-0) [68,](#page-26-0) [69\]](#page-26-0)), but the scope of the existing results in this direction has been somewhat limited, and many of these available results have been derived based on case-by-case analyses. Roughly speaking, there are three approaches so far in the literature:

1. Breslow and Wellner [\[16,](#page-24-0) [17\]](#page-24-0) developed theory in the context of two-phase sampling with phase II a simple sampling without replacement sampling design. The key idea therein is to view the Horvitz–Thompson empirical process conditionally as an exchangeably weighted bootstrap empirical process [\[60\]](#page-25-0). This idea is further exploited in [\[69\]](#page-26-0) in the context of calibrated Horvitz–Thompson empirical processes. A similar bootstrap approach is adopted in [\[68\]](#page-26-0) in the setting of stratified sampling with potential overlaps.

Received April 2019; revised December 2019.

*[MSC2020 subject classifications.](https://mathscinet.ams.org/mathscinet/msc/msc2020.html)* Primary 60E15; secondary 62G05.

*Key words and phrases.* Complex sampling design, empirical process, uniform limit theorems.

2. Bertail et al. [\[7\]](#page-23-0) derived a Donsker theorem for the Bernoulli sampling design and other sampling designs that are close enough to the rejective sampling design (= high entropy designs) under a uniform entropy condition on the indexing function class. Their techniques heavily rely on the conditional independence of the inclusion indicators.

3. Conti [\[26\]](#page-24-0) and Boistard et al. [\[11\]](#page-24-0) established Donsker theorems over one class {**1***(*−∞*,t*] : *t* ∈ R} under sampling designs with increasing level of generality, by explicit calculations that verify the one-dimensional tightness condition.

The apparent case-by-case complication here is that complex sampling designs typically induce complicated dependence structure between the samples, so in order to use existing techniques from empirical process theory, certain latent independence or exchangeability structure needs to be identified in a case-by-case routine.

On the other hand, some structural commonality is indeed hinted at by the results proved in the above cited papers: uniform laws of large numbers (i.e., Glivenko–Cantelli theorems) and uniform central limit theorems (i.e., Donsker theorems) hold under rather minimal conditions on the indexing function classes. The intriguing question naturally arises:

QUESTION 1.1. Does there exist any general approach to proving uniform limit theorems for the Horvitz–Thompson empirical process under natural conditions, without being confined to a particular form of the sampling design?

A possible solution to this very natural question, however, appears far from obvious from the previously described approaches. The challenges involved here were already noted in Lin [\[50\]](#page-25-0) as ". . . *To our knowledge there does not exist a general theory on conditions required for the tightness and weak convergence of Horvitz–Thompson processes*. . . ," dating back to as early as 2000. One of the goals of this paper is to address Question 1.1 in an appropriate general framework that includes a wide variety of sampling designs. Part of the philosophical difficulty in such a general approach is that there is an easily believable impression that any general attempt at establishing global uniform limit theorems for the Horvitz–Thompson empirical process, must necessarily give general recipes for establishing finite-dimensional convergence of the Horvitz–Thompson empirical process. In the specific context of Donsker theorems, this impression pushes one to think about the "right conditions" under which at least central limit theorems hold for a single function under various different sampling designs—a task that usually already requires a case-by-case study.

In this paper, we show that this easily believable impression need not be the rule in the context of uniform limit theorems for Horvitz–Thompson empirical processes, at least in the superpopulation framework adopted in [\[11,](#page-24-0) [65\]](#page-26-0) with uniformly positive first-order inclusion probabilities. The major "change of thinking" adopted in the current paper, interestingly, indicates that *the problem of deriving uniform limit theorems for Horvitz–Thompson empirical processes is not really more difficult than that of establishing the corresponding finite-dimensional limit theorems, once the usual complexity conditions on the function class are satisfied*. In the context of Donsker theorems, this amounts to saying that, as long as the Horvitz–Thompson empirical process converges finite-dimensionally, weak convergence at the process level follows almost automatically. Since finite-dimensional convergence is necessary for weak convergence of the process to hold, the real point here is to separate the problem of establishing finite-dimensional convergence of the Horvitz–Thompson empirical process from that of establishing a uniform limit theorem. The approach here is in part inspired by a multiplier inequality developed in a recent work of the authors [\[40\]](#page-25-0), which holds regardless of the dependence structure among the multipliers, given sufficient independence structure between the multipliers and the samples.

Establishing global uniform limit theorems serves as a first step in understanding the behavior of these Horvitz–Thompson empirical processes. In typical semi/nonparametric applications, it is also of crucial importance to understand the local behavior of these empirical processes. To this end, we further study the local behavior of the Horvitz–Thompson empirical process by characterizing its local asymptotic modulus and proving several ratio-type limit theorems. These local uniform limit theorems show that the Horvitz–Thompson empirical process typically has similar local behavior compared to its empirical process counterpart. Similar global and local uniform limit theorems are established for the calibrated version of the Horvitz–Thompson empirical processes. Some other variants of Horvitz–Thompson empirical processes are discussed. Conditional versions of the uniform limit theorems are also established.

As an illustration and a proof of concept of the utility of our global and local uniform limit theorems (and related techniques), we apply these new tools to a variety of important statistical problems, including (i) *M*-estimation, or *empirical risk minimization*, in a general nonparametric model, (ii) *Z*-estimation in a general semiparametric model and (iii) frequentist theory of pseudo-Bayesian procedures (i.e., theory of posterior contraction rates and Bernstein–von Mises type theorems), all based on weighted likelihood. Several concrete examples are illustrated to further demonstrate the applicability of these general results.

The rest of the paper is organized as follows. Section 2 is devoted to a general probabilistic framework for complex sampling designs and detailed illustrations of the theory in the context of a number of examples. Section [3](#page-6-0) studies the global and local uniform limit theorems for the Horvitz–Thompson empirical process. Section [4](#page-15-0) gives applications of the theory developed in Section [3](#page-6-0) to the statistical problems listed above. Proofs are collected in the Appendix (see Supplementary Material [\[41\]](#page-25-0)).

1.2. *Notation*. For a real-valued measurable function *f* defined on *(*X *,*A*,P)* and *p* ≥ 1, *f Lp(P)* ≡ *(P*|*f* | *p)*1*/p* denotes the usual *Lp*-norm under *P* , and *f* ∞ ≡ *f L*<sup>∞</sup> ≡ sup*x*∈<sup>X</sup> <sup>|</sup>*f (x)*|. *<sup>f</sup>* is said to be *<sup>P</sup>* -centered if *Pf* <sup>=</sup> 0. *Lp(g,B)* denotes the *Lp(P)*-ball centered at *g* with radius *B*. For simplicity, we write *Lp(B)* ≡ *Lp(*0*,B)*.

Let *(*F*,* ·*)* be a subset of the normed space of real functions *<sup>f</sup>* : <sup>X</sup> <sup>→</sup> <sup>R</sup>. Let N *(ε,*F*,* ·*)* be the *ε*-covering number, and let N[]*(ε,*F*,* ·*)* be the *ε*-bracketing number; see page 83 of [\[81\]](#page-26-0) for more details. To avoid unnecessary measurability digressions, we assume that <sup>F</sup> is countable throughout the article. As usual, for any *<sup>φ</sup>* : <sup>F</sup> <sup>→</sup> <sup>R</sup>, we write *φ(f )*<sup>F</sup> for sup*<sup>f</sup>* <sup>∈</sup><sup>F</sup> <sup>|</sup>*φ(f )*|.

Throughout the article, *ε*1*,...,εn* will be i.i.d. Rademacher random variables independent of all other random variables. *Cx* will denote a generic constant that depends only on *x*, whose numeric value may change from line to line unless otherwise specified. *a <sup>x</sup> b* and *a <sup>x</sup> b* mean *a* ≤ *Cxb* and *a* ≥ *Cxb*, respectively, and *a <sup>x</sup> b* means *a <sup>x</sup> b* and *a <sup>x</sup> b* [*a b* means *a* ≤ *Cb* for some absolute constant *C*]. For two real numbers, *a*, *b*, *a* ∨ *b* ≡ max{*a,b*} and *a* ∧ *b* ≡ min{*a,b*}. For two sequence of nonnegative real numbers {*an*}, {*bn*}, *an ( )bn* means lim*<sup>n</sup> an/bn* = 0*(*∞*)*. We slightly abuse notation by defining log*(x)* ≡ log*(x* ∨ *e)* (and similarly for log log*(x)*).

## **2. Sampling designs.**

2.1. *Setup*. Let *UN* ≡ {1*,...,N*}, and S*<sup>N</sup>* ≡ {{*s*1*,...,sn*} : 0 ≤ *n* ≤ *N,si* ∈ *UN ,si* = *sj ,* ∀*i* = *j* } (*n* = 0 corresponds to the empty set) be the collection of subsets of *UN* . We adopt the superpopulation framework as in [\[65\]](#page-26-0): Let *(*Y*,*B<sup>Y</sup> *)*, *(*Z*,*BZ*)* be measurable spaces, and {*(Yi,Zi)* ∈ Y × Z} *N <sup>i</sup>*=<sup>1</sup> be i.i.d. superpopulation samples defined on the probability space *(*<sup>X</sup> *,*A*,P(Y,Z))* <sup>≡</sup> *(*<sup>Y</sup> <sup>×</sup> <sup>Z</sup>*,*B<sup>Y</sup> <sup>⊗</sup> <sup>B</sup>Z*,P(Y,Z))*. Here, *<sup>Y</sup>(N)* <sup>≡</sup> *(Y*1*,...,YN )* is the vector of interest, and *Z(N)* ≡ *(Z*1*,...,ZN )* is an auxiliary vector. A sampling design is a function <sup>p</sup> : <sup>S</sup>*<sup>N</sup>* <sup>×</sup> <sup>Z</sup>⊗*<sup>N</sup>* → [0*,* <sup>1</sup>] such that:

- <span id="page-3-0"></span>1. for all *<sup>s</sup>* <sup>∈</sup> <sup>S</sup>*<sup>N</sup>* , *<sup>z</sup>(N)* → <sup>p</sup>*(s,z(N))* is measurable,
- 2. for all *<sup>z</sup>(N)* <sup>∈</sup> <sup>Z</sup>⊗*<sup>N</sup>* , *<sup>s</sup>* → <sup>p</sup>*(s,z(N))* is a probability measure.

The probability space we work with that includes both the superpopulation and the designspace is the same product space *(*S*<sup>N</sup>* <sup>×</sup> <sup>X</sup> *,σ(*S*<sup>N</sup> )* <sup>×</sup> <sup>A</sup>*,*P*)* as constructed in [\[11\]](#page-24-0). We include the construction here for convenience of the reader: the probability measure P is uniquely defined through its restriction on all rectangles: for any *(s,E)* ∈ S*<sup>N</sup>* × A (note that S*<sup>N</sup>* is a finite set),

(2.1) 
$$
\mathbb{P}(s \times E) \equiv \int_{E} \mathfrak{p}(s, z^{(N)}(\omega)) dP_{(Y,Z)}(\omega) \equiv \int_{E} \mathbb{P}_d(s, \omega) dP_{(Y,Z)}(\omega),
$$

where P*<sup>d</sup> (s,ω)* ≡ p*(s,z(N)(ω))*. We also use *P* to denote the marginal law of *Y* for notational convenience.

Given *(Y (N),Z(N))* and a sampling design p, let {*ξi*} *N <sup>i</sup>*=<sup>1</sup> ⊂ [0*,* 1] be random variables defined on *(*S*<sup>N</sup>* <sup>×</sup> <sup>X</sup> *,σ(*S*<sup>N</sup> )* <sup>×</sup> <sup>A</sup>*,*P*)* with *πi* <sup>≡</sup> *πi(Z(N))* <sup>≡</sup> <sup>E</sup>[*ξi*|*Z(N)*]. We further assume that {*ξi*} *N <sup>i</sup>*=<sup>1</sup> are independent of *<sup>Y</sup> (N)* conditionally on *<sup>Z</sup>(N)*. Typically, we take *ξi* <sup>≡</sup> **<sup>1</sup>***i*∈*s*, where *s* ∼ p*(*·*,Z(N))*, to be the indicator of whether or not the *i*th sample *Yi* is observed (and in this case *πi(Z(N))* = *<sup>s</sup>*∈S*<sup>N</sup>* :*i*∈*<sup>s</sup>* <sup>p</sup>*(s,Z(N))*), but we do not require this structure a priori. The *πi*'s are often referred to as the first-order inclusion probabilities, and *πij* ≡ *πij (Z(N))* ≡ E[*ξiξj* |*Z(N)*] are the second-order inclusion probabilities.

We define the Horvitz–Thompson empirical measure and empirical process as follows: for {*πi*}, {*ξi*}, {*Yi*} as above

$$
\mathbb{P}_N^{\pi}(f) \equiv \frac{1}{N} \sum_{i=1}^N \frac{\xi_i}{\pi_i} f(Y_i), \quad f \in \mathcal{F},
$$

and the associated Horvitz–Thompson empirical process

$$
\mathbb{G}_N^{\pi}(f) \equiv \sqrt{N} \big( \mathbb{P}_N^{\pi} - P \big)(f), \quad f \in \mathcal{F}.
$$

The name of such an empirical process goes back to [\[44\]](#page-25-0), in which P*<sup>π</sup> <sup>N</sup> (Y)* <sup>≡</sup> *<sup>N</sup>*−<sup>1</sup> *<sup>N</sup> <sup>i</sup>*=1*(ξi/ πi)Yi* is used as an estimator for the population mean *P(Y)* ≡ E*Y*∼*<sup>P</sup> Y* . The usual empirical measure and empirical process (i.e., with *ξi/πi* ≡ 1 for all *i* = 1*,...,N*) will be denoted by P*<sup>N</sup>* , G*<sup>N</sup>* , respectively.

ASSUMPTION A. Consider the following conditions on the sampling design p:

(A1) min1<sup>≤</sup>*i*≤*<sup>N</sup> πi* ≥ *π*<sup>0</sup> holds for some nonrandom *π*<sup>0</sup> *>* 0. (A2-LLN) <sup>1</sup> *N <sup>N</sup> <sup>i</sup>*=1*( ξi πi* <sup>−</sup> <sup>1</sup>*)* <sup>=</sup> <sup>o</sup>P*(*1*)*. (A2-CLT) √ 1 *N <sup>N</sup> <sup>i</sup>*=1*( ξi πi* <sup>−</sup> <sup>1</sup>*)* <sup>=</sup> <sup>O</sup>P*(*1*)*.

(A1) is a common assumption in the literature. (A2-LLN) says that the weights {*ξi/πi*} satisfy a law of large numbers; while (A2-CLT) says that the weights {*ξi/πi*} have a <sup>√</sup>*<sup>N</sup>* rate of convergence (so that a uniform central limit theorem for the more complicated Horvitz– Thompson empirical process G*<sup>π</sup> <sup>N</sup>* can be possible). As we will see below in the examples, a generic way of verifying these conditions is to obtain a good estimate on the correlations {*πij* − *πiπj* }*i*=*<sup>j</sup>* . Conditions on (even higher order) correlations are very common in the literature; cf. [\[10–12,](#page-23-0) [18\]](#page-24-0).

## <span id="page-4-0"></span>2.2. *Examples of sampling designs*.

EXAMPLE 2.1 (Sampling without replacement). A simple random sampling without replacement (SWOR) design <sup>p</sup> is such that for all *<sup>z</sup>(N)* <sup>∈</sup> <sup>Z</sup>⊗*<sup>N</sup>* , <sup>p</sup>*(*·*,z(N))* is the sampling without replacement design with cardinality *n(z(N))*. In this case, the parameter in this sampling design is *n(z(N))* and *(ξ*1*,...,ξN )* is a random permutation of the vector that contains 1 in the first *n(z(N))* components and 0 otherwise. Then

$$
\pi_i(z^{(N)}) = \mathbb{E}[\xi_i|z^{(N)}] = \frac{n(z^{(N)})}{N}.
$$

Condition (A1) holds if *n(z(N))/N* ≥ *c* for some constant *c >* 0. Condition (A2) is trivially satisfied since *<sup>N</sup> <sup>i</sup>*=<sup>1</sup> *ξi* = *n(z(N))*, and hence

$$
\sum_{i=1}^{N} \left( \frac{\xi_i}{\pi_i} - 1 \right) = \left( \frac{1}{n(z^{(N)})/N} \cdot \sum_{i=1}^{N} \xi_i \right) - N = 0.
$$

EXAMPLE 2.2 (Bernoulli sampling). A Bernoulli sampling design p is such that for all *<sup>z</sup>(N)* <sup>∈</sup> <sup>Z</sup>⊗*<sup>N</sup>* and *<sup>s</sup>* <sup>∈</sup> <sup>S</sup>*<sup>N</sup>* ,

$$
\mathfrak{p}(s, z^{(N)}) = \prod_{i \in s} \pi_i(z^{(N)}) \prod_{i \notin s} (1 - \pi_i(z^{(N)})).
$$

In other words, conditionally on auxiliary random variables *Z(N)*, the *ξi*'s are independent Bernoulli random variables with success probability *πi(Z(N))*, so the parameters in this sampling design are {*πi(Z(N))* : 1 ≤ *i* ≤ *N*}. Note that we allow {*πi(Z(N))*} to be unequal. Condition (A1) holds if *πi(Z(N))* ≥ *c* for some constant *c >* 0. Since

$$
\mathbb{E}\bigg(\frac{1}{\sqrt{N}}\sum_{i=1}^{N}\bigg(\frac{\xi_{i}}{\pi_{i}}-1\bigg)\bigg)^{2} = \mathbb{E}_{(Y^{(N)}, Z^{(N)})}\bigg[\mathbb{E}_{\xi^{(N)}}\frac{1}{N}\sum_{i=1}^{N}\bigg(\frac{\xi_{i}}{\pi_{i}}-1\bigg)^{2}\bigg] = \mathcal{O}(1),
$$

condition (A2) is satisfied.

EXAMPLE 2.3 (Rejective sampling and high entropy sampling). A rejective sampling design r maximizes the entropy functional p → *<sup>s</sup>*∈S*<sup>N</sup>* p*(s)*log*(*p*(s))* over all sampling designs of fixed size *n* with the constraint that the first-order inclusion probabilities equal *(π*1*,...,πN )* (cf. [\[38\]](#page-25-0)). The parameters in this sampling design are *n* and {*πi(z(N))* : 1 ≤ *i* ≤ *N*}. r can also be realized as a conditional Bernoulli sampling design with appropriate success probabilities *(p*1*,...,pN )*: for all *<sup>z</sup>(N)* <sup>∈</sup> <sup>Z</sup>⊗*<sup>N</sup>* and *<sup>s</sup>* <sup>∈</sup> <sup>S</sup>*<sup>N</sup>* ,

$$
\mathfrak{r}(s, z^{(N)}) \propto \prod_{i \in s} p_i(z^{(N)}) \prod_{i \notin s} (1 - p_i(z^{(N)})) \mathbf{1}_{|s|=n},
$$

where *<sup>N</sup> <sup>i</sup>*=<sup>1</sup> *pi(z(N))* = *n*. The relationship between *pi* and *πi* is given in, for example, the statement and proof of Theorem 5.1 of [\[37\]](#page-25-0).

Condition (A1) holds if *πi(Z(N))* ≥ *c* for some constant *c >* 0. Let *dN* ≡ *<sup>N</sup> <sup>i</sup>*=<sup>1</sup> *πi(z(N))* × *(*1 − *πi(z(N)))*, and suppose that there exists some constant *K >* 0 such that for *N* large enough

$$
\frac{N}{d_N} \leq K.
$$

Then we have

$$
\mathbb{E}\left(\frac{1}{\sqrt{N}}\sum_{i=1}^{N}\left(\frac{\xi_{i}}{\pi_{i}}-1\right)\right)^{2}
$$
\n
$$
=\mathbb{E}_{Y^{(N)},Z^{(N)}}\left[\mathbb{E}_{\xi^{(N)}}\frac{1}{N}\left(\sum_{i=1}^{N}\left(\frac{\xi_{i}}{\pi_{i}}-1\right)^{2}+\sum_{i\neq j}\left(\frac{\xi_{i}}{\pi_{i}}-1\right)\left(\frac{\xi_{j}}{\pi_{j}}-1\right)\right)\right]
$$
\n
$$
\lesssim 1+\mathbb{E}_{Y^{(N)},Z^{(N)}}\left[N^{-1}\sum_{i\neq j}|\pi_{ij}-\pi_{i}\pi_{j}|\right]=\mathcal{O}(1),
$$

where in the last inequality we used an old result due to Hajék (cf. Theorem 5.2 of [\[37\]](#page-25-0)). Hence condition (A2) is satisfied under [\(2.2\)](#page-4-0).

Assuming (for simplicity) that 0 *<* inf*<sup>i</sup> πi* ≤ sup*<sup>i</sup> πi <* 1. Then Theorems 1 and 2 in [\[6\]](#page-23-0) showed that high entropy designs satisfy a central limit theorem. More precisely, any sampling design p with first-order inclusion probabilities *(π*1*,...,πN )* and the property that

$$
D_{\mathrm{KL}}(\mathfrak{p} \parallel \mathfrak{r}) = \sum_{s \in \mathcal{S}_N} \mathfrak{p}(s) \log \frac{\mathfrak{p}(s)}{\mathfrak{r}(s)} \to 0
$$

satisfies a CLT. An alternative argument can be found in the discussions after Proposition [3.4](#page-8-0) below. In particular, all such high entropy designs satisfy conditions (A1)–(A2-CLT) under 0 *<* inf*<sup>i</sup> πi* ≤ sup*<sup>i</sup> πi <* 1. The examples in this regard examined in [\[6\]](#page-23-0) include Rao–Sampford sampling and successive sampling (under some scaling conditions).

EXAMPLE 2.4 (Stratified sampling). Suppose that *UN* is partitioned into {*UN (*1*),..., UN (k)*} according to the auxiliary variables *Z(N)* (we omit such dependence for simplicity). In other words, *<sup>k</sup>* =<sup>1</sup> *UN ()* <sup>=</sup> *UN* , *UN ()* <sup>∩</sup> *UN ( )* = ∅ for = and |*UN ()*| = *N* with *<sup>k</sup>* =<sup>1</sup> *N* <sup>=</sup> *<sup>N</sup>*. Let *<sup>n</sup>*1*,...,nk* be such that *<sup>k</sup>* =<sup>1</sup> *n* = *n*. Within each stratum *UN ()*, we draw *n* ≤ *N* samples *s* without replacement. The overall sample is *s* = *<sup>k</sup>* =<sup>1</sup> *s*. The parameters in this sampling design are the partition {*UN ()* : 1 ≤ ≤ *k*} and {*n(Z(N))* : 1 ≤ ≤ *k*}. Similar to the calculations in Example [2.1,](#page-4-0) since *<sup>i</sup>*∈*s ξi* = *n*, we have

$$
\sum_{i=1}^{N} \left( \frac{\xi_i}{\pi_i} - 1 \right) = \sum_{\ell=1}^{k} \left( \frac{1}{n_{\ell}/N_{\ell}} \sum_{i \in s_{\ell}} \xi_i \right) - N = \left( \sum_{\ell=1}^{k} N_{\ell} \right) - N = 0.
$$

Hence (A2) is satisfied. (A1) holds if *n/N* ≥ *c* for some constant *c >* 0.

EXAMPLE 2.5 (Stratified sampling with overlap). Recently, [\[68\]](#page-26-0) studied an interesting extension of the stratified sampling design as follows: suppose that {*UN (*1*),...,UN (k)*} ⊂ *UN* are *k* potentially overlapping "data sources" determined by the auxiliary variables *Z(N)*, where *k* is a fixed integer. Let *N* ≡ |*UN ()*|. For each source *UN ()*, we draw *n* ≤ *N* samples *s* without replacement. The overall sample is *s* = *<sup>k</sup>* =<sup>1</sup> *s*, which may include duplicate samples due to the overlapping nature of the data sources. The parameters in this sampling design are the same as the above example. This sampling scheme is also known as multipleframe surveys; cf. [\[42,](#page-25-0) [43,](#page-25-0) [51\]](#page-25-0).

Let *π*¯ *() <sup>i</sup>* ≡ *n/N* if *i* ∈ *UN ()* be the sampling probability of unit *i* in the data source *UN ()*, and let *ξ*¯*() <sup>i</sup>* be the indicator of whether or not unit *i* is sampled in *UN ()*. Following [\[68\]](#page-26-0), we consider the following variant of the Horvitz–Thompson empirical measure (or *Hartley empirical measure* as it is named in [\[68\]](#page-26-0)):

$$
\mathbb{P}_{N}^{H}(f) \equiv \frac{1}{N} \sum_{i=1}^{N} \sum_{\ell=1}^{k} \frac{\bar{\xi}_{i}^{(\ell)} \rho_{i}^{(\ell)}}{\bar{\pi}_{i}^{(\ell)}} \mathbf{1}_{i \in U_{N}(\ell)} f(Y_{i}),
$$

<span id="page-6-0"></span>and the associated (Hartley) empirical process

$$
\mathbb{G}_N^H(f) \equiv \sqrt{N} (\mathbb{P}_N^H - P)(f).
$$

Here, the weights {*ρ() <sup>i</sup>* <sup>≡</sup> *<sup>ρ</sup>() <sup>i</sup> (z(N))* ∈ [0*,* <sup>1</sup>]} are such that *<sup>k</sup>* =<sup>1</sup> *<sup>ρ</sup>() <sup>i</sup> (z(N))* = 1 and that *ρ() <sup>i</sup>* = 0 if *i /*∈ *UN ()*. Now letting

(2.3) 
$$
\pi_i \equiv \prod_{\ell=1}^k \bar{\pi}_i^{(\ell)}, \qquad \xi_i \equiv \sum_{\ell=1}^k \left( \mathbf{1}_{i \in U_N(\ell)} \bar{\xi}_i^{(\ell)} \rho_i^{(\ell)} \prod_{\ell' \neq \ell} \bar{\pi}_i^{(\ell')} \right) \in [0, 1],
$$

we see that the Hartley empirical measure P*<sup>H</sup> <sup>N</sup>* and the associated empirical process <sup>G</sup>*<sup>H</sup> <sup>N</sup>* reduces to the Horvitz–Thompson empirical measure and empirical process with {*πi,ξi*} specified in (2.3).

Condition (A1) holds if *n/N* ≥ *c* for some constant *c >* 0 (by noting that *k* is a fixed constant that does not depend on *Z(N)*). Now we verify (A2). Note that

$$
\frac{1}{\sqrt{N}} \sum_{i=1}^{N} \left( \frac{\xi_i}{\pi_i} - 1 \right) = \frac{1}{\sqrt{N}} \left[ \sum_{i=1}^{N} \sum_{\ell=1}^{k} \frac{\bar{\xi}_i^{(\ell)} \rho_i^{(\ell)}}{\bar{\pi}_i^{(\ell)}} \mathbf{1}_{i \in U_N(\ell)} - N \right]
$$

$$
= \sum_{\ell=1}^{k} \frac{1}{\sqrt{N}} \sum_{i=1}^{N} \left( \frac{\bar{\xi}_i^{(\ell)}}{\bar{\pi}_i^{(\ell)}} - 1 \right) \rho_i^{(\ell)} \mathbf{1}_{i \in U_N(\ell)} = \mathcal{O}_{\mathbb{P}}(1),
$$

where the last line follows by computing the second moment:

$$
\mathbb{E}\Bigg[\frac{1}{\sqrt{N}}\sum_{i=1}^{N}\Big(\frac{\bar{\xi}_{i}^{(\ell)}}{\bar{\pi}_{i}^{(\ell)}}-1\Big)\rho_{i}^{(\ell)}\mathbf{1}_{i\in U_{N}(\ell)}\Bigg]^{2}\n\times 1+\frac{1}{N}\sum_{i\neq j\in U_{N}(\ell)}\mathbb{E}_{(Y^{(N)},Z^{(N)})}\Bigg[\Big|\mathbb{E}_{\xi^{(N)}}\Big(\frac{\bar{\xi}_{i}^{(\ell)}}{\bar{\pi}_{i}^{(\ell)}}-1\Big)\Big(\frac{\bar{\xi}_{j}^{(\ell)}}{\bar{\pi}_{j}^{(\ell)}}-1\Big)\Big|\Bigg]=\mathcal{O}(1).
$$

This verifies (A2-CLT).

From the above derivation, it is easy to see that (A1)–(A2-CLT) hold with the sampling without replacement design replaced by Bernoulli sampling and rejective sampling designs.

We also note that different choices of the weights {*ρ() <sup>i</sup>* <sup>≡</sup> *<sup>ρ</sup>() <sup>i</sup> (z(N))* ∈ [0*,* 1]} lead to different asymptotic variances. Since this issue is not the main concern of this paper, we refer the readers to [\[68\]](#page-26-0) for the optimal choice of weights in the context of Bernoulli sampling and sampling without replacement designs.

**3. Theory.** In this section, we will be mainly interested in the global and local behavior of the Horvitz–Thompson empirical process. In particular, we prove a Glivenko–Cantelli theorem and a Donsker theorem that provide global information concerning the Horvitz– Thompson empirical process in the limit. As will be seen, our formulation requires almost minimal conditions. We further study local behavior of the Horvitz–Thompson empirical process by characterizing its local asymptotic modulus and several ratio limit theorems. Understanding the local behavior of the Horvitz–Thompson empirical process plays a key role in applications to statistical problems as will be demonstrated in Section [4.](#page-15-0) Corresponding results for the calibrated version of the Horvitz–Thompson empirical process are also included. We also discuss uniform limit theorems for some variants of the Horvitz–Thompson empirical process and their conditional versions thereof. Finally, we present some positive and negative results on CLTs when the condition (A1) fails.

<span id="page-7-0"></span>3.1. *Global and local limit theorems*. First we study the Glivenko–Cantelli theorem. We say that <sup>F</sup> is *<sup>P</sup>* -Glivenko–Cantelli if and only if sup*<sup>f</sup>* <sup>∈</sup><sup>F</sup> <sup>|</sup>*(*P*<sup>N</sup>* <sup>−</sup> *P)(f )*| = <sup>o</sup>P*(*1*)*.

THEOREM 3.1 (Glivenko–Cantelli theorem). *Suppose that* (*A*1) *and* (*A*2*-LLN*) *hold*. *If* F *is P -Glivenko–Cantelli with PF <* ∞ *for some measurable envelope F*, *then*

$$
\sup_{f \in \mathcal{F}} \left| \left( \mathbb{P}_N^{\pi} - P \right) (f) \right| = \mathfrak{o}_{\mathbb{P}}(1).
$$

Recall the notion of weak convergence in the Hoffmann–Jørgensen sense: Let {*X(f )*}*<sup>f</sup>* <sup>∈</sup><sup>F</sup> be a bounded process whose finite-dimensional laws correspond to the finite dimensional projections of a tight Borel law on ∞*(*F*)*. Let {*XN (f )*}*<sup>f</sup>* <sup>∈</sup><sup>F</sup> be bounded processes. We say that *XN <sup>X</sup>* in ∞*(*F*)* if and only if <sup>E</sup><sup>∗</sup>*H(XN )* <sup>→</sup> <sup>E</sup>*H(X)*˜ for all *<sup>H</sup>* <sup>∈</sup> *Cb(*<sup>∞</sup>*(*F*))*, where *Cb(*<sup>∞</sup>*(*F*))* denotes all bounded continuous functions on ∞*(*F*)*, and *<sup>X</sup>*˜ is a measurable version of *X* with separable range (so *H(X)*˜ is measurable). Equivalently, *d*BL*(XN ,X)*˜ → 0, where *d*BL is the dual bounded Lipschitz metric (cf. p. 246 of [\[35\]](#page-25-0)). It is also well known that *XN <sup>X</sup>* in ∞*(*F*)* if and only if *XN* converges to *<sup>X</sup>* finite-dimensionally, and there exists a pseudo-metric *d* on F such that for any *δN* → 0,

$$
\sup_{d(f,g)\leq\delta_N}|X_N(f)-X_N(g)|=\mathfrak{o}_{\mathbb{P}}(1).
$$

We refer the readers to [\[35,](#page-25-0) [81\]](#page-26-0) for more details. We say that F is *P* -Donsker if and only if <sup>G</sup>*<sup>N</sup>* <sup>G</sup> in ∞*(*F*)* where <sup>G</sup> is a *<sup>P</sup>* -Brownian bridge process.

THEOREM 3.2 (Donsker theorem). *Suppose that* (*A*1) *and* (*A*2*-CLT*) *hold*. *Further assume that*:

(D1) G*<sup>π</sup> <sup>N</sup> converges finite-dimensionally to a Gaussian process* <sup>G</sup>*<sup>π</sup>* .

(D2) F *is P -Donsker*.

*Then* <sup>G</sup>*<sup>π</sup> admits a tight measurable version in* ∞*(*F*) for which*, *using the same notation*,

$$
\mathbb{G}_N^{\pi} \rightsquigarrow \mathbb{G}^{\pi} \quad in \ \ell^{\infty}(\mathcal{F}).
$$

Clearly, the finite-dimensional convergence condition (D1) above is necessary for a uniform central limit theorem in ∞*(*F*)*. (D2) is also minimal. One intriguing feature of Theorem 3.2 is that a uniform central limit theorem follows essentially automatically as long as the *finite-dimensional convergence property of the Horvitz–Thompson empirical process is verified*. A similar phenomenon was also observed in [\[72\]](#page-26-0) in a univariate non-i.i.d. case.

REMARK 3.3. F is assumed to be countable for simplicity for Theorems 3.1 and 3.2. For the general uncountable case, we may use outer probability for the statement and proofs of these theorems.

Although being necessary, establishing a finite-dimensional CLT for G*<sup>π</sup> <sup>N</sup>* and identifying the covariance structure of G*<sup>π</sup>* can be a nontrivial problem for general sampling designs; see, for example, [\[5,](#page-23-0) [6,](#page-23-0) [23,](#page-24-0) [30,](#page-24-0) [37,](#page-25-0) [61–64,](#page-25-0) [82\]](#page-26-0). Below we exploit one possible strategy, inspired by [\[11\]](#page-24-0), for identifying the covariance structure of G*<sup>π</sup>* .

<span id="page-8-0"></span>PROPOSITION 3.4. *Suppose* (*A*1) *and the following conditions hold*:

(F1) *There exists q* ∈ [4*,*∞] *such that for any i*.*i*.*d*. *random variables* {*Vi*} *defined on (*X *,*A*,P(Y,Z)) with V*1*Lq (P(Y,Z)) <* ∞,

$$
\frac{1}{S_N} \left( \frac{1}{N} \sum_{i=1}^N \frac{\xi_i}{\pi_i} V_i - \frac{1}{N} \sum_{i=1}^N V_i \right) \rightsquigarrow \mathcal{N}(0, 1)
$$

*holds under* <sup>P</sup>*<sup>d</sup> (*·*,ω)* (*notation defined in* [\(2.1\)](#page-3-0)) *for P(Y,Z)-a*.*s*. *<sup>ω</sup>* <sup>∈</sup> <sup>X</sup> . *Here*, *SN is the designbased variance given by*

$$
S_N^2 \equiv \frac{1}{N^2} \sum_{1 \le i, j \le N} \frac{\pi_{ij} - \pi_i \pi_j}{\pi_i \pi_j} V_i V_j.
$$

(F2) *The* (*essentially*) *first-order inclusion probabilities satisfy*

$$
\frac{1}{N} \sum_{i=1}^{N} \frac{\pi_{ii} - \pi_i^2}{\pi_i^2} \to_{P(Y,Z)} \mu_{\pi 1},
$$

*for some nonrandom μπ*<sup>1</sup> ∈ R.

(F3) *The second-order inclusion probabilities satisfy*

$$
\sup_{N \in \mathbb{N}} \sup_{1 \le i \ne j \le N} N |\pi_{ij} - \pi_i \pi_j| \le K, \qquad \frac{1}{N} \sum_{i \ne j} \frac{\pi_{ij} - \pi_i \pi_j}{\pi_i \pi_j} \to_{P(Y, Z)} \mu_{\pi 2},
$$

*where K >* 0 *is an absolute constant*, *and μπ*<sup>2</sup> ∈ R *is nonrandom*.

*If* <sup>F</sup> *is such that FLq (P) <sup>&</sup>lt;* <sup>∞</sup> *for <sup>q</sup>* ∈ [4*,*∞] *that verifies* (*F*1), *then* <sup>G</sup>*<sup>π</sup> <sup>N</sup> converges finitedimensionally to a Gaussian process* G*<sup>π</sup> whose covariance structure is given by the following*: *for any f,g* ∈ F,

$$
Cov(\mathbb{G}^{\pi}(f), \mathbb{G}^{\pi}(g)) = (1 + \mu_{\pi 1})P(fg) - (1 - \mu_{\pi 2})(Pf)(Pg)
$$
  
=  $P(fg) - (Pf)(Pg) + \mu_{\pi 1}P(fg) + \mu_{\pi 2}(Pf)(Pg).$ 

The above covariance formula can be inferred from the decomposition

$$
\mathbb{G}_N^{\pi} = \sqrt{N}(\mathbb{P}_N^{\pi} - P) = \sqrt{N}(\mathbb{P}_N - P) + \sqrt{N}(\mathbb{P}_N^{\pi} - \mathbb{P}_N),
$$

where the covariance structure of the second term <sup>√</sup>*N(*P*<sup>π</sup> <sup>N</sup>* <sup>−</sup> <sup>P</sup>*<sup>N</sup> )* can be deduced from conditions (F1)–(F3). These conditions are also used in [\[11\]](#page-24-0): (F1) corresponds to (HT1) in [\[11\]](#page-24-0), (F2) corresponds to condition (i) in Proposition 3.1 in [\[11\]](#page-24-0), and (F3) corresponds to (C2) and condition (ii) in Proposition 3.1 in [\[11\]](#page-24-0). Combined with Proposition 3.4, we see that Theorem [3.2](#page-7-0) extends Proposition 3.2 of [\[11\]](#page-24-0) in at least the following directions: (i) we work with a general *P* -Donsker class F with *FLq (P) <* ∞ instead of one particular class {**1***(*−∞*,t*] : *t* ∈ R}, and (ii) we weaken conditions for the sampling designs, that is, (C3)–(C4) in [\[11\]](#page-24-0) are no longer required. We should, however, remind readers that Proposition 3.4 is not exhaustive for identifying the covariance structure of G*<sup>π</sup>* and, therefore, it is possible that the current conditions in Proposition 3.4 can be further weakened via other approaches.

The conditions in Proposition 3.4 are verified in [\[11\]](#page-24-0) under a slightly different setting, but for the convenience of the reader, we provide some details for various sampling designs. Below we take *q* = 4; see also Table [1](#page-9-0) for a summary.

• For sampling without replacement, *πii* = *πi* = *n/N* and *πij* = *n(n* − 1*)/N(N* − 1*)* for *i* = *j* . If *n/N* → *λ* ∈ *(*0*,* 1*)*, (F1) can be verified using Hajék's rank central limit theorem (cf. [\[36\]](#page-25-0), or Proposition A.5.3 of [\[81\]](#page-26-0)), and (F2)–(F3) are satisfied with *μπ*<sup>1</sup> = *λ*−<sup>1</sup> −1 and *μπ*<sup>2</sup> = 1 − *λ*<sup>−</sup>1. The cases for stratified sampling with/without overlaps can be considered analogously.

<span id="page-9-0"></span>TABLE 1 *Values of μπ*1, *μπ*2 *for different sampling designs*. *Here*, *<sup>λ</sup>* <sup>=</sup> lim*<sup>N</sup> n/N*, *<sup>A</sup>* <sup>=</sup> lim*<sup>N</sup> <sup>N</sup>*−<sup>1</sup> *<sup>N</sup> <sup>i</sup>*=<sup>1</sup> *<sup>π</sup>*−<sup>1</sup> *<sup>i</sup>* , *<sup>d</sup>* <sup>=</sup> lim*<sup>N</sup> <sup>N</sup>*−<sup>1</sup> *<sup>N</sup> <sup>i</sup>*=<sup>1</sup> *πi(*1 − *πi)*

- For Bernoulli sampling, *πii* = *πi* and *πij* = *πiπj* for *i* = *j* . If {*πi*} *N <sup>i</sup>*=<sup>1</sup> ⊂ [*ε,* 1 − *ε*]*(ε >* 0*)*, (F1) can be verified using the Lindeberg–Feller central limit theorem, and (F2)–(F3) are satisfied with *μπ*<sup>1</sup> = lim*<sup>N</sup> N*−<sup>1</sup> *<sup>N</sup> i*=1*(π*−<sup>1</sup> *<sup>i</sup>* − 1*)* and *μπ*<sup>2</sup> = 0.
- For rejective sampling with first-order inclusion probabilities {*πi*} *N <sup>i</sup>*=<sup>1</sup> ⊂ [*ε,* 1 − *ε*]*(ε >* 0*)*, let *dN* = *<sup>N</sup> <sup>i</sup>*=<sup>1</sup> *πi(*1 − *πi)*. (F1) can be verified by Theorem 1 of [\[6\]](#page-23-0). Using Theorem 1 of [\[10\]](#page-23-0), (F2)–(F3) are satisfied with *μπ*<sup>1</sup> = lim*<sup>N</sup> N*−<sup>1</sup> *<sup>N</sup> i*=1*(π*−<sup>1</sup> *<sup>i</sup>* − 1*)* and

$$
\mu_{\pi 2} = \lim_{N} \left[ -\frac{1}{N} \sum_{i \neq j} \frac{(1 - \pi_i)(1 - \pi_j)}{d_N} + \mathcal{O}(N d_N^{-2}) \right] = -d^{-1} (1 - \lambda)^2,
$$

provided *n/N* → *λ* ∈ *(*0*,* 1*)* and *dN /N* → *d*. The covariance structure of G*<sup>π</sup>* with high entropy sampling designs is the same as the rejective sampling design, which can be verified using the same arguments in pages 1754–1755 of [\[11\]](#page-24-0).

Hence, under the assumptions of Proposition [3.4,](#page-8-0) the covariance formula for G*<sup>π</sup>* can be written more explicitly: for any *f,g* ∈ F,

$$
Cov(\mathbb{G}^{\pi}(f), \mathbb{G}^{\pi}(g))
$$
  
= 
$$
\begin{cases} \lambda^{-1}(P(fg) - (Pf)(Pg)) & \text{under SWOR}, \\ A \cdot P(fg) - (Pf)(Pg) & \text{under Bernoulli}, \\ A \cdot P(fg) - [1 + d^{-1}(1 - \lambda)^{2}](Pf)(Pg) & \text{under Rejective}. \end{cases}
$$

Here, *λ* = lim*<sup>N</sup> n/N*, *A* = lim*<sup>N</sup> N*−<sup>1</sup> *<sup>N</sup> <sup>i</sup>*=<sup>1</sup> *<sup>π</sup>*−<sup>1</sup> *<sup>i</sup>* , *<sup>d</sup>* <sup>=</sup> lim*<sup>N</sup> <sup>N</sup>*−<sup>1</sup> *<sup>N</sup> <sup>i</sup>*=<sup>1</sup> *πi(*1 − *πi)* (the convergence is all in probability sense).

Our next goal is to study the local behavior of the Horvitz–Thompson empirical process. Although being of crucial importance in applications to semi/nonparametric statistics, to the best knowledge of the authors, this issue has not been addressed in the literature.

We first study *local asymptotic modulus* of the Horvitz–Thompson empirical process, which has been considered historically for VC-type classes of sets and function classes in [\[3,](#page-23-0) [33,](#page-24-0) [34\]](#page-24-0) in the context of usual empirical processes. As will be clear below, one of the strengths of the formulation of our theorems is that finite-dimensional convergence of G*<sup>π</sup> <sup>N</sup>* is not required for studying the local behavior of G*<sup>π</sup> <sup>N</sup>*—we only require that the weights have a <sup>√</sup>*<sup>N</sup>* convergence rate as in (A2-CLT).

Before formally stating the results on the local behavior of the Horvitz–Thompson empirical process, we need some definitions.

DEFINITION 3.5. A *local asymptotic modulus* of the Horvitz–Thompson empirical process indexed by a class of functions F is an increasing function *φ(*·*)* for which there exist <span id="page-10-0"></span>some *rN δN* <sup>≤</sup> <sup>1</sup>*/*2, both nonincreasing with *<sup>N</sup>* → <sup>√</sup>*NδN* nondecreasing, such that

(3.1) 
$$
\sup_{f \in \mathcal{F}: r_N^2 < P f^2 \leq \delta_N^2} \frac{|\mathbb{G}_N^\pi(f)|}{\phi(\sigma_P f)} = \mathcal{O}_{\mathbb{P}}(1).
$$

Here, *σ*<sup>2</sup> *<sup>P</sup> (f )* = Var*<sup>P</sup> (f )*.

DEFINITION 3.6. We say that F satisfies an entropy condition with exponent *α* ∈ *(*0*,* 2*)* if either

$$
\sup_{Q} \log \mathcal{N}(\varepsilon || F ||_{L_2(Q)}, \mathcal{F}, L_2(Q)) \lesssim \varepsilon^{-\alpha},
$$

where the supremum is over all finitely discrete measures *Q* on *(*X *,*A*)*; or

$$
\log \mathcal{N}_{[]}(\varepsilon, \mathcal{F}, L_2(P)) \lesssim \varepsilon^{-\alpha}.
$$

The entropy condition is well understood in the literature; we only refer the readers to [\[35,](#page-25-0) [77,](#page-26-0) [81\]](#page-26-0) for various examples in this regard.

THEOREM 3.7. *Suppose that* (*A*1) *and* (*A*2*-CLT*) *hold and* F *is a uniformly bounded class satisfying an entropy condition with exponent α* ∈ *(*0*,* 2*)*. *Then ωα(t)* = *t* 1−*<sup>α</sup>* <sup>2</sup> *is a local asymptotic modulus for the Horvitz–Thompson empirical process indexed by* F, *that is*, (3.1) *holds with φ* = *ωα*.

The local asymptotic modulus is a key step in understanding the behavior of the Horvitz– Thompson empirical process at a local level. This will be useful in applications in the next section. The local asymptotic modulus *ωα* cannot be improved in general; this can be shown for the usual empirical process indexed by *α-full class* (which essentially requires a lower bound for the entropy number in a more local sense; cf. [\[33\]](#page-24-0)).

One may also invert the above viewpoint by fixing one particular weight function *φ* and asking for the rate of convergence of the corresponding weighted Horvitz–Thompson empirical process. Below are two particular choices: the first one (3.2) uses *φ(x)* = *x*, and the second one (3.3) uses (essentially) *φ(x)* = *x*2.

THEOREM 3.8. *Suppose that* (*A*1) *and* (*A*2*-CLT*) *hold and* F *is a uniformly bounded class satisfying an entropy condition with exponent α* ∈ *(*0*,* 2*)*. *Let rN N*−1*/(α*+2*)* . *Then*

(3.2) 
$$
N^{1/(\alpha+2)} \sup_{f \in \mathcal{F}: \sigma_P f \ge r_N} \frac{|(\mathbb{P}_N^{\pi} - P)(f)|}{\sigma_P f} = \mathcal{O}_{\mathbb{P}}(1).
$$

*If furthermore* F *takes value in* [0*,* 1], *then for any LN* → ∞,

(3.3) 
$$
\sup_{f \in \mathcal{F}: Pf \geq L_N \cdot r_N} \left| \frac{\mathbb{P}_N^{\pi} f}{Pf} - 1 \right| = \mathfrak{o}_{\mathbb{P}}(1).
$$

Results analogous to (3.2)–(3.3) have been derived in the case of i.i.d. sampling in [\[55,](#page-25-0) [73–](#page-26-0) [75,](#page-26-0) [83\]](#page-26-0) for uniform empirical processes on (subsets of) R (or R*<sup>d</sup>* ), and are further investigated in [\[3\]](#page-23-0) for VC classes of sets, and extended by [\[33,](#page-24-0) [34\]](#page-24-0) who studied more general VC-subgraph classes.

Note that (3.3) can be viewed as a uniform law of large numbers for the weighted Horvitz– Thompson empirical process. We can also establish a central limit theorem for the weighted Horvitz–Thompson empirical process, analogous to the development in [\[1–3,](#page-23-0) [33\]](#page-24-0) for the usual empirical process.

<span id="page-11-0"></span>THEOREM 3.9. *Suppose that* (*A*1) *and* (*A*2*-CLT*) *hold*, *and that* F *is a uniformly bounded class satisfying an entropy condition with exponent α* ∈ *(*0*,* 2*)*. *Let φ* : R≥<sup>0</sup> → R≥<sup>0</sup> *be such that φ(*0*)* = 0 *and that*

$$
\frac{\phi(t)}{t^{1-\frac{\alpha}{2}}(\log\log(1/t))^{1/2}} \to \infty
$$

*as t* → 0. *If rN N*−1*/(α*+2*) and* G*<sup>π</sup> <sup>N</sup> converges finite-dimensionally to a Gaussian process* <sup>G</sup>*<sup>π</sup>* , *then* <sup>G</sup>*<sup>π</sup> admits a tight measurable version in* ∞*(*F*) for which*, *using the same notation*,

$$
\frac{\mathbb{G}_N^{\pi}(f)}{\phi(\sigma_P f)} \mathbf{1}_{\sigma_P f > r_N} \leadsto \frac{\mathbb{G}^{\pi}(f)}{\phi(\sigma_P f)} \quad \text{in } \ell^{\infty}(\mathcal{F}).
$$

The weight function in the above theorem is required to be only slightly stronger than the local asymptotic modulus by an iterated logarithmic factor. This is very natural: the weight function cannot beat the local asymptotic modulus for a weighted CLT to hold, so the condition (3.4) is optimal up to an iterated logarithmic factor.

REMARK 3.10. The countability assumption on F in Theorems [3.7–](#page-10-0)3.9 is used at a technical level via the one-sided Talagrand-type concentration inequality (cf. Proposition C.2). One may assume, for instance, pointwise measurability (cf. Example 2.3.4 of [\[81\]](#page-26-0)) for F to handle the uncountable class.

3.2. *Calibration*. In practice, since the Horvitz–Thompson estimator may be severely inefficient, calibration of the weights is often used to improve efficiency [\[28,](#page-24-0) [52\]](#page-25-0). The main purpose of this section, instead of proposing new calibration methods or addressing efficiency issues, rests in demonstrating that our theoretical results are still valid for the Horvitz– Thompson empirical process with calibrated weights.

To illustrate this, we consider one popular calibration method that aims at matching the population mean for the Horvitz–Thompson estimator [\[28\]](#page-24-0). Let <sup>Z</sup> <sup>⊂</sup> <sup>R</sup>*<sup>d</sup>* be a compact set, and *<sup>G</sup>* : <sup>R</sup> <sup>→</sup> <sup>R</sup>≥0. Let *<sup>α</sup>*ˆ*<sup>N</sup>* <sup>∈</sup> <sup>A</sup>*c*, where <sup>A</sup>*<sup>c</sup>* is a compact set of <sup>R</sup>*<sup>d</sup>* , be defined via

$$
\frac{1}{N} \sum_{i=1}^{N} \frac{\xi_i G(Z_i^\top \hat{\alpha}_N)}{\pi_i} Z_i = \frac{1}{N} \sum_{i=1}^{N} Z_i.
$$

Then the *calibrated Horvitz–Thompson empirical measure* and *calibrated Horvitz–Thompson empirical process* are defined by

$$
\mathbb{P}_N^{\pi,c}(f) \equiv \frac{1}{N} \sum_{i=1}^N \frac{\xi_i G(Z_i^\top \hat{\alpha}_N)}{\pi_i} f(Y_i), \quad f \in \mathcal{F},
$$

and

$$
\mathbb{G}_N^{\pi,c}(f) \equiv \sqrt{N} \big( \mathbb{P}_N^{\pi,c} - P \big)(f), \quad f \in \mathcal{F},
$$

respectively.

Our next theorem asserts that as long as *α*ˆ*<sup>N</sup>* converges to the "truth" 0 (which can be defined to be another value, but we use 0 for notational convenience) sufficiently fast, the global and local theorems also hold for the calibrated Horvitz–Thompson empirical process.

THEOREM 3.11. *Suppose G(*0*)* = 1, *G (*0*) >* 0. *Let* F *be a class of measurable functions with a measurable envelope F*.

(1) *Let the assumptions in Theorem* [3.1](#page-7-0) *hold with PF <* ∞. *If α*ˆ*<sup>N</sup>* = oP*(*1*)*, *then the conclusion of Theorem* [3.1](#page-7-0) *holds with* P*<sup>π</sup> <sup>N</sup> replaced by* <sup>P</sup>*π,c <sup>N</sup>* .

<span id="page-12-0"></span>(2) *Let the assumptions in Theorem* [3.2](#page-7-0) *hold with PF*<sup>2</sup> *<* ∞ (*but the finite-dimensional convergence condition is replaced by* G*π,c <sup>N</sup> converges finite-dimensionally to some Gaussian process* <sup>G</sup>*π,c*). *If* <sup>√</sup>*Nα*ˆ*<sup>N</sup>* <sup>=</sup> <sup>O</sup>P*(*1*)*, *then* <sup>G</sup>*π,c admits a tight measurable version in* ∞*(*F*) for which*, *using the same notation*,

$$
\mathbb{G}_N^{\pi,c} \rightsquigarrow \mathbb{G}^{\pi,c} \quad in \ \ell^{\infty}(\mathcal{F}).
$$

(3) *If* <sup>√</sup>*Nα*ˆ*<sup>N</sup>* <sup>=</sup> <sup>O</sup>P*(*1*)*, *then under the same conditions as in Theorems* [3.7,](#page-10-0) [3.8](#page-10-0) *and* [3.9](#page-11-0) (*but the finite-dimensional convergence condition is replaced by* G*π,c <sup>N</sup> converges finitedimensionally to some Gaussian process* G*π,c*), *the respective conclusions hold for the calibrated Horvitz–Thompson empirical process*.

The structural commonality in the above theorem is characterized by the <sup>√</sup>*N*-rate of the estimate *<sup>α</sup>*ˆ*<sup>N</sup>* . Establishing a <sup>√</sup>*N*-rate for *<sup>α</sup>*ˆ*<sup>N</sup>* is not hard: in fact we can use Theorem 3.3.1 of [\[81\]](#page-26-0) for such a purpose by verifying the asymptotic equicontinuity of the Horvitz–Thompson empirical process.

Below we exploit one possible strategy for this via the method of Proposition [3.4.](#page-8-0)

PROPOSITION 3.12. *Assume the conditions of Proposition* [3.4](#page-8-0) *and Theorem* [3.11](#page-11-0) *hold*, *and that πi* ≡ *πi(Zi) for i* = 1*,...,N*. *Further assume that G is continuous with its derivative G continuous in a neighborhood of* 0, *and the map α* → *PZ*[*G(Zα*−1*)Z*] *has a unique zero at* 0, *and PZ(ZZ)* ∈ R*d*×*<sup>d</sup> is invertible*. *Then*

(3.5) 
$$
\sqrt{N}\hat{\alpha}_N = -(G'(0))^{-1} (P_Z(ZZ^{\top}))^{-1} (\mathbb{G}_N^{\pi} - \mathbb{G}_N)Z + \mathfrak{o}_{\mathbb{P}}(1).
$$

*Furthermore*, G*π,c <sup>N</sup> converges finite-dimensionally to a tight Gaussian process* <sup>G</sup>*π,c whose covariance structure is given by the following*: *for any f,g* ∈ F,

$$
Cov(\mathbb{G}^{\pi,c}(f), \mathbb{G}^{\pi,c}(g))
$$
  
=  $P(fg) - (Pf)(Pg) + \mu_{\pi 1}P_{(Y,Z)}(\mathcal{T}(f)\mathcal{T}(g)) + \mu_{\pi 2}(P_{(Y,Z)}\mathcal{T}(f))(P_{(Y,Z)}\mathcal{T}(g)).$ 

*Here*, *the operator* <sup>T</sup> : <sup>R</sup>Y×<sup>Z</sup> <sup>→</sup> <sup>R</sup>Y×<sup>Z</sup> *is defined by*

$$
\mathcal{T}(f)(y, z) = f(y) - P_{(Y, Z)}((\xi/\pi) f(Y) Z^{\top}) (P_Z(ZZ^{\top}))^{-1} z.
$$

As we will see in the proofs, the asymptotic expansion for <sup>√</sup>*Nα*ˆ*<sup>N</sup>* in (3.5) plays a crucial role in identifying the covariance structure of G*π,c*. Although here we only study one particular calibration method that matches the population mean, other calibration methods are also possible. Typically different calibration methods only differ in terms of the exact form of the corresponding operators T ; see, for example, [\[69\]](#page-26-0) for various calibration methods under the (two-phase) stratified sampling design.

3.3. *Other variants*. Our global limit theorems in Theorems [3.1](#page-7-0) and [3.2](#page-7-0) can be used for several other variants of the Horvitz–Thompson empirical processes studied in [\[11\]](#page-24-0). We illustrate this by considering Donsker theorems for the variants as detailed below.

First, consider <sup>√</sup>*n(*P*<sup>π</sup> <sup>N</sup>* <sup>−</sup> <sup>P</sup>*<sup>N</sup> )*. We have the following.

COROLLARY 3.13. *Suppose that* (*A*1) *and* (*A*2*-CLT*) *hold*, *and that* F *is P -Donsker*. *Further suppose that the conditions in Proposition* [3.4](#page-8-0) *hold*, *and that n/N* → *λ* ∈ *(*0*,* 1*)*. *Then* <sup>√</sup>*n(*P*<sup>π</sup> <sup>N</sup>* <sup>−</sup>P*<sup>N</sup> ) converges weakly in* ∞*(*F*) to a Gaussian process* <sup>G</sup>¯ *<sup>π</sup> whose covariance structure is given by the following*: *for any f,g* ∈ F,

$$
Cov(\bar{\mathbb{G}}^{\pi}(f), \bar{\mathbb{G}}^{\pi}(g)) = \lambda(\mu_{\pi 1} P(fg) + \mu_{\pi 2}(Pf)(Pg))
$$
  
= 
$$
\begin{cases} (1 - \lambda)(P(fg) - (Pf)(Pg)) & \text{under SWOR}, \\ \lambda(A - 1) \cdot P(fg) & \text{under Bernoulli}, \\ \lambda((A - 1) \cdot P(fg) - d^{-1}(1 - \lambda)^{2}(Pf)(Pg)) & \text{under Rejective}. \end{cases}
$$

*Here*, *λ* = lim*<sup>N</sup> n/N*, *A* = lim*<sup>N</sup> N*−<sup>1</sup> *<sup>N</sup> <sup>i</sup>*=<sup>1</sup> *<sup>π</sup>*−<sup>1</sup> *<sup>i</sup>* , *<sup>d</sup>* <sup>=</sup> lim*<sup>N</sup> <sup>N</sup>*−<sup>1</sup> *<sup>N</sup> <sup>i</sup>*=<sup>1</sup> *πi(*1 − *πi)*.

The covariance formula above is a direct consequence of the assumptions in Proposition [3.4.](#page-8-0) Furthermore, the above corollary extends Theorem 3.1 of [\[11\]](#page-24-0) from the onedimensional case <sup>F</sup> = {**1***(*−∞*,t*] : *<sup>t</sup>* <sup>∈</sup> <sup>R</sup>} to a general setting.

Next, consider the Hájek empirical process. Let

$$
\mathbb{P}_N^{\pi,H}(f) \equiv \frac{1}{\hat{N}} \sum_{i=1}^N \frac{\xi_i}{\pi_i} f(Y_i), \quad \hat{N} \equiv \sum_{i=1}^N \frac{\xi_i}{\pi_i}
$$

be the Hájek empirical measure. We have the following.

COROLLARY 3.14. *Suppose that* (*A*1) *and* (*A*2*-CLT*) *hold*, *and that* F *is P -Donsker*. *Further suppose that the conditions in Proposition* [3.4](#page-8-0) *hold*, *and that n/N* → *λ* ∈ *(*0*,* 1*)*. *Then* <sup>√</sup>*n(*P*π,H <sup>N</sup>* <sup>−</sup> <sup>P</sup>*<sup>N</sup> ) converges weakly to a Gaussian process* <sup>G</sup>¯ *π,H whose covariance structure is given by the following*: *for any f,g* ∈ F,

$$
Cov(\bar{\mathbb{G}}^{\pi,H}(f), \bar{\mathbb{G}}^{\pi,H}(g)) = \lambda \mu_{\pi 1} (P(fg) - (Pf)(Pg))
$$
  
= 
$$
\begin{cases} (1 - \lambda)(P(fg) - (Pf)(Pg)) & \text{under SWOR}, \\ \lambda(A - 1) \cdot (P(fg) - (Pf)(Pg)) & \text{under Bernoulli}, \\ \lambda(A - 1) \cdot (P(fg) - (Pf)(Pg)) & \text{under Rejective}. \end{cases}
$$

*Here*, *λ* = lim*<sup>N</sup> n/N*, *A* = lim*<sup>N</sup> N*−<sup>1</sup> *<sup>N</sup> <sup>i</sup>*=<sup>1</sup> *<sup>π</sup>*−<sup>1</sup> *<sup>i</sup>* .

As we will see in the proofs, the covariance structure of the limit of <sup>√</sup>*n(*P*π,H <sup>N</sup>* <sup>−</sup>P*<sup>N</sup> )* is the same as that of

$$
f \mapsto \frac{1}{\sqrt{N}} \sum_{i=1}^{N} \left( \frac{\xi_i}{\pi_i} - 1 \right) (f(Y_i) - Pf)
$$

up to a factor of *λ*, which can be determined by the conditions of Proposition [3.4.](#page-8-0) Furthermore, the above corollary extends Theorem 4.2 of [\[11\]](#page-24-0), again from the one-dimensional case to a general setting.

REMARK 3.15. Under (F3), since the harmonic mean is less than the arithmetic mean, we have *A*−<sup>1</sup> = lim*<sup>N</sup> (N*−<sup>1</sup> *<sup>N</sup> <sup>i</sup>*=<sup>1</sup> *<sup>π</sup>*−<sup>1</sup> *<sup>i</sup> )*−<sup>1</sup> <sup>≤</sup> lim*<sup>N</sup> (N*−<sup>1</sup> *<sup>N</sup> <sup>i</sup>*=<sup>1</sup> *πi)* <sup>=</sup> lim*<sup>N</sup> <sup>n</sup> <sup>N</sup>* = *λ*, where the next to last equality follows by computing the second moment and using (F3). It then follows that *λ(A* − 1*)* ≥ 1 − *λ* under (F3).

3.4. *Conditional limit theorems*. In this section, we consider conditional versions of the (global) uniform limit theorems. For clarity of presentation, following [\[24\]](#page-24-0) and [\[84\]](#page-26-0), we introduce the following notion.

DEFINITION 3.16. Let {*N* }*N*∈<sup>N</sup> be a sequence of random variables defined on *(*S*<sup>N</sup>* × <sup>X</sup> *,σ(*S*<sup>N</sup> )*×A*,*P*)*. We say that *N* is of order <sup>o</sup>P*<sup>d</sup> (*1*)* in *P(Y,Z)*-probability if for any *ε,δ >* 0, we have *P(Y,Z)(*P*d*<sup>|</sup>*(Y,Z)(*|*N* | *> ε) > δ)* → 0 as *N* → ∞.

Below we establish conditional versions of Glivenko–Cantelli and Donsker theorems for P*π <sup>N</sup>* <sup>−</sup> <sup>P</sup>*<sup>N</sup>* .

COROLLARY 3.17 (Conditional Glivenko–Cantelli theorem). *Suppose that* (*A*1) *and* (*A*2*-LLN*) *hold*. *If* F *is P -Glivenko–Cantelli*, *then*

$$
\sup_{f \in \mathcal{F}} \left| \left( \mathbb{P}_N^{\pi} - \mathbb{P}_N \right) (f) \right| = \mathfrak{o}_{\mathbb{P}_d}(1) \quad \text{in } P_{(Y,Z)}\text{-}probability.
$$

COROLLARY 3.18 (Conditional Donsker theorem). *Suppose that* (*A*1) *and* (*A*2*-CLT*) *hold*, *and that* F *is P -Donsker*. *Further suppose that the conditions in Proposition* [3.4](#page-8-0) *hold*, *and that n/N* → *λ* ∈ *(*0*,* 1*)*. *Then*

$$
\sqrt{n}(\mathbb{P}_N^{\pi}-\mathbb{P}_N)\rightsquigarrow\bar{\mathbb{G}}^{\pi}
$$
 in  $\ell^{\infty}(\mathcal{F})$  in  $P_{(Y,Z)}$ -probability.

*Here*, G¯ *<sup>π</sup> is a Gaussian process whose covariance structure is given in Corollary* [3.13.](#page-12-0)

The precise meaning of the above conditional Donsker theorem is that *d*BL*,d (* <sup>√</sup>*n(*P*<sup>π</sup> <sup>N</sup>* − <sup>P</sup>*<sup>N</sup> ),*G¯ *<sup>π</sup> )* <sup>≡</sup> sup*H*∈BL1*(*∞*(*F*))*|E<sup>∗</sup> *d*|*(Y,Z)H(*√*n(*P*<sup>π</sup> <sup>N</sup>* <sup>−</sup> <sup>P</sup>*<sup>N</sup> ))* <sup>−</sup> <sup>E</sup>*H(*G¯ *<sup>π</sup> )*| → 0 in *P(Y,Z)* probability.

3.5. *Positive and negative results for CLTs when* (*A*1) *fails*. Let

$$
Z_N \equiv \frac{1}{\sqrt{N}} \sum_{i=1}^N \left(\frac{\xi_i}{\pi_i} - 1\right) = \mathbb{G}_N^{\pi}(1).
$$

We present a negative result concerning CLTs for *ZN* when (A1) fails.

PROPOSITION 3.19. *Let eN be such that eN* 0 *and NeN* ∞. *There exists a Bernoulli sampling scheme with equal first-order including probabilities satisfying* min1<sup>≤</sup>*i*≤*<sup>N</sup> πi* = *eN* , *such that a CLT fails for ZN* .

PROOF. Let the *ξi*'s be i.i.d. Bern*(eN )* random variables independent of *Yi*'s. Then *πi* = *eN* for 1 ≤ *i* ≤ *N*. First, note that

$$
\text{Var}\bigg[\frac{1}{\sqrt{N}}\sum_{i=1}^{N}\bigg(\frac{\xi_i}{\pi_i}-1\bigg)\bigg]=\frac{1}{N}\sum_{i=1}^{N}\pi_i^{-2}\text{Var}(\xi_i)=\frac{1-e_N}{e_N}.
$$

Hence for *Z* ∼ N *(*0*,* 1*)*, we have by Liapunov's CLT and uniform integrability

$$
\mathbb{E}\left|e_N^{1/2}\cdot\frac{1}{\sqrt{N}}\sum_{i=1}^N\left(\frac{\xi_i}{\pi_i}-1\right)\right|\to\mathbb{E}|Z|.
$$

This shows that *ZN* is not bounded in probability, and hence the CLT fails.

<span id="page-15-0"></span>Therefore, when (A1) fails, the conclusion of Theorem [3.2](#page-7-0) with <sup>√</sup>*<sup>N</sup>* normalization does not hold without further conditions.

It is easy to note from the proof above that the problem can be fixed if we change normalization from <sup>√</sup>*<sup>N</sup>* to <sup>√</sup>*n*. Specifically, if *ξi*'s are i.i.d. Bern*(eN )* with *eN* 0, *NeN* <sup>∞</sup>, then

$$
\sqrt{\frac{n}{N}}Z_N = \frac{\sqrt{n}}{N} \sum_{i=1}^N \left(\frac{\xi_i}{\pi_i} - 1\right) = \sqrt{\frac{n}{N}} \mathbb{G}_N^{\pi}(1) \quad \text{or} \quad e_N^{1/2} \mathbb{G}_N^{\pi}(1)
$$

converges to a normal random variable. This phenomenon can be generalized much further to a uniform central limit theorem as follows.

PROPOSITION 3.20. *Let* {*eN* } ⊂ *(*0*,* 1] *be such that NeN* ∞. *Suppose that the firstorder probabilities are equal in that π*<sup>1</sup> =···= *πN* = *eN with sampling indicators ξi's independent from Yi's*, *and that* <sup>√</sup> 1 *NeN <sup>N</sup> <sup>i</sup>*=1*(ξi* − *eN )* = OP*(*1*)*. *Further assume that*:

1. *e* 1*/*2 *<sup>N</sup>* <sup>G</sup>*<sup>π</sup> <sup>N</sup> converges finite-dimensionally to a Gaussian process* <sup>G</sup>*<sup>π</sup>* 0 . 2. F *is P -Donsker*.

*Then* G*<sup>π</sup>* <sup>0</sup> *admits a tight measurable version in* ∞*(*F*) for which*, *using the same notation*,

$$
e_N^{1/2} \mathbb{G}_N^{\pi} \rightsquigarrow \mathbb{G}_0^{\pi}
$$
 in  $\ell^{\infty}(\mathcal{F})$ .

One may wonder to what extent the idea above can be further generalized to the situation where the first-order inclusion probabilities can be unequal. However, as the following further counterexample shows, in such situations a CLT becomes impossible in general.

PROPOSITION 3.21. *Let α* ∈ *(*0*,* 1*)*. *Then there exists a Bernoulli sampling scheme with unequal first-order inclusion probabilities such that n/N* 0, *n/N<sup>α</sup>* ∞ *and* min1<sup>≤</sup>*i*≤*<sup>N</sup> πi* 0, *and the random variable <sup>n</sup> <sup>N</sup>* · *ZN is not bounded in probability*.

Consequently, in sharp contrast to Proposition 3.20, there is no hope for a general Donsker theory with <sup>√</sup>*<sup>n</sup>* normalization if min1<sup>≤</sup>*i*≤*<sup>N</sup> πi* 0 as long as the first-order inclusion probabilities are unequal. The proposition above is actually proving a much more negative phenomenon: although a CLT is possible under equal *πi*'s for Bernoulli sampling in the whole regime *n/N* 0, *n* ∞, such CLTs are not possible for *any convergence regime of n/N* 0, as soon as one allows unequal *πi*'s. The failure of CLTs with <sup>√</sup>*<sup>n</sup>* normalization here is particularly striking, as one would heuristically imagine that *n* is the effective sample size. The main trouble here, however, is that when (A1) fails with unequal first-order inclusion probabilities, the variance pattern of *ξi*'s can be arbitrarily complicated.

**4. Applications.** In this section, we apply the new tools developed in Section [3](#page-6-0) in statistical problems including:

- 1. *M*-estimation (or *empirical risk minimization*) in a general nonparametric model;
- 2. *Z*-estimation in a general semiparametric model;

3. frequentist theory for pseudo-Bayes procedures, namely, theory of pseudo-posterior contraction rates and Bernstein–von Mises type theorems,

where the usual likelihood is replaced by the Horvitz–Thompson weighted likelihood. We will not consider the calibrated version of these problems for simplicity of exposition, given that the corresponding theory has been fully developed in Section [3.](#page-6-0) These problems are not meant to be exhaustive; they are demonstrated as an illustration and a proof of concept of the new tools.

<span id="page-16-0"></span>4.1. *M-Estimation*. Consider the canonical *empirical risk minimization* problem (or "*M*estimation") based on weighted likelihood:

(4.1) 
$$
\hat{f}_N^{\pi} \equiv \arg \min_{f \in \mathcal{F}} \mathbb{P}_N^{\pi} f.
$$

The quality of the estimator defined in (4.1) is evaluated through the *excess risk* of *f*ˆ*<sup>π</sup> N* , denoted <sup>E</sup>*<sup>P</sup> (f*ˆ*<sup>π</sup> <sup>N</sup> )*, where

$$
\mathcal{E}_P(f) \equiv Pf - \inf_{g \in \mathcal{F}} Pg, \quad f \in \mathcal{F}.
$$

The problem of studying excess risk of empirical risk minimizers under the usual empirical measure has been extensively studied in the 2000s; we only refer the reader to [\[33,](#page-24-0) [46,](#page-25-0) [47\]](#page-25-0) and references therein. Under the Horvitz–Thompson empirical measure, [\[25\]](#page-24-0) studied risk bounds for the binary classification problem under sampling designs that are close to the rejective sampling design. Our goal here will be a study of the excess risk for the *M*-estimator based on weighted likelihood as defined in (4.1) for the general empirical risk minimization problem under general sampling designs.

To this end, let <sup>F</sup><sup>E</sup> *(δ)* ≡ {*<sup>f</sup>* <sup>∈</sup> <sup>F</sup> : <sup>E</sup>*<sup>P</sup> (f ) < δ*2}, let *ρP* : <sup>F</sup> <sup>×</sup> <sup>F</sup> <sup>→</sup> <sup>R</sup>≥<sup>0</sup> be such that *ρ*2 *<sup>P</sup> (f,g)* <sup>≥</sup> *P(f* <sup>−</sup>*g)*2−*(P(f* <sup>−</sup>*g))*2, and *D(δ)* <sup>≡</sup> sup*f,g*∈F<sup>E</sup> *(δ) ρP (f,g)*. Now we may prove the following theorem.

THEOREM 4.1. *Suppose* (*A*1) *holds*. *Suppose that there exist some L >* 0, *κ* ≥ 1 *such that D(δ)* <sup>≤</sup> *<sup>L</sup>*· *<sup>δ</sup>*1*/κ* , *and that* <sup>F</sup> *is uniformly bounded and satisfies an entropy condition with exponent α* ∈ *(*0*,* 2*)*. *Then there exist some constants* {*Ci*} 3 *<sup>i</sup>*=<sup>1</sup> *only depending on π*0, *L*, *κ*, *α such that for any s,t* ≥ 0, *with*

$$
r_N \geq C_1 N^{-\frac{\kappa}{4\kappa - 2 + \alpha}} + C_2 \bigg[ \left( \frac{s \vee t^2}{N} \right)^{\frac{\kappa}{4\kappa - 2}} \vee \frac{s}{N} \bigg],
$$

*it holds that*

$$
\mathbb{P}(\mathcal{E}_P(\hat{f}_N^{\pi}) \ge r_N^2) \le \frac{C_3}{s} e^{-s/C_3} + \mathbb{P}\Bigg(\bigg|\frac{1}{\sqrt{N}}\sum_{i=1}^N\bigg(\frac{\xi_i}{\pi_i}-1\bigg)\bigg|>t\Bigg).
$$

As an illustration of Theorem 4.1, we consider below two standard settings, regression and classification, similar to the development in [\[33\]](#page-24-0). For simplicity of exposition, we also assume that (A2-CLT) holds.

EXAMPLE 4.2 (Bounded regression). Let {*(Xi,Yi)* <sup>∈</sup> <sup>X</sup> × [−1*,* <sup>1</sup>]}*<sup>N</sup> <sup>i</sup>*=<sup>1</sup> denote the i.i.d. copies of the pairs consisting of covariates *Xi* and responses *Yi*. Our goal is to estimate the regression function *g*0*(x)* ≡ E[*Y* |*X* = *x*] using the weighted least squares method:

$$
\hat{g}_N^{\pi} \equiv \arg\min_{g \in \mathcal{G}} \sum_{i=1}^N \frac{\xi_i}{\pi_i} (Y_i - g(X_i))^2,
$$

where G is a function class containing functions taking values in [−1*,* 1], and the weights {*ξi,πi*} may depend on auxiliary information *<sup>Z</sup>(N)*. To apply Theorem 4.1, let <sup>F</sup> <sup>≡</sup> {*fg(x,y)* <sup>≡</sup> *(y* <sup>−</sup> *g(x))*<sup>2</sup> : *<sup>g</sup>* <sup>∈</sup> <sup>G</sup>}. Then following the arguments in page 1208 of [\[33\]](#page-24-0), we have <sup>E</sup>*P(X,Y)(fg)* = *<sup>g</sup>* <sup>−</sup> *<sup>g</sup>*0<sup>2</sup> *<sup>L</sup>*2*(PX)* and we may take *κ* = 1. If <sup>G</sup> satisfies an entropy condition with exponent *α* ∈ *(*0*,* 2*)*, it is easy to verify that the same holds for F*,* and hence Theorem 4.1 yields

$$
\|\hat{g}_N^{\pi} - g_0\|_{L_2(P)}^2 = \mathcal{O}_{\mathbb{P}}(N^{-\frac{2}{2+\alpha}}),
$$

a very typical rate in the regression problem.

<span id="page-17-0"></span>EXAMPLE 4.3 (Classification). Let {*(Xi,Yi)* <sup>∈</sup> <sup>X</sup> × {0*,* <sup>1</sup>}}*<sup>N</sup> <sup>i</sup>*=<sup>1</sup> denote the i.i.d. copies of the pairs consisting of covariates *Xi* and responses *Yi*. A classifier *g* : X → {0*,* 1} over a class G has a generalization error *P(X,Y)(Y* = *g(X))*. The excess risk for a classifier *g* over G under law *P(X,Y)* is given by

$$
\mathcal{E}_{P(X,Y)}(g) \equiv P_{(X,Y)}(Y \neq g(X)) - \inf_{g' \in \mathcal{G}} P_{(X,Y)}(Y \neq g'(X)).
$$

It is known that for a given law *P(X,Y)* on *(X,Y)*, the minimal generalized error is attained by a Bayes classifier *g*0*(x)* ≡ **1***η(x)*≥1*/*<sup>2</sup> where *η(x)* ≡ E[*Y* |*X* = *x*]; cf. [\[29\]](#page-24-0). In the setting of complex sampling design, it is natural to estimate *g*<sup>0</sup> by minimizing the weighted training error:

$$
\hat{g}_N^{\pi} \equiv \arg\min_{g \in \mathcal{G}} \sum_{i=1}^N \frac{\xi_i}{\pi_i} \mathbf{1}_{Y_i \neq g(X_i)},
$$

where *g*<sup>0</sup> ∈ G is a collection of classifiers. To apply Theorem [4.1,](#page-16-0) let F ≡ {*fg* ≡ **1***y*=*g(x)* : *g* ∈ G}. Suppose the following margin condition (cf. [\[54,](#page-25-0) [76\]](#page-26-0)) holds for some *c >* 0, *κ* ≥ 1: for all *g* ∈ G*,*

$$
(4.2) \t\t\t\mathcal{E}_{P_{(X,Y)}}(g) \ge c \Pi^{\kappa}\big(g(X) \ne g_0(X)\big),
$$

where is the marginal law of *X* under *P* . Following page 1212 of [\[33\]](#page-24-0), we may choose *D(δ) <sup>δ</sup>*1*/κ* , and hence if the collection of classifiers <sup>G</sup> satisfies an entropy condition with exponent *<sup>α</sup>* <sup>∈</sup> *(*0*,* <sup>2</sup>*)*, using *(fg*<sup>1</sup> <sup>−</sup> *fg*<sup>2</sup> *)*<sup>2</sup> <sup>≤</sup> *(g*<sup>1</sup> <sup>−</sup> *<sup>g</sup>*2*)*2, we see that <sup>F</sup> also satisfies the same entropy condition, and hence

$$
P_{(X,Y)}(Y \neq \hat{g}_N^{\pi}(X)) - \inf_{g' \in \mathcal{G}} P_{(X,Y)}(Y \neq g'(X)) = \mathcal{O}_{\mathbb{P}}(N^{-\frac{\kappa}{2\kappa - 1 + \alpha/2}}),
$$

a very typical rate in the classification problem.

4.2. *Z-Estimation*. The method of *Z*-estimation that produces estimators by finding those values of the parameters which zero out a set of estimating equations is well understood by now under the usual empirical measure; see [\[78,](#page-26-0) [81\]](#page-26-0) for a comprehensive treatment. With the Horvitz–Thompson empirical measure, [\[16,](#page-24-0) [17,](#page-24-0) [68,](#page-26-0) [69\]](#page-26-0) considered weighted likelihood estimation under stratified sampling designs, both with and without overlaps. The goal of this section is to give a unified theoretical treatment for the *Z*-estimation problem under general sampling designs.

Let *θ*ˆ*<sup>π</sup> <sup>N</sup>* ∈ solve the (possibly infinite-dimensional) estimating equations based on weighted likelihood:

$$
\mathbb{P}_N^{\pi} \psi_{\hat{\theta}_N^{\pi}, h} = 0 \quad \text{for all } h \in \mathcal{H},
$$

while the "truth" *θ*<sup>0</sup> ∈ solves the population equations

$$
P\psi_{\theta_0,h}=0 \quad \text{for all } h \in \mathcal{H}.
$$

Let *N ,* : <sup>→</sup> ∞*(*H*)* be given by *N (θ)(h)* <sup>≡</sup> <sup>P</sup>*<sup>π</sup> <sup>N</sup> ψθ,h* and *(θ)(h)* ≡ *Pψθ,h*. We assume that H is countable without loss of generality.

THEOREM 4.4. *Suppose that* (*A*1) *and* (*A*2*-CLT*) *hold*, *and that the following conditions hold*:

(Z1) *The map is Fréchet differentiable at θ*<sup>0</sup> *with a continuously invertible derivative* ˙ *<sup>θ</sup>*<sup>0</sup> .

<span id="page-18-0"></span>(Z2) *The stochastic equicontinuity condition holds*:

$$
\|\mathbb{G}_N(\psi_{\hat{\theta}_N^{\pi},h} - \psi_{\theta_0,h})\|_{\mathcal{H}} = \mathfrak{o}_{\mathbb{P}}(1+\sqrt{N}\|\hat{\theta}_N^{\pi} - \theta_0\|)
$$

*and* {*ψθ*0*,h* : *h* ∈ H} *is a P -Glivenko–Cantelli class*. *If θ*ˆ*<sup>π</sup> <sup>N</sup>* →**<sup>P</sup>** *θ*0, *then* √

$$
\sqrt{N}(\hat{\theta}_N^{\pi}-\theta_0)=-\dot{\Psi}_{\theta_0}^{-1}\mathbb{G}_N^{\pi}\psi_{\theta_0,+}+\mathfrak{o}_{\mathbb{P}}(1).
$$

This theorem is comparable to the standard *Z*-Theorem 3.3.1 in [\[81\]](#page-26-0), but here we work in the context of *Z*-estimation under weighted likelihood. Note that our conditions are are almost identical to the standard *Z*-Theorem, many examples for which Theorem [4.4](#page-17-0) applies can be found in Section 3.3 of [\[81\]](#page-26-0) (see also [\[78,](#page-26-0) [79\]](#page-26-0)). In particular, (Z2) is imposed for the usual empirical process G*<sup>N</sup>* , and can be easily checked if a Donsker property for the class {*ψθ,h* − *ψθ*0*,h* : *θ* − *θ*0 ≤ *δ,h* ∈ H} holds. We omit these details here.

Now consider estimation of a finite-dimensional parameter in the presence of an infinitedimensional nuisance parameter, that is, estimation in a semiparametric model. Following [\[24,](#page-24-0) [53\]](#page-25-0), we use the following general semiparametric framework: Consider a model {*Pθ,η* : *(θ,η)* <sup>∈</sup> <sup>R</sup>*<sup>d</sup>* <sup>×</sup> <sup>H</sup>}, where <sup>H</sup> is an infinite dimensional Hilbert space with norm ·H. Suppose that the true parameter is *(θ*0*,η*0*)*. An estimator *(θ*ˆ*<sup>π</sup> <sup>N</sup> , <sup>η</sup>*ˆ*<sup>π</sup> <sup>N</sup> )* of *(θ*0*,η*0*)* usually takes the form

(4.3) 
$$
(\hat{\theta}_N^{\pi}, \hat{\eta}_N^{\pi}) := \arg \sup \mathbb{P}_N^{\pi} m_{\theta, \eta},
$$

where *mθ,η* is often the log likelihood function (for *n* = 1). However, here we will work with a more general *Z*-estimation framework.

For any fixed *η* ∈ H, let *η(t)* be a smooth curve at *t* = 0 with *η(*0*)* = *η* and *a* = *(∂/∂t)η(t)*|*t*=<sup>0</sup> for some *a* ∈ H. Denote A ⊂ H the collection for all such admissible *a*'s. Now let *mθ (θ,η)* = *∂θm(θ,η)* ∈ R*<sup>d</sup>* , *mη(θ,η)*[*a*] = *(∂/∂t)m(θ,η(t))*|*t*=<sup>0</sup> with *∂tη(t)*|*t*=<sup>0</sup> = *a* ∈ A. The second derivatives can be defined in a similar fashion. Suppose further the following orthogonality condition hold: there exists *A*<sup>∗</sup> = *(a*<sup>∗</sup> <sup>1</sup> *,...,a*<sup>∗</sup> *<sup>d</sup> )* <sup>∈</sup> <sup>A</sup>*<sup>d</sup>* so that for any *<sup>A</sup>* <sup>∈</sup> <sup>A</sup>*<sup>d</sup>* , it holds that

(4.4) 
$$
P_{\theta_0,\eta_0}(m_{\theta\eta}(\theta_0,\eta_0)[A]-m_{\eta\eta}[A^*][A])=0.
$$

This condition is commonly adopted in semiparametric literature to handle the case when nuisance parameter is not <sup>√</sup>*n*-estimable; see, for example, Condition 2, page 555 in [\[45\]](#page-25-0).<sup>1</sup>

Define the *efficient score function m(θ,η)* ˜ = *mθ (θ,η)* − *mη(θ,η)*[*A*∗] (since if *m* is the log likelihood function, *m*˜ typically becomes the efficient score function). Then (4.4) can be rewritten as following: for any *<sup>A</sup>* <sup>∈</sup> <sup>A</sup>*<sup>d</sup>* ,

(4.5) 
$$
P_{\theta_0, \eta_0} \tilde{m}_{\eta}(\theta_0, \eta_0)[A] = 0.
$$

We assume that the true parameter *(θ*0*,η*0*)* zeros out the population estimating equation:

$$
(4.6) \t\t P_{\theta_0,\eta_0} \tilde{m}(\theta_0,\eta_0) = 0.
$$

To allow some flexibility in the framework, the estimators *(θ*ˆ*<sup>π</sup> <sup>N</sup> , <sup>η</sup>*ˆ*<sup>π</sup> <sup>N</sup> )* are assumed to approximately zero out the Horvitz–Thompson empirical estimating equation:

(4.7) 
$$
\mathbb{P}_N^{\pi} \tilde{m}(\hat{\theta}_N^{\pi}, \hat{\eta}_N^{\pi}) = \mathfrak{o}_{\mathbb{P}}(N^{-1/2}).
$$

It is easy to see that the above condition is satisfied if (4.3) holds. Note here our general condition also includes the case where *η*ˆ*<sup>π</sup> <sup>N</sup>* may depend on *<sup>θ</sup>*ˆ*<sup>π</sup> <sup>N</sup>* , for example, profile likelihood estimation.

<sup>1</sup>See also condition A3 in [\[85\]](#page-26-0), page 2138; condition (4) in [\[53\]](#page-25-0), page 196; condition (4) in [\[24\]](#page-24-0), page 2887.

<span id="page-19-0"></span>THEOREM 4.5. *Suppose that* (*A*1) *holds*, *and that* [\(4.5\)](#page-18-0)*–*[\(4.7\)](#page-18-0) *hold*. *Further assume the following conditions*:

(S1) *The matrix Iθ*0*,η*<sup>0</sup> ≡ −*Pθ*0*,η*0*m*˜ *<sup>θ</sup> (θ*0*,η*0*)* ∈ R*d*×*<sup>d</sup> is nonsingular*.

(S2) *θ*ˆ*<sup>π</sup> <sup>N</sup>* <sup>−</sup> *<sup>θ</sup>*0∨ˆ*η<sup>π</sup> <sup>N</sup>* <sup>−</sup> *<sup>η</sup>*0<sup>H</sup> <sup>=</sup> <sup>O</sup>P*(N*<sup>−</sup>*β) holds for some β >* <sup>1</sup>*/*4.

(S3) *The model is smooth in the sense that*

$$
\|P_{\theta_0,\eta_0}(\tilde{m}(\theta,\eta)-\tilde{m}(\theta_0,\eta_0)-\tilde{m}_{\theta}(\theta_0,\eta_0)(\theta-\theta_0))\|
$$
  
=  $\mathcal{O}(\|\theta-\theta_0\|^2 \vee \|\eta-\eta_0\|_{\mathcal{H}}^2)$ 

*holds for (θ,η) close enough to (θ*0*,η*0*)*.

(S4) *For any C >* 0,

$$
\sup_{\|\theta-\theta_0\| \vee \|\eta-\eta_0\|_{\mathcal{H}} \le CN^{-\beta}} \left| \mathbb{G}_N(\tilde{m}(\theta,\eta)-\tilde{m}(\theta_0,\eta_0))\right| = \mathfrak{o}_{\mathbb{P}}(1).
$$

*Then*

$$
\sqrt{N}(\hat{\theta}_N^{\pi} - \theta_0) = I_{\theta_0, \eta_0}^{-1} \mathbb{G}_N^{\pi} \tilde{m}(\theta_0, \eta_0) + \mathfrak{o}_{\mathbb{P}}(1).
$$

Conditions (S1)–(S4) are all standard assumptions in semiparametric literature, and can be verified in numerous models, including the Cox model with right censored/current status data, partially linear model, panel count data (with covariates), etc. Here, we only consider the partially linear model; detailed verifications for other models can be found in, for example, [\[24,](#page-24-0) [53,](#page-25-0) [68,](#page-26-0) [69,](#page-26-0) [85\]](#page-26-0).

EXAMPLE 4.6 (Partially linear model). Consider the following model:

$$
Y_i = X_i^{\top} \theta_0 + f_0(W_i) + e_i, \quad i = 1, ..., N,
$$

where *Yi*'s are the responses, {*(Xi,Wi)* ∈ [−1*,* 1] *<sup>d</sup>* × [0*,* 1]}'s are i.i.d. covariates, and *ei*'s are i.i.d. normal errors independent of the covariates. The "true signal" *θ*<sup>0</sup> ∈ R*<sup>d</sup>* and *f*<sup>0</sup> : [0*,* 1] → R is a "smooth" function. For ease of exposition, we will consider the parameter space ≡ {*(θ,f )* : *θ*<sup>1</sup> ≤ 1*, f* ∞ ≤ 1*,J(f )* ≤ *M*} for some *M >* 0, and here *J* <sup>2</sup>*(f )* := 1 <sup>0</sup> *(f (t))*<sup>2</sup> <sup>d</sup>*t*. Now with *<sup>λ</sup>*¯ *<sup>N</sup> <sup>N</sup>*−2*/*5, let

(4.8) 
$$
(\hat{\theta}_N^{\pi}, \hat{f}_N^{\pi}) := \arg \min_{(\theta, f) \in \Xi} \left[ \mathbb{P}_N^{\pi} (Y - X^{\top} \theta - f(W))^{2} + \bar{\lambda}_N^{2} J^{2}(f) \right].
$$

To put the model into our framework, let *m(θ,f )* := −*(y* − *xθ* − *f (w))*2. Then for any admissible *a*, *b*, we have

$$
m_{\theta}(\theta, f) = 2x(y - x^{\top}\theta - f(w)), \qquad m_f(\theta, f)[a] = 2a(w)(y - x^{\top}\theta - f(w)),
$$
  
\n
$$
m_{\theta f}(\theta, f)[b] = -2xb(w), \qquad m_{ff}(\theta, f)[a][b] = -2a(w)b(w).
$$

Now let *A*∗*(W)* = E[*X*|*W*] ∈ R*<sup>d</sup>* . Then a direct calculation verifies [\(4.4\)](#page-18-0). Thus we can take

(4.9) 
$$
\tilde{m}(\theta, f) = 2(y - x^{\top} \theta - f(w))(x - \mathbb{E}[X|W = w]).
$$

[\(4.6\)](#page-18-0) is immediately verified; [\(4.7\)](#page-18-0) can also be verified by taking partial derivatives in the definition (4.8) and noting that *λ*¯ <sup>2</sup> *<sup>N</sup>* <sup>=</sup> <sup>o</sup>*(N*−1*/*2*)*. Now we verify (S1)–(S4). (S1) will be satisfied if the matrix *Iθ*0*,η*<sup>0</sup> ≡ 2E[*(X* − E[*X*|*W*]*)X*] = 2E[*(X* − E[*X*|*W*]*)*⊗2] is nonsingular. (S2) can be verified with *β* = 2*/*5 along the lines of Lemma 25.88 in [\[80\]](#page-26-0) with the tools developed in Section [3.](#page-6-0) (S3) is trivially satisfied since *m*˜ is linear in *θ* and *f* . (S4) is also easy to verify. Hence we have shown that under the same conditions as in Lemma 25.88 of [\[80\]](#page-26-0),

$$
\sqrt{N}(\hat{\theta}_N^{\pi} - \theta_0) = I_{\theta_0, \eta_0}^{-1} \mathbb{G}_N^{\pi} \tilde{m}(\theta_0, \eta_0) + \mathfrak{o}_{\mathbb{P}}(1).
$$

4.3. *Frequentist theory for pseudo-Bayesian procedures*. Suppose the i.i.d. superpopulation variables of interest {*Yi*} *N <sup>i</sup>*=<sup>1</sup> have law *Pf*<sup>0</sup> where *f*<sup>0</sup> belongs to a statistical model F and {*Pf* }*<sup>f</sup>* <sup>∈</sup><sup>F</sup> is dominated by a *σ* -finite measure *μ*. A Bayesian approach assigns a prior *N* on the model F and makes estimation/inference based on the posterior distribution. In the case where all the superpopulation {*Yi*} *N <sup>i</sup>*=<sup>1</sup> are available, by Bayes' formula, the posterior distribution, that is, a random measure on F, is defined as follows: for a measurable subset *B* ⊂ F,

$$
(4.10) \quad \Pi_N(B|Y^{(N)}) \equiv \frac{\int_B \prod_{i=1}^N p_f(Y_i) \, d\Pi_N(f)}{\int \prod_{i=1}^N p_f(Y_i) \, d\Pi_N(f)} = \frac{\int_B \exp(N \mathbb{P}_N \log p_f) \, d\Pi_N(f)}{\int \exp(N \mathbb{P}_N \log p_f) \, d\Pi_N(f)},
$$

where *pf (*·*)* denotes the probability density function of *Pf* with respect to the dominating measure *μ*.

In the current superpopulation setup with complex sampling designs, we may naturally replace the usual empirical measure P*<sup>N</sup>* in (4.10) by the Horvitz–Thompson empirical measure P*<sup>π</sup> <sup>N</sup>* to define the *pseudo-posterior distribution with weighted likelihood* as follows: for a measurable subset *B* ⊂ F,

$$
(4.11) \quad \Pi_N^{\pi}(B|D^{(N)}) \equiv \frac{\int_B \prod_{i=1}^N p_f(Y_i)^{\xi_i/\pi_i} d\Pi_N(f)}{\int \prod_{i=1}^N p_f(Y_i)^{\xi_i/\pi_i} d\Pi_N(f)} = \frac{\int_B \exp(N \mathbb{P}_N^{\pi} \log p_f) d\Pi_N(f)}{\int \exp(N \mathbb{P}_N^{\pi} \log p_f) d\Pi_N(f)}.
$$

Recall here *D(N)* ≡ *(Y(N),Z(N),ξ (N),π(N))*. Note that since *<sup>N</sup> <sup>i</sup>*=<sup>1</sup> *pf (Yi)ξi/πi* is not a true likelihood because of the power *ξi/πi*, the resulting expression is not a posterior distribution in the usual sense, and hence we call it a *pseudo-posterior based on weighted likelihood*. Bayesian inference based on (4.11) in the complex sampling setting is initiated in [\[71\]](#page-26-0), and is further developed in [\[49\]](#page-25-0). As we will see below, one particular advantage of the pseudoposterior distribution with weighted likelihood defined above is that we may obtain a complete frequentist theory for pseudo-Bayes procedures analogous to that based on observing the whole superpopulation {*Yi*} *N <sup>i</sup>*=1.

4.3.1. *Pseudo-posterior contraction rate theory*. We say that the pseudo-posterior distribution with weighted likelihood, namely *π <sup>N</sup> (*·|*D(N))*, contracts at a rate *δN* with respect to a metric *d* if

$$
P_{f_0} \Pi_N^{\pi} \big( f \in \mathcal{F} : d^2(f, f_0) > L_N \delta_N^2 | D^{(N)} \big) \to 0
$$

for any *LN* → ∞.

Our first goal in this section is to develop some useful results in deriving such pseudoposterior contraction rates for the pseudo-posterior distribution using weighted likelihood. We will use (essentially the same) machinery developed in [\[39\]](#page-25-0) (which we find easier to adapt in the current context than the standard machinery [\[31,](#page-24-0) [32\]](#page-24-0)). For some *v >* 0, *c* ∈ [0*,*∞*)* let

$$
\psi_{v,c}(\lambda) = v\lambda^2 \cdot \mathbf{1}_{|\lambda| \le 1/c} + \infty \cdot \mathbf{1}_{|\lambda| > 1/c}
$$

denote the local quadratic function.

THEOREM 4.7. *Suppose* (*A*1) *holds and the following conditions hold*:

(B1) (*Local Gaussianity condition*) *There exist some constants c*<sup>1</sup> *>* 0 *and κ* = *(κg,κ)* ∈ *(*0*,*∞*)* × [0*,*∞*) such that for all <sup>n</sup>* <sup>∈</sup> <sup>N</sup>, *and <sup>f</sup>*0*,f*<sup>1</sup> <sup>∈</sup> <sup>F</sup>,

$$
P_{f_0} \exp\biggl[\lambda\biggl(\log\frac{p_{f_0}}{p_{f_1}} - P_{f_0}\log\frac{p_{f_0}}{p_{f_1}}\biggr)\biggr] \leq c_1 \exp\bigl[\psi_{\kappa_g d^2(f_0, f_1), \kappa_{\Gamma}}(\lambda)\bigr].
$$

<span id="page-21-0"></span>*Here*, *<sup>d</sup>* : <sup>F</sup> <sup>×</sup> <sup>F</sup> <sup>→</sup> <sup>R</sup>≥<sup>0</sup> *is a symmetric function satisfying*

$$
c_2 \cdot d^2(f_0, f_1) \le P_{f_0} \log \frac{p_{f_0}}{p_{f_1}} \le c_3 \cdot d^2(f_0, f_1)
$$

*for some constants c*2*,c*<sup>3</sup> *>* 0.

(B2) (*Local entropy condition*) *There exist some* {*δN* }*N*∈<sup>N</sup> *such that*

$$
1+\sup_{\varepsilon>\delta_N}\log\mathcal{N}(c_5\varepsilon,\{f\in\mathcal{F}:d(f,f_0)\leq 2\varepsilon\},d)\leq c_4N\delta_N^2,
$$

*where c*<sup>4</sup> ∈ *(*0*,* 1*)*, *c*<sup>5</sup> ∈ *(*0*,* 1*/*4*) depend on* {*ci*} 3 *<sup>i</sup>*=1. (B3) (*Prior mass condition*) *For all j* ∈ N,

$$
\frac{\Pi_N(\lbrace f \in \mathcal{F} : j\delta_N < d(f, f_0) \le (j+1)\delta_N \rbrace)}{\Pi_N(d(f, f_0) \le \delta_N)} \le \exp\left(c_6 j^2 N \delta_N^2\right),
$$

*where c*<sup>6</sup> *>* 0 *is a small enough constant depending on* {*ci*} 3 *<sup>i</sup>*=1.

*Then*

$$
P_{f_0} \Pi_N^{\pi}(f \in \mathcal{F} : d^2(f, f_0) > C_1 \delta_N^2 | D^{(N)} \rangle \le C_2 \exp(-N \delta_N^2 / C_2).
$$

*Here*, *C*1*,C*<sup>2</sup> *>* 0 *only depend on* {*ci*} 3 *<sup>i</sup>*=<sup>1</sup> *and κ*.

The local Gaussianity condition (B1) can be easily verified in a wide range of experiments including regression/density estimation/Gaussian autoregression/Gaussian time series/covariance matrix estimation, etc. (B2)–(B3) are standard conditions in the literature. In particular, (B3) allows the exact <sup>√</sup>*<sup>N</sup>* parametric pseudo-posterior contraction rate, which will be useful below. It is also possible to consider hierarchical priors to formulate a similar theorem as in [\[39\]](#page-25-0)—in essence all examples therein can be considered here (except for regression where random design instead of fixed design is needed to maintain the i.i.d. property of the superpopulation {*Yi*} *N <sup>i</sup>*=1).

4.3.2. *Bernstein–von Mises theorem*. Next, we will be interested in a more precise limiting distribution of the pseudo-posterior distribution with weighted likelihood, that is, a Bernstein–von Mises type theorem. To this end, we work with a finite-dimensional model, where is a compact subset of R*<sup>d</sup>* . Let *θ*<sup>0</sup> ∈ , an interior point of , be the true parameter. Let N*μ,* denote the *d*-dimensional normal distribution with mean *μ* and covariance matrix .

THEOREM 4.8. *Suppose that* (*A*1) *and* (*A*2*-CLT*) *hold*. *Further assume the following conditions*:

(Bv1) (*Experiment*) *The map θ* → log*pθ (x)* = *θ (x) is differentiable at θ*<sup>0</sup> *for all x with derivative* ˙ *<sup>θ</sup>*<sup>0</sup> *(x)*, *and for θ*1, *θ*<sup>2</sup> *close enough to θ*,

$$
\left|\ell_{\theta_1}(x) - \ell_{\theta_2}(x)\right| \le m(x) \|\theta_1 - \theta_2\|
$$

*holds for some Pθ*<sup>0</sup> *-square integrable function m*. *Furthermore*, *the log-likelihood ratio* {log *pθ pθ*<sup>0</sup> }*θ*∈ *satisfies the local Gaussianity condition*, *and is twice differentiable under Pθ*<sup>0</sup> *with a nonsingular Hessian Iθ*<sup>0</sup> : *for θ close enough to θ*0,

$$
P_{\theta_0} \log \frac{p_{\theta_0}}{p_{\theta}} = \frac{1}{2} (\theta - \theta_0) I_{\theta_0} (\theta - \theta_0) + \mathfrak{o} (\|\theta - \theta_0\|^2).
$$

(Bv2) (*Prior*) *The prior has a Lebesgue density bounded away from* 0 *and* ∞ *on* .

<span id="page-22-0"></span>*Then the pseudo-posterior distribution with weighted likelihood <sup>π</sup> <sup>N</sup> converges to a sequence of normal distributions in the total variational distance*:

$$
\sup_{B}\left|\Pi_{N}^{\pi}(\sqrt{N}(\theta-\theta_{0})\in B|D^{(N)})-\mathcal{N}_{I_{\theta_{0}}^{-1}\mathbb{G}_{N}^{\pi}\ell_{\theta_{0}},I_{\theta_{0}}^{-1}}(B)\right|=\mathfrak{o}_{\mathbb{P}}(1).
$$

Note that in finite-dimensional problems, the efficient score *m*˜ in Theorem [4.5](#page-19-0) can usually be taken as ˙ *<sup>θ</sup>*<sup>0</sup> . Then under the regularity conditions as in Theorem [4.5,](#page-19-0) we have the usual interpretation of the Bernstein–von Mises theorem in our context of weighted likelihood estimation: the sequence of pseudo-posterior distributions with weighted likelihood resembles that of progressively sharpened normal distributions centered at the maximum weighted likelihood estimator *θ*ˆ*<sup>π</sup> N* :

$$
\sup_{B} \left| \Pi_{N}^{\pi}(\theta \in B|D^{(N)}) - \mathcal{N}_{\hat{\theta}_{N}^{\pi}, N^{-1}I_{\theta_{0}}^{-1}}(B) \right| = o_{\mathbb{P}}(1).
$$

4.3.3. *Inference using the Bernstein–von Mises theorem*. The Bernstein–von Mises theorem in the i.i.d. sampling models justifies the frequentist validity of the credible sets of the posterior distribution for the purpose of inference. The situation in the complex sampling setting is however more subtle. As will be clear from the discussion below, the structure of the Bernstein–von Mises theorem in Theorem [4.8](#page-21-0) shows that: (1) vanilla credible sets may not lead to valid frequentist inference procedure; (2) but at the same time suggests the construction of a corrected credible set with asymptotically valid coverage.

To see (1), suppose *CN* = *CN (D(N))* ⊂ is a (vanilla) *(*1 − *α)* credible set, that is, *<sup>π</sup> <sup>N</sup> (θ* <sup>∈</sup> *CN* <sup>|</sup>*D(N))* <sup>=</sup> <sup>1</sup> <sup>−</sup> *<sup>α</sup>*. Then by the Bernstein–von Mises theorem in Theorem [4.8,](#page-21-0) <sup>N</sup>0*,I ((NIθ*<sup>0</sup> *)*1*/*<sup>2</sup>*(CN* <sup>−</sup> *<sup>θ</sup>*¯*<sup>π</sup> <sup>N</sup> ))* <sup>→</sup> <sup>1</sup> <sup>−</sup> *<sup>α</sup>* in <sup>P</sup>-probability. Here, *<sup>θ</sup>*¯*<sup>π</sup> <sup>N</sup>* <sup>≡</sup> *<sup>θ</sup>*<sup>0</sup> <sup>+</sup> *<sup>N</sup>*−1*/*2*<sup>I</sup>* <sup>−</sup><sup>1</sup> *θ*0 G*π N* ˙ *<sup>θ</sup>*<sup>0</sup> . In other words, *CN* = *θ*¯*<sup>π</sup> <sup>N</sup>* <sup>+</sup> *<sup>I</sup>* <sup>−</sup>1*/*<sup>2</sup> *θ*0 *BN /* <sup>√</sup>*<sup>N</sup>* for some random *BN* such that <sup>N</sup>0*,I (BN )* <sup>→</sup> <sup>1</sup> <sup>−</sup> *<sup>α</sup>* in P-probability. By Proposition [3.4,](#page-8-0) we have G*<sup>π</sup> N* ˙ *<sup>θ</sup>*<sup>0</sup> →*<sup>d</sup>* N *(*0*,(*1 + *μπ*1*)Iθ*<sup>0</sup> *)*, and hence the frequentist coverage for the credible set *CN* is

$$
\mathbb{P}_{\theta_0}(\theta_0 \in C_N) = \mathbb{P}(\theta_0 \in \bar{\theta}_N^{\pi} + I_{\theta_0}^{-1/2} B_N / \sqrt{N}) = \mathbb{P}(I_{\theta_0}^{1/2} \sqrt{N} (\theta_0 - \bar{\theta}_N^{\pi}) \in B_N)
$$
  
=  $\mathbb{P}(-I_{\theta_0}^{-1/2} \mathbb{G}_N^{\pi} \dot{\ell}_{\theta_0} \in B_N) = \mathbb{E} \mathcal{N}_{0,I} (B_N / (1 + \mu_{\pi 1})^{1/2}) + o(1),$ 

which does not converge to 1 − *α* in general as *N* → ∞ as long as *μπ*<sup>1</sup> = 0.

Fortunately, the vanilla credible set *CN* can be corrected as follows. Let *θ*ˆ*<sup>π</sup> <sup>N</sup>* = arg sup*θ*∈ P*<sup>π</sup> <sup>N</sup>* log*pθ* be the maximum weighted likelihood estimator as in Section [4.2.](#page-17-0) Then under regularity conditions, *θ*ˆ*<sup>π</sup> <sup>N</sup>* <sup>=</sup> *<sup>θ</sup>*<sup>0</sup> <sup>+</sup> *<sup>N</sup>*−1*/*2*<sup>I</sup>* <sup>−</sup><sup>1</sup> *θ*0 G*π N* ˙ *<sup>θ</sup>*<sup>0</sup> + oP*(*1*)*. Now for any *CN* = *CN (D(N))* such that *<sup>π</sup> <sup>N</sup> (θ* <sup>∈</sup> *CN* <sup>|</sup>*D(N))* <sup>=</sup> <sup>1</sup> <sup>−</sup> *<sup>α</sup>*, let

(4.13) 
$$
C_N^* = \hat{\theta}_N^{\pi} + (1 + \mu_{\pi 1})^{1/2} (C_N - \hat{\theta}_N^{\pi}).
$$

Note again that *CN* = *θ*ˆ*<sup>π</sup> <sup>N</sup>* <sup>+</sup>*<sup>I</sup>* <sup>−</sup>1*/*<sup>2</sup> *θ*0 *BN /* <sup>√</sup>*<sup>N</sup>* for some random *BN* such that <sup>N</sup>0*,I (BN )* <sup>→</sup> <sup>1</sup>−*<sup>α</sup>* in P-probability, so we have

$$
\mathbb{P}_{\theta_0}(\theta_0 \in C_N^*) = \mathbb{P}_{\theta_0}(\theta_0 \in \hat{\theta}_N^{\pi} + (1 + \mu_{\pi 1})^{1/2} (C_N - \hat{\theta}_N^{\pi}))
$$
  
\n
$$
= \mathbb{P}_{\theta_0} (I_{\theta_0}^{1/2} \sqrt{N} (\theta_0 - \hat{\theta}_N^{\pi}) \in (1 + \mu_{\pi 1})^{1/2} \cdot I_{\theta_0}^{1/2} \sqrt{N} (C_N - \hat{\theta}_N^{\pi}))
$$
  
\n
$$
= \mathbb{P}_{\theta_0} (-(1 + \mu_{\pi 1})^{-1/2} I_{\theta_0}^{-1/2} \mathbb{G}_N^{\pi} \dot{\ell}_{\theta_0} \in B_N) + o(1)
$$
  
\n
$$
\to 1 - \alpha.
$$

Hence the corrected credible set *C*<sup>∗</sup> *<sup>N</sup>* (4.13) has the correct coverage.

<span id="page-23-0"></span>The construction of the corrected credible set in [\(4.13\)](#page-22-0) is generic, regardless of different sampling schemes as long as *μπ*<sup>1</sup> is known (cf. Table [1\)](#page-9-0). Such a unified construction of corrected credible sets based on the pseudo-posterior distributions could bring significant advantage for the purpose of inference in the complex sampling setting, in that it alleviates complicated bootstrap methods whose design architectures must adapt to the dependence structure in each and every different sampling schemes (e.g., [\[66,](#page-26-0) [67\]](#page-26-0) in the context of twophase sampling). In Appendix D, we present an illustrative example and simulation results in the context of one-dimensional Gaussian location model with a Gaussian prior for the phenomenon discussed above.

Finally, we remark that there are interesting recent developments in semiparametric and nonparametric Bernstein–von Mises theorems; cf. [\[19–22,](#page-24-0) [58,](#page-25-0) [59\]](#page-25-0). It is an interesting open question to extend the Bernstein–von Mises theorem and the correction scheme [\(4.13\)](#page-22-0) to these more complicated settings with complex sampling designs.

**Acknowledgments.** The authors would like to thank Thomas Lumley for several useful suggestions. Helpful and constructive comments from an Associate Editor and two referees are greatly appreciated.

The research of Q. Han was supported in part by NSF Grant DMS-1916221. The research of J. A. Wellner was supported in part by NSF Grant DMS-1566514, NI-AID grant 2R01 AI291968-04, a Simons Fellowship via the Newton Institute (INI-program STS 2018), Cambridge University and the Saw Swee Hock Visiting Professorship of Statistics at the National University of Singapore (in 2019).

## SUPPLEMENTARY MATERIAL

**Supplement: Proofs** (DOI: [10.1214/20-AOS1964SUPP;](https://doi.org/10.1214/20-AOS1964SUPP) .pdf). In the supplement, we provide proofs for the results in the main paper.

## REFERENCES

- [1] ALEXANDER, K. S. (1985). Rates of growth for weighted empirical processes. In *Proceedings of the Berkeley Conference in Honor of Jerzy Neyman and Jack Kiefer*, *Vol*. *II* (*Berkeley*, *Calif*., 1983). *Wadsworth Statist*.*/Probab*. *Ser*. 475–493. Wadsworth, Belmont, CA. [MR0822047](http://www.ams.org/mathscinet-getitem?mr=0822047)
- [2] ALEXANDER, K. S. (1987). The central limit theorem for weighted empirical processes indexed by sets. *J*. *Multivariate Anal*. **22** 313–339. [MR0899666](http://www.ams.org/mathscinet-getitem?mr=0899666) [https://doi.org/10.1016/0047-259X\(87\)90093-5](https://doi.org/10.1016/0047-259X(87)90093-5)
- [3] ALEXANDER, K. S. (1987). Rates of growth and sample moduli for weighted empirical processes indexed by sets. *Probab*. *Theory Related Fields* **75** 379–423. [MR0890285](http://www.ams.org/mathscinet-getitem?mr=0890285)<https://doi.org/10.1007/BF00318708>
- [4] BARRETT, G. F. and DONALD, S. G. (2009). Statistical inference with generalized Gini indices of inequality, poverty, and welfare. *J*. *Bus*. *Econom*. *Statist*. **27** 1–17. [MR2484980](http://www.ams.org/mathscinet-getitem?mr=2484980) [https://doi.org/10.1198/jbes.](https://doi.org/10.1198/jbes.2009.0001) [2009.0001](https://doi.org/10.1198/jbes.2009.0001)
- [5] BERGER, Y. G. (1998). Rate of convergence for asymptotic variance of the Horvitz–Thompson estimator. *J*. *Statist*. *Plann*. *Inference* **74** 149–168. [MR1665125](http://www.ams.org/mathscinet-getitem?mr=1665125) [https://doi.org/10.1016/S0378-3758\(98\)00107-4](https://doi.org/10.1016/S0378-3758(98)00107-4)
- [6] BERGER, Y. G. (1998). Rate of convergence to normal distribution for the Horvitz–Thompson estimator. *J*. *Statist*. *Plann*. *Inference* **67** 209–226. [MR1624693](http://www.ams.org/mathscinet-getitem?mr=1624693) [https://doi.org/10.1016/S0378-3758\(97\)00107-9](https://doi.org/10.1016/S0378-3758(97)00107-9)
- [7] BERTAIL, P., CHAUTRU, E. and CLÉMENÇON, S. (2017). Empirical processes in survey sampling with (conditional) Poisson designs. *Scand*. *J*. *Stat*. **44** 97–111. [MR3619696](http://www.ams.org/mathscinet-getitem?mr=3619696) [https://doi.org/10.1111/sjos.](https://doi.org/10.1111/sjos.12243) [12243](https://doi.org/10.1111/sjos.12243)
- [8] BHATTACHARYA, D. (2007). Inference on inequality from household survey data. *J*. *Econometrics* **137** 674–707. [MR2354960](http://www.ams.org/mathscinet-getitem?mr=2354960)<https://doi.org/10.1016/j.jeconom.2005.09.003>
- [9] BHATTACHARYA, D. and MAZUMDER, B. (2011). A nonparametric analysis of black–white differences in intergenerational income mobility in the United States. *Quant*. *Econ*. **2** 335–379.
- [10] BOISTARD, H., LOPUHAÄ, H. P. and RUIZ-GAZEN, A. (2012). Approximation of rejective sampling inclusion probabilities and application to high order correlations. *Electron*. *J*. *Stat*. **6** 1967–1983. [MR3020253](http://www.ams.org/mathscinet-getitem?mr=3020253)<https://doi.org/10.1214/12-EJS736>
- <span id="page-24-0"></span>[11] BOISTARD, H., LOPUHAÄ, H. P. and RUIZ-GAZEN, A. (2017). Functional central limit theorems for single-stage sampling designs. *Ann*. *Statist*. **45** 1728–1758. [MR3670194](http://www.ams.org/mathscinet-getitem?mr=3670194) [https://doi.org/10.1214/](https://doi.org/10.1214/16-AOS1507) [16-AOS1507](https://doi.org/10.1214/16-AOS1507)
- [12] BREIDT, F. J. and OPSOMER, J. D. (2000). Local polynomial regresssion estimators in survey sampling. *Ann*. *Statist*. **28** 1026–1053. [MR1810918](http://www.ams.org/mathscinet-getitem?mr=1810918)<https://doi.org/10.1214/aos/1015956706>
- [13] BRESLOW, N., MCNENEY, B. and WELLNER, J. A. (2003). Large sample theory for semiparametric regression models with two-phase, outcome dependent sampling. *Ann*. *Statist*. **31** 1110–1139. [MR2001644](http://www.ams.org/mathscinet-getitem?mr=2001644) <https://doi.org/10.1214/aos/1059655907>
- [14] BRESLOW, N. E., LUMLEY, T., BALLANTYNE, C. M., CHAMBLESS, L. E. and KULICH, M. (2009). Improved Horvitz–Thompson estimation of model parameters from two-phase stratified samples: Applications in epidemiology. *Stat*. *Biosci*. **1** 32–49.
- [15] BRESLOW, N. E., LUMLEY, T., BALLANTYNE, C. M., CHAMBLESS, L. E. and KULICH, M. (2009). Using the whole cohort in the analysis of case-cohort data. *Am*. *J*. *Epidemiol*. **169** 1398–1405.
- [16] BRESLOW, N. E. and WELLNER, J. A. (2007). Weighted likelihood for semiparametric models and twophase stratified samples, with application to Cox regression. *Scand*. *J*. *Stat*. **34** 86–102. [MR2325244](http://www.ams.org/mathscinet-getitem?mr=2325244) <https://doi.org/10.1111/j.1467-9469.2006.00523.x>
- [17] BRESLOW, N. E. and WELLNER, J. A. (2008). A *Z*-theorem with estimated nuisance parameters and correction note for: "Weighted likelihood for semiparametric models and two-phase stratified samples, with application to Cox regression" [Scand. J. Statist. **34** (2007), no. 1, 86–102; MR2325244]. *Scand*. *J*. *Stat*. **35** 186–192. [MR2391566](http://www.ams.org/mathscinet-getitem?mr=2391566)<https://doi.org/10.1111/j.1467-9469.2007.00574.x>
- [18] CARDOT, H., CHAOUCH, M., GOGA, C. and LABRUÈRE, C. (2010). Properties of design-based functional principal components analysis. *J*. *Statist*. *Plann*. *Inference* **140** 75–91. [MR2568123](http://www.ams.org/mathscinet-getitem?mr=2568123) [https://doi.org/10.](https://doi.org/10.1016/j.jspi.2009.06.012) [1016/j.jspi.2009.06.012](https://doi.org/10.1016/j.jspi.2009.06.012)
- [19] CASTILLO, I. (2012). A semiparametric Bernstein–von Mises theorem for Gaussian process priors. *Probab*. *Theory Related Fields* **152** 53–99. [MR2875753](http://www.ams.org/mathscinet-getitem?mr=2875753)<https://doi.org/10.1007/s00440-010-0316-5>
- [20] CASTILLO, I. and NICKL, R. (2013). Nonparametric Bernstein–von Mises theorems in Gaussian white noise. *Ann*. *Statist*. **41** 1999–2028. [MR3127856](http://www.ams.org/mathscinet-getitem?mr=3127856)<https://doi.org/10.1214/13-AOS1133>
- [21] CASTILLO, I. and NICKL, R. (2014). On the Bernstein–von Mises phenomenon for nonparametric Bayes procedures. *Ann*. *Statist*. **42** 1941–1969. [MR3262473](http://www.ams.org/mathscinet-getitem?mr=3262473)<https://doi.org/10.1214/14-AOS1246>
- [22] CASTILLO, I. and ROUSSEAU, J. (2015). A Bernstein–von Mises theorem for smooth functionals in semiparametric models. *Ann*. *Statist*. **43** 2353–2383. [MR3405597](http://www.ams.org/mathscinet-getitem?mr=3405597)<https://doi.org/10.1214/15-AOS1336>
- [23] CHAUVET, G. (2015). Coupling methods for multistage sampling. *Ann*. *Statist*. **43** 2484–2506. [MR3405601](http://www.ams.org/mathscinet-getitem?mr=3405601) <https://doi.org/10.1214/15-AOS1348>
- [24] CHENG, G. and HUANG, J. Z. (2010). Bootstrap consistency for general semiparametric *M*-estimation. *Ann*. *Statist*. **38** 2884–2915. [MR2722459](http://www.ams.org/mathscinet-getitem?mr=2722459)<https://doi.org/10.1214/10-AOS809>
- [25] CLÉMENÇON, S., BERTAIL, P. and PAPA, G. (2016). Learning from survey training samples: Rate bounds for Horvitz–Thompson risk minimizers. In *Asian Conference on Machine Learning* 142–157.
- [26] CONTI, P. L. (2014). On the estimation of the distribution function of a finite population under high entropy sampling designs, with applications. *Sankhya B* **76** 234–259. [MR3302272](http://www.ams.org/mathscinet-getitem?mr=3302272) [https://doi.org/10.](https://doi.org/10.1007/s13571-014-0083-x) [1007/s13571-014-0083-x](https://doi.org/10.1007/s13571-014-0083-x)
- [27] DAVIDSON, R. (2009). Reliable inference for the Gini index. *J*. *Econometrics* **150** 30–40. [MR2525992](http://www.ams.org/mathscinet-getitem?mr=2525992) <https://doi.org/10.1016/j.jeconom.2008.11.004>
- [28] DEVILLE, J.-C. and SÄRNDAL, C.-E. (1992). Calibration estimators in survey sampling. *J*. *Amer*. *Statist*. *Assoc*. **87** 376–382. [MR1173804](http://www.ams.org/mathscinet-getitem?mr=1173804)
- [29] DEVROYE, L., GYÖRFI, L. and LUGOSI, G. (1996). *A Probabilistic Theory of Pattern Recognition*. *Applications of Mathematics* (*New York*) **31**. Springer, New York. [MR1383093](http://www.ams.org/mathscinet-getitem?mr=1383093) [https://doi.org/10.1007/](https://doi.org/10.1007/978-1-4612-0711-5) [978-1-4612-0711-5](https://doi.org/10.1007/978-1-4612-0711-5)
- [30] FULLER, W. A. (2011). *Sampling Statistics* **560**. Wiley.
- [31] GHOSAL, S., GHOSH, J. K. and VAN DER VAART, A. W. (2000). Convergence rates of posterior distributions. *Ann*. *Statist*. **28** 500–531. [MR1790007](http://www.ams.org/mathscinet-getitem?mr=1790007)<https://doi.org/10.1214/aos/1016218228>
- [32] GHOSAL, S. and VAN DER VAART, A. (2007). Convergence rates of posterior distributions for non-i.i.d. observations. *Ann*. *Statist*. **35** 192–223. [MR2332274](http://www.ams.org/mathscinet-getitem?mr=2332274)<https://doi.org/10.1214/009053606000001172>
- [33] GINÉ, E. and KOLTCHINSKII, V. (2006). Concentration inequalities and asymptotic results for ratio type empirical processes. *Ann*. *Probab*. **34** 1143–1216. [MR2243881](http://www.ams.org/mathscinet-getitem?mr=2243881) [https://doi.org/10.1214/](https://doi.org/10.1214/009117906000000070) [009117906000000070](https://doi.org/10.1214/009117906000000070)
- [34] GINÉ, E., KOLTCHINSKII, V. and WELLNER, J. A. (2003). Ratio limit theorems for empirical processes. In *Stochastic Inequalities and Applications*. *Progress in Probability* **56** 249–278. Birkhäuser, Basel. [MR2073436](http://www.ams.org/mathscinet-getitem?mr=2073436)
- <span id="page-25-0"></span>[35] GINÉ, E. and NICKL, R. (2016). *Mathematical Foundations of Infinite-Dimensional Statistical Models*. *Cambridge Series in Statistical and Probabilistic Mathematics* **40**. Cambridge Univ. Press, New York. [MR3588285](http://www.ams.org/mathscinet-getitem?mr=3588285)<https://doi.org/10.1017/CBO9781107337862>
- [36] HÁJEK, J. (1961). Some extensions of the Wald–Wolfowitz–Noether theorem. *Ann*. *Math*. *Stat*. **32** 506–523. [MR0130707](http://www.ams.org/mathscinet-getitem?mr=0130707)<https://doi.org/10.1214/aoms/1177705057>
- [37] HÁJEK, J. (1964). Asymptotic theory of rejective sampling with varying probabilities from a finite population. *Ann*. *Math*. *Stat*. **35** 1491–1523. [MR0178555](http://www.ams.org/mathscinet-getitem?mr=0178555)<https://doi.org/10.1214/aoms/1177700375>
- [38] HÁJEK, J. (1981). *Sampling from a Finite Population*. *Statistics*: *Textbooks and Monographs* **37**. Dekker, New York. Edited by Václav Dupac, With a foreword by P. K. Sen. ˇ [MR0627744](http://www.ams.org/mathscinet-getitem?mr=0627744)
- [39] HAN, Q. (2017). Bayes model selection. arXiv preprint, [arXiv:1704.07513](http://arxiv.org/abs/arXiv:1704.07513).
- [40] HAN, Q. and WELLNER, J. A. (2019). Convergence rates of least squares regression estimators with heavytailed errors. *Ann*. *Statist*. **47** 2286–2319. [MR3953452](http://www.ams.org/mathscinet-getitem?mr=3953452)<https://doi.org/10.1214/18-AOS1748>
- [41] HAN, Q. and WELLNER, J. A. (2021). Supplement to "Complex sampling designs: Uniform limit theorems and applications." <https://doi.org/10.1214/20-AOS1964SUPP>
- [42] HARTLEY, H. O. (1962). Multiple frame surveys. In *Proceedings of the Social Statistics Section* **19** 203– 206. American Statistical Association, Washington, DC.
- [43] HARTLEY, H. O. (1974). Multiple frame methodology and selected applications. *Sankhya¯* **36** 118.
- [44] HORVITZ, D. G. and THOMPSON, D. J. (1952). A generalization of sampling without replacement from a finite universe. *J*. *Amer*. *Statist*. *Assoc*. **47** 663–685. [MR0053460](http://www.ams.org/mathscinet-getitem?mr=0053460)
- [45] HUANG, J. (1996). Efficient estimation for the proportional hazards model with interval censoring. *Ann*. *Statist*. **24** 540–568. [MR1394975](http://www.ams.org/mathscinet-getitem?mr=1394975)<https://doi.org/10.1214/aos/1032894452>
- [46] KOLTCHINSKII, V. (2006). Local Rademacher complexities and oracle inequalities in risk minimization. *Ann*. *Statist*. **34** 2593–2656. [MR2329442](http://www.ams.org/mathscinet-getitem?mr=2329442)<https://doi.org/10.1214/009053606000001019>
- [47] KOLTCHINSKII, V. (2011). *Oracle Inequalities in Empirical Risk Minimization and Sparse Recovery Problems*. *Lecture Notes in Math*. **2033**. Springer, Heidelberg. Lectures from the 38th Probability Summer School held in Saint-Flour, 2008, École d'Été de Probabilités de Saint-Flour. [Saint-Flour Probability Summer School]. [MR2829871](http://www.ams.org/mathscinet-getitem?mr=2829871)<https://doi.org/10.1007/978-3-642-22147-7>
- [48] KOSOROK, M. R. (2008). *Introduction to Empirical Processes and Semiparametric Inference*. *Springer Series in Statistics*. Springer, New York. [MR2724368](http://www.ams.org/mathscinet-getitem?mr=2724368)<https://doi.org/10.1007/978-0-387-74978-5>
- [49] LEÓN-NOVELO, L. G. and SAVITSKY, T. D. (2019). Fully Bayesian estimation under informative sampling. *Electron*. *J*. *Stat*. **13** 1608–1645. [MR3939589](http://www.ams.org/mathscinet-getitem?mr=3939589)<https://doi.org/10.1214/19-ejs1538>
- [50] LIN, D. Y. (2000). On fitting Cox's proportional hazards models to survey data. *Biometrika* **87** 37–47. [MR1766826](http://www.ams.org/mathscinet-getitem?mr=1766826)<https://doi.org/10.1093/biomet/87.1.37>
- [51] LOHR, S. and RAO, J. N. K. (2006). Estimation in multiple-frame surveys. *J*. *Amer*. *Statist*. *Assoc*. **101** 1019–1030. [MR2324141](http://www.ams.org/mathscinet-getitem?mr=2324141)<https://doi.org/10.1198/016214506000000195>
- [52] LUMLEY, T., SHAW, P. A. and DAI, J. Y. (2011). Connections between survey calibration estimators and semiparametric models for incomplete data. *Int*. *Stat*. *Rev*. **79** 200–220. [https://doi.org/10.1111/](https://doi.org/10.1111/j.1751-5823.2011.00138.x) [j.1751-5823.2011.00138.x](https://doi.org/10.1111/j.1751-5823.2011.00138.x)
- [53] MA, S. and KOSOROK, M. R. (2005). Robust semiparametric M-estimation and the weighted bootstrap. *J*. *Multivariate Anal*. **96** 190–217. [MR2202406](http://www.ams.org/mathscinet-getitem?mr=2202406)<https://doi.org/10.1016/j.jmva.2004.09.008>
- [54] MAMMEN, E. and TSYBAKOV, A. B. (1999). Smooth discrimination analysis. *Ann*. *Statist*. **27** 1808–1829. [MR1765618](http://www.ams.org/mathscinet-getitem?mr=1765618)<https://doi.org/10.1214/aos/1017939240>
- [55] MASON, D. M., SHORACK, G. R. and WELLNER, J. A. (1983). Strong limit theorems for oscillation moduli of the uniform empirical process. *Z*. *Wahrsch*. *Verw*. *Gebiete* **65** 83–97. [MR0717935](http://www.ams.org/mathscinet-getitem?mr=0717935) <https://doi.org/10.1007/BF00534996>
- [56] NAN, B., KALBFLEISCH, J. D. and YU, M. (2009). Asymptotic theory for the semiparametric accelerated failure time model with missing data. *Ann*. *Statist*. **37** 2351–2376. [MR2543695](http://www.ams.org/mathscinet-getitem?mr=2543695) [https://doi.org/10.1214/](https://doi.org/10.1214/08-AOS657) [08-AOS657](https://doi.org/10.1214/08-AOS657)
- [57] NAN, B. and WELLNER, J. A. (2013). A general semiparametric Z-estimation approach for case-cohort studies. *Statist*. *Sinica* **23** 1155–1180. [MR3114709](http://www.ams.org/mathscinet-getitem?mr=3114709)
- [58] NICKL, R. (2017). Bernstein–von Mises theorems for statistical inverse problems I: Schrödinger equation, arXiv preprint [arXiv:1707.01764](http://arxiv.org/abs/arXiv:1707.01764).
- [59] NICKL, R. and SÖHL, J. (2019). Bernstein–von Mises theorems for statistical inverse problems II: Compound Poisson processes. *Electron*. *J*. *Stat*. **13** 3513–3571. [MR4013745](http://www.ams.org/mathscinet-getitem?mr=4013745) [https://doi.org/10.1214/](https://doi.org/10.1214/19-ejs1609) [19-ejs1609](https://doi.org/10.1214/19-ejs1609)
- [60] PRÆSTGAARD, J. and WELLNER, J. A. (1993). Exchangeably weighted bootstraps of the general empirical process. *Ann*. *Probab*. **21** 2053–2086. [MR1245301](http://www.ams.org/mathscinet-getitem?mr=1245301)
- [61] ROSÉN, B. (1965). Limit theorems for sampling from finite populations. *Ark*. *Mat*. **5** 383–424. [MR0177437](http://www.ams.org/mathscinet-getitem?mr=0177437) <https://doi.org/10.1007/BF02591138>
- <span id="page-26-0"></span>[62] ROSÉN, B. (1967). On the central limit theorem for a class of sampling procedures. *Z*. *Wahrsch*. *Verw*. *Gebiete* **7** 103–115. [MR0210181](http://www.ams.org/mathscinet-getitem?mr=0210181)<https://doi.org/10.1007/BF00536324>
- [63] ROSÉN, B. (1972). Asymptotic theory for successive sampling with varying probabilities without replacement. I, II. *Ann*. *Math*. *Stat*. **43** 373–397; ibid. 43 (1972), 748–776. [MR0321223](http://www.ams.org/mathscinet-getitem?mr=0321223) [https://doi.org/10.](https://doi.org/10.1214/aoms/1177692620) [1214/aoms/1177692620](https://doi.org/10.1214/aoms/1177692620)
- [64] ROSÉN, B. (1974). Asymptotic theory for Des Raj's estimator. I, II. *Scand*. *J*. *Stat*. **1** 71–83; ibid. 1 (1974), no. 3, 135–144. [MR0375560](http://www.ams.org/mathscinet-getitem?mr=0375560)
- [65] RUBIN-BLEUER, S. and SCHIOPU KRATINA, I. (2005). On the two-phase framework for joint model and design-based inference. *Ann*. *Statist*. **33** 2789–2810. [MR2253102](http://www.ams.org/mathscinet-getitem?mr=2253102) [https://doi.org/10.1214/](https://doi.org/10.1214/009053605000000651) [009053605000000651](https://doi.org/10.1214/009053605000000651)
- [66] SAEGUSA, T. (2014). Bootstrapping two-phase sampling. arXiv preprint, [arXiv:1406.5580](http://arxiv.org/abs/arXiv:1406.5580).
- [67] SAEGUSA, T. (2015). Variance estimation under two-phase sampling. *Scand*. *J*. *Stat*. **42** 1078–1091. [MR3426311](http://www.ams.org/mathscinet-getitem?mr=3426311)<https://doi.org/10.1111/sjos.12152>
- [68] SAEGUSA, T. (2019). Large sample theory for merged data from multiple sources. *Ann*. *Statist*. **47** 1585– 1615. [MR3911123](http://www.ams.org/mathscinet-getitem?mr=3911123)<https://doi.org/10.1214/18-AOS1727>
- [69] SAEGUSA, T. and WELLNER, J. A. (2013). Weighted likelihood estimation under two-phase sampling. *Ann*. *Statist*. **41** 269–295. [MR3059418](http://www.ams.org/mathscinet-getitem?mr=3059418)<https://doi.org/10.1214/12-AOS1073>
- [70] SÄRNDAL, C.-E., SWENSSON, B. and WRETMAN, J. (1992). *Model Assisted Survey Sampling*. *Springer Series in Statistics*. Springer, New York. [MR1140409](http://www.ams.org/mathscinet-getitem?mr=1140409)<https://doi.org/10.1007/978-1-4612-4378-6>
- [71] SAVITSKY, T. D. and TOTH, D. (2016). Bayesian estimation under informative sampling. *Electron*. *J*. *Stat*. **10** 1677–1708. [MR3522657](http://www.ams.org/mathscinet-getitem?mr=3522657)<https://doi.org/10.1214/16-EJS1153>
- [72] SHORACK, G. R. (1973). Convergence of reduced empirical and quantile processes with application to functions of order statistics in the non-I.I.D. case. *Ann*. *Statist*. **1** 146–152. [MR0336776](http://www.ams.org/mathscinet-getitem?mr=0336776)
- [73] SHORACK, G. R. and WELLNER, J. A. (1982). Limit theorems and inequalities for the uniform empirical process indexed by intervals. *Ann*. *Probab*. **10** 639–652. [MR0659534](http://www.ams.org/mathscinet-getitem?mr=0659534)
- [74] STUTE, W. (1982). The oscillation behavior of empirical processes. *Ann*. *Probab*. **10** 86–107. [MR0637378](http://www.ams.org/mathscinet-getitem?mr=0637378)
- [75] STUTE, W. (1984). The oscillation behavior of empirical processes: The multivariate case. *Ann*. *Probab*. **12** 361–379. [MR0735843](http://www.ams.org/mathscinet-getitem?mr=0735843)
- [76] TSYBAKOV, A. B. (2004). Optimal aggregation of classifiers in statistical learning. *Ann*. *Statist*. **32** 135–166. [MR2051002](http://www.ams.org/mathscinet-getitem?mr=2051002)<https://doi.org/10.1214/aos/1079120131>
- [77] VAN DE GEER, S. A. (2000). *Applications of Empirical Process Theory*. *Cambridge Series in Statistical and Probabilistic Mathematics* **6**. Cambridge Univ. Press, Cambridge. [MR1739079](http://www.ams.org/mathscinet-getitem?mr=1739079)
- [78] VAN DER VAART, A. (2002). Semiparametric statistics. In *Lectures on Probability Theory and Statistics* (*Saint-Flour*, 1999). *Lecture Notes in Math*. **1781** 331–457. Springer, Berlin. [MR1915446](http://www.ams.org/mathscinet-getitem?mr=1915446)
- [79] VAN DER VAART, A. W. (1995). Efficiency of infinite-dimensional *M*-estimators. *Stat*. *Neerl*. **49** 9–30. [MR1333176](http://www.ams.org/mathscinet-getitem?mr=1333176)<https://doi.org/10.1111/j.1467-9574.1995.tb01452.x>
- [80] VAN DER VAART, A. W. (1998). *Asymptotic Statistics*. *Cambridge Series in Statistical and Probabilistic Mathematics* **3**. Cambridge Univ. Press, Cambridge. [MR1652247](http://www.ams.org/mathscinet-getitem?mr=1652247) [https://doi.org/10.1017/](https://doi.org/10.1017/CBO9780511802256) [CBO9780511802256](https://doi.org/10.1017/CBO9780511802256)
- [81] VAN DER VAART, A. W. and WELLNER, J. A. (1996). *Weak Convergence and Empirical Processes*: *With Applications to Statistics*. *Springer Series in Statistics*. Springer, New York. [MR1385671](http://www.ams.org/mathscinet-getitem?mr=1385671) <https://doi.org/10.1007/978-1-4757-2545-2>
- [82] VÍŠEK, J. Á. (1979). Asymptotic distribution of simple estimate for rejective, Sampford and successive sampling. In *Contributions to Statistics* 263–275. Reidel, Dordrecht. [MR0561274](http://www.ams.org/mathscinet-getitem?mr=0561274)
- [83] WELLNER, J. A. (1978). Limit theorems for the ratio of the empirical distribution function to the true distribution function. *Z*. *Wahrsch*. *Verw*. *Gebiete* **45** 73–88. [MR0651392](http://www.ams.org/mathscinet-getitem?mr=0651392) [https://doi.org/10.1007/](https://doi.org/10.1007/BF00635964) [BF00635964](https://doi.org/10.1007/BF00635964)
- [84] WELLNER, J. A. and ZHAN, Y. (1996). Bootstrapping Z-estimators. Technical Report 308, Univ. Washington Dept. Statistics.
- [85] WELLNER, J. A. and ZHANG, Y. (2007). Two likelihood-based semiparametric estimation methods for panel count data with covariates. *Ann*. *Statist*. **35** 2106–2142. [MR2363965](http://www.ams.org/mathscinet-getitem?mr=2363965) [https://doi.org/10.1214/](https://doi.org/10.1214/009053607000000181) [009053607000000181](https://doi.org/10.1214/009053607000000181)