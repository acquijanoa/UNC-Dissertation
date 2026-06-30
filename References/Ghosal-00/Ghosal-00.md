## CONVERGENCE RATES OF POSTERIOR DISTRIBUTIONS

By Subhashis Ghosal, Jayanta K. Ghosh and Aad W. van der Vaart

Free University Amsterdam, Indian Statistical Institute and Free University Amsterdam

We consider the asymptotic behavior of posterior distributions and Bayes estimators for infinite-dimensional statistical models. We give general results on the rate of convergence of the posterior measure. These are applied to several examples, including priors on finite sieves, log-spline models, Dirichlet processes and interval censoring.

1. Introduction. Suppose that we observe a random sample X1--Xn from a distribution P with density p relative to some reference measure on the sample space - - -. The unknown distribution is known to belong to some model , a set of probability measures on the sample space. Given some prior distribution n on the set , the posterior distribution is the random measure given by

$$
\Pi_n(B|X_1,\ldots,X_n)=\frac{\int_B \prod_{i=1}^n p(X_i) d\Pi_n(P)}{\int \prod_{i=1}^n p(X_i) d\Pi_n(P)}.
$$

If the distribution P is considered random and distributed according to , as it is in Bayesian inference, then the posterior distribution is the conditional distribution of P given the observations. The prior is, of course, a measure on some σ-field on and we must assume that the expressions in the display are well defined. In particular, we assume that the map x p- → px is measurable for the product σ-field on × . It will be silently understood that the sets of which we compute prior or posterior measures are measurable.

In this paper we study the frequentist properties of the posterior distribution as n → ∞, assuming that the observations are a random sample from some fixed measure P0. In particular, we study the rate at which this random distribution converges to P0. The posterior is said to be consistent if, as a random measure, it concentrates on arbitrarily small neighborhoods of P0, with probability tending to 1 or almost surely, as n → ∞. We study the rate at which such neighborhoods may decrease to zero meanwhile still capturing most of the posterior mass.

If = Pθ θ ∈  is parametrized by a parameter θ, then usually the prior is constructed by putting a measure on the parameter set . If is a

Received July 1998; revised December 1999.

AMS 1991 subject classifications. 62G15, 62G20, 62F25

Key words and phrases. Infinite dimensional model, posterior distribution, rate of convergence, sieves, splines

subset of a finite-dimensional Euclidean space and the dependence of θ → Pθ is sufficiently regular, then it is well known that the posterior distribution of θ achieves the optimal rate of convergence, as n → ∞ [see, for example, Le Cam (1973) and Ibragimov and Has'minskii (1981)]. In particular, if the model θ → Pθ is suitably differentiable, then the rate of the posterior mean for <sup>θ</sup> is <sup>√</sup><sup>n</sup> and the posterior distribution, when rescaled, tends to a normal distribution with covariance the inverse Fisher information, according to the Bernstein–von Mises theorem. In that case the posterior expectation is an asymptotically efficient estimator for the parameter under some integrability conditions.

Much less is known on the behavior of posterior distributions for infinitedimensional models. Most of the known results in this area address consistency issues. A famous theorem by Doob (1949) shows that consistency obtains on a set of prior measure 1, but his result concludes nothing on consistency at a particular true distribution of interest. Schwartz (1965) gives results that do apply to a particular true distribution. She shows that the posterior distribution is consistent if the true distribution P<sup>0</sup> can be suitably tested versus the complements of neighborhoods of P<sup>0</sup> and Kullback–Leibler neighborhoods of P<sup>0</sup> receive positive probabilities under the prior. Examples by, among others, Freedman (1963, 1965) and Diaconis and Freedman (1986) show that the situation is more complicated, even though, perhaps, these examples put too much emphasis on the situations where Bayes estimation does not work. A number of recent papers consider the consistency with a particular interest in the infinite-dimensional situation. Barron, in an unpublished paper, refines Schwartz's theorem [see Proposition 2 of Barron, Schervish and Wasserman (1999)] in a way that is particularly suitable for prior measures on infinitedimensional spaces of densities. Ghosal, Ghosh and Ramamoorthi (1999a) use this extension to study consistency in the variation distance for Dirichlet mixture priors. For reviews on posterior consistency in infinite dimensions, see Ghosal, Ghosh and Ramamoorthi (1999b) and Wasserman (1998).

Le Cam [(1986), pages 509–529], addresses rates of convergence of Bayes estimators in an abstract setting. Our methods are clearly related to the methods used by Le Cam. A crucial distinction appears to be that Le Cam appears to base his argument on the prior mass present in fairly small balls (the sets V in his Lemma 1 on page 510, later chosen such that PV is close to P<sup>n</sup> 0 ), whereas our result is based on having sufficient prior mass in balls of radius equal to the rate of convergence that we wish to obtain. The behavior of product densities in these bigger balls appears not to be determined by the metric distance of the marginal components alone. Instead we use a combination of Kullback–Leibler numbers and distances on the log likelihood ratio ratios. Another distinction is that we consider rates of posterior measures, whereas Le Cam considers the rate of "formal Bayes estimators." For these reasons our results appear not to be covered by Le Cam's Theorem 1 on page 513 (in which distances and other quantities are in terms of the product measures P<sup>n</sup>). However, we would like to acknowledge the great importance of Le Cam's work to the present paper. In particular, we are indebted to the part of Le Cam's work as it was extended by Birge (1983) to general metric spaces. ´

Shortly after completing this paper, we learned of independent work by Shen and Wasserman (1998), who also address rates of convergence.

The construction of prior measures on infinite-dimensional models is not a trivial matter and has also received recent attention. This development started with the introduction of Dirichlet processes by Ferguson (1973, 1974). Given computing algorithms such as Markov chain Monte Carlo methods and powerful computing machines, implementation of Bayesian methods has now become feasible even for many complicated priors and infinite dimensional models.

In Section 2 we present a main result and several variations concerning the rate of convergence of the posterior relative to the total variation, Hellinger and L2-metrics. Every time the two main elements characterizing the rate of convergence are the size of the model (measured by covering numbers or existence of certain tests) and the amount of prior mass given to a shrinking ball around the true measure. Actually, the size of the model comes in only to guarantee the existence of certain tests of the true measure versus the complement of a shrinking ball around it, and conditions can be put in terms of such tests instead. Conditions of this form go back to Schwartz (1965) and Le Cam (1973). We discuss testing in Section 7, and reformulate our main result in terms of tests in this section. The proofs of the main results are contained in Section 8 following the discussion of the existence of tests. In Section 2 we also note that a rate of convergence for the posterior automatically entails the existence of point estimators with the same rate.

We apply the general result to several examples. In Section 3 we consider discrete priors constructed on ε-nets over the model. In Section 4 we discuss Bayes estimators based on the log spline models for density estimation discussed by Stone (1986). In Section 6 we consider finite-dimensional models. In Section 6 we discuss applications to Dirichlet priors.

The notation is used to denote inequality up to a universal multiplicative constant, or up to a constant that is fixed throughout. We define the Hellinger distance hp q or hP- Q between two probability densities or measures by the <sup>L</sup>2-distance between the root densities <sup>√</sup><sup>p</sup> and <sup>√</sup>q. The total variation distance is the L1-distance. (Some authors define these distances with an additional factor 1/2.)

2. Main results. Let X1--Xn be distributed according to some distribution P<sup>0</sup> and let n be a sequence of prior probability measures supported on some set of probability measures . Let d be either the variation or the Hellinger metric on . If the set of densities is uniformly bounded, then we may also choose d equal to the L2-distance. This metric is used in condition (2.2) of the following theorem and also in its assertion.

Let Dε- d denote the ε-packing number of . This is the maximal number of points in such that the distance between every pair is at least ε. It is easy to see that this is related to the ε-covering number Nε- d-, which is the minimal number of balls of radius ε needed to cover , by the inequalities

(2.1) 
$$
N(\varepsilon, \mathscr{P}, d) \le D(\varepsilon, \mathscr{P}, d) \le N(\varepsilon/2, \mathscr{P}, d).
$$

Because we are only interested in rates of convergence, the additional constant 2 is of no real importance in the following, and covering numbers may replace packing numbers throughout. The set of centers of a minimal set of balls of radius ε covering is called an ε-net.

A good early reference on entropy numbers is the paper Kolmogorov and Tikhomirov (1961). Alternative references are Dudley (1984) and van der Vaart and Wellner (1996).

We use the notation Pf to abbreviate f dP, and, later on, nf for n<sup>−</sup><sup>1</sup> <sup>n</sup> <sup>i</sup>=<sup>1</sup> fXi-.

Theorem 2.1. Suppose that for a sequence εn with εn <sup>→</sup> <sup>0</sup> and nε<sup>2</sup> <sup>n</sup> → ∞, a constant C > 0 and sets <sup>n</sup> ⊂ , we have

(2.2) 
$$
\log D(\varepsilon_n, \mathscr{P}_n, d) \leq n \varepsilon_n^2,
$$

(2.3) 
$$
\Pi_n(\mathscr{P}\setminus\mathscr{P}_n)\leq \exp\bigl(-n\,\varepsilon_n^2(C+4)\bigr),
$$

$$
(2.4) \qquad \Pi_n\Big(P\colon -P_0\Big(\log \frac{p}{p_0}\Big)\leq \varepsilon_n^2,\ P_0\Big(\log \frac{p}{p_0}\Big)^2\leq \varepsilon_n^2\Big)\geq \exp\big(\neg n\varepsilon_n^2 C\big).
$$

Then for sufficiently large M, we have that nP dP- P0- ≥ MεnX1--Xn-<sup>→</sup> <sup>0</sup> in <sup>P</sup><sup>n</sup> <sup>0</sup> -probability.

The first and third conditions of the theorem are the essential ones. Condition (2.3) allows some additional flexibility, but should first be understood as expressing that <sup>n</sup> is almost the support of the prior (in which case its left side is zero and the condition is trivially satisfied).

Condition (2.2) requires that the "model" <sup>n</sup> be not too big. It is true for every ε <sup>n</sup> ≥ εn as soon as it is true for εn and can thus be seen as defining a minimal possible value of εn. Condition (2.2) ensures the existence of certain tests, as discussed in Section 7, and could be replaced by a testing condition. Note that the metric d used here reappears in the assertion of the theorem. Since the total variation metric is bounded above by twice the Hellinger metric, the assertion of the theorem using the Hellinger metric is stronger, but also condition (2.2) will be more restrictive, so that we really have two theorems. In the case that the densities are uniformly bounded, we even have a third theorem, when using the L2-distance, which in that case will be bounded above by a multiple of the Hellinger distance. If the densities are also uniformly bounded and uniformly bounded away from zero, then these three distances are equivalent and are also equivalent to the Kullback–Leibler number and L2-norm appearing in condition (2.4). See, for example, Lemmas 8.2 and 8.3 and (8.6).

A rate εn satisfying (2.2) for = <sup>n</sup> and d the Hellinger metric is often viewed as giving the "optimal" rate of convergence for estimators of P relative to the Hellinger metric, given the model . Under certain conditions, such as likelihood ratios bounded away from zero and infinity, this is proved as a theorem by Birge (1983) and Le Cam (1973, 1986). From Birg ´ e's work it is clear ´ that condition (2.2) is the correct expression of the complexity of the model, as relating to estimating the true density relative to the Hellinger distance, if this is to be given in terms of metric entropy. A weaker, but more involved, condition is in terms of the existence of certain tests. We give a generalization of the theorem using tests in Section 7.

Condition (2.4) is the other main determinant of the posterior rate given by the theorem. It requires that the prior measures put a sufficient amount of mass near the true measure P0. Here "near" is measured through a combination of the Kullback–Leibler divergence of p and p<sup>0</sup> and the L2P0--norm of logp/p0-. Again this condition is satisfied for ε <sup>n</sup> ≥ εn if it is satisfied for εn and thus is another restriction on a minimal value of εn. The form of this condition can be motivated from entropy considerations. Suppose that we wish to satisfy (2.4) for the minimal εn satisfying (2.2) with <sup>n</sup> = , that is, for the optimal rate of convergence for the model. Furthermore, for the sake of the argument assume that all distances used are equivalent. Then a minimal εn-cover of consists of expnε<sup>2</sup> n balls. If the prior n would spread its mass uniformly over , then every ball would obtain mass approximately exp−Cnε<sup>2</sup> n-. (The constant C expresses the constants in comparing the distances and the fact that the balls of radius εn may overlap.) On the other hand, if n is not "uniform," then we should expect (2.4) to fail for some P<sup>0</sup> ∈ . Here we must admit that "uniform" priors do not exist in infinite-dimensional models and actually condition (2.4) is stronger than needed and will be improved ahead in Theorem 2.4. However, a rough implication of the condition is that n should be "uniformly spread" in order for the posterior distribution to attain the optimal rate of convergence.

Condition (2.3), combined with (2.2), can be interpreted as saying that a part of that barely receives prior mass need not be small. The sets <sup>n</sup> may be thought of as "sieves" approximating the parameter space, which capture most of the prior probability. This type of condition has received much attention in the discussion of consistency issues [see Barron, Schervish and Wasserman (1998)], but plays a smaller role in the present paper. Of course, condition (2.3) is trivially satisfied for <sup>n</sup> = ; we can make this choice if condition (2.2) holds with <sup>n</sup> = itself.

The assertion of the theorem is an in-probability statement that the posterior mass outside a large ball of radius proportional to εn is approximately zero. The in-probability statement can be improved to an almost sure assertion, but under stronger conditions. We present two results.

Let h be the Hellinger distance and write log<sup>+</sup> x for log x-∨ 0.

Theorem 2.2. Suppose that conditions (2.2) and (2.3) hold as in the preceding theorem and in addition <sup>n</sup> exp−Bnε<sup>2</sup> n-< ∞ for every B > 0 and

(2.5) 
$$
\Pi_n\Big(P\colon h^2(P, P_0)\Big\|\frac{p_0}{p}\Big\|_{\infty}\leq \varepsilon_n^2\Big)\geq \exp\big(\frac{-n\varepsilon_n^2C}{\varepsilon_n^2}\big).
$$

Then for sufficiently large M, we have that nP dP- P0- ≥ MεnX1--Xn-<sup>→</sup> <sup>0</sup> in <sup>P</sup><sup>n</sup> <sup>0</sup> -almost surely.

Theorem 2.3. Suppose that conditions (2.2) and (2.3) hold as in the preceding theorem and in addition <sup>n</sup> exp <sup>−</sup> Bnε<sup>2</sup> n < ∞ for every B > 0 and for a given function m with P0m < ∞,

$$
\Pi_n(P: 18h^2(P, P_0) \Big( \log_+ (\sqrt{P_0 m}/h(P, P_0)) + \Phi^{-1}(h^2(P, P_0)) \Big) \le \varepsilon_n^2, \frac{p_0}{p} \le m \Big) \ge \exp(-n\varepsilon_n^2 C),
$$

where <sup>−</sup>1ε- = sup M M- ≥ ε is the inverse of the function M- = P0m1 m ≥ M/M. Then for sufficiently large M, we have that nP dP- P0-≥ MεnX1--Xn-<sup>→</sup> <sup>0</sup> in <sup>P</sup><sup>n</sup> <sup>0</sup> -almost surely.

If the quotients p0/p are uniformly bounded, then condition (2.5) simply requires that shrinking Hellinger balls possess a sufficient amount of prior mass. Then a fairly symmetric statement is obtained when combined with condition (2.2) for the Hellinger metric d: if we can cover the model with not too many Hellinger balls and the Hellinger ball around P<sup>0</sup> contains a sufficient amount of mass, then the rate of convergence relative to the Hellinger distance is εn.

Lemmas 8.2 and 8.3 in Section 8 relate the Kullback–Leibler divergence and L2-norm of logp/p0 to <sup>h</sup>2P- P0- p0/p <sup>∞</sup> and imply that the conditions of Theorem 2.2 are essentially stronger than those of Theorem 2.1.

Condition (2.6) is milder in its control of p/p<sup>0</sup> than (2.5) by allowing a general bound m that need only satisfy a moment condition. However, in comparison with (2.5) it will be satisfied only for somewhat bigger εn, due to the presence of the term involving log and <sup>−</sup>1.

In general, good control on the quotients p/p<sup>0</sup> is needed next to the closeness of p to p<sup>0</sup> relative to, for example, the Hellinger metric, because the product measures P<sup>n</sup> and P<sup>n</sup> <sup>0</sup> can be arbitrarily far apart as n → ∞ within balls of radii εn, for the values of εn bigger than n<sup>−</sup>1/<sup>2</sup> that we are considering here. The bound on p/p<sup>0</sup> together with the distance ensures that P<sup>n</sup> and P<sup>n</sup> 0 are still "close" enough on an exponential scale. Only prior mass on such close alternatives helps to increase the rate of convergence of the posterior.

One deficit of the theorems as presented so far is that they do not satisfactorily cover finite-dimensional models. When applied to such models, they would yield the rate 1/ <sup>√</sup><sup>n</sup> times a logarithmic factor rather than 1/ <sup>√</sup><sup>n</sup> itself. Similarly, the theorems may also yield unnecessary logarithmic factors when applied to priors constructed on a sequence of finite-dimensional sieves. To improve this situation we must refine both the entropy condition (2.2) and the prior mass condition (2.4). The following generalization of Theorem 2.1 is more complicated but does yield the right result in the finite-dimensional situation. It is essential for our examples using spline approximations in Section 4. Theorems 2.2 and 2.3 can be generalized similarly. Let

$$
B_n(\varepsilon) = \Big\{ P: \ -P_0\Big(\log \frac{p}{p_0}\Big) \leq \varepsilon^2, \ \ P_0\Big(\log \frac{p}{p_0}\Big)^2 \leq \varepsilon^2 \Big\}.
$$

Theorem 2.4. Suppose that for a sequence εn with εn <sup>→</sup> <sup>0</sup> and such that nε<sup>2</sup> <sup>n</sup> is bounded away from zero, every sufficiently large j and sets <sup>n</sup> ⊂ , we have

$$
(2.7) \ \log D\Big(\frac{\varepsilon}{2}, \{P \in \mathscr{P}_n : \varepsilon \leq d(P, P_0) \leq 2\varepsilon\}, d\Big) \leq n\varepsilon_n^2 \quad \text{for every } \varepsilon \geq \varepsilon_n,
$$

(2.8) 
$$
\frac{\Pi_n(\mathscr{P} - \mathscr{P}_n)}{\Pi_n(B_n(\varepsilon_n))} = o(\exp(-2n\varepsilon_n^2)),
$$

(2.9) 
$$
\frac{\Pi_n(P: j\varepsilon_n \lt d(P, P_0) \leq 2j\varepsilon_n)}{\Pi_n(B_n(\varepsilon_n))} \leq \exp(Kn\varepsilon_n^2 j^2/2).
$$

Here K is the universal testing constant appearing in (7.1) and (7.2). Then for every Mn → ∞, we have that nP dP- P0- ≥ MnεnX1--Xn- → 0 in P<sup>n</sup> <sup>0</sup> -probability.

Convergence of the posterior distribution at the rate εn implies the existence of point estimators, which are Bayes in that they are based on the posterior distribution, that converge at least as fast as εn in the frequentist sense. One possible construction is to define Pˆ <sup>n</sup> as the (near) maximizer of

$$
Q \mapsto \Pi_n(P: d(P, Q) < \varepsilon_n | X_1, \ldots, X_n).
$$

Theorem 2.5. Suppose that n<sup>P</sup> dP- P0- ≥ εnX1--Xn converges to 0, almost surely (respectively, in probability) under P<sup>n</sup> <sup>0</sup> and let Pˆ <sup>n</sup> maximize, up to o1-, the function Q → n P dP- Q- < εnX1--Xn . Then dPˆ n- P0- ≤ 2εn eventually almost surely (respectively, in probability) under P<sup>n</sup> 0 .

Proof. By definition, the εn-ball around <sup>P</sup><sup>ˆ</sup> <sup>n</sup> contains at least as much posterior probability as the εn-ball around P0. The latter, by posterior convergence at rate εn, has posterior probability close to unity. Therefore, these two balls cannot be disjoint, for otherwise, the total posterior mass would exceed unity. Now apply the triangle inequality. ✷

The preceding construction actually applies to general statistical models and posterior distributions, and the theorem is well-known. [See, e.g., Le Cam (1986) or Le Cam and Yang (1990).] If we use the Hellinger or total variation metric (or some other bounded metric whose square is convex), then an alternative is to use the posterior expectation, which typically has a similar property. By Jensen's inequality and convexity of <sup>P</sup> → <sup>d</sup><sup>2</sup>P- P0-,

$$
d^2\Big(\int P d\Pi_n(P|X_1,\ldots,X_n),P_0\Big) \leq \int d^2(P,P_0) d\Pi_n(P|X_1,\ldots,X_n)
$$
  

$$
\leq \varepsilon_n^2 + d_\infty^2 \Pi_n(P: d(P,P_0) > \varepsilon_n | X_1,\ldots,X_n),
$$

where <sup>d</sup><sup>∞</sup> is a bound on the maximal distance (<sup>√</sup> 2 and 2, respectively, for Hellinger and variation distance). To obtain the desired result, we now need that the posterior probability of the complement of the εn-ball around p<sup>0</sup> converges to zero at least at the order ε<sup>2</sup> <sup>n</sup>. This is usually the case, in particular under the conditions of Theorems 2.2 and 2.3, whose proofs yield the exponential order exp−Bnε<sup>2</sup> n-. (We use the square of the distance, because the Hellinger distance is not convex. With the total variation distance the argument would work also with the distance itself.)

More generally, we could use the minimizer Pˆ <sup>n</sup> of

$$
Q \mapsto \int \ell_n(d(Q, P)) d\Pi_n(P|X_1, \ldots, X_n),
$$

for appropriate loss functions "n. Such estimators are called formal Bayes estimators in Le Cam (1986).

On the one hand, Theorem 2.5 shows that we can construct good estimators from the posterior if the posterior converges at a good rate. On the other hand, it shows that the posterior cannot converge at a rate faster than the optimal rate of convergence for point estimators. We use this argument in a number of examples to show that the posterior converges at the best possible rate. Of course, our arguments have nothing to say about the best possible constants. Furthermore, for many priors the rate may be suboptimal.

3. Priors based on finite approximating sets. In this section, we construct, under bracketing entropy conditions, priors based on uniform distributions on carefully chosen finite sets for which the posterior converges at the best possible rate. Priors based on uniform distributions on finite subsets are introduced by Ghosal, Ghosh and Ramamoorthi (1997) as the Bayesian default priors for nonparametric problems. They establish posterior consistency for such priors under mild entropy conditions. In the present case, the prior is constructed more carefully to achieve the optimal rate of convergence as well.

Given two functions l u → the bracket l u is defined as the set of all functions f → such that l ≤ f ≤ u everywhere. The bracket is said to be of size ε relative to the distance d if dl u- < ε. In this section we use the Hellinger distance h as the distance d and restrict the brackets to consisting of nonnegative functions, which are assumed to be integrable relative to a reference measure µ. Let N ε- h be the minimal number of brackets of size ε needed to cover . Because a bracket of size ε is contained in the ball of radius ε/2 around its midpoint, it follows that Nε/2- h- ≤ N ε- hand hence the present bracketing numbers are bigger than the packing numbers Dε- h defined previously [see (2.1)]. However, in many examples there is also an equality in the other direction, up to a constant, and bracketing and packing numbers give equivalent results. The corresponding bracketing entropy is defined as the logarithm of the bracketing number N ε- h-.

We shall construct a discrete prior supported on densities constructed from minimal sets of brackets for the Hellinger distance. For a given number εn > 0, let n be the uniform discrete measure on the N εn- h densities obtained by covering with a minimal set of εn-brackets and next renormalizing the upper bounds of the brackets to integrate to 1. Thus if l1 u1-lN uN are the N = N εn- h brackets, then n is the uniform measure on the N functions uj/ uj dµ. Next set

$$
\Pi = \sum_{n \in \mathbb{N}} \lambda_n \Pi_n
$$

for a given sequence λn with λn ≥ 0 and <sup>n</sup> λn = 1.

Theorem 3.1. Suppose that εn are numbers decreasing in <sup>n</sup> such that log N εn- h-<sup>≤</sup> nε<sup>2</sup> <sup>n</sup> for every n and nε<sup>2</sup> n/ log n → ∞. Construct the prior as given previously for a sequence λn such that λn > 0 for all n and log λ−<sup>1</sup> <sup>n</sup> = Olog n-. Then the conditions of Theorem 2.2 are satisfied for εn a sufficiently large multiple of the present εn and hence the corresponding posterior converges at the rate εn almost surely, for every P<sup>0</sup> ∈ , relative to the Hellinger distance.

Proof. The prior gives probability 1 to the set <sup>=</sup> <sup>∞</sup> <sup>j</sup>=<sup>1</sup> <sup>j</sup> for <sup>j</sup> the N εj- hfunctions in the support of j. We claim that

(3.1) 
$$
D(8\varepsilon_n, \mathscr{Q}, h) \leq \exp(2n\varepsilon_n^2).
$$

To see this, we first note that, given an ε-bracket l u that contains a probability density p, with ·<sup>2</sup> the norm in L2µ-,

$$
1 \leq \left(\int u \, d\mu\right)^{1/2} = \|\sqrt{u}\|_2 \leq \|\sqrt{u} - \sqrt{p}\|_2 + \|\sqrt{p}\|_2 = h(u, p) + 1 \leq 1 + \varepsilon,
$$
  

$$
h\left(p, \frac{u}{\int u \, d\mu}\right) \leq h(p, u) + h\left(u, \frac{u}{\int u \, d\mu}\right) \leq 2\varepsilon.
$$

Therefore, <sup>n</sup> is a 2εn-net over : every point of is within distance 2εn of some point in n. Since for j>n every point of <sup>j</sup> is within distance 2εj ≤ 2εn of , it follows that <sup>n</sup> is also a 4εn-net over j. This being true for every j>n it follows that <sup>n</sup> is a 4εn-net over <sup>j</sup>≥<sup>n</sup> <sup>j</sup> and hence, trivially, <sup>j</sup>≤<sup>n</sup> <sup>j</sup> is a 4εn-net over . The cardinality of the latter net is at most nN εn- h-<sup>≤</sup> expnε<sup>2</sup> <sup>n</sup> + log n-<sup>≤</sup> exp2nε<sup>2</sup> nfor sufficiently large n.By virtue of the relationship (2.1) between covering numbers and packing numbers, we obtain (3.1). This verifies condition (2.2) with <sup>n</sup> taken equal to and εn taken equal to eight times the present εn.

If u is the upper limit of the εn-bracket containing p0, then

$$
\frac{p_0}{(u/\int u\,d\mu)}\leq \int u\,d\mu\leq (1+\varepsilon_n)^2.
$$

It follows that for large <sup>n</sup> the set of points <sup>p</sup> such that <sup>h</sup><sup>2</sup>p p0- p0/p <sup>∞</sup> <sup>≤</sup> <sup>8</sup>ε<sup>2</sup> n contains at least the function u/ u dµ and hence has prior mass at least

$$
\lambda_n \frac{1}{N_{[\,]}(\varepsilon_n, \mathscr{P}, h)} \ge \exp[-n \varepsilon_n^2 - O(\log n)] \ge \exp(-2n \varepsilon_n^2),
$$

for large n. This verifies condition (2.5) for εn a multiple of the present εn.

Since condition (2.3) is trivially satisfied for <sup>n</sup> = , the proof is complete. ✷

There are many specific examples in which the preceding theorem applies. The situation here is similar to that in recent papers on rates of convergence of (sieved) maximum likelihood estimators, as in Birge and Massart (1993, ´ 1997, 1998), Wong and Shen (1995) or Chapter 3.4 of van der Vaart and Wellner (1996). It is interesting to note that these authors also use brackets, whereas Birge (1983), in his study of the metric entropy of statistical ´ models, uses ε-nets. This is because the cited papers are concerned with a particular type of estimator (namely, minimum contrast estimators), whereas Birge (1983) uses special constructs, called ´ d-estimators. It appears that for good behavior of Bayes estimators on nets we also need some special property of the nets, such as available from nets obtained from brackets.

We include two concrete examples.

Example 3.1 (Smooth densities). Suppose that consists of all measures with densities whose roots √p belong to a fixed multiple of the unit ball of the Holder class ¨ Cα0- 1, for some fixed α > 0. [See, e.g., van der Vaart and Wellner (1996) for a precise definition of this space of functions.] By results of Kolmogorov and Tihomirov (1961), the ε-entropy numbers of this unit ball relative to the uniform norm are bounded by a multiple of 1/ε-<sup>1</sup>/α. [Their result is reproduced in Theorem 2.7.1 of van der Vaart and Wellner (1996).] Because we can construct upper and lower brackets from uniform approximations, this shows that the bracketing Hellinger entropies grow like ε<sup>−</sup>1/α, so that we can take εn of the order n<sup>−</sup>α/2α+1 to satisfy the relation log N εn- h-<sup>≤</sup> nε<sup>2</sup> n. This rate is known to be the frequentist optimal rate for estimators. From Theorem 2.5, we therefore conclude that the prior constructed above achieves the optimal rate of convergence for the posterior.

Upper brackets are, in principle, available from the classical proof of Kolmogorov and Tihomirov (1961). Alternatively, we may use more modern classes of approximating functions, such as wavelets or splines.

Example 3.2 (Monotone densities). Suppose that consists of all monotone decreasing densities on a compact interval in , bounded above by a fixed constant. The root of a monotone density is monotone and hence the bracketing entropy of for the Hellinger distance is bounded by the L2-entropy for the set of monotone functions. This is of the order 1/ε [e.g., van der Vaart and Wellner (1996), Theorem 2.7.5], whence we obtain a n<sup>−</sup>1/3-rate of convergence of the posterior. Again this rate cannot be improved.

Inspection of the proof of the theorem shows that the lower bounds of the brackets are not really needed. The theorem can be generalized by defining upper bracketing numbers Nε- h as the minimal number of functions u1-um such that for every p ∈ there exist a function ui such that both p ≤ ui and hui p- < ε. Next we construct a prior as before. These upper bracketing numbers are clearly smaller than the bracketing numbers N ε- h-. We have formulated the theorem using the better known bracketing numbers, because we do not know any example where this generalization could be useful.

The preceding theorem implicitly requires that the model be totally bounded for the Hellinger metric. A simple modification works for countable unions of totally bounded models, provided that we use a sequence of priors. Suppose that the bracketing numbers of are infinite, but there exist subsets <sup>n</sup> ↑ with finite bracketing numbers. Let εn be numbers such that log N εn n h-<sup>≤</sup> nε<sup>2</sup> <sup>n</sup>. Then we construct n as the discrete uniform distribution on renormalized upper brackets of a minimal set of εn-brackets over n, as before. Then the posterior relative to prior n achieves the convergence rate εn. (Note that this time we do not construct a fixed prior = <sup>n</sup> λnn, but use the prior n when n observations are available.)

In the preceding we start with a condition on the entropies with bracketing even though we apply Theorem 2.2, which demands control over metric entropies only. This is because Theorem 2.2 also requires control over the likelihood ratios. If, for instance, the densities would be uniformly bounded away from zero and infinity, so that the quotients p0/p are uniformly bounded, then we can replace the bracketing entropy in Theorem 3.1 by ordinary entropy. Alternatively, if the set of densities possesses an integrable envelope function, then we can construct priors achieving the rate εn determined by the covering numbers up to logarithmic factors. Here we define εn as the minimal solution of the equation log Nε- h-<sup>≤</sup> nε<sup>2</sup> and <sup>N</sup>ε- h denotes the Hellinger covering number (without bracketing). The construction, described briefly below, parallels Theorem 6 of Wong and Shen (1995) for sieved maximum likelihood estimators.

We assume that the set of densities has a µ-integrable envelope function: a measurable function m with m dµ < ∞ such that p ≤ m for every p ∈ . Given εn > 0 let s<sup>1</sup>n-sNn<sup>n</sup> be a minimal εn-net over [hence Nn = Nεn- h-] and put

$$
g_{j,n} = (s_{j,n}^{1/2} + \varepsilon_n m^{1/2})^2 / c_{j,n},
$$

where cj<sup>n</sup> is a constant ensuring that gj<sup>n</sup> is a probability density. Finally, let n be the uniform discrete measure on g<sup>1</sup> n-gNn<sup>n</sup> and let = <sup>∞</sup> <sup>n</sup>=<sup>1</sup> λnn be a convex combination of the n as before.

Theorem 3.2. Suppose that εn are numbers decreasing in <sup>n</sup> such that log Nεn- h-<sup>≤</sup> nε<sup>2</sup> <sup>n</sup> for every n and nε<sup>2</sup> n/ log n → ∞. Construct the prior as given previously for a sequence λn such that λn > 0 for all n and log λ<sup>−</sup><sup>1</sup> <sup>n</sup> = Olog n-. Assume m is a µ-integrable envelope. Then the corresponding posterior converges at the rate εn log1/εn in probability, relative to the Hellinger distance.

Proof. The proof follows as before, but this time we apply Theorem 2.1, using the observation of Wong and Shen (1995) that for any p ∈ such that hp sj n- ≤ εn we have that hp gj n- = Oεn and that p/gj<sup>n</sup> is bounded above by a multiple of ε−<sup>2</sup> <sup>n</sup> . This verifies (2.4) with εn replaced by a multiple of εn log1/εn through a use of Theorem 5 of Wong and Shen (1995), the relevant part of which is reproduced below as Lemma 8.6. ✷

4. Log spline models. In this section we apply the general results to prior distributions on log spline models for densities. Log spline models for density estimation have been used, among others, by Stone (1990), who shows that the sieved maximum likelihood estimator attains the optimal rate of convergence for estimating a smooth density. As shown by Stone (1994) they can be extended to higher dimensions by using tensor splines, but following Stone (1990), we restrict ourselves to the one-dimensional case.

We assume that the observations are sampled from a density p<sup>0</sup> on the unit interval 0- 1 in the real line that is bounded away from zero and infinity. Our choice of priors will yield the optimal rate of convergence of the posterior if the density p<sup>0</sup> belongs to the Holder space ¨ Cα0- 1. (This is the set of all functions that have α<sup>0</sup> derivatives, for α<sup>0</sup> the largest integer strictly smaller than α, with the α0th derivative being Lipschitz of order α − α0.)

Our prior measures will not be supported on the set of smooth functions, but on exponential families constructed from a spline basis. Fix some "order" q, a natural number, throughout this section. Let K be another natural number, which will increase with n, and partition the half-open unit interval 0- 1 into K subintervals k−1-/K k/K for k = 1--K. Consider the linear space of splines of order q relative to this partition, that is, all functions f 0- 1- → whose restriction to every of the partioning intervals k − 1-/K k/K is a polynomial of degree strictly less than q and, in the case that q ≥ 2, that are q − 2 times continuously differentiable on 0- 1-. It can be shown that this is a J = q + K − 1-dimensional vector space. A convenient basis is the set of B-splines B1--BJ, defined, for example, in de Boor (1978). More precisely, let B1--BJ be the B-splines of order q for the known sequence

$$
\overbrace{0,0,\ldots,0}^{q \text{ times}}, \frac{1}{K}, \frac{2}{K},\ldots, \frac{K-1}{K}, \overbrace{1,1,\ldots,1}^{q \text{ times}},
$$

as defined on page 108 of de Boor (1978). The exact nature of these functions does not matter to us here, except for the following properties [cf. de Boor (1978), pages 109 and 110]:

1. 
$$
B_j \geq 0, \qquad j=1,\ldots,J
$$

2. 
$$
\sum_{j=1}^{J} B_j \equiv 1
$$

- 3. Bj is supported inside an interval of length q/K
- 4. at most q functions Bj are nonzero at every given x.

The first two properties express that the basis elements form a partition of unity, and the third and fourth properties mean that their supports are close to being disjoint if K is very large relative to q.

For θ ∈ <sup>J</sup> let θTB = <sup>j</sup> θjBj and define

$$
p_{\theta}(x) = \exp(\theta^T B(x) - c(\theta)),
$$
  $e^{c(\theta)} = \int_0^1 \exp(\theta^T B(x)) dx.$ 

Thus pθ belongs to a J-dimensional exponential family, with the B-spline functions as sufficient statistics. Since the B-splines add up to unity, the family is actually of dimension J − 1 and we could restrict θ to the subset <sup>0</sup> = θ ∈ <sup>J</sup> θT1 = 0. The true density p<sup>0</sup> of the observations need not be of the form pθ for some θ. (Hence we make a difference between p<sup>0</sup> and p<sup>0</sup> for 0 ∈ <sup>J</sup>; this should not lead to confusion as p<sup>0</sup> does not play a role.) In the following we construct a prior measure n on the set of probability densities on 0- 1 by choosing a prior on 0, which next induces a prior on the probability densities pθ through the map θ → pθ.

For q = 1 the linear space of splines consists of histograms with cell boundaries k/K for k = 0- 1--K. Since exponentials of histograms are histograms, our construction therefore contains priors constructed on histograms as a special case.

Since the true density p<sup>0</sup> need not belong to this "log spline model," we must ensure that it is approximated sufficiently closely by some pθ. To approximate sufficiently many p<sup>0</sup> it is necessary to let the dimension J−1 of the log spline models tend to infinity with n. Here we fix the order q and let the number K of partioning sets tend to infinity. If we focus on α-smooth densities p0, then the minimal rate at which J = Jn must grow is determined by the following lemma, taken from de Boor [(1978), page 170]. Let f∞ = sup0≤x≤<sup>1</sup> fx- be the supremum norm, and let f<sup>α</sup> be the seminorm

$$
||f||_{\alpha}=\sup_{x\neq y}\frac{|f^{(\alpha_0)}(x)-f^{(\alpha_0)}(y)|}{|x-y|^{\alpha-\alpha_0}}.
$$

Because we assume that p<sup>0</sup> is bounded away from zero (and infinity) the function p<sup>0</sup> is in C<sup>α</sup>0- 1 if and only if log p<sup>0</sup> ∈ Cα0-1.

Lemma 4.1. Let <sup>q</sup> <sup>≥</sup> α > <sup>0</sup>. There exists a constant <sup>C</sup> depending only on <sup>q</sup> and α such that, for every p<sup>0</sup> ∈ Cα0-1 that is bounded away from zero,

$$
\inf_{\theta\in\mathbb{R}^J} \|\theta^TB-\log p_0\|_\infty \leq C J^{-\alpha}\|\log p_0\|_\alpha.
$$

It is easy to see from this, as we show in part below, that the root of the Kullback–Leibler divergence and the Hellinger distance between p<sup>0</sup> and the closest pθ are of the order J<sup>−</sup><sup>α</sup> as well. Since a ball of radius εn around p<sup>0</sup> must contain prior mass in order to satisfy (2.9), the rate of convergence εn of the posterior can certainly not be faster than J<sup>−</sup><sup>α</sup>. The minimum distance of alternatives to allow appropriate tests, determined by (2.7), will be shown to satisfy nε<sup>2</sup> <sup>n</sup> ≥ Jn. Together with the previous restriction on εn this will yield a rate of convergence of n<sup>−</sup>α/2α+1-, for Jn <sup>∼</sup> <sup>n</sup><sup>1</sup>/2α+1-. This is also the rate of convergence of the sieved maximum likelihood estimator, found by Stone (1990). It is well known that this rate is optimal for α-smooth densities.

To make this precise we start with stating some lemmas that connect distances and norms on the densities pθ with the J-dimensional Euclidean norm θ and infinity norm θ∞ = max<sup>j</sup> θj. Let f be the L20- 1 norm of f and write a b if a ≤ Cb for a constant C that is universal or depends only on q (which is fixed throughout) and not on K. Most of these are known from or implicit in Stone (1986, 1990) or the literature on approximation theory.

Lemma 4.2. For any <sup>θ</sup> <sup>∈</sup> <sup>J</sup>,

$$
\|\theta\|_{\infty} \lesssim \|\theta^T B\|_{\infty} \le \|\theta\|_{\infty},
$$
  

$$
\|\theta\| \lesssim \sqrt{J} \|\theta^T B\| \lesssim \|\theta\|.
$$

Proof. The first inequality is proved by de Boor [(1978), page 156, Corollary 3]. The second is immediate from the fact that the B-spline basis forms a partition of unity. The third and fourth inequalities are stated in Stone [(1986), equation (12)]. As their full proofs are not in one place, we sketch the argument for completeness.

Let Ii be the interval i − q-/K ∨ 0 i/K ∧ 1 . By (2) on page 155 of de Boor (1978), we have

$$
\sum_i \theta_i^2 \lesssim \sum_i \|\theta^T B_{|I_i}\|_\infty^2 \lesssim \sum_i K \|\theta^T B_{|I_i}\|^2.
$$

The last inequality follows, because <sup>θ</sup>TBIi consists of at most <sup>q</sup> polynomial pieces, each on an interval of length 1/K, and the supremum norm of a polynomial of order q on an interval of length L is bounded by 1/ √ L times the L2-norm, up to a constant depending on q. [To see the last: the squared L20- <sup>1</sup>-norm of the polynomial <sup>x</sup> → <sup>q</sup>−<sup>1</sup> <sup>j</sup>=<sup>0</sup> αjx<sup>j</sup> on 0- 1 is the quadratic form αTEUqU<sup>T</sup> <sup>q</sup> α for Uq = 1- U- - Uq−<sup>1</sup> and U a uniform 0- 1 variable. The second moment matrix EUqU<sup>T</sup> <sup>q</sup> is nonsingular and hence the quadratic form is bounded below by a constant times α<sup>2</sup> ≥ α<sup>2</sup> <sup>∞</sup>.) This yields the third inequality.

By property (3) of the B-spline basis at most q elements Bjx are nonzero for every given x, say for j ∈ Jx-. Therefore,

$$
\left(\theta^T B(x)\right)^2 = \left(\sum_{j\in J(x)} \theta_j B_j(x)\right)^2 \leq \sum_{j\in J(x)} \theta_j^2 B_j^2(x) q,
$$

by the Cauchy–Schwarz inequality. Since each Bj is supported on an interval of length proportional to 1/J and takes its values in 0- 1, its L20- 1-norm is of the order 1/ √ J. Combined with the preceding display this yields

$$
\int_0^1 (\theta^T B(x))^2 dx \lesssim \frac{q}{J} \|\theta\|^2.
$$

This yields the fourth inequality. ✷

Lemma 4.3. For any <sup>θ</sup> <sup>∈</sup> <sup>J</sup> such that <sup>θ</sup>T<sup>1</sup> <sup>=</sup> <sup>0</sup>,

$$
\|\theta\|_{\infty} \lesssim \|\log p_{\theta}\|_{\infty} \leq 2\|\theta\|_{\infty}.
$$

Proof. By the second inequality in Lemma 4.2 we have that θTB∞ <sup>≤</sup> θ∞, whence <sup>e</sup>cθ is contained in the interval e<sup>−</sup>M eM for <sup>M</sup> = θ∞, by its definition, so that cθ- ≤ θ∞. Consequently, by the triangle inequality

$$
\|\log p_{\theta}\|_{\infty}=\|\theta^T B-c(\theta)\|_{\infty}\leq 2\|\theta\|_{\infty}.
$$

This yields the inequality on the right.

For the inequality on the left, we note that, since θT1 = 0,

$$
\begin{aligned} \left| c(\theta) \right| \, &= \, \left| (\theta - c(\theta)1)^T 1 \right| \frac{1}{J} \leq \left\| \theta - c(\theta)1 \right\|_{\infty} \| 1 \|_1 \frac{1}{J} \\ &\lesssim \, \left\| (\theta - c(\theta)1)^T B \right\|_{\infty} = \| \log p_\theta \|_{\infty}, \end{aligned}
$$

by Lemma 4.2. Consequently, by Lemma 4.2 and the triangle inequality,

θ∞ θTB∞ <sup>≤</sup> θTB <sup>−</sup> <sup>c</sup>θ- <sup>∞</sup> <sup>+</sup> cθ- <sup>2</sup> log pθ∞

This concludes the proof. ✷

As a consequence of the preceding lemma, a set of densities pθ is uniformly bounded away from 0 and ∞ if and only if the norms θ∞ of the corresponding set of θ are bounded. This is true uniformly in J ∈ .

Lemma 4.4. For every θ1 θ<sup>2</sup> such that 1<sup>T</sup>θ<sup>1</sup> − θ2-= 0,

$$
\inf_{x, \; \theta} \, p_{\theta}(x) \Big( \frac{\|\theta_1-\theta_2\|^2}{J} \wedge 1 \Big) \! \lesssim \! h^2(\text{$P_{\theta_1}$}, \text{$P_{\theta_2}$}) \! \lesssim \! \sup_{x, \; \theta} p_{\theta}(x) \Big( \frac{\|\theta_1-\theta_2\|^2}{J} \Big),
$$

where the infimum and supremum are taken over all θ on the line segment between θ<sup>1</sup> and θ<sup>2</sup> and all x ∈ 0-1.

Proof. By direct calculation and Taylor's theorem, we have

$$
h^2(P_{\theta_1}, P_{\theta_2}) = 2\left(1 - \exp\left[c\left(\frac{1}{2}\theta_1 + \frac{1}{2}\theta_2\right) - \frac{1}{2}c(\theta_1) - \frac{1}{2}c(\theta_2)\right]\right)
$$
  
= 2\left(1 - \exp\left[-\left(\frac{1}{16}\right)(\theta\_1 - \theta\_2)^T\left(\ddot{c}(\tilde{\theta}) + \ddot{c}(\tilde{\theta})\right)(\theta\_1 - \theta\_2)\right]\right),

for θ˜ and ˜ θ˜ vectors on the line segment between θ<sup>1</sup> and θ<sup>2</sup> and c¨θ the Hessian of c. By the well-known properties of exponential families, we have

$$
\tau^T \ddot{c}(\theta)\tau = \operatorname{var}_{\theta} \tau^T B = \inf_{\mu \in \mathbb{R}} \int_0^1 ((\tau(x) - \mu \mathbf{1})^T B(x))^2 p_{\theta}(x) dx,
$$

since 1TB ≡ 1. Up to bounds below and above on pθ the right side is equivalent to the infimum over µ of the squared L20- 1-norm of τ − µ1-TB.By Lemma 4.2 the latter is comparable to the infimum over <sup>µ</sup> of <sup>τ</sup> <sup>−</sup> <sup>µ</sup>12/J, which is equal to τ2/J if 1Tτ <sup>=</sup> 0.

We can finish the proof by applying this in the first display, together with the inequalities 1 <sup>−</sup> <sup>e</sup>−<sup>x</sup> <sup>≤</sup> <sup>x</sup> for <sup>x</sup> <sup>≥</sup> 0 and 1 <sup>−</sup> <sup>e</sup>−<sup>x</sup> <sup>≥</sup> <sup>1</sup> <sup>2</sup> x ∧ 1 for x ≥ 0 and cx- ∧ 1 ≥ cx ∧ 1for x ≥ 0 and c ≤ 1. ✷

By combining these lemmas we see that the Hellinger distance hPθ<sup>1</sup> - Pθ<sup>2</sup> and 1/ √ J times the J-dimensional Euclidean distance θ<sup>1</sup> − θ2 are proportional, uniformly in J and in θ1 θ<sup>2</sup> having uniformly bounded coordinates. This combined with the estimate on the distance of p<sup>0</sup> to the set of pθ given by Lemma 4.1 reduces the verification of (2.7) and (2.9) to calculations in the Euclidean setting.

We are now ready to prove the following theorem. By Lemma 4.3 there exists a constant <sup>d</sup> such that <sup>d</sup>θ∞ ≤ log pθ∞ for every <sup>θ</sup> <sup>∈</sup> <sup>J</sup> with <sup>θ</sup>T<sup>1</sup> <sup>=</sup> <sup>0</sup> and every J ∈ . We shall assume that the prior is chosen as roughly uniform on a large box −M-M <sup>J</sup>. This corresponds to densities pθ that are bounded and bounded away from zero by at least a small constant.

Theorem 4.5. Suppose that n has a density with respect to Lebesgue measure on θ ∈ Jn <sup>θ</sup>T<sup>1</sup> <sup>=</sup> <sup>0</sup> for Jn <sup>∼</sup> <sup>n</sup>1/2α+1 whose minimum and maximal values on −M-M Jn are bounded below and above by terms of the orders cJn and CJn , respectively, for positive constants c- C, and which vanishes outside −M-M Jn . Let M ≥ 1. Then for every p<sup>0</sup> ∈ Cα0- 1 for q ≥ α ≥ 1/2 such that log <sup>p</sup>0∞ <sup>≤</sup> <sup>1</sup> <sup>2</sup>dM the conditions of Theorem 2.4 are satisfied for εn a large multiple of n<sup>−</sup>α/2α+1 and <sup>n</sup> the support of n, and hence the posterior rate of convergence is n<sup>−</sup>α/2α+1-.

Proof. Let <sup>θ</sup><sup>0</sup> minimize <sup>θ</sup> → log pθ <sup>−</sup> log <sup>p</sup>0∞ over <sup>θ</sup> <sup>∈</sup> <sup>J</sup> such that θ<sup>T</sup>1 = 0. We first show that, for constants C1- C depending on p0, α and q only,

(4.1) 
$$
h(p_{\theta_0}, p_0) \leq C_1 \| \log p_{\theta_0} - \log p_0 \|_{\infty} \leq C J^{-\alpha}.
$$

By Lemma 4.1 there exists θ<sup>∗</sup> such that θ∗-TB <sup>−</sup> log <sup>p</sup>0∞ <sup>J</sup><sup>−</sup><sup>α</sup>. Taking exponentials we see that this implies that expθ∗-TB−p0∞ <sup>J</sup><sup>−</sup><sup>α</sup>, and next, by integrating this inequality, that exp <sup>c</sup>θ∗- −1 <sup>J</sup><sup>−</sup><sup>α</sup>, whence cθ∗- <sup>J</sup><sup>−</sup><sup>α</sup>. Consequently, log pθ<sup>∗</sup> <sup>−</sup> log <sup>p</sup>0∞ <sup>J</sup><sup>−</sup><sup>α</sup>, whence <sup>θ</sup><sup>∗</sup> minimizes <sup>θ</sup> → log pθ <sup>−</sup> log <sup>p</sup>0∞ up to a multiple of <sup>J</sup><sup>−</sup><sup>α</sup>. Since the set of pθ is the same whether <sup>θ</sup> is restricted to satisfy θ<sup>T</sup>1 = 0 or not, the second inequality in the display follows by the definition of θ0. The first now follows easily, since p<sup>0</sup> and hence pθ<sup>0</sup> is bounded away from zero and infinity.

Thus the Hellinger ball of radius ε around P<sup>0</sup> is contained in a multiple of the Hellinger ball of radius ε + J<sup>−</sup><sup>α</sup> around Pθ<sup>0</sup> , whence by Lemma 4.4, for any ε > 0 and suitable constants A- B-<sup>C</sup>, since θ0∞ <sup>≤</sup> <sup>1</sup> <sup>2</sup>M + o1 by the assumption that log <sup>p</sup>0∞ <sup>≤</sup> <sup>1</sup> <sup>2</sup>dM,

$$
\{P_{\theta}: h(P_{\theta}, P_{0}) \leq \varepsilon, \|\theta\|_{\infty} \leq M\}
$$
  

$$
\subset \{P_{\theta}: Ah(P_{\theta}, P_{\theta_{0}}) \leq \varepsilon + J^{-\alpha}, \|\theta\|_{\infty} \leq M\}
$$
  

$$
\subset \{P_{\theta}: B\frac{\|\theta - \theta_{0}\|}{\sqrt{J}} \wedge 1 \leq \varepsilon + J^{-\alpha}, \|\theta - \theta_{0}\| \leq 2\sqrt{J}M, \|\theta\|_{\infty} \leq M\}
$$
  

$$
\subset \{P_{\theta}: \|\theta - \theta_{0}\| \leq C\sqrt{J}(\varepsilon + J^{-\alpha})M\},
$$

since x x ∧ 1 ≤ ε x ≤ M⊂ x x ≤ εM for M ≥ 1. Hence, in view of Example 7.1 [or Pollard (1990), Lemma 4.1] and Lemma 4.4, for constants E-F,

$$
\begin{aligned} D&\Big(\frac{\varepsilon}{2},\big\{P_{\theta}:\,h(P_{\theta},\,P_{0})\leq 2\varepsilon,\,\,\|\theta\|_{\infty}\leq M\big\},\,h\Big)\\ &\leq D\Big(E\,\varepsilon \sqrt{J},\,\big\{\theta;\,\|\theta-\theta_{0}\|\leq 2C\sqrt{J}(\varepsilon+J^{-\alpha})M\big\},\|\cdot\|\Big)\\ &\leq \Big(\frac{F\sqrt{J}(\varepsilon+J^{-\alpha})M}{\varepsilon\sqrt{J}}\Big)^J.\end{aligned}
$$

Therefore, we can verify (2.7) for <sup>n</sup> <sup>=</sup> Pθ θ∞ <sup>≤</sup> <sup>M</sup> and every εn such that

$$
J_n \log \Bigl( 1 + \frac{J_n^{-\alpha}}{\varepsilon_n} \Bigr) \lesssim n \varepsilon_n^2.
$$

Next, we have, with vol<sup>J</sup> the volume of the J − 1--dimensional unit ball,

$$
\Pi_n(P_\theta; h(P_\theta, P_0) \leq 2j\varepsilon, \|\theta\|_\infty \leq M)
$$
  
\$\leq \sup\_{\|\theta\|\_\infty \leq M} \pi\_n(\theta) (2C\sqrt{J}(j\varepsilon + J^{-\alpha})M)^J \text{vol}\_J\$.

By Lemma 4.3 and the assumption that p<sup>0</sup> is bounded, the norms p0/pθ∞ are uniformly bounded over θ ranging over a set of bounded θ∞. Therefore, in view of Lemma 4.4 and (4.1), uniformly in θ∞ ≤ M,

$$
h^2(\,p_\theta,\,p_0)\Big\|\frac{p_0}{p_\theta}\Big\|_\infty \, \lesssim h^2(\,p_\theta,\,p_{\theta_0}) + h^2(\,p_{\theta_0},\,p_0) \, \lesssim \frac{\|\theta-\theta_0\|^2}{J} + J^{-2\alpha}.
$$

We conclude that for εn bigger than a sufficiently large multiple of J<sup>−</sup><sup>α</sup>,

$$
\begin{aligned} \Pi_n&\Big(P_\theta\text{: }h^2(P_\theta,\,p_0)\Big\|\frac{p_0}{p_\theta}\Big\|_\infty\lesssim \varepsilon_n^2\Big)\\ &\geq \Pi_n\Big(\theta\text{: } \|\theta\|_\infty\leq M,\ \|\theta-\theta_0\|/\sqrt{J}\leq \varepsilon_n-J^{-\alpha}\Big)\\ &\geq \inf_{\|\theta\|_\infty\leq M} \pi_n(\theta)\text{vol}\big\{\theta\text{: } \|\theta\|_\infty\leq M,\|\theta-\theta_0\|\leq \frac{1}{2}\varepsilon_n\sqrt{J}\big\}\\ &=\inf_{\|\theta\|_\infty\leq M} \pi_n(\theta)\Big(\frac{1}{2}\varepsilon_n\sqrt{J}\Big)^\prime\text{vol}_J,\end{aligned}
$$

since θ∞ ≤ θ0∞+θ−θ0 √ <sup>J</sup> ≤ θ0∞+<sup>1</sup> <sup>2</sup> εn √ J ≤ M eventually, if θ−θ0 ≤ 1 <sup>2</sup> εn √ J. By assumption, the first term is of the order cJ. Thus condition (2.9) is satisfied if, for all sufficiently large j,

$$
J\log j \lesssim n\epsilon_n^2 j^2 \quad \text{and} \quad \epsilon_n \gtrsim J^{-\alpha}.
$$

This gives εn of the order 1/nα/2α+1 for Jn of the order ε <sup>−</sup>1/α <sup>n</sup> . ✷

5. Finite-dimensionalmodels. Although in this paper we are primarily interested in infinite-dimensional models, it is desirable to have a unified theory applicable to both finite- and infinite-dimensional models. In this section we show that Theorem 2.4 yields the right rate of convergence for finitedimensional models.

Let pθ θ ∈  be a family of densities parametrized by a Euclidean parameter θ running through a set ⊂ <sup>d</sup>. Assume that for every θ θ1 θ<sup>2</sup> ∈ and some α > 0,

$$
\begin{array}{l} \displaystyle - P_{\theta_0} \log \frac{p_\theta}{p_{\theta_0}} \, \lesssim \, \|\theta-\theta_0\|^{2\alpha}, \\\\ \displaystyle P_{\theta_0} \Bigl( \log \frac{p_\theta}{p_{\theta_0}} \Bigr)^2 \, \lesssim \, \|\theta-\theta_0\|^{2\alpha}, \\\\ \displaystyle \|\theta_1-\theta_2\|^\alpha \, \lesssim \, h(P_{\theta_1},\,P_{\theta_2}) \lesssim \|\theta_1-\theta_2\|^\alpha. \end{array}
$$

Assume that the prior measure possesses a density that is uniformly bounded away from zero and infinity on . In this situation the posterior rate of convergence is 1/ <sup>√</sup><sup>n</sup> relative to the Hellinger distance <sup>h</sup>. Under the assumptions, this translates into a n<sup>1</sup>/2α--rate of convergence of the posterior for θ in the Euclidean distance.

Theorem 5.1. Under the conditions listed previously and θ<sup>0</sup> interior to , the conditions of Theorem 2.4 are satisfied for = <sup>n</sup> = Pθ θ ∈ , the Hellinger distance d and εn a sufficiently large multiple of 1/ <sup>√</sup>n.

Proof. The left side of condition (2.7) is seen to be bounded by a constant in Example 7.1 in the case that α = 1. The case of general α is not different. It follows that (2.7) is satisfied for εn <sup>=</sup> M/√<sup>n</sup> and sufficiently large <sup>M</sup>.

In order to verify (2.9) we calculate

$$
\frac{\Pi_n(P_\theta\colon h(P_\theta,P_{\theta_0})\leq j\varepsilon_n)}{\Pi_n(P_\theta\colon -P_{\theta_0}\log(p_\theta/p_{\theta_0})\leq \varepsilon_n^2,\ P_{\theta_0}(\log(p_\theta/p_{\theta_0}))^2\leq \varepsilon_n^2)}\\\leq \frac{\Pi_n(\theta\colon \|\theta-\theta_0\|\leq A(j\varepsilon_n)^{1/\alpha})}{\Pi_n(\theta\colon \|\theta-\theta_0\|\leq B\varepsilon_n^{1/\alpha})}\leq C\Big(\frac{A}{B}\Big)^d\,j^{d/\alpha},
$$

for constants A- B defined by the conditions preceding the theorem, and a constant C depending on the prior density only. It follows that (2.9) is satisfied easily for εn <sup>=</sup> M/√<sup>n</sup> and sufficiently large <sup>M</sup>. ✷

It may be noted that our conditions preclude unbounded parameter spaces : we cannot have that the Hellinger distance is bounded below by a multiple of the Euclidean distance unless the latter is bounded, since the Hellinger distance is uniformly bounded above. This could be improved by replacing condition (2.7) by a testing condition. The lower bound on the Hellinger distance is used only to verify (2.7), which in turn is used only to ensure the existence of tests of θ<sup>0</sup> versus the complements of balls of ε around θ0. For most classical parametric models such tests exist. In fact, existence of uniformly consistent tests of the outside of a compact neigborhood of θ<sup>0</sup> already implies existence of tests with exponential error probabilities (see Lemma 7.2), and this would be sufficient to reduce the problem to a bounded parameter set, to which the preceding theorem applies. Note that the conditions are very reasonable for equal to a small neighborhood of θ0. See also Le Cam (1973) and Le Cam and Yang (1990).

6. Priors based on Dirichlet processes. In this section we apply the general theorems to priors based on Dirichlet processes. A major difficulty is the computation of the prior mass, as in conditions (2.4) or (2.5). We present one such computation and expect that future papers will address more problems of this sort. We shall need an estimate of the probability of an L1-ball under a Dirichlet distribution given by the following lemma.

Lemma 6.1. Let X1--XN be distributed according to the Dirichlet distribution on the N-simplex with parameters m# α1-αN-, where Aε ≤ αi ≤ 1 and <sup>N</sup> <sup>i</sup>=<sup>1</sup> αi <sup>=</sup> <sup>m</sup> for some constant <sup>A</sup>. Let x10-xN0 be any point on the N-simplex. There exist positive constants c and C depending only on A such that, for ε ≤ 1/N,

(6.1) 
$$
\Pr\bigg(\sum_{i=1}^N |X_i - x_{i0}| \leq 2\varepsilon\bigg) \geq C \exp\bigg(-cN \log \frac{1}{\varepsilon}\bigg).
$$

Proof. Find an index <sup>i</sup> such that xi<sup>0</sup> <sup>≥</sup> <sup>1</sup>/N. By relabelling, we can assume that <sup>i</sup> <sup>=</sup> <sup>N</sup>. If xi <sup>−</sup> xi0 ≤ <sup>ε</sup><sup>2</sup> for <sup>i</sup> <sup>=</sup> <sup>1</sup>--N − 1, then

$$
\sum_{i=1}^{N-1} x_i \le 1 - x_{N0} + (N-1)\varepsilon^2 \le (N-1)(\varepsilon^2 + 1/N) \le 1 - \varepsilon^2 < 1.
$$

Hence there exists x = x1-xN in the simplex with these first N − 1 coordinates. Furthermore, <sup>N</sup> <sup>i</sup>=<sup>1</sup> xi−xi0 ≤ <sup>2</sup> N−<sup>1</sup> <sup>i</sup>=<sup>1</sup> xi−xi0 ≤ <sup>2</sup>ε<sup>2</sup>N−1- ≤ 2ε. Therefore the probability on the left-hand side of (6.1) is bounded below by

$$
P(|X_i - x_{i0}| \le \varepsilon^2, i = 1, ..., N - 1)
$$
  
\n
$$
\ge \frac{\Gamma(m)}{\prod_{i=1}^N \Gamma(\alpha_i)} \prod_{i=1}^{N-1} \int_{\max((x_{i0} - \varepsilon^2), 0)}^{\min((x_{i0} + \varepsilon^2), 1)} x_i^{\alpha_i - 1} dx_i.
$$

We use here that <sup>1</sup> <sup>−</sup> <sup>N</sup>−<sup>1</sup> <sup>i</sup>=<sup>1</sup> xiαN−<sup>1</sup> <sup>≥</sup> 1, since αN <sup>≤</sup> 1. Similarly, since αi <sup>≤</sup> <sup>1</sup> for every i, we can lower bound the integrand by 1 and note that the interval of integration contains at least an interval of length <sup>ε</sup>2. Since α:α- = :α+1- ≤ 1 for 0 < α ≤ 1 we can bound the last display from below by

$$
\Gamma(m)\varepsilon^{2(N-1)}\prod_{i=1}^N\alpha_i\geq \Gamma(A)\varepsilon^{2(N-1)}(A\varepsilon)^N\geq C\exp\biggl(-cN\log\frac{1}{\varepsilon}\biggr).
$$

This concludes the proof. ✷

Example 6.1 (Current status censoring). Let Y1--Yn be an i.i.d. sample from a distribution F and C1--Cn be an independent i.i.d. sample from a distribution G, both on 0-∞-. Suppose that we observe Xi = =i- Ci for i = 1-n, where =i = 1 Yi ≤ Ci and would like to estimate F. The density function pF of X with respect to the product of counting measure on 0- 1 and a dominating measure for G at δ cis given by

$$
p_F(\delta, c) = F(c)^{\delta} (1 - F(c))^{1-\delta} g(c).
$$

Since this factorizes in parts depending on F and G only, if we put a product prior on the pair F- G and next compute the posterior for F only, then the part involving G will cancel out. Therefore, it is equivalent to treat g as a known density and put no prior on g.

We assume that G is supported on some compact interval a b and that the true distribution F<sup>0</sup> is continuous and has support which extends to the left and the right of a b. [Hence F0a−- > 0 and F0b- < 1.] As a prior measure on F we consider a Dirichlet prior with base measure α that has a positive, continuous density on a compact interval containing a b. We shall show that the conditions of Theorem 2.2 are satisfied for εn a large multiple of <sup>n</sup><sup>−</sup>1/<sup>3</sup>log <sup>n</sup>-<sup>1</sup>/3. This is very close to the optimal rate of convergence in this model, which is n<sup>−</sup>1/3. We do not exclude the possibility that this small discrepancy is due to suboptimal estimates of the prior mass in the following, and not a deficit of Dirichlet priors. We note that the priors based on ε-nets given in Section 3 do lead to a posterior rate of convergence of n<sup>−</sup>1/3, as the bracketing entropy for this model is of the order 1/ε.

Since the roots of the densities pF are essentially pairs of two bounded monotone functions and the Hellinger distance is the L2-distance between the root densities, the Hellinger entropy of the model PF F ∈ , where is the set of all distribution functions on 0-∞-, can be estimated by the estimate of the entropy of the space of uniformly bounded monotone functions. Thus it is of the order 1/ε [see Theorem 2.7.5 of van der Vaart and Wellner (1996)]. Therefore condition (2.2) is verified for εn equal or bigger than n<sup>−</sup>1/3.

Under our conditions, F<sup>0</sup> is bounded away from 0 and 1 on the interval a b that contains all observation times Ci. Consequently, the quotients pF<sup>0</sup> /pFxare uniformly bounded away from zero and infinity, uniformly in F that are uniformly close to F<sup>0</sup> on the interval a b. The squared Hellinger distance is equal to

$$
h^{2}(P_{F}, P_{F_{0}}) = \int |F^{1/2}(c) - F_{0}^{1/2}(c)|^{2} dG(c)
$$
  
+ 
$$
\int |(1 - F(c))^{1/2} - (1 - F_{0}(c))^{1/2}|^{2} dG(c)
$$
  

$$
\leq C \sup_{c \in [a, b]} |F(c) - F_{0}(c)|^{2},
$$

for a constant depending on F0. Thus, to verify (2.5) it suffices to estimate the prior mass of a Kolmogorov–Smirnov ball of radius εn around F0. Given ε > 0, partition the positive half line in intervals E1--EN such that F0Ei- ≤ ε and Aε ≤ αEi- ≤ 1 for every i and some fixed constant A. We can achieve this with N = O1/ε intervals. By Lemma 6.1, the set F <sup>N</sup> <sup>i</sup>=<sup>1</sup> FEi-−F0Ei- ≤ ε has probability of the order exp−c1/ε log1/ε--. For every F in this set, the Kolmogorov–Smirnov distance to F<sup>0</sup> is of the order ε. We conclude that the prior mass in a Hellinger ball of radius a large multiple of ε is of the order exp−c1/ε log1/ε--. Thus condition (2.5) is verified for εn a large multiple of <sup>n</sup>−1/3log <sup>n</sup>-<sup>1</sup>/3.

7. Existence of tests. In this section we consider some results on the existence of tests of P<sup>0</sup> versus the complement P dP- P0- > ε of the ball of radius ε around P0. The existence of certain tests is a main element in the proofs of Theorems 2.1–2.3 and is guaranteed by entropy bounds. At the end of this section we state a theorem on the rate of convergence of posteriors directly in terms of tests.

Appropriate tests can be built up from tests of P<sup>0</sup> versus balls P dP- P1-≤ η for given P1. Throughout this section we use a distance d such that for every pair P<sup>0</sup> and P<sup>1</sup> in the model there exist tests φn such that, for some universal constant K,

(7.1) 
$$
P_0^n \phi_n \le \exp(-Knd^2(P_0, P_1)),
$$

(7.2) 
$$
\sup_{d(P, P_1) < d(P_0, P_1)/2} P^n(1 - \phi_n) \leq \exp(-Knd^2(P_0, P_1)).
$$

This is true both for d equal to the total variation distance and for d equal to the Hellinger distance. (The constant 2 has no particular interest and is not optimal; any constant bigger than 1 is possible and would do for our purposes.)

More generally, it is known from Birge (1984) and Le Cam (1986) (see ´ Lemma 4 on page 478) that given any two convex sets <sup>0</sup> and <sup>1</sup> of probability measures, there exist tests φn such that

(7.3) 
$$
\sup_{P \in \mathscr{P}_0} P^n \phi_n \leq \exp(n \log \rho(\mathscr{P}_0, \mathscr{P}_1)),
$$

(7.4) 
$$
\sup_{P \in \mathscr{P}_1} P^n(1 - \phi_n) \leq \exp(n \log \rho(\mathscr{P}_0, \mathscr{P}_1)),
$$

where ρ0- 1-<sup>=</sup> <sup>1</sup> <sup>−</sup> <sup>1</sup> 2h<sup>2</sup>0- 1 is the Hellinger affinity, and h0- 1 is the minimum of hP0- P1 over <sup>P</sup><sup>0</sup> <sup>∈</sup> <sup>0</sup> and <sup>P</sup><sup>1</sup> <sup>∈</sup> 1. Because log <sup>ρ</sup> ≤ −<sup>1</sup> 2h<sup>2</sup> this gives exponential decrease of the error probabilities, with the exponent proportional to <sup>−</sup>nh20- 1-. This general result brings out the special role of the Hellinger distance (even though in some situations it may be preferable to work with the log Hellinger affinity directly).

If a distance d is bounded above by the Hellinger distance, then the ball P dP- P1- < dP0- P1-/2 is at Hellinger distance at least dP0- P1-/2 from P0. Thus if this ball is a convex set of probability measures, then (7.1) and (7.2) is satisfied for <sup>d</sup> (with <sup>K</sup> <sup>=</sup> <sup>1</sup> <sup>2</sup> ), by the general results of Birge and Le Cam. ´ This argument immediately gives (7.1) and (7.2) for the Hellinger distance itself and the total variation distance, which satisfies dP−dQ ≤ 2hP- Q-, by the Cauchy–Schwarz inequality. If the set of probability densities under consideration is uniformly bounded, then it also gives (7.1) and (7.2) for the L2-distance, because this is then also bounded by a multiple of the Hellinger distance.

The next step is to combine the tests for balls (which are convex) into a test for the complements of balls, which are nonconvex. The following result is related to Lemma 2.1 in Birge (1983). The number ´ Dε in its first condition is related to the measure of metric dimension used by Birge and Le Cam. ´ The number supε≥εn Dε is almost identical to what Le Cam (1986) calls the dimension of for the pair d εn-.

Theorem 7.1. Suppose that for some nonincreasing function <sup>D</sup>ε-, some εn ≥ 0 and every ε>εn,

$$
D\Big(\frac{\varepsilon}{2},\big\{P;\,\varepsilon\leq d(P,\,P_0)\leq 2\varepsilon\big\},\,d\Big)\leq D(\varepsilon).
$$

Then for every ε>εn there exist tests φn (depending on ε > 0) such that, for a universal constant K and every j ∈ ,

(7.5) 
$$
P_0^n \phi_n \leq D(\varepsilon) \exp(-Kn\varepsilon^2) \frac{1}{1 - \exp(-Kn\varepsilon^2)},
$$

sup dP- P0->jε <sup>P</sup><sup>n</sup><sup>1</sup> <sup>−</sup> φn-<sup>≤</sup> exp <sup>−</sup>Knε<sup>2</sup>j<sup>2</sup> (7.6)

Proof. For a given <sup>j</sup> <sup>∈</sup> choose a maximal jε/2 separated set of points in Sj <sup>=</sup> P jε < dP- P0-≤j + 1ε . This yields a set S <sup>j</sup> of at most Djεpoints and every P ∈ Sj is within distance jε/2 of at least one of these points. (Take S <sup>j</sup> empty and adapt the following in the obvious way if Sj is empty.) For every such point P<sup>1</sup> ∈ S <sup>j</sup> there exists a test ωn with the properties as in (7.1) and (7.2). Let φn be the maximum of all tests attached in this way to some point P<sup>1</sup> ∈ S <sup>j</sup> for some j ∈ . Then

$$
P_0^n \phi_n \leq \sum_j \sum_{P_1 \in S'_j} \exp(-Knj^2 \varepsilon^2) \leq \sum_{j \in \mathbb{N}} D(j\varepsilon) \exp(-Knj^2 \varepsilon^2),
$$
  

$$
\sup_{P \in \bigcup_{i \geq j} S_i} P^n(1 - \phi_n) \leq \sup_{i \geq j} \exp(-Kni^2 \varepsilon^2).
$$

The right sides can be further bounded as desired. [Note that Djε- ≤ Dεfor every j ∈ , by assumption.] ✷

One possible choice for Dε is the ε-packing number Dε/2- d-. This is a bigger number, but in many infinite-dimensional situations this does not appear to yield a real loss. On the other hand, the theorem is needed as stated if is finite-dimensional.

Example 7.1. Suppose that <sup>=</sup>Pθ θ ∈ ⊂ <sup>m</sup> and, for given constants A and B and · the m-dimensional Euclidean norm,

$$
d(P_{\theta}, P_{\theta_0}) \ge A \|\theta - \theta_0\|,
$$
  

$$
d(P_{\theta_1}, P_{\theta_2}) \le B \|\theta_1 - \theta_2\|.
$$

(Since both the Hellinger and total variation metric are bounded, the first can be true with d one of these distances only if is bounded.) The ε-packing number of the m-dimensional unit ball is bounded above by 6/ε<sup>m</sup> [e.g., Pollard (1990), Lemma 4.1]. Thus

$$
D\Big(k\varepsilon,\big\{\theta\in\mathbb{R}^m\colon \|\theta-\theta_0\|\leq l\varepsilon\big\},\|\cdot\|\Big)\leq \bigg(\frac{6l}{k}\bigg)^m.
$$

It follows that

$$
D(\varepsilon,\{P_{\theta}: d(P_{\theta}, P_{\theta_0}) \leq 2\varepsilon\}, d) \leq D\bigg(\frac{\varepsilon}{B},\bigg\{\theta: \|\theta-\theta_0\| \leq \frac{2\varepsilon}{A}\bigg\}, \|\cdot\|\bigg) \leq \bigg(\frac{12B}{A}\bigg)^m.
$$

Thus we can take Dε in Theorem 7.1 independent of ε, but increasing exponentially with the dimension (if A/B is fixed). In comparison, the numbers Dε- hare of the order ε<sup>−</sup><sup>m</sup>.

It is known from Le Cam (1973) that even for a fixed ε there need not exist a consistent sequence of tests of P<sup>0</sup> versus P ∈ dP- P0- > ε . The preceding theorem shows that total boundedness of [which is equivalent to Dε- d being finite for every ε > 0] is sufficient for the existence of such a test. However, this is not necessary. One example showing this is given by (7.1) and (7.2) applied with = P<sup>0</sup> P dP- P1- < dP0- P1-/2, because a total variation or Hellinger ball is usually not totally bounded. A classical example is as follows.

Example 7.2. The collection of all normal distributions <sup>N</sup>θ- 1 on is not totally bounded for the Hellinger or total variation distances, but certainly there are good tests of H<sup>0</sup> θ = 0 versus H<sup>1</sup> θ > ε. [Actually, in this case the affinity satisfies log ρ Nθ- 1-- Nθ - 1- = −θ−θ -<sup>2</sup>/8 and hence we could apply the general form of (7.3) and (7.4) in combination with the Euclidean distance to obtain good tests through the approach of Theorem 7.1, even for unbounded alternatives. For other parametric models the log affinity is typically not nicely related to the Euclidean distance and this approach would fail.]

On the other hand, if a fixed part of can be uniformly consistently tested versus P0, then it can also be tested with exponentially small error probabilities. This implies that such a fixed part can be ignored for our purposes, in that it is not a loss of generality in the main result to assume that the prior only charges the remaining part of (and P0). The error probabilities of the tests φn given in the following lemma are of smaller order than the error probabilities of the tests in Theorem 7.1 if ε = εn → 0 in the latter theorem. This can be useful to reduce the model to a totally bounded submodel, by trimming away parts that can easily be tested by ad hoc arguments. The following lemma is a consequence of results of Le Cam (1973).

Lemma 7.2. Suppose that there exist tests ωn such that for fixed sets <sup>0</sup> and <sup>1</sup> of probability measures

$$
\sup_{P_0\in\mathscr{P}_0}P_0^n\omega_n\to 0,\qquad \sup_{P\in\mathscr{P}_1}P^n(1-\omega_n)\to 0.
$$

Then there exist tests φn and constants K > 0 such that

$$
\sup_{P_0 \in \mathscr{P}_0} P_0^n \phi_n \le e^{-Kn}, \qquad \sup_{P \in \mathscr{P}_1} P^n (1 - \omega_n) \le e^{-Kn}.
$$

In view of the fact that, apparently, entropy conditions are not always appropriate to ensure the existence of tests, it is fruitful to formulate a theorem on rates of convergence directly in terms of existence of tests. The following is a result of this type.

Theorem 7.3. Suppose that (2.8) and (2.9) hold for a sequence εn with εn <sup>→</sup> <sup>0</sup> and nε<sup>2</sup> <sup>n</sup> bounded away from zero and sets <sup>n</sup> ⊂ , and in addition suppose that their exists a sequence of tests φn such that for some constant K > 0 and for every sufficiently large j,

$$
(7.7) \t\t P_0^n \phi_n \to 0,
$$

(7.8) 
$$
\sup_{P \in \mathscr{P}_n: \ \varepsilon_n j < d(P, P_0) \leq 2j\varepsilon_n} P^n(1 - \phi_n) \leq \exp(-Kn\varepsilon_n^2 j^2).
$$

Then for any Mn → ∞, we have that nP dP- P0- ≥ MnεnX1--Xn- → 0 in P<sup>n</sup> <sup>0</sup> -probability.

8. Proof of Theorems 2.1–2.3. In the proof of Theorem 2.1 we use the following simple lemma, which will need to be replaced by more complicated results for the proofs of Theorems 2.2 and 2.3.

Lemma 8.1. For every ε > 0 and probability measure on the set

$$
\{P: -P_0 \log (p/p_0) \le \varepsilon^2, P_0(\log (p/p_0))^{2} \le \varepsilon^2\}
$$

we have, for every C > 0,

$$
P_0^n\bigg(\int \prod_{i=1}^n \frac{p}{p_0}(X_i) d\Pi(P) \le \exp(-(1+C)n\varepsilon^2)\bigg) \le \frac{1}{C^2n\varepsilon^2}.
$$

Proof. By Jensen's inequality applied to the logarithm,

$$
\log \int \prod_{i=1}^n \frac{p}{p_0}(X_i) d\Pi(P) \geq \sum_{i=1}^n \int \log \frac{p}{p_0}(X_i) d\Pi(P).
$$

Thus the probability is bounded by, with <sup>n</sup> <sup>=</sup> <sup>√</sup>n<sup>n</sup> <sup>−</sup> <sup>P</sup>0 the empirical process,

$$
P_0^n\bigg(\mathbb{G}_n\int\log\frac{p}{p_0}\,d\Pi(P)\leq-\sqrt{n}(1+C)\varepsilon^2-\sqrt{n}P_0\int\log\frac{p}{p_0}\,d\Pi(P)\bigg).
$$

By Fubini's theorem and the assumption on the expression on the right of the inequality sign is bounded by −√nε2C. An application of Chebyshev's inequality yields the upper bound

$$
\frac{\mathrm{var}\int \log(p/p_0)(X_1)\,d\Pi(P)}{C^2n\,\varepsilon^4}\leq \frac{P_0\int (\log(p/p_0))^2\,d\Pi(P)}{C^2n\,\varepsilon^4}.
$$

by another application of Jensen's inequality. The right side is bounded by C2nε2-<sup>−</sup><sup>1</sup> by the assumption on . This concludes the proof. ✷

Proof of Theorem 2.1. For every ε > <sup>2</sup>εn we have by (2.2),

$$
\log D\bigg(\frac{\varepsilon}{2},\mathscr{P}_n,d\bigg) \leq \log D(\varepsilon_n,\mathscr{P}_n,d) \leq n \varepsilon_n^2.
$$

Therefore, by Theorem 7.1, applied with Dε-<sup>=</sup> expnε<sup>2</sup> n- (constant in ε) and ε = Mεn and j = 1 in its assertion, where M ≥ 2 is a large constant to be chosen later, there exist tests φn that satisfy

(8.1) 
$$
P_0^n \phi_n \le \exp(n \varepsilon_n^2)
$$

$$
\times \exp(-KnM^2 \varepsilon_n^2) \frac{1}{1 - \exp(-KnM^2 \varepsilon_n^2)},
$$

$$
\sup_{P \in \mathscr{P}_n: d(P, P_0) > M \varepsilon_n} P^n (1 - \phi_n) \le \exp(-KnM^2 \varepsilon_n^2).
$$

By the first condition (8.1) it follows that, if KM<sup>2</sup> <sup>−</sup> <sup>1</sup> > K, as <sup>n</sup> → ∞,

E<sup>P</sup><sup>0</sup> n P dP- P0- ≥ εnX1--Xn φn <sup>≤</sup> <sup>P</sup><sup>n</sup> <sup>0</sup>φn <sup>≤</sup> <sup>2</sup>e<sup>−</sup>Knε<sup>2</sup> <sup>n</sup> (8.3) By Fubini's theorem and the fact that P0p/p0-≤ 1,

$$
\mathbf{E}_{P_0}\int_{\mathscr{P}(\mathscr{P}_n)}\prod_{i=1}^n\frac{p}{p_0}(X_i)\,d\Pi_n(P)\leq\Pi_n(\mathscr{P}-\mathscr{P}_n).
$$

Combining the above assertion with (8.2) we see that

$$
\begin{aligned} \mathcal{E}_{P_0} \int_{P: d(P, P_0) > M_{\varepsilon_n}} & \prod_{i=1}^n \frac{p}{p_0} (X_i) \, d\Pi_n(P) (1 - \phi_n) \\ &\leq \Pi_n(\mathscr{P} \mathscr{P}_n) + \int_{P \in \mathscr{P}_n: d(P, P_0) > M_{\varepsilon_n}} P^n (1 - \phi_n) \, d\Pi_n(P) \\ &\leq \Pi_n(\mathscr{P} \mathscr{P}_n) + \exp(-KnM^2 \varepsilon_n^2) \leq 2\exp(-n\varepsilon_n^2(C+4)), \end{aligned}
$$

for <sup>M</sup> <sup>≥</sup> <sup>C</sup> <sup>+</sup> <sup>4</sup>-/K. By Lemma 8.1, we have with probability tending to 1, with Bn <sup>=</sup> P − P<sup>0</sup> logp/p0-<sup>≤</sup> <sup>ε</sup><sup>2</sup> n- P0logp/p0--<sup>2</sup> <sup>≤</sup> <sup>ε</sup><sup>2</sup> n ,

$$
\int \prod_{i=1}^n \frac{p}{p_0}(X_i) d\Pi_n(P) \ge \exp(-2n\epsilon_n^2) \Pi_n(B_n) \ge \exp(-n\epsilon_n^2(2+C)),
$$

by assumption (2.4). If An is the event that this inequality is true, so that P<sup>n</sup> <sup>0</sup> An-→ 1, then it follows that

$$
\mathbf{E}_{P_0}\Pi_n(P: d(P, P_0) > M\varepsilon_n|X_1, \ldots, X_n)(1-\phi_n)\mathbf{1}_{A_n}
$$
  
\n
$$
\leq \exp(n\varepsilon_n^2(2+C))2\exp(-n\varepsilon_n^2(C+4)) \to 0.
$$

This concludes the proof. ✷

For the proof of Theorem 2.2 we need a replacement of Lemma 8.1 that gives a faster rate of convergence in its statement. We can achieve this by controlling the quotients p/p0. First, if one has uniform control from below, then the Hellinger distance and the Kullback–Leibler information are comparable. The following lemma can be found in Birge and Massart (1998) [see their (7.6)]. ´

Lemma 8.2. For any pair of probability measures P and P0,

$$
h^2(P, P_0) \le -P_0 \log \frac{p}{p_0} \le 2h^2(P, P_0) \left[1 + \log \left\|\frac{p_0}{p}\right\|_{\infty}\right] \le 2h^2(P, P_0) \left\|\frac{p_0}{p}\right\|_{\infty}
$$

A second lemma is a comparison of a certain exponential moment and the Hellinger distance. This exponential moment [called the "Bernstein norm" in van der Vaart and Wellner (1996) even though it is not a norm] is essential in Bernstein's inequality. Birge and Massart (1993) used this "norm" to derive ´ results on rates of convergence of minimum contrast estimators.

Lemma 8.3. For any pair of probability measures P and P0,

$$
(8.4) \tP_0(\exp(|\log(p/p_0)|)-1-|\log(p/p_0)|)\leq 2h^2(P,P_0)\left\|\frac{p_0}{p}\right\|_{\infty}.
$$

Proof. For every <sup>c</sup> <sup>≤</sup> 0 and <sup>x</sup> <sup>≥</sup> <sup>c</sup> we have the inequality

(8.5) 
$$
(e^{|x|}-1-|x|)\leq 2e^{|c|}(e^{x/2}-1)^2.
$$

If <sup>e</sup><sup>−</sup><sup>c</sup> = p0/p∞, then logp/p0- ≥ c and hence the integrand on the left side of (8.4) is bounded above by

$$
2e^{-c}\left(\exp\left(\frac{1}{2}\log(p/p_0)\right)-1\right)^2=2e^{-c}\left(\sqrt{\frac{p}{p_0}}-1\right)^2.
$$

The integral of the right side with respect to P<sup>0</sup> is equal to 2e−<sup>c</sup> times the squared Hellinger distance. ✷

The "Bernstein norm" of logp/p0 dominates all moments of order greater than or equal to 2 of logp/p0 up to constants, including the second moment up to a factor 2. Therefore, when combined the preceding two lemmas show that

(8.6) 
$$
\left\{P: h^2(P, P_0)\middle\|\frac{p_0}{p}\middle\|_{\infty} \leq \varepsilon^2\right\} \subset \left\{P: P_0 \log \frac{p_0}{p} \leq 2\varepsilon^2, P_0\left(\log \frac{p}{p_0}\right)^2 \leq 4\varepsilon^2\right\}.
$$

This shows that condition (2.4) is weaker than condition (2.5), up to constants. Actually controlling all moments is more than what is needed. Another possible extension of Lemma 8.1 would be to replace the second moment of logp/p0 by a higher moment (and use Markov's inequality at the end of the proof). This would give a result good enough for the proof of Theorem 2.2 provided the higher moment is chosen "high enough" (dependent on the order εn, faster convergence to zero needing a higher moment). We have chosen here to forego such refinements and obtain an exponential inequality under a somewhat stronger assumption.

We are ready for an adaptation of Lemma 8.1

Lemma 8.4. For every ε > 0 and probability measure on the set

(8.7) 
$$
\{P: h^2(P, P_0) \|p_0/p\|_{\infty} \leq \varepsilon^2\},\
$$

we have, for a universal constant B > 0,

(8.8) 
$$
P_0^n\bigg(\int \prod_{i=1}^n \frac{p}{p_0}(X_i) d\Pi(P) \le \exp(-3n\varepsilon^2)\bigg) \le \exp(-Bn\varepsilon^2).
$$

Proof. Lemma 8.2 gives that <sup>−</sup>P<sup>0</sup> log p/p<sup>0</sup> <sup>≤</sup> <sup>2</sup>ε<sup>2</sup> for every <sup>P</sup> in the set (8.7), which has -probability 1. Furthermore, by Lemma 8.3,

> P0exp logp/p0- − 1 − logp/p0--<sup>≤</sup> <sup>2</sup>ε<sup>2</sup>

By monotonicity and convexity of the function y → e<sup>y</sup> − 1 − y on 0-∞ and Jensen's inequality,

$$
P_0\bigg(\exp\bigg(\bigg|\int \log(p/p_0) d\Pi(P)\bigg|\bigg) - 1 - \bigg|\int \log(p/p_0) d\Pi(P)\bigg|\bigg)
$$
  

$$
\leq P_0\int (\exp(|\log(p/p_0)|) - 1 - |\log(p/p_0)|) d\Pi(P) \leq 2\varepsilon^2,
$$

by Fubini's theorem. By the lemma below the same bound is true for <sup>1</sup> <sup>2</sup> times the variable logp/p0dP centered at its expectation. Therefore, rewriting the probability on the left side of (8.8) as in the proof of Lemma 8.1, we see that it is bounded above by

$$
P_0^n\bigg(\mathbb{G}_n\int\bigg(\log\frac{p}{p_0}\bigg)d\Pi(P)\leq-3\sqrt{n}\varepsilon^2+\sqrt{n}2\varepsilon^2\bigg)\leq\exp\bigg(-D\frac{n\varepsilon^4}{\varepsilon^2+\sqrt{n}\varepsilon^2/\sqrt{n}}\bigg),
$$

by (the refined version of) Bernstein's inequality. [see, e.g., Lemma 2.2.11 of van der Vaart and Wellner (1996).] ✷

Lemma 8.5. If <sup>ψ</sup> 0-∞- → is convex and nondecreasing, then EψX − EX- ≤ Eψ2Xfor every random variable X.

Proof. The map <sup>y</sup> → <sup>ψ</sup>y is convex on . If X is an independent copy of X, then the left side is equal to EψX−EX - ≤ EψX−X -, by Jensen's inequality. Next bound X − X ≤X+X and use the monotonicity and convexity of ψ again to bound the expectation by E<sup>1</sup> <sup>2</sup> ψ2X-+ψ2X --, which is the right side. ✷

Proof of Theorem 2.2. The proof of Theorem 2.2 follows the same lines as the proof of Theorem 2.1. The difference is that we use Lemma 8.4 instead of Lemma 8.1 to ensure that the probability of the events An converges to 1 at an exponential rate. By inspecting the proof, we conclude that for some B1-B<sup>2</sup> > 0 and M chosen as before,

$$
P_0(\Pi_n(P: d(P, P_0) > M\varepsilon_n | X_1, \ldots, X_n) \geq \exp(-B_1 n \varepsilon_n^2))
$$

converges to zero at the rate exp−B2nε<sup>2</sup> n-. Since <sup>n</sup> exp−B2nε<sup>2</sup> n- < ∞, almost sure convergence follows by the Borel–Cantelli lemma. ✷

For the proof of Theorem 2.3 we need other variations on the preceding lemmas. The following lemma follows from Theorem 5 of Wong and Shen (1995). Let log<sup>+</sup> x = log x-∨ 0.

Lemma 8.6. For any pair of probability measures P and P<sup>0</sup> such that hP- P0- ≤ 044 and P0p0/p-< ∞,

$$
\begin{aligned} -P_0\log\frac{p}{p_0} &\leq 18h^2(P,P_0)\bigg(1+\log_{+}\frac{\sqrt{P_0(p_0/p)}}{h(P,P_0)}\bigg),\\ P_0\bigg(\log\frac{p}{p_0}\bigg)^2 &\leq 5h^2(P,P_0)\bigg(1+\log_{+}\frac{\sqrt{P_0(p_0/p)}}{h(P,P_0)}\bigg)^2. \end{aligned}
$$

Lemma 8.7. For any pair of probability measures P and P<sup>0</sup> such that P0p0/p-< ∞,

$$
P_0(\exp(|\log(p/p_0)|)-1-|\log(p/p_0)|)\leq 4h^2(P,P_0)(1+\Phi^{-1}(h^2(P,P_0))).
$$

for <sup>−</sup>1ε- = sup M M- ≥ ε the inverse of the function M- = P0p0/p-1 p0/p ≥ M/M.

Proof. Set <sup>m</sup> <sup>=</sup> <sup>p</sup>0/p. By inequality (8.5) in the proof of Lemma 8.3, the left side is bounded above by

$$
\begin{aligned} 2P_0\bigg(\sqrt{\frac{p}{p_0}}-1\bigg)^21\{p\geq p_0\}+2P_0\frac{p_0}{p}\bigg(\sqrt{\frac{p}{p_0}}-1\bigg)^21\{pM\}, \end{aligned}
$$

for every M > 0. The function is left continuous and strictly decreasing from infinity at 0 to 0 at a point <sup>τ</sup> ≤ ∞. If we choose <sup>M</sup> <sup>=</sup> <sup>−</sup>1h2P- P0--, then MM-<sup>≥</sup> Mh2P- P0- ≥ MM+- = P0m1 m>M. The right side of the last display can now be bounded by an expression as in the lemma. ✷

Lemma 8.8. For a given function <sup>m</sup> let <sup>−</sup><sup>1</sup>ε- = sup M M- ≥ ε be the inverse function of M- = P0m1 m ≥ M/M. For every ε ∈ 0- 044 and probability measure on the set

$$
\left\{P: p_0/p \le m, 18h^2(P, P_0)\bigg(1 + \log_+\frac{\sqrt{P_0 m}}{h(P, P_0)} + \Phi^{-1}(h^2(P, P_0))\bigg) \le \varepsilon^2\right\}
$$

we have, for a universal constant B > 0,

(8.9) 
$$
P_0^n\bigg(\int \prod_{i=1}^n \frac{p}{p_0}(X_i) d\Pi(P) \le \exp(-2n\varepsilon^2)\bigg) \le \exp(-Bn\varepsilon^2).
$$

Proof. This follows the same lines as the proof of Lemma 8.4, now substituting Lemmas 8.6 and 8.7 for Lemmas 8.2–8.3.

Proof of Theorem 2.3. This is identical to the proof of Theorem 2.2, except that we use Lemma 8.8 instead of Lemma 8.4.

Proof of Theorem 2.4. The first part of the proof is identical to the first part of the proof of Theorem 2.1, except that we choose the tests φn to satisfy (8.1) and [instead of 8.2] for every j ∈ ,

(8.10) 
$$
\sup_{P \in \mathscr{P}_n: d(P, P_0) > M \varepsilon_n j} P^n (1 - \phi_n) \le \exp(-KnM^2 \varepsilon_n^2 j^2).
$$

We also choose M large enough to ensure that the right side of (8.1) and hence the left side of (8.3) converges to zero. Defining Sn<sup>j</sup> <sup>=</sup> P ∈ <sup>n</sup> Mεnj < dP- P0- ≤ Mεnj + 1- and using (8.10), we obtain

$$
\mathbf{E}_{P_0}\int_{S_{n,j}}\prod_{i=1}^n\frac{p}{p_0}(X_i)\,d\Pi_n(P)(1-\phi_n)\leq \exp(-KnM^2\varepsilon_n^2j^2)\Pi_n(S_{n,j}).
$$

Fix some C<sup>0</sup> ≥ 1. By Lemma 8.1, we have on an event An with probability at least 1 − nε<sup>2</sup> nC<sup>2</sup> 0-<sup>−</sup>1,

$$
\int \prod_{i=1}^n \frac{p}{p_0}(X_i) d\Pi_n(P) \ge \exp(-2C_0 n \varepsilon_n^2) \Pi_n(B_n(\varepsilon_n)).
$$

Hence, by assumption (2.9), for every sufficiently large J,

$$
\begin{split} &\mathbf{E}_{P_0}\Pi_n\big(P\in\mathscr{P}_n: d(P, P_0) > J\varepsilon_n M|X_1, \ldots, X_n\big)(1-\phi_n)1_{A_n} \\ &\leq \sum_{j\geq J} \frac{\exp(-KnM^2\varepsilon_n^2 j^2)\Pi_n(S_{n,j})}{\exp(-2C_0n\varepsilon_n^2)\Pi_n(B_n(\varepsilon_n))} \\ &\leq \sum_{j\geq J} \exp(-n\varepsilon_n^2(KM^2 j^2 - 2C_0 - \frac{1}{2}KM^2 j^2)). \end{split}
$$

This converges to zero as <sup>J</sup> → ∞ if nε<sup>2</sup> <sup>n</sup> is bounded away from zero. Next

$$
\mathbf{E}_{P_0}\Pi_n(P\notin\mathscr{P}_n|X_1,\ldots,X_n)(1-\phi_n)\mathbf{1}_{A_n}\leq \frac{\Pi_n(\mathscr{P}\mathscr{P}_n)}{\exp(-2C_0n\epsilon_n^2\Pi_n)(B_n(\epsilon_n))}.
$$

We may assume that either nε<sup>2</sup> <sup>n</sup> is bounded or nε<sup>2</sup> <sup>n</sup> → ∞; otherwise we argue along subsequences. If nε<sup>2</sup> <sup>n</sup> is bounded, then we first choose C<sup>0</sup> large but fixed so as to make P<sup>n</sup> <sup>0</sup> An as large as desired. Then the right side of the preceding display converges to zero by assumption (2.8). If nε<sup>2</sup> <sup>n</sup> → ∞, then we choose <sup>C</sup><sup>0</sup> <sup>=</sup> 1, in which case <sup>P</sup><sup>n</sup> <sup>0</sup> An- → 1 and again the right side of the preceding display converges to zero. ✷

Proof of Theorem 7.3. This is essentially contained in the proof of Theorem 2.4 (take M = 1).

Acknowledgment. We thank Lucien Birge for insightful discussions that ´ have led to an improved presentation (and some corrections), in particular relating to Section 7.

## REFERENCES

- Barron, A., Schervish, M. J. and Wasserman, L. (1999). The consistency of posterior distributions in nonparametric problems. Ann. Statist. 27 536–561.
- Birge, L. ´ (1983). Approximation dans les espaces metriques et th ´ eorie de l'estimation. ´ Z. Wahrsch. Verw. Gebiete 65 181–238.
- Birge, L. ´ (1984). Sur un theor ´ eme de minimax et son application aux tests. ` Probab. Math. Statist. 3 259–282.
- Birge, L. ´ and Massart, P. (1993). Rates of convergence for minimum contrast estimators. Probab. Theory Related Fields 97 113–150.
- Birge, L. ´ and Massart, P. (1997). From model selection to adaptive estimation. In Festschrift for Lucien Le Cam (G. Yang and D. Pollard, eds.) 55–87. Springer, New York.
- Birge, L. ´ and Massart, P. (1998). Minimum contrast estimators on sieves: exponential bounds and rates of convergence. Bernoulli 4 329–375.
- de Boor, C. (1978). A Practical Guide to Splines. Springer, New York.
- Diaconis, P. and Freedman, D. (1986). On the consistency of Bayes estimates (with discussion). Ann. Statist. 14 1–67.
- Doob, J. L. (1949). Le Calcul des Probabilites et ses Applications. ´ Coll. Int. du CNRS 13 23–27.
- Dudley, R. M. (1984). A course on empirical processes. Lectures Notes in Math. 1097 2–141. Springer, Berlin.
- Ferguson, T. S. (1973). A Bayesian analysis of some nonparametric problems. Ann. Statist. 1 209–230.
- Ferguson, T. S. (1974). Prior distribution on the spaces of probability measures. Ann. Statist. 2 615–629.
- Freedman, D. A. (1963). On the asymptotic behavior of Bayes' estimates in the discrete case. Ann. Math. Statist. 34 1194–1216.
- Freedman, D. A. (1965). On the asymptotic behavior of Bayes' estimates in the discrete case II. Ann. Math. Statist. 36 454–456.
- Ghosal, S., Ghosh, J. K. and Ramamoorthi, R. V. (1997). Non-informative priors via sieves and packing numbers. In Advances in Statistical Decision Theory and Applications (S. Panchapakeshan and N. Balakrishnan eds.) 129–140. Birkhauser, Boston. ¨
- Ghosal, S., Ghosh, J. K. and Ramamoorthi R. V. (1999a). Posterior consistency of Dirichlet mixtures in density estimation. Ann. Statist. 27 143–158.
- Ghosal, S., Ghosh, J. K. and Ramamoorthi, R. V. (1999b). Consistency issues in Bayesian nonparametrics. In Asymptotics, Nonparametrics and Time Series: A Tribute to Madan Lal Puri (Subir Ghosh, ed.) 639–667. Dekker, New York.
- Ibragimov, I. A. and Has'minskii, R. Z. (1981). Statistical Estimation: Asymptotic Theory. Springer, New York.
- Kolmogorov, A. N. and Tikhomirov, V. M. (1961). Epsilon-entropy and epsilon-capacity of sets in function spaces. Amer. Math. Soc. Trans. Ser. 2 17 277–364.
- Le Cam, L. M. (1973). Convergence of estimates under dimensionality restrictions. Ann. Statist. 1 38–53.
- Le Cam, L. M. (1986). Asymptotic Methods in Statistical Decision Theory. Springer, New York.
- Le Cam, L. M. and Yang, G. (1990). Asymptotics in Statistics: Some Basic Concepts. Springer, New York.
- Pollard, D. (1990). Empirical Processes: Theory and Applications. IMS, Hayward, CA and Amer. Statist. Assoc., Alexandria, VA.
- Schwartz, L. (1965). On Bayes procedures. Z. Wahrsch. Verw. Gebiete 4 10–26.
- Shen, X. and Wasserman, L. (1999). Rates of convergence of posterior distributions. Preprint.
- Stone, C. J. (1986). The dimensionality reduction principle for generalized additive models. Ann. Statist. 14 590–606.
- Stone, C. J. (1990). Large-sample inference for log-spline models. Ann. Statist. 18 717–741.

Stone, C. J. (1994). The use of polynomial splines and their tensor products in multivariate function estimation (with discussion). Ann. Statist. 22 118–184.

- van der Vaart, A. W. and Wellner, J. A. (1996). Weak Convergence and Empirical Processes. Springer, New York.
- Wasserman, L. (1998). Asymptotic properties of nonparametric Bayesian procedures. Practical Nonparametric and Semiparametric Bayesian Statistics. Lecture Notes in Statist. 133 293–304. Springer, New York.
- Wong, W. H. and Shen, X. (1995). Probability inequalities for likelihood ratios and convergence rates of sieve MLEs. Ann. Statist. 23 339–362.

S. Ghosal A. W. van der Vaart Department of Mathematics Free University De Boelelaan 1081a 1081 HV Amsterdam Netherlands E-mail: aad@cs.vu.nl

J. K. Ghosh Statistics and Mathematics Unit Indian Statistical Institute 203 B.T. Road Calcutta 700 035 India