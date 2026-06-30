# BASS Travel Stipend Application: Research Summary & Abstract

## Applicant Rationale for BASS Attendance
I am highly interested in attending the Biopharmaceutical Applied Statistics Symposium (BASS XXXIII). As a doctoral candidate in Biostatistics, my career goal is to work at the intersection of biostatistics, causal inference, and pharmaceutical development, specifically focused on utilizing Real-World Evidence (RWE) and generalizability methods. Participating in BASS will allow me to present my research on causal transportability to leading experts from academia, the pharmaceutical industry, and regulatory agencies (such as the FDA), obtaining vital feedback on the regulatory and practical viability of my methods. Furthermore, the symposium’s focus on timely statistical applications in clinical trials and post-market safety will provide me with a deeper understanding of the regulatory landscape, directly informing the final stages of my dissertation.

---

## Dissertation/Doctoral Research Abstract

### Presentation Details
**Title:** Random Weight BART for Causal Transportability: Generalizing GLP-1 + SGLT2 Efficacy to the HCHS/SOL Hispanic Target Population  
**Presenter:** Alvaro Quijano, Doctoral Candidate in Biostatistics  

### Biography
Alvaro Quijano is a doctoral candidate in Biostatistics. His dissertation research focuses on developing robust non-parametric Bayesian methods for complex survey designs, with applications in causal inference, precision medicine, and health economics. His current work addresses the theoretical and computational challenges of generalizability and transportability, developing machine learning methods to transport treatment effects from clinical trials to target populations represented in large-scale databases like HCHS/SOL and NHANES.

### Abstract
Randomized Controlled Trials (RCTs) are the gold standard for establishing treatment efficacy due to their internal validity. However, strict inclusion and exclusion criteria often limit their external validity, making clinical trial results difficult to generalize to the broader, heterogeneous real-world patient population. Causal transportability methods address this by using participation or propensity weighting to transport estimated treatment effects from a clinical trial to a target population represented by large-scale complex health surveys. A key clinical application of interest is transporting the treatment effects from the **DURATION-8 clinical trial (evaluating exenatide plus dapagliflozin combination therapy)** to estimate the real-world population impact of adding GLP-1 and SGLT2 agents to patients already receiving standard **RAAS (Renin-Angiotensin-Aldosterone System) medication** for cardiovascular and renal protection.

An unresolved challenge in this setting is the dual source of uncertainty: the causal model estimation uncertainty from the trial, and the sampling design variance (stratification, clustering, and weights) from the target survey database. Standard pseudo-Bayesian approaches treat the survey weights as fixed, precise information, which systematically underestimates the posterior variance of the transported treatment effect and leads to overconfident regulatory or clinical decisions in target minority populations.

To resolve this, we present the **Random Weight BART** framework for causal transportability. This methodology combines Bayesian Additive Regression Trees (BART)—allowing for flexible, non-parametric estimation of heterogeneous treatment effects and complex interactions with baseline RAAS medication use—with the unified bootstrap resampling weights of Kim et al. (2024). By drawing a fresh set of bootstrap weights from the target survey's design at each MCMC iteration, our method stochastically perturbs the transition probabilities during the Bayesian backfitting steps. This natively propagates both causal model uncertainty and complex survey design variance into the posterior draws, achieving frequentist sandwich covariance coverage. 

We demonstrate through simulations and an application transporting treatment effects from DURATION-8 to the **Hispanic Community Health Study / Study of Latinos (HCHS/SOL)** cohort that our framework provides robust, population-representative treatment effect estimates with correctly calibrated credible intervals, preventing false-positive findings in real-world evidence translation.

---

## Academic Background & Experience Summary
* **Academic Background:** Doctoral Candidate in Biostatistics. Coursework includes advanced statistical inference, causal inference, Bayesian data analysis, survival analysis, machine learning, and survey sampling theory.
* **Research Focus:** Non-parametric Bayesian statistics, causal inference, generalizability and transportability, complex survey designs, and machine learning (BART, DPMM).
* **Relevant Experience:** Extensive experience working with large-scale national health databases (HCHS/SOL, MEPS, NHANES) and clinical datasets. Skilled in statistical programming in R, Python, and C++, with a focus on implementing scalable MCMC samplers for high-dimensional clinical and epidemiological data.
