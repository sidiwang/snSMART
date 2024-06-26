---
title: "snSMART: An R Package for Small Sample, Sequential, Multiple Assignment, Randomized Trial Data Analysis"
tags:
  - rare disease
  - clinical trial
  - snSMART
  - Bayesian
authors:  
  - name: Sidi Wang
    orcid: 0000-0003-4838-0842
    affiliation: 1
  - name: Fang Fang
    orcid: 0000-0002-7089-3591
    affiliation: 1
  - name: Roy Tamura
    affiliation: 2
  - name: Thomas Braun
    orcid: 0000-0002-7113-6998
    affiliation: 1
  - name: Kelley M Kidwell
    orcid: 0000-0002-1717-4483
    affiliation: 1
affiliations:
  - name: University of Michigan
    index: 1
  - name: University of South Florida
    index: 2
date: 22 May 2024
bibliography: paper.bib

---

# Summary

Small sample, sequential, multiple assignment, randomized trials (snSMARTs) are multistage trials with the primary goal of estimating treatment effects. Participants are first randomized to one of the first stage treatments, then participants may be re-randomized in the second stage depending on their outcome or response to first stage treatment. In an snSMART, treatment effects in the first stage are estimated using data from both stages, which results in more precise estimates. The design of the snSMART may also help to improve participant recruitment and retention over standard rare disease designs. In the past few years, substantial progress has been made regarding the development of statistical methods for snSMART designs. To better facilitate the application of these statistical methods, we introduce the R package snSMART in this paper, which can be used to analyze data from an existing snSMART design or aid in the design and sample size of an snSMART.

# Statement of need

The design and methods of snSMARTs generally apply to any disorder or disease that affects a small group of people and remains stable over the duration of the trial. Recent progress includes methods for snSMARTs with three active treatments [@wei2018bayesian; @wei2020sample; @chao2020dynamic], group sequential designs [@chao2020bayesian], placebo and two dose levels [@fang2021bayesian], and continuous outcomes [@hartman2021design]. However, there is a lack of software for these methods. The snSMART R package fills this gap by providing sample size calculations and trial data analysis using Bayesian and frequentist approaches. To our knowledge, no other R packages offer similar functions.

## snSMART comparing two dose levels with placebo

Here we explain one of the snSMART designs in detail. This snSMART design investigates the response rate of one experimental treatment at low and high doses compared with placebo [@fang2021bayesian]. In such an snSMART design (Figure \ref{fig:snSMART-dose}), participants are equally assigned to either receive placebo, low dose, or high dose in the first stage. Participants receive their initial treatment for a pre-specified amount of time until the measurement of their responses at the end of stage 1. In the second stage, all participants who received placebo or low dose in the first stage are re-randomized to either low or high dose regardless of their first stage response. Participants who responded to high dose are re-randomized between low and high dose, whereas those who did not respond to high dose receive high dose again in the second stage. The main goal of this snSMART is to estimate first stage low and high dose treatment response rates and compare them to placebo by modeling the pooled data from two stages.

![\label{fig:snSMART-dose}](dose_snSMART.png)

@fang2021bayesian used the Bayesian joint stage model (BJSM) from @wei2018bayesian in Equations \ref{eq:1} and \ref{eq:2}. For this design $m = P, L, H$ and $m' = L, H$ where P = placebo, L = low dose and H = high dose. The prior distribution for the response rate of placebo may be informed by natural history studies or previous trials and specified such that $\pi_P \sim Beta(\zeta_n, \eta_n)$. A weak tendency for the doses of the drug to have greater response rates than the effect of placebo can be assumed through $\log(\pi_L/\pi_P) \sim N(\mu, \sigma^2)$ and $\log(\pi_H/\pi_P) \sim N(\mu, \sigma^2)$. The prior distributions for the linkage parameters may vary with guidance of specifying $\beta_{0m}, \beta_{1m} \sim Gamma(\omega, \psi)$. 

The BJSM is specified as follows:

\begin{equation}
Y_{i1m}|\pi_m \sim Bernoulli(\pi_m) \label{eq:1}
\end{equation}
\begin{equation}
Y_{i2m'}|Y_{i1m}, \pi_{m'}, \beta_{1m}, \beta_{0m} \sim Bernoulli((\beta_{1m}\pi_{m})^{Y_{i1m}}(\beta_{0m}\pi_{m'})^{1-Y_{i1m}}) \label{eq:2}
\end{equation}
for $i = 1,...,N$; and $j = 1, 2$; where 
$Y_{ijm}$ is the outcome for participant $i$ at stage $j$ for treatment $m$ and takes the value 1 for response to treatment and 0 for no response; 
$N$ is the total sample size;
$\beta_{0m}$ and $\beta_{1m}$ are the linkage parameters for non-responders and responders, respectively;
$\pi_m$ is the first stage response rate for treatment $m$;
$\beta_{1m}\pi_{m}$ is the second stage response rate for first stage responders; and 
$\beta_{0m}\pi_{m'}$ is the second stage response rate for non-responders to treatment $m$ in the first stage who receive treatment $m'$ in the second stage. 

Table: Summary of the functionality of the snSMART package.

+----------------------------------------------+----------------------------------------------------------------------------+
| Function                                     | Description                                                                | 
|                                              |                                                                            |       
+:=============================================+:===========================================================================+
| *BJSM functions*                             |                                                                            | 
+----------------------------------------------+----------------------------------------------------------------------------+
| BJSM_binary                                  | BJSM binary (3AT or P2D snSMART)                                           |
+----------------------------------------------+----------------------------------------------------------------------------+
| BJSM_c                                       | BJSM (3AT snSMART with a mapping function and continuous outcome)          |
+----------------------------------------------+----------------------------------------------------------------------------+
| group_seq                                    | BJSM (interim analysis and final analysis of group sequential 3AT snSMART) |
+----------------------------------------------+----------------------------------------------------------------------------+
+----------------------------------------------+----------------------------------------------------------------------------+
| *Frequentist functions*                      |                                                                            | 
+----------------------------------------------+----------------------------------------------------------------------------+
| LPJSM_binary                                 | LPJSM (3AT or P2D snSMART)                                                 |
+----------------------------------------------+----------------------------------------------------------------------------+
+----------------------------------------------+----------------------------------------------------------------------------+
| *Sample size calculation*                    |                                                                            | 
+----------------------------------------------+----------------------------------------------------------------------------+
| sample_size                                  | 3AT snSMART sample size calculation                                        |
+----------------------------------------------+----------------------------------------------------------------------------+
+----------------------------------------------+----------------------------------------------------------------------------+
| *S3 summary and print methods*               |                                                                            |
+----------------------------------------------+----------------------------------------------------------------------------+
| for class ’BJSM_binary’                      | Summarize and print ’BJSM_binary’ object                                   |
+----------------------------------------------+----------------------------------------------------------------------------+
| for class ’BJSM_binary_dose’                 | Summarize and print ’BJSM_binary_dose’ object                              |
+----------------------------------------------+----------------------------------------------------------------------------+
| for class ’BJSM_c’                           | Summarize and print ’BJSM_c’ object                                        |
+----------------------------------------------+----------------------------------------------------------------------------+
| for class ’group_seq’                        | Summarize and print ’group_seq’ object                                     |
+----------------------------------------------+----------------------------------------------------------------------------+
| for class ’LPJSM_binary’                     | Summarize and print ’LPJSM_binary’ object                                  |
+----------------------------------------------+----------------------------------------------------------------------------+
| for class ’sim_group_seq’                    | Summarize and print ’sim_group_seq’ object                                 |
+==============================================+============================================================================+

# The `snSMART` package

The `snSMART` package can be downloaded from the Comprehensive R Archive Network (CRAN) at https://cran.r-project.org/web/packages/snSMART/index.html. Please install the [JAGS library](https://sourceforge.net/projects/mcmc-jags/) [@plummer2003jags] before using this package. This package provides a sample size calculation and multiple methods of trial data analysis for various snSMART designs. We also summarized the functionality of all the `snSMART` functions included in this package in Table 1. 

After installing package `snSMART`, load the package:

```r
library("snSMART")
```
We assume that we have data from an snSMART in a matrix with 4 columns: ``treatment_stageI, response_stageI, treatment_stageII``, and ``response_stageII``. The rows of the dataset correspond to each participant in the trial. For example, here we use a dataset of 30 total individuals. 

```r
head(data_dose)
```
```
  response_stageI treatment_stageI response_stageII treatment_stageII
1               1                1                1                 3
2               0                1                1                 2
3               0                1                1                 2
4               0                1                1                 2
5               0                1                1                 2
6               0                1                1                 2
```
## Bayesian joint stage model (BJSM)

The ``BJSM_binary`` function uses the BJSM algorithm to estimate first-stage response rates by borrowing information from both stages. Users specify priors, MCMC details, and BJSM model type (six beta or two beta). The function creates a ``BJSM_binary`` object with posterior samples, linkage parameter estimates, and treatment response rates.

To analyze trial data comparing dose levels with placebo, assume the prior distribution of \(\pi_P\) as \(Beta(3, 17)\), \(\beta_{jm}\) as \(Gamma(2, 2)\), and the treatment effect ratio as \(Normal(0.2, 100)\). Label placebo as 1, low dose as 2, and high dose as 3 in the dataset. The output is a ``BJSM_dose_binary`` object with posterior samples and estimates of linkage parameters, treatment response rates, and pairwise response rate differences.

```r
BJSM_dose_result <- BJSM_binary(
  data = data_dose, prior_dist = c("beta", "gamma"),
  pi_prior = c(3, 17), normal.par = c(0.2, 100), beta_prior = c(2, 2),
  n_MCMC_chain = 2, n.adapt = 1000, BURN.IN = 10000, MCMC_SAMPLE = 60000, ci = 0.95
)
summary(BJSM_dose_result)
```

The response rates for placebo, low dose and high dose are estimated to be 0.08 (95\% credible interval (CI): 0.02 - 0.15), 0.4 (95\% CI: 0.28 - 0.52), and 0.74 (95\% CI: 0.59 - 0.87) respectively. Other estimated outcomes are also clearly presented in the R output.

# Discussion

We introduced and demonstrated the `snSMART` package for analyzing snSMART data using Bayesian methods. BJSM is often more efficient, but frequentist methods are recommended for sensitivity analysis. The package will be updated with new designs and methods to aid in finding effective treatments for rare diseases.

# Acknowledgments
This work was supported by a Patient-Centered Outcomes Research Institute (PCORI) award (ME-1507-31108). We thank Boxian Wei, Yan-Cheng Chao, and Holly Hartman for contributing their original R code used in the creation of this package. 

# References
