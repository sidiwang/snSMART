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

As a novel approach to study treatments in small samples, snSMART, or small sample, sequential, multiple assignment, randomized trial, designs and methods have been developed in recent years. The design and methods of snSMARTs generally apply to any disorder or disease that affects a small group of people and remains stable over the duration of the trial. A two-stage snSMART design first randomizes all participants to one of the first-stage treatments, then conducts second stage randomization based on its first stage treatment outcome. The goal of an snSMART is to identify the superior first stage treatment using data from both stages. An snSMART design repeats the same first stage treatments in the second stage of the trial, and may allow responders to the first stage treatment to stay on that same treatment in the second stage. Thus, the snSMART design allows for the estimation of first stage treatment effects by utilizing data from both stages, which is helpful for achieving more precise estimates.

In the past few years, substantial progress has been made regarding the development of statistical methods for analyzing trial data from the following snSMART designs: an snSMART with three active treatments [@wei2018bayesian; @wei2020sample; @chao2020dynamic], a group sequential snSMART with three active treatments [@chao2020bayesian], an snSMART with placebo and two dose levels of one treatment [@fang2021bayesian], and an snSMART with continuous outcomes [@hartman2021design]. However, there is a lack of statistical software to disseminate the methods. Our R package `snSMART` fills this gap by providing functions for calculating the required sample size for an snSMART with three active treatments, and analyzing trial data using both Bayesian and frequentist approaches. To the extent of our knowledge, there are no existing R packages that provide similar functions. 


## snSMART comparing two dose levels with placebo (P2D)
Here we explain one of the snSMART designs in detail. This snSMART design investigates the response rate of one experimental treatment at low and high doses compared with placebo [@fang2021bayesian]. In such an snSMART design (Figure \ref{fig:snSMART-dose}), participants are equally assigned to either receive placebo, low dose, or high dose in the first stage. Participants receive their initial treatment for a pre-specified amount of time until the measurement of their responses at the end of stage 1. In the second stage, all participants who received placebo or low dose in the first stage are re-randomized to either low or high dose regardless of their first stage response. Participants who responded to high dose are re-randomized between low and high dose, whereas those who did not respond to high dose receive high dose again in the second stage. The main goal of this snSMART is to estimate first stage low and high dose treatment response rates and compare them to placebo by modeling the pooled data from two stages.

![\label{fig:snSMART-dose}](dose_snSMART.png)

@fang2021bayesian used the Bayesian joint stage model (BJSM) from @wei2018bayesian in Equations \ref{eq:1} and \ref{eq:2} with six linkage parameters that varied by first stage treatments. For this design $m = P, L, H$ and $m' = L, H$ where P = placebo, L = low dose and H = high dose. The prior distribution for the response rate of placebo may be informed by natural history studies or previous trials and specified such that $\pi_P \sim Beta(\zeta_n, \eta_n)$. A weak tendency for the doses of the drug to have greater response rates than the effect of placebo can be assumed through $\log(\pi_L/\pi_P) \sim N(\mu, \sigma^2)$ and $\log(\pi_H/\pi_P) \sim N(\mu, \sigma^2)$. The prior distributions for the linkage parameters may vary with guidance of specifying $\beta_{0m}, \beta_{1m} \sim Gamma(\omega, \psi)$. As with the other designs and analytic methods, posterior samples are drawn through MCMC sampling.

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

Table 1: Summary of the functionality of the snSMART package.
+----------------------------------------------+---------------------------------------------------------------------------------------------------------+
| Function                                     | Description                                                                                             | 
|                                              |                                                                                                         |       
+:=============================================+:========================================================================================================+
| *BJSM functions*                             |                                                                                                         | 
+----------------------------------------------+---------------------------------------------------------------------------------------------------------+
| BJSM_binary                                  | BJSM binary (3AT or P2D snSMART)                                                                        |
| BJSM_c                                       | BJSM (3AT snSMART with a mapping function and continuous outcome)                                       |
| group_seq                                    | BJSM (interim analysis and final analysis of group sequential 3AT snSMART)                              |
+----------------------------------------------+---------------------------------------------------------------------------------------------------------+
| *Frequentist functions*                      |                                                                                                         | 
+----------------------------------------------+---------------------------------------------------------------------------------------------------------+
| LPJSM_binary                                 | LPJSM (3AT or P2D snSMART)                                                                              |
+----------------------------------------------+---------------------------------------------------------------------------------------------------------+
| *Sample size calculation*                    |                                                                                                         | 
+----------------------------------------------+---------------------------------------------------------------------------------------------------------+
| sample_size                                  | 3AT snSMART sample size calculation                                                                     |
+----------------------------------------------+---------------------------------------------------------------------------------------------------------+
| *S3 summary and print methods*               |                                                                                                         | 
+----------------------------------------------+---------------------------------------------------------------------------------------------------------+
| for class ’BJSM_binary’                      | Summarize and print ’BJSM_binary’ object                                                                |
| for class ’BJSM_binary_dose’                 | Summarize and print ’BJSM_binary_dose’ object                                                           |
| for class ’BJSM_c’                           | Summarize and print ’BJSM_c’ object                                                                     |
| for class ’group_seq’                        | Summarize and print ’group_seq’ object                                                                  |
| for class ’LPJSM_binary’                     | Summarize and print ’LPJSM_binary’ object                                                               |
| for class ’sim_group_seq’                    | Summarize and print ’sim_group_seq’ object                                                              |
+==============================================+=========================================================================================================+

# The `snSMART` package

The `snSMART` package can be downloaded from the Comprehensive R Archive Network (CRAN) at https://cran.r-project.org/web/packages/snSMART/index.html. Please install the [JAGS library](https://sourceforge.net/projects/mcmc-jags/) [@plummer2003jags] before using this package. This package provides a sample size calculation and multiple methods of trial data analysis for various snSMART designs. A brief introduction of the snSMART design comparing two dose levels with placebo is provided in Section 2. We also summarized the functionality of all the `snSMART` functions included in this package in Table 1. The `snSMART` package includes the option to set the number of "adaptation iterations" in the functions that perform MCMC computation. This option aims to provide users with the complete functionality of the  ``jags.model`` function available in the `rjags` package [@rjags].

After installing package `snSMART`, load the package:

```r
library("snSMART")
```
We assume that we have data from an snSMART with three active treatments in a matrix with 4 columns: ``treatment_stageI, response_stageI, treatment_stageII``, and ``response_stageII``. The rows of the dataset correspond to each participant in the trial. For example, here we use a dataset of 30 total individuals. 

```r
head(data_binary)
```

## Bayesian joint stage model (BJSM)

The function ``BJSM_binary`` implements the BJSM algorithm that borrows information across both stages to estimate the first stage response rates of each treatment based on trial data. Users specify the prior distributions for all treatment response rates and linkage parameters; the details of MCMC simulations; and the type of BJSM model they want to implement, i.e., six beta model or two beta model. This function creates an object of class '``BJSM_binary``' that contains a list of components: posterior samples of linkage parameters and treatment response rates; estimates of linkage parameters, treatment response rates, and pairwise response rate differences of all treatments.

To analyze trial data collected from an snSMART that compares dose levels with placebo. Here we call the function assuming the prior distribution of $\pi_P$ being $Beta(3, 17)$, the prior distribution of $\beta_{jm}$ being $Gamma(2, 2)$ and the normal priors of the logarithm of treatment effect ratio being $Normal(0.2, 100)$. Under the dose level design scenario, users should label placebo treatment as 1, low dose treatment as 2 and high dose treatment as 3 in the trial dataset. Output under this setting is an object of class '``BJSM_dose_binary``' that contains a list of components: posterior samples of linkage parameters and treatment response rates; estimates of linkage parameters, treatment response rates, and pairwise response rate differences of all treatments. 

```r
BJSM_dose_result <- BJSM_binary(
  data = data_dose, prior_dist = c("beta", "gamma"),
  pi_prior = c(3, 17), normal.par = c(0.2, 100), beta_prior = c(2, 2),
  n_MCMC_chain = 2, n.adapt = 1000, BURN.IN = 10000, MCMC_SAMPLE = 60000, ci = 0.95
)
summary(BJSM_dose_result)
```

The response rates for placebo, low dose and high dose are estimated to be 0.08 (95\% credible interval (CI): 0.02 - 0.15), 0.4 (95\% CI: 0.28 - 0.52), and 0.74 (95\% CI: 0.59 - 0.87) respectively. Other estimated outcomes are also clearly presented in the R output above.

# Discussion
In this paper, we introduced and demonstrated the `snSMART` package to analyze data from one of the snSMART designs using Bayesian methods. All statistical methods covered in the package are listed in the table above. While the authors of these methods generally found the BJSM to produce more efficient estimates than the frequentist method `LPJSM`, it is recommended to use the frequentist approach as a sensitivity analysis to the Bayesian approach. This package predominantly presents analytic functions, however, these functions may also be useful in the design phase. To conduct simulation studies, users can simulate snSMART data and implement the analysis functions included in this package to assess the operating characteristics of the design.

This R package will continue to be updated as more snSMART designs and methods are developed. We hope that this package translates snSMART design and methods into finding more effective treatments for rare diseases. 

# Acknowledgments
This work was supported by a Patient-Centered Outcomes Research Institute (PCORI) award (ME-1507-31108). We thank Boxian Wei, Yan-Cheng Chao, and Holly Hartman for contributing their original R code used in the creation of this package. 

# References
