
# snSMART

<!-- badges: start -->

![cran-badge example for snSMART
package](http://www.r-pkg.org/badges/version/snSMART)
[![R-CMD-check](https://github.com/sidiwang/snSMART/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sidiwang/snSMART/actions/workflows/R-CMD-check.yaml)
![RStudio CRAN
downloads](http://cranlogs.r-pkg.org/badges/grand-total/snSMART)
![RStudio CRAN monthly
downloads](http://cranlogs.r-pkg.org/badges/snSMART) [![License: GPL
v2](https://img.shields.io/badge/License-GPL_v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/)
[![Codecov test
coverage](https://codecov.io/gh/sidiwang/snSMART/branch/main/graph/badge.svg)](https://app.codecov.io/gh/sidiwang/snSMART?branch=main)
<!-- badges: end -->

The aim of the **snSMART** R package is to consolidate data simulation,
sample size calculation and analysis functions for several snSMART
(small sample sequential, multiple assignment, randomized trial) designs
under one library.

An snSMART is a multi-stage trial design where for a two-stage design,
randomization in the second stage depends on the outcome to first stage
treatment. snSMART designs require that the same outcome is measured at
the end of the first stage and at the end of the second stage.
Additionally, the length of the first stage of the trial must be the
same amount of time as the length for the second stage. snSMARTs are
motivated by obtaining more information from a small sample of
individuals with the primary goal to identify the superior first stage
treatment or dosage level using both stages of data. Data are shared
across the two stages of the snSMART design to more precisely estimate
the effect of the treatments given in the first stage.

## Installation

You can install the package from CRAN:

``` r
install.packages("snSMART", repos = "http://cran.us.r-project.org")
```

    ## 
    ## The downloaded binary packages are in
    ##  /var/folders/5_/sp2r7r_s5snf4hq634xmk_0r0000gn/T//RtmpNRlYCm/downloaded_packages

``` r
library(snSMART)
```

    ## Loading required package: EnvStats

    ## 
    ## Attaching package: 'EnvStats'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     predict, predict.lm

Or get the development version from GitHub:

``` r
# Install devtools first if you haven't done so
library(devtools)
# install snSMART
devtools::install_github("sidiwang/snSMART")
library(snSMART)
```

## snSMART designs and functions covered in this package

- snSMART with binary outcome (3 active treatments or placebo and 2 dose
  level, non-responders re-randomized)
  - Bayesian Joint Stage Model (BJSM) analysis function - `BJSM_binary`
  - Log Poisson Joint Stage Model (JSRM) analysis function -
    `JSRM_binary`
  - sample size calculation function - `sample_size`
  - Group Sequential snSMART BJSM analysis function - `group_seq`
- snSMART with mapping function (3 active treatments, re-randomization
  depends on continuous outcome at stage 1; continuous outcomes)
  - BJSM analysis function - `BJSM_c`

## Example

`BJSM_binary`: We call the `BJSM_binary` function using data from an
snSMART with 30 total individuals. We assumed a six beta model with the
priors for $\pi_A, \pi_B$ and $\pi_C$ being $Beta(0.4, 1.6)$, the prior
for $\beta_{0m}$ being $Beta(1.6, 0.4)$, and the prior for $\beta_{1m}$
being $Pareto(3, 1)$. The coverage probability for credible intervals is
set to 0.95, and the expected response rate of DTR will also be
calculated.

``` r
mydata <- data_binary
BJSM_result <- BJSM_binary(
  data = mydata, prior_dist = c("beta", "beta", "pareto"),
  pi_prior = c(0.4, 1.6, 0.4, 1.6, 0.4, 1.6), beta_prior = c(1.6, 0.4, 3, 1),
  n_MCMC_chain = 1, n.adapt = 1000, MCMC_SAMPLE = 20, ci = 0.95,
  six = TRUE, DTR = TRUE, verbose = FALSE
)
summary(BJSM_result)
```

    ## 
    ## Treatment Effects Estimate:
    ##       Estimate Std. Error C.I.    CI low   CI high
    ## trtA 0.4046939 0.10427533 0.95 0.2103586 0.6583262
    ## trtB 0.3900296 0.07865828 0.95 0.2501780 0.5374523
    ## trtC 0.6180690 0.08709180 0.95 0.4925770 0.7772426
    ## 
    ## Differences between Treatments:
    ##           Estimate Std.Error C.I.     CI low     CI high
    ## diffAB  0.01466433 0.1287876 0.95 -0.2202233  0.32554853
    ## diffBC -0.22803936 0.1271915 0.95 -0.5270646 -0.06082518
    ## diffAC -0.21337503 0.1493975 0.95 -0.4840569  0.04220399
    ## 
    ## Linkage Parameter Estimate:
    ##         Estimate Std. Error C.I.    CI low   CI high
    ## beta0A 0.9595907 0.04719401 0.95 0.8544802 0.9998446
    ## beta0B 0.9575869 0.05057664 0.95 0.8615737 0.9997691
    ## beta0C 0.6930520 0.24176704 0.95 0.3444717 0.9917792
    ## beta1A 1.5171964 0.38056743 0.95 1.0808245 2.3459093
    ## beta1B 1.7896323 0.56996358 0.95 1.1208255 2.9052050
    ## beta1C 1.3295262 0.23142426 0.95 1.0005293 1.8001938
    ## 
    ## Expected Response Rate of Dynamic Treatment Regimens (DTR):
    ##         Estimate Std. Error C.I.    CI low   CI high
    ## rep_AB 0.4697516 0.09275624 0.95 0.2917702 0.6296532
    ## rep_AC 0.6015571 0.07162453 0.95 0.4451989 0.7184334
    ## rep_BA 0.5012073 0.08332682 0.95 0.3709581 0.6217467
    ## rep_BC 0.6294671 0.09964662 0.95 0.4624197 0.7828289
    ## rep_CA 0.6175090 0.11258901 0.95 0.3954040 0.8142629
    ## rep_CB 0.6065809 0.10694339 0.95 0.4134393 0.7979918

`LPJSM_binary`: Here, we call the LPJSM_binary mirroring our example for
the BJSM_binary above.

``` r
LPJSM_result <- LPJSM_binary(data = data_binary, six = TRUE, DTR = TRUE)
summary(LPJSM_result)
```

    ## 
    ## GEE output:
    ## 
    ## Call:
    ## geepack::geeglm(formula = Y ~ alphaA + alphaB + alphaC + gamma1A + 
    ##     gamma2A + gamma1B + gamma2B + gamma1C + gamma2C - 1, family = poisson(link = "log"), 
    ##     data = geedata, id = ptid, corstr = "independence")
    ## 
    ##  Coefficients:
    ##         Estimate  Std.err     Wald Pr(>|W|)    
    ## alphaA   -1.2155   0.4030    9.096  0.00256 ** 
    ## alphaB   -0.9845   0.3403    8.370  0.00381 ** 
    ## alphaC   -0.8444   0.3036    7.737  0.00541 ** 
    ## gamma1A   1.2155   0.4030    9.096  0.00256 ** 
    ## gamma2A   0.7029   0.3664    3.680  0.05508 .  
    ## gamma1B   0.9845   0.3403    8.370  0.00381 ** 
    ## gamma2B   0.8271   0.3497    5.594  0.01803 *  
    ## gamma1C   0.8444   0.3036    7.737  0.00541 ** 
    ## gamma2C -42.2449   0.6085 4820.409  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation structure = independence 
    ## Estimated Scale Parameters:
    ## 
    ##             Estimate Std.err
    ## (Intercept)   0.4096  0.1259
    ## Number of clusters:   30  Maximum cluster size: 2 
    ## 
    ## Treatment Effect Estimate:
    ##      Estimate Std. Error
    ## trtA   0.2966     0.1195
    ## trtB   0.3736     0.1271
    ## trtC   0.4298     0.1305
    ## 
    ## Expected Response Rate of Dynamic Treatment Regimens (DTR):
    ##        Estimate Std. Error
    ## rep_AB   0.8274     0.1328
    ## rep_AC   0.9072     0.2341
    ## rep_BA   0.7984     0.1403
    ## rep_BC   0.9892     0.1427
    ## rep_CA   0.4298     0.1305
    ## rep_CB   0.4298     0.1305

`sample_size`: In this example, we call the function to request the
sample size needed per arm with the following assumptions: the response
rates for treatments A, B, and C are 0.7, 0.5 and 0.25, respectively;
$\beta_1$ is assumed to be 1.4; $\beta_0$ is assumed to be 0.5; the
coverage rate for the posterior difference of top two treatments is set
to 0.9; the ‘power’ is set to 0.8; the prior sample size is 4 for
treatment A, 2 for treatment B and 3 for treatment C; the prior means
are 0.65 for treatment A, 0.55 for treatment B and 0.25 for treatment C

``` r
sampleSize <- sample_size(
  pi = c(0.7, 0.5, 0.25), beta1 = 1.4, beta0 = 0.5,
  coverage = 0.9, power = 0.8, mu = c(0.65, 0.55, 0.25),
  n = c(4, 2, 3)
)
```

    ## With given settings, the estimated sample size per arm for an snSMART is: 34
    ## This implies that for an snSMART with sample size of 34 per arm (102 in total for three treatments):
    ## The probability of successfully identifying the best treatment is 0.8 when the difference of response rates between the best and second best treatment is at least 0.2, and the response rate of the best treatment is 0.7

`group_seq`: This function either outputs which treatment arm should be
dropped if or provides a full BJSM analysis based on the complete
dataset if . The dataset and are provided in the package for
illustration purpose.

``` r
result1 <- group_seq(
  data = groupseqDATA_look1, interim = TRUE, drop_threshold_pair = c(0.5, 0.4),
  prior_dist = c("beta", "beta", "pareto"), pi_prior = c(0.4, 1.6, 0.4, 1.6, 0.4, 1.6),
  beta_prior = c(1.6, 0.4, 3, 1), MCMC_SAMPLE = 6000, n.adapt = 1000, n_MCMC_chain = 1
)
```

    ## 
    ## Interim Analysis Outcome:

    ## Threshold tau_l is set to:

    ## 0.5

    ## 
    ## Threshold psi_l is set to:

    ## 0.4

    ## 
    ## Step 1: Arm C's interim posterior probability of having the greatest response is bigger than threshold

    ## 0.5

    ## 

    ## Step 2: Arm A's interim posterior probability of having the lowest response is higher

    ## Arm A is dropped

    ## 

``` r
result2 <- group_seq(
  data = groupseqDATA_full, interim = FALSE, prior_dist = c("beta", "beta", "pareto"),
  pi_prior = c(0.4, 1.6, 0.4, 1.6, 0.4, 1.6), beta_prior = c(1.6, 0.4, 3, 1),
  MCMC_SAMPLE = 60000, BURN.IN = 10000, n_MCMC_chain = 1, ci = 0.95, DTR = TRUE
)
summary(result2)
```

    ## 
    ## Treatment Effects Estimate:
    ##      Estimate Std. Error C.I. CI low CI high
    ## trtA   0.3009    0.04858 0.95 0.2073  0.3965
    ## trtB   0.4734    0.03911 0.95 0.3974  0.5503
    ## trtC   0.6769    0.04219 0.95 0.5959  0.7610
    ## 
    ## Differences between Treatments:
    ##        Estimate Std.Error C.I.  CI low  CI high
    ## diffAB  -0.1724   0.06193 0.95 -0.2945 -0.05246
    ## diffBC  -0.2035   0.05663 0.95 -0.3137 -0.09219
    ## diffAC  -0.3759   0.06264 0.95 -0.4936 -0.24832
    ## 
    ## Linkage Parameter Estimate:
    ##        Estimate Std. Error C.I. CI low CI high
    ## beta0A   0.8728    0.10923 0.95 0.6662   1.000
    ## beta0B   0.7207    0.11468 0.95 0.5401   1.000
    ## beta0C   0.8729    0.12545 0.95 0.6231   1.000
    ## beta1A   1.4711    0.36005 0.95 1.0000   2.164
    ## beta1B   1.3608    0.16853 0.95 1.0301   1.673
    ## beta1C   1.4633    0.09484 0.95 1.2800   1.641

`BJSM_c`: Below, we call the function assuming the mean and standard
deviation of the normal prior being 50 and 50 for all three treatments,
and the standard deviation of the prior distribution of $\phi_3$ being
20. The number of MCMC chain is set to 1 with 1,000 adaptation
iterations and 5,000 total iterations.

``` r
BJSM_result <- BJSM_c(
  data = trialDataMF, xi_prior.mean = c(50, 50, 50), xi_prior.sd = c(50, 50, 50),
  phi3_prior.sd = 20, n_MCMC_chain = 1, n.adapt = 1000, MCMC_SAMPLE = 5000, BURN.IN = 1000,
  ci = 0.95, n.digits = 5
)
summary(BJSM_result)
```

    ## 
    ## Parameter Estimation:
    ##         Estimate   CI     CI_low  CI_high
    ## V1[1,1]  110.526 0.95  6.351e+01 167.1243
    ## V1[2,1]   84.797 0.95  4.640e+01 130.1543
    ## V1[1,2]   84.797 0.95  4.640e+01 130.1543
    ## V1[2,2]   77.082 0.95  4.386e+01 116.8413
    ## V2[1,1]  122.695 0.95  6.342e+01 195.2136
    ## V2[2,1]   28.201 0.95 -2.200e+01  74.1924
    ## V2[1,2]   28.201 0.95 -2.200e+01  74.1924
    ## V2[2,2]  107.332 0.95  5.682e+01 173.0125
    ## phi1       0.184 0.95  3.294e-05   0.3795
    ## phi3       4.004 0.95  2.738e+00   5.3468
    ## xi_[A]    51.126 0.95  4.724e+01  54.9047
    ## xi_[B]    62.030 0.95  5.836e+01  65.7607
    ## xi_[C]    69.035 0.95  6.526e+01  72.7477

This R package will continue to be updated as more snSMART designs and
methods are developed. We hope that this package translates snSMART
design and methods into finding more effective treatments for rare
disease.

## References

Chao, Y.C., Trachtman, H., Gipson, D.S., Spino, C., Braun, T.M. and
Kidwell, K.M., 2020. Dynamic treatment regimens in small n, sequential,
multiple assignment, randomized trials: An application in focal
segmental glomerulosclerosis. Contemporary clinical trials, 92,
p.105989.

Chao, Y.C., Braun, T.M., Tamura, R.N. and Kidwell, K.M., 2020. A
Bayesian group sequential small n sequential multiple‐assignment
randomized trial. Journal of the Royal Statistical Society: Series C
(Applied Statistics), 69(3), pp.663-680.

Fang, F., Hochstedler, K.A., Tamura, R.N., Braun, T.M. and Kidwell,
K.M., 2021. Bayesian methods to compare dose levels with placebo in a
small n, sequential, multiple assignment, randomized trial. Statistics
in Medicine, 40(4), pp.963-977.

Hartman, H., Tamura, R.N., Schipper, M.J. and Kidwell, K.M., 2021.
Design and analysis considerations for utilizing a mapping function in a
small sample, sequential, multiple assignment, randomized trials with
continuous outcomes. Statistics in Medicine, 40(2), pp.312-326.

Wei, B., Braun, T.M., Tamura, R.N. and Kidwell, K., 2020. Sample size
determination for Bayesian analysis of small n sequential, multiple
assignment, randomized trials (snSMARTs) with three agents. Journal of
Biopharmaceutical Statistics, 30(6), pp.1109-1120.

Wei, B., Braun, T.M., Tamura, R.N. and Kidwell, K.M., 2018. A Bayesian
analysis of small n sequential multiple assignment randomized trials
(snSMARTs). Statistics in medicine, 37(26), pp.3723-3732.
