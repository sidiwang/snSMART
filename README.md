
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
    ##  /var/folders/5_/sp2r7r_s5snf4hq634xmk_0r0000gn/T//RtmpxqLzkg/downloaded_packages

``` r
library(snSMART)
```

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
    ## trtA 0.3724607 0.09596326 0.95 0.2529876 0.6262884
    ## trtB 0.3937709 0.12073746 0.95 0.1868479 0.6406974
    ## trtC 0.4490142 0.08668038 0.95 0.3508578 0.6316868
    ## 
    ## Differences between Treatments:
    ##           Estimate Std.Error C.I.     CI low    CI high
    ## diffAB -0.02131011 0.1671619 0.95 -0.2897015 0.32646743
    ## diffBC -0.05524338 0.1671081 0.95 -0.3362856 0.21362700
    ## diffAC -0.07655348 0.1023327 0.95 -0.3605940 0.08820396
    ## 
    ## Linkage Parameter Estimate:
    ##         Estimate Std. Error C.I.    CI low   CI high
    ## beta0A 0.9675508 0.03408395 0.95 0.8698091 0.9970421
    ## beta0B 0.9906452 0.01737068 0.95 0.9412510 0.9996143
    ## beta0C 0.7724702 0.25633534 0.95 0.2095624 0.9925097
    ## beta1A 1.6468066 0.46080230 0.95 1.0075555 2.7302687
    ## beta1B 1.6115237 0.49898527 0.95 1.0037997 3.1514827
    ## beta1C 1.7958893 0.38568874 0.95 1.2192022 2.4783575
    ## 
    ## Expected Response Rate of Dynamic Treatment Regimens (DTR):
    ##         Estimate Std. Error C.I.    CI low   CI high
    ## rep_AB 0.4767614 0.09769684 0.95 0.3093855 0.6403936
    ## rep_AC 0.5061297 0.10840983 0.95 0.3193806 0.7170460
    ## rep_BA 0.4885528 0.12068662 0.95 0.3309846 0.7400024
    ## rep_BC 0.5353517 0.11543082 0.95 0.3868979 0.7559936
    ## rep_CA 0.5127254 0.09818464 0.95 0.3216308 0.6994843
    ## rep_CB 0.5291309 0.10477125 0.95 0.3362417 0.7004188

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
library("EnvStats")
```

    ## 
    ## Attaching package: 'EnvStats'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     predict, predict.lm

    ## The following object is masked from 'package:base':
    ## 
    ##     print.default

``` r
sampleSize <- sample_size(
  pi = c(0.7, 0.5, 0.25), beta1 = 1.4, beta0 = 0.5,
  coverage = 0.9, power = 0.8, mu = c(0.65, 0.55, 0.25),
  n = c(4, 2, 3)
)
```

    ## With given settings, the estimated sample size per arm for an snSMART is: 34
    ## This implies that for an snSMART with sample size of 34 per arm (102 in total for three agents):
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
    ## trtA   0.3012    0.04875 0.95 0.2047  0.3950
    ## trtB   0.4729    0.03947 0.95 0.3960  0.5503
    ## trtC   0.6751    0.04099 0.95 0.5925  0.7542
    ## 
    ## Differences between Treatments:
    ##        Estimate Std.Error C.I.  CI low  CI high
    ## diffAB  -0.1717   0.06239 0.95 -0.2929 -0.04934
    ## diffBC  -0.2022   0.05592 0.95 -0.3094 -0.08953
    ## diffAC  -0.3739   0.06275 0.95 -0.4968 -0.25160
    ## 
    ## Linkage Parameter Estimate:
    ##        Estimate Std. Error C.I. CI low CI high
    ## beta0A   0.8739    0.10955 0.95 0.6669   1.000
    ## beta0B   0.7238    0.11415 0.95 0.5443   1.000
    ## beta0C   0.8739    0.12504 0.95 0.6241   1.000
    ## beta1A   1.4726    0.36688 0.95 1.0000   2.175
    ## beta1B   1.3637    0.17166 0.95 1.0262   1.678
    ## beta1C   1.4671    0.09288 0.95 1.2916   1.652

`BJSM_c`: Below, we call the function assuming the mean and standard
deviation of the normal prior being 50 and 50 for all three treatments,
and the standard deviation of the prior distribution of $\phi_3$ being
20. The number of MCMC chain is set to 1 with 1,000 adaptation
iterations and 5,000 total iterations.

``` r
BJSM_result <- BJSM_c(
  data = trialDataMF, xi_prior.mean = c(50, 50, 50), xi_prior.sd = c(50, 50, 50),
  phi3_prior.sd = 20, n_MCMC_chain = 1, n.adapt = 1000, MCMC_SAMPLE = 5000,
  ci = 0.95, n.digits = 5
)
summary(BJSM_result)
```

    ## 
    ## Parameter Estimation:
    ##            Estimate   CI     CI_low   CI_high
    ## phi1        0.17753 0.95  5.414e-05  0.376034
    ## phi3        4.01432 0.95  2.706e+00  5.339355
    ## rho[1,1,1]  0.00976 0.95  4.846e-03  0.014915
    ## rho[2,1,1]  0.06318 0.95  3.533e-02  0.094202
    ## rho[1,2,1] -0.00261 0.95 -6.442e-03  0.001507
    ## rho[2,2,1] -0.06957 0.95 -1.062e-01 -0.037752
    ## rho[1,1,2] -0.00261 0.95 -6.442e-03  0.001507
    ## rho[2,1,2] -0.06957 0.95 -1.062e-01 -0.037752
    ## rho[1,2,2]  0.01123 0.95  5.720e-03  0.017143
    ## rho[2,2,2]  0.09072 0.95  4.980e-02  0.134680
    ## xi_[1]     51.09928 0.95  4.718e+01 54.731318
    ## xi_[2]     62.04339 0.95  5.837e+01 65.903243
    ## xi_[3]     68.98975 0.95  6.531e+01 72.827499

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
