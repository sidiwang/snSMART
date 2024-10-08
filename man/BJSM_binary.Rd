% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BJSM_binary.R
\name{BJSM_binary}
\alias{BJSM_binary}
\alias{print.summary.BJSM_binary}
\alias{print.BJSM_binary}
\alias{print.summary.BJSM_dose_binary}
\alias{print.BJSM_dose_binary}
\title{BJSM for snSMART (3 active treatments/placebo and 2 dose level) with binary outcome}
\usage{
BJSM_binary(
  data,
  prior_dist,
  pi_prior,
  normal.par,
  beta_prior,
  n_MCMC_chain = 1,
  n.adapt,
  BURN.IN = 100,
  thin = 1,
  MCMC_SAMPLE,
  ci = 0.95,
  six = TRUE,
  DTR = TRUE,
  jags.model_options = NULL,
  coda.samples_options = NULL,
  verbose = FALSE,
  ...
)

\method{print}{summary.BJSM_binary}(x, ...)

\method{print}{BJSM_binary}(x, ...)

\method{print}{summary.BJSM_dose_binary}(x, ...)

\method{print}{BJSM_dose_binary}(x, ...)
}
\arguments{
\item{data}{trial data with 4 columns: \code{treatment_stageI, response_stageI,
treatment_stageII} and \code{response_stageII}. Missing data is allowed in stage 2.}

\item{prior_dist}{for 3 active treatment design: vector of three values
("prior distribution for \code{pi}", "prior distribution for \code{beta0}",
"prior distribution for \code{beta1}"). User can choose from "gamma", "beta", "pareto".
e.g. prior_dist = c("beta", "beta", "pareto"); for dose level design: vector of two
values ("prior distribution for \code{pi_P}", "prior distribution for \code{beta}")}

\item{pi_prior}{for 3 active treatment design: vector of six values (a, b, c, d, e, f),
where a and b are the parameter \code{a} and parameter \code{b} of the prior distribution
for \code{pi_1A}, c and d are the parameter \code{a} and parameter \code{b} of the prior
distribution for \code{pi_1B}, and e and f are the parameter \code{a} and parameter
\code{b} of the prior distribution for \code{pi_1C}. for dose level design: vector
of two values (a, b). \code{a} is the parameter \code{a} of the prior distribution for
\code{pi} (response rate) of placebo. \code{b} is the parameter \code{b} of the prior
distribution for \code{pi} of placebo. Please check the \code{Details} section for
more explanation}

\item{normal.par}{for dose level design: vector of two values (normal.mean, normal.var).
our function assumes that the logarithm of treatment effect ratio follows a Gaussian
prior distribution \eqn{N(\mu, \sigma^2)}, that is \eqn{log(\pi_L/\pi_P)~N(normal.mean, normal.var)},
and \eqn{log(\pi_H/\pi_P)~N(normal.mean, normal.var)}. \code{normal.mean} is the mean of
this Gaussian prior. \code{normal.var} is the variance of this Gaussian prior distribution}

\item{beta_prior}{for 3 active treatment design: vector of four values (a, b, c, d).
\code{a} is the value of parameter \code{a} of the prior distribution for linkage parameter
\code{beta_0} or \code{beta_0m}, \code{b} is the value of parameter \code{b} of the
prior distribution for linkage parameter \code{beta_0} or \code{beta_0m}. \code{c}
is the value of parameter \code{a} of the prior distribution for linkage parameter
\code{beta_1} or \code{beta_1m}. \code{d} is the value of parameter \code{b} of the
prior distribution for linkage parameter \code{beta_1} or \code{beta_1m}. for
dose level design: vector of two values (a, b). \code{a} is the parameter \code{a} of the
prior distribution for linkage parameter \code{beta}. \code{b} is the parameter b of
the prior distribution  for linkage parameter \code{beta}. Please check the \code{Details}
section for more explanation}

\item{n_MCMC_chain}{number of MCMC chains, default to 1.}

\item{n.adapt}{the number of iterations for adaptation}

\item{BURN.IN}{number of burn-in iterations for MCMC}

\item{thin}{thinning interval for monitors}

\item{MCMC_SAMPLE}{number of iterations for MCMC}

\item{ci}{coverage probability for credible intervals, default = 0.95}

\item{six}{TRUE or FALSE. If TRUE, will run the six beta model (allow for estimating
\code{beta_0m} and \code{beta_1m} values that differ among different treatments m), if FALSE
will run the two beta model. default = TRUE. Only need to specify this for 3 active
treatment design.}

\item{DTR}{TRUE or FALSE. If TRUE, will also return the expected response rate of
dynamic treatment regimens. default = TRUE. Only need to specify this for 3 active
treatment design.}

\item{jags.model_options}{a list of optional arguments that are passed to \code{jags.model()} function.}

\item{coda.samples_options}{a list of optional arguments that are passed to \code{coda.samples()} function.}

\item{verbose}{TRUE or FALSE. If FALSE, no function message and progress bar will be
printed.}

\item{...}{further arguments. Not currently used.}

\item{x}{object to summarize.}
}
\value{
\item{posterior_sample}{an \code{mcmc.list} object generated through the \code{coda.samples()} function,
which includes posterior samples of the link parameters and response rates generated through the MCMC
process}
\item{pi_hat_bjsm}{estimate of response rate/treatment effect}

\item{se_hat_bjsm}{standard error of the response rate}

\item{ci_pi_A(P), ci_pi_B(L), ci_pi_C(H)}{x\% credible intervals for treatment A(P), B(L), C(H)}

\item{diff_AB(PL), diff_BC(LH). diff_AC(PH)}{estimate of differences between
treatments A(P) and B(L), B(L) and C(H), A(P) and C(H)}

\item{ci_diff_AB(PL), ci_diff_BC(LH), ci_diff_AC(PH)}{x\% credible intervals
for the estimated differences between treatments A(P) and B(L), B(L) and C(H), A(P) and C(H)}

\item{se_AB(PL), se_BC(LH), se_AC(PH)}{standard error for the estimated differences
between treatments A(P) and B(L), B(L) and C(H), A(P) and C(H)}

\item{beta0_hat, beta1_hat}{linkage parameter \code{beta0} and \code{beta1} estimates}

\item{se_beta0_hat, se_beta1_hat}{standard error of the estimated value of linkage
parameter \code{beta0} and \code{beta1}}

\item{ci_beta0_hat, ci_beta1_hat}{linkage parameter \code{beta0} and \code{beta1}
credible interval}

\item{pi_DTR_est}{expected response rate of dynamic treatment regimens (DTRs)}

\item{pi_DTR_se}{standard error for the estimated DTR response rate}

\item{ci_pi_AB, ci_pi_AC, ci_pi_BA, ci_pi_BC, ci_pi_CA, ci_pi_CB}{x\% credible intervals for the estimated DTR response rate}
}
\description{
This function implements the BJSM (Bayesian Joint Stage Modeling) method which
borrows information across both stages to estimate the individual response rate
of each treatment/dose level in a snSMART design with binary outcomes.
}
\details{
For \code{gamma} distribution, \code{prior.a} is the shape parameter \code{r}, \code{prior.b} is the rate parameter \code{lambda}. For \code{beta} distribution, \code{prior.a} is the shape parameter \code{a}, \code{prior.b} is the shape parameter \code{b}.
For \code{pareto} distribution, \code{prior.a} is the scale parameter \code{alpha}, \code{prior.b} is the shape parameter \code{c} (see jags user manual).

The individual response rate is regarded as a permanent feature of the treatment. The second stage outcome is modeled conditionally on the first stage results linking the first and
second stage response probabilities through linkage parameters. The first stage response rate is denoted as \eqn{\pi_m} for treatment \eqn{m}. In the two \eqn{\beta} model, the second stage response rate for first stage responders is equal to \eqn{\beta_1\pi_m}. For nonresponders to treatment \eqn{m} in the first stage who
receive treatment \eqn{m'} in the second the stage, the second stage response rate in the second stage is equal to \eqn{\beta_0\pi_{m'}}. In the six \eqn{\beta} model, the second stage response rate of the first stage responders to treatment m is denoted by \eqn{\beta_{1m}\pi_m}, and the second stage response rate of the non-responders
to first stage treatment $m$ who receive treatment \eqn{m'} in the second stage is denoted by \eqn{\beta_{0m}\pi_{m'}}. All the \eqn{\beta}s are linkage parameters.

Please refer to the paper listed under \code{reference} section for standard snSMART trial design and detailed definition of parameters.

Note that this package does not include the JAGS library, users need to install JAGS separately. Please check this page for more details: \url{https://sourceforge.net/projects/mcmc-jags/}
}
\examples{
mydata <- data_binary

BJSM_result <- BJSM_binary(
  data = mydata, prior_dist = c("beta", "beta", "pareto"),
  pi_prior = c(0.4, 1.6, 0.4, 1.6, 0.4, 1.6), beta_prior = c(1.6, 0.4, 3, 1),
  n_MCMC_chain = 1, n.adapt = 1000, MCMC_SAMPLE = 2000, ci = 0.95,
  six = TRUE, DTR = TRUE, verbose = FALSE
)

\donttest{
BJSM_result2 <- BJSM_binary(
  data = mydata, prior_dist = c("beta", "beta", "pareto"),
  pi_prior = c(0.4, 1.6, 0.4, 1.6, 0.4, 1.6), beta_prior = c(1.6, 0.4, 3, 1),
  n_MCMC_chain = 1, n.adapt = 10000, MCMC_SAMPLE = 60000, ci = 0.95,
  six = FALSE, DTR = FALSE, verbose = FALSE
)

summary(BJSM_result)
summary(BJSM_result2)
}

data <- data_dose
BJSM_dose_result <- BJSM_binary(
  data = data_dose, prior_dist = c("beta", "gamma"),
  pi_prior = c(3, 17), normal.par = c(0.2, 100), beta_prior = c(2, 2),
  n_MCMC_chain = 2, n.adapt = 1000, MCMC_SAMPLE = 6000, ci = 0.95, verbose = FALSE
)

summary(BJSM_dose_result)

}
\references{
Wei, B., Braun, T.M., Tamura, R.N. and Kidwell, K.M., 2018. A Bayesian analysis of small n sequential multiple assignment randomized trials (snSMARTs).
Statistics in medicine, 37(26), pp.3723-3732. \doi{10.1002/sim.7900}

Chao, Y.C., Trachtman, H., Gipson, D.S., Spino, C., Braun, T.M. and Kidwell, K.M., 2020. Dynamic treatment regimens in small n, sequential, multiple assignment, randomized trials: An application in focal segmental glomerulosclerosis. Contemporary clinical trials, 92, p.105989. \doi{10.1016/j.cct.2020.105989}

Fang, F., Hochstedler, K.A., Tamura, R.N., Braun, T.M. and Kidwell, K.M., 2021. Bayesian methods to compare dose levels with placebo in a small n,
sequential, multiple assignment, randomized trial. Statistics in Medicine, 40(4), pp.963-977. \doi{10.1002/sim.8813}
}
\seealso{
\code{\link{LPJSM_binary}} \cr
\code{\link{sample_size}}
}
