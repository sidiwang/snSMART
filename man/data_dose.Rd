% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{data_dose}
\alias{data_dose}
\title{Data Dose Level}
\description{
sample dataset of snSMART (dose level treatment) with binary outcomes
}
\examples{
mydata = data_dose
BJSM_dose_result = BJSM_binary(data = data_dose, prior_dist = c("beta", "gamma"),
    pi_prior = c(3, 17), normal.par = c(0.2, 100), beta_prior = c(2, 2),
    n_MCMC_chain = 2, BURN.IN = 10000, MCMC_SAMPLE = 60000, ci = 0.95)

}
\keyword{data}
