% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{groupseqDATA_look1}
\alias{groupseqDATA_look1}
\title{Group sequential data look 1}
\format{
This data frame contains the following columns:
\describe{
\item{time.1st.trt}{first treatment time}
\item{time.1st.resp}{first response time}
\item{time.2nd.trt}{second treatment time}
\item{time.2nd.resp}{second response time}
\item{trt.1st}{treatment arm for first treatment}
\item{resp.1st}{response for first treatment}
\item{trt.2nd}{treatment arm for second treatment}
\item{resp.2nd}{response for second treatment}
}
}
\usage{
groupseqDATA_look1
}
\description{
sample synthetic dataset of group sequential trial design snSMART, can be used for interim analysis
}
\examples{
mydata <- groupseqDATA_look1

result1 <- group_seq(
  data = mydata, interim = TRUE, drop_threshold_pair = c(0.5, 0.4),
  prior_dist = c("beta", "beta", "pareto"), pi_prior = c(0.4, 1.6, 0.4, 1.6, 0.4, 1.6),
  beta_prior = c(1.6, 0.4, 3, 1), MCMC_SAMPLE = 6000, n.adapt = 1000, n_MCMC_chain = 1
)

}
\keyword{data}
