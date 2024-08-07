% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/group_seq.R
\name{summary.group_seq}
\alias{summary.group_seq}
\title{Summarizing BJSM fits}
\usage{
\method{summary}{group_seq}(object, ...)
}
\arguments{
\item{object}{an object of class "\code{group_seq}", usually, a result of a call to \code{\link{group_seq}}}

\item{...}{further arguments. Not currently used.}
}
\value{
\describe{
\item{Treatment Effects Estimate}{a 3 x 5 matrix with columns for the estimated treatment effects, its standard error, coverage probability of its credible interval, lower bound for its credible interval and higher bound for its credible interval}
\item{Differences between Treatments}{a 3 x 5 matrix with columns for the estimated differences in treatment effects between two treatments, its standard error, coverage probability of its credible interval, lower bound and higher bound of the credible interval}
\item{Linkage Parameter Estimate}{a 2 x 5 matrix, if the two beta model is fitted, or a 6 x 5 matrix, if the six beta model is fitted, with columns for the estimated linkage parameters}
\item{Expected Response Rate of Dynamic Treatment Regimens (DTR)}{only when \code{DTR = TRUE}}
}
}
\description{
\code{summary} method for class "\code{group_seq}"
}
