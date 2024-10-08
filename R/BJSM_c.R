#' BJSM continuous (snSMART with three active treatments and a continuous outcome design)
#'
#' BJSM (Bayesian Joint Stage Modeling) method that borrows information across both stages
#' to estimate the individual response rate of each treatment (with continuous
#' outcome and a mapping function).
#'
#' @param data trial ddatset with columns: \code{id, trt1} (treatment 1), \code{stage1outcome, stay}
#'  (stay = 1 if patient stay on the same treatment in stage 2, otherwise stay = 0),
#'  \code{trt2} (treatment 2), \code{stage2outcome}
#' @param xi_prior.mean a 3-element vector of mean of the prior distributions
#'  (normal distribution) for \code{xi}s (treatment effect). Please check the `Details`
#'  section for more explaination
#' @param xi_prior.sd a 3-element vector of standard deviation of the prior distributions
#'  (normal distribution) for \code{xi}s (treatment effect). Please check the `Details`
#'  section for more explaination
#' @param phi3_prior.sd standard deviation of the prior distribution (folded normal
#'  distribution) of \code{phi3} (if the patient stays on the same treatment, \code{phi3}
#'  is the cumulative effect of stage 1 that occurs on the treatment longer term).
#'  Please check the `Details` section for more explaination
#' @param n_MCMC_chain number of MCMC chains, default to 1
#' @param n.adapt the number of iterations for adaptation
#' @param BURN.IN number of burn-in iterations for MCMC
#' @param thin thinning interval for monitors
#' @param MCMC_SAMPLE number of iterations for MCMC
#' @param ci coverage probability for credible intervals, default = 0.95
#' @param n.digits number of digits to keep in the final estimation of treatment effect
#' @param x an object of class "`BJSM_c`", usually, a result of a call to \code{\link{BJSM_c}}
#' @param verbose TRUE or FALSE. If FALSE, no function message and progress bar will be
#'  printed.
#' @param jags.model_options a list of optional arguments that are passed to \code{jags.model()} function.
#' @param coda.samples_options a list of optional arguments that are passed to \code{coda.samples()} function.

#'
#' @details
#' section 2.2.1 and 2.2.2 of the paper listed under `reference` provides a detailed
#' description of the assumptions and prior distributions of the model.
#'
#' Note that this package does not include the JAGS library, users need to install JAGS separately. Please check this page for more details: \url{https://sourceforge.net/projects/mcmc-jags/}
#' @return
#' \item{posterior_sample}{an \code{mcmc.list} object generated through the \code{coda.samples()} function,
#'    which includes posterior samples of the link parameters and response rates generated through the MCMC
#'    process}
#' \item{mean_estimate}{BJSM estimate of each parameter:
#'     \enumerate{
#'          \item \code{phi1} - lingering effect of the first treatment
#'          \item `phi3` - if the patient stays on the same treatment, \code{phi3} is the cumulative effect of stage 1 that occurs on the treatment longer term
#'          \item `xi_j` - the expected effect of treatment j, j = 1, 2, 3 in the first stage
#'          \item \code{V1},\code{V2} are the variance-covariance matrix of the multivariate distribution. \code{V1} is for patients who stay on the same treatment,
#'          and \code{V2} is for patients who switch treatments. This allows those who stay on the same treatment to have a different correlation between stage one
#'          stage two outcomes than those who switch treatments.
#'    }}
#' \item{ci_estimate}{x% credible interval for each parameter. By default round to
#'     2 decimal places, if more decimals are needed, please access the results by
#'     `[YourResultName]$ci_estimates$CI_low` or `[YourResultName]$ci_estimates$CI_high` }
#'
#' @examples
#' trialData <- trialDataMF
#'
#' BJSM_result <- BJSM_c(
#'   data = trialData, xi_prior.mean = c(50, 50, 50),
#'   xi_prior.sd = c(50, 50, 50), phi3_prior.sd = 20, n_MCMC_chain = 1,
#'   n.adapt = 1000, MCMC_SAMPLE = 5000, BURIN.IN = 1000, ci = 0.95, n.digits = 5, verbose = FALSE
#' )
#'
#' summary(BJSM_result)
#' print(BJSM_result)
#' @references
#' Hartman, H., Tamura, R.N., Schipper, M.J. and Kidwell, K.M., 2021. Design and analysis considerations for utilizing a mapping function in a small sample,
#' sequential, multiple assignment, randomized trials with continuous outcomes. Statistics in Medicine, 40(2), pp.312-326. \doi{10.1002/sim.8776}
#'
#' @rdname BJSM_c
#' @export

BJSM_c <- function(data, xi_prior.mean, xi_prior.sd, phi3_prior.sd, n_MCMC_chain, n.adapt,
                   MCMC_SAMPLE, ci = 0.95, n.digits, thin = 1, BURN.IN = 100,
                   jags.model_options = NULL, coda.samples_options = NULL, verbose = FALSE, ...) {
  quiet <- FALSE
  progress.bar <- "text"

  if (verbose == FALSE) {
    quiet <- TRUE
    progress.bar <- "none"
  }

  # bug files written to temporary directory on function call to satisfy CRAN
  # requirements of not accessing user's system files

  # "csnSMART.bugs"
  csnSMART_file <- tempfile(fileext = ".bug")
  writeLines(csnSMART_text(), con = csnSMART_file)

  trialData <- data

  NUM_ARMS <- length(unique(trialData$trt1[!is.na(trialData$trt1)]))

  jag <- do.call(rjags::jags.model, c(list(
    file = csnSMART_file,
    data = list(
      n = length(unique(trialData$id)),
      ntrt = NUM_ARMS,
      Y = cbind(trialData$stage1outcome, trialData$stage2outcome),
      trt1 = trialData$trt1,
      trt2 = trialData$trt2,
      stay = trialData$stay,
      stay1 = trialData$stay + 1,
      xi_prior.mean = xi_prior.mean,
      xi_prior.sd = 1 / (xi_prior.sd^2),
      phi3_prior.sd = 1 / (phi3_prior.sd^2)
    ),
    n.chains = n_MCMC_chain, n.adapt = n.adapt, quiet = quiet, jags.model_options
  )))
  update(jag, BURN.IN, progress.bar = progress.bar)
  posterior_sample <- do.call(rjags::coda.samples, c(list(
    model = jag,
    variable.names = c("xi_", "phi1", "phi3", "V1", "V2"),
    n.iter = MCMC_SAMPLE,
    thin = thin, progress.bar = progress.bar, coda.samples_options
  )))


  out_post <- as.data.frame(posterior_sample[[1]])


  result <- list(
    "posterior_sample" = posterior_sample, # posterior samples of the link parameters and response rates generated through the MCMC process
    "mean_estimate" = round(apply(out_post, 2, mean), n.digits), # estimate of each parameter
    "ci_estimate" = bayestestR::ci(out_post, ci = ci, method = "HDI")
  ) # x% credible intervals for each parameter

  class(result) <- "BJSM_c"

  return(result)
}

#' @rdname BJSM_c
#' @param object object to summarize.
#' @param ... further arguments. Not currently used.
#' @export
#'
summary.BJSM_c <- function(object, ...) {
  out <- as.data.frame(cbind(object$mean_estimate, object$ci_estimate))
  colnames(out)[1] <- "Estimate"
  rownames(out)[11:13] <- c("xi_[A]", "xi_[B]", "xi_[C]")
  obj <- list(out = out)
  class(obj) <- "summary.BJSM_c"
  obj
}

#' @rdname BJSM_c
#' @param x object to print
#' @param ... further arguments. Not currently used.
#' @export
#' @export print.summary.BJSM_c
#'
print.summary.BJSM_c <- function(x, ...) {
  cat("\nParameter Estimation:\n")
  print(x$out[-2])
}

#' @rdname BJSM_c
#' @param x object to print
#' @param ... further arguments. Not currently used.
#' @export
#' @export print.BJSM_c
#'
print.BJSM_c <- function(x, ...) {
  cat("\nTreatment Effect Estimation:\n")
  out <- as.data.frame(cbind(x$mean_estimate, x$ci_estimate))
  colnames(out)[1] <- "Estimate"
  rownames(out)[11:13] <- c("trtA", "trtB", "trtC")
  print(out[-c(seq(1, 10, 1)), -2])
}
