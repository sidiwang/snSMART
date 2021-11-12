#' Trial dataset (dose level snSMART design)
#'
#' Generate trial dataset based on the exact number of responders/non-responders for dose level snSMART design (placebo, low, high dose; binary outcome) without missing value.
#' Useful for recording real snSMART trial result or generating specific simulation scenario. Not to be confused with \code{\link{data_simulation_dose}}.
#'
#' @param trt vector of 3 values - c(number of people who receive placebo in stage 1, number of people who receive low dose treatment in stage 1, number of people who receive high dose treatment in stage 1)
#' @param resp vector of 3 values - c(number of people who respond to placebo in stage 1, number of people who respond to low dose treatment in stage 1, number of people who respond to high dose treatment in stage 1)
#' @param trt_rep_P vector of 2 values - c(number of 1st stage responders to placebo who receive low dose treatment in stage 2, number of 1st stage responders to placebo who receive high dose treatment in stage 2)
#' @param trt_rep_L vector of 2 values - c(number of 1st stage responders to low dose treatment who receive low dose treatment in stage 2, number of 1st stage responders to low dose treatment who receive high dose treatment in stage 2)
#' @param trt_rep_H vector of 2 values - c(number of 1st stage responders to high dose treatment who receive low dose treatment in stage 2, number of 1st stage responders to high dose treatment who receive high dose treatment in stage 2)
#' @param resp_rep_P vector of 2 values - c(number of 1st stage responders to placebo who also respond to low dose treatment in stage 2, number of 1st stage responders to placebo who also respond to high dose treatment in stage 2)
#' @param resp_rep_L vector of 2 values - c(number of 1st stage responders to low dose treatment who also respond to low dose treatment  in stage 2, number of 1st stage responders to low dose treatment who also respond to high dose treatment in stage 2)
#' @param resp_rep_H vector of 2 values - c(number of 1st stage responders to high dose treatment who also respond to low dose treatment in stage 2, number of 1st stage responders to high dose treatment who also respond to high dose treatment in stage 2)
#' @param trt_nrep_P vector of 2 values - c(number of 1st stage non-responders to placebo who receive low dose treatment in stage 2, number of 1st stage non-responders to placebo who receive high dose treatment in stage 2)
#' @param trt_nrep_L vector of 2 values -c(number of 1st stage non-responders to low dose treatment who receive low dose treatment again in stage 2, number of 1st stage non-responders to low dose treatment who receive high dose treatment in stage 2)
#' @param trt_nrep_H number of 1st stage non-responders to high dose treatment who receive high dose treatment again in stage 2
#' @param resp_nrep_P vector of 2 values - c(number of 1st stage non-responders to placebo who respond to low dose treatment in stage 2, number of 1st stage non-responders to placebo who respond to high dose treatment in stage 2)
#' @param resp_nrep_L vector of 2 values - c(number of 1st stage non-responders to low dose treatment who respond to low dose treatment in stage 2, number of 1st stage non-responders to low dose treatment who respond to high dose treatment in stage 2)
#' @param resp_nrep_H number of 1st stage non-responders to high dose treatment who respond to high dose treatment in stage 2
#'
#' @return a `matrix` of the trial dataset with 4 columns: treatment_stageI, response_stageI, treatment_stageII, response_stageII

#'
#' @examples
#' mydata = trial_dataset_dose(trt = c(30, 30, 30), resp = c(5, 10, 15), trt_rep_P = c(3, 2), trt_rep_L = c(5, 5),
#'      trt_rep_H = c(8, 7), resp_rep_P = c(1, 2), resp_rep_L = c(2, 3), resp_rep_H = c(4, 6), trt_nrep_P = c(10, 15),
#'      trt_nrep_L = c(10, 10), trt_nrep_H = 15, resp_nrep_P = c(7, 8), resp_nrep_L = c(7, 6), resp_nrep_H = 10)
#'
#' @references
#' Fang, F., Hochstedler, K.A., Tamura, R.N., Braun, T.M. and Kidwell, K.M., 2021. Bayesian methods to compare dose levels with placebo in a small n,
#' sequential, multiple assignment, randomized trial. Statistics in Medicine, 40(4), pp.963-977.
#'
#' @seealso
#' \code{\link{data_simulation_dose}} \cr
#' \code{\link{BJSM_binary_dose}} \cr
#' \code{\link{JSRM_binary_dose}}
#'
#' @export
#'

trial_dataset_dose <- function(trt, resp,
                               trt_rep_P, trt_rep_L, trt_rep_H,
                               resp_rep_P, resp_rep_L, resp_rep_H,
                               trt_nrep_P, trt_nrep_L, trt_nrep_H,
                               resp_nrep_P, resp_nrep_L, resp_nrep_H){

  trtP_I = trt[1]
  trtL_I = trt[2]
  trtH_I = trt[3]

  respP_I = resp[1]
  respL_I = resp[2]
  respH_I = resp[3]

  trtPL_II = trt_rep_P[1]
  trtPH_II = trt_rep_P[2]

  trtLL_II = trt_rep_L[1]
  trtLH_II = trt_rep_L[2]

  trtHL_II = trt_rep_H[1]
  trtHH_II = trt_rep_H[2]

  respPL_II = resp_rep_P[1]
  respPH_II = resp_rep_P[2]

  respLL_II = resp_rep_L[1]
  respLH_II = resp_rep_L[2]

  respHL_II = resp_rep_H[1]
  respHH_II = resp_rep_H[2]

  trtNPL_II = trt_nrep_P[1]
  trtNPH_II = trt_nrep_P[2]

  trtNLH_II = trt_nrep_L[2]
  trtNLL_II = trt_nrep_L[1]

  trtNHH_II = trt_nrep_H

  respNPL_II = resp_nrep_P[1]
  respNPH_II = resp_nrep_P[2]

  respNLH_II = resp_nrep_L[2]
  respNLL_II = resp_nrep_L[1]

  respNHH_II = resp_nrep_H


  data_P.L.Y <- data.frame(treatment_stageI = c(rep(1, trtPL_II + trtNPL_II)),
                           response_stageI =  c(rep(1, trtPL_II), rep(0, trtNPL_II)),
                           treatment_stageII = rep(2, trtPL_II + trtNPL_II),
                           response_stageII = c(rep(1, respPL_II), rep(0, trtPL_II - respPL_II), rep(1, respNPL_II), rep(0, trtNPL_II - respNPL_II)))
  data_P.H.Y <- data.frame(treatment_stageI = c(rep(1, trtPH_II + trtNPH_II)),
                           response_stageI =  c(rep(1, trtPH_II), rep(0, trtNPH_II)),
                           treatment_stageII = rep(3, trtPH_II + trtNPH_II),
                           response_stageII = c(rep(1, respPH_II), rep(0, trtPH_II - respPH_II), rep(1, respNPH_II), rep(0, trtNPH_II - respNPH_II)))

  if (trtP_I != sum(trtPL_II, trtNPL_II, trtPH_II, trtNPH_II)){
    data_P <- data.frame(treatment_stageI = rep(1, trtP_I - trtPL_II - trtNPL_II - trtPH_II - trtNPH_II),
                         response_stageI = c(rep(1, respP_I - trtPL_II - trtPH_II),
                                             rep(0, trtP_I - respP_I - trtNPL_II - trtNPH_II)),
                         treatment_stageII = NA,
                         response_stageII = NA)
  }

  data_L.L.Y <- data.frame(treatment_stageI = c(rep(2, trtLL_II + trtNLL_II)),
                           response_stageI =  c(rep(1, trtLL_II), rep(0, trtNLL_II)),
                           treatment_stageII = rep(2, trtLL_II + trtNLL_II),
                           response_stageII = c(rep(1, respLL_II), rep(0, trtLL_II - respLL_II), rep(1, respNLL_II), rep(0, trtNLL_II - respNLL_II)))
  data_L.H.Y <- data.frame(treatment_stageI = c(rep(2, trtLH_II + trtNLH_II)),
                           response_stageI =  c(rep(1, trtLH_II), rep(0, trtNLH_II)),
                           treatment_stageII = rep(3, trtLH_II + trtNLH_II),
                           response_stageII = c(rep(1, respLH_II), rep(0, trtLH_II - respLH_II), rep(1, respNLH_II), rep(0, trtNLH_II - respNLH_II)))

  if (trtL_I != sum(trtLL_II, trtNLL_II, trtLH_II, trtNLH_II)){
    data_L <- data.frame(treatment_stageI = rep(2, trtL_I - trtLL_II - trtNLL_II - trtLH_II - trtNLH_II),
                         response_stageI = c(rep(1, respL_I - trtLL_II - trtLH_II),
                                             rep(0, trtL_I - respL_I - trtNLL_II - trtNLH_II)),
                         treatment_stageII = NA,
                         response_stageII = NA)
  }

  trtNHL_II = respNHL_II = 0
  data_H.L.Y <- data.frame(treatment_stageI = c(rep(3, trtHL_II + trtNHL_II)),
                           response_stageI =  c(rep(1, trtHL_II), rep(0, trtNHL_II)),
                           treatment_stageII = rep(2, trtHL_II + trtNHL_II),
                           response_stageII = c(rep(1, respHL_II), rep(0, trtHL_II - respHL_II), rep(1, respNHL_II), rep(0, trtNHL_II - respNHL_II)))
  data_H.H.Y <- data.frame(treatment_stageI = c(rep(3, trtHH_II + trtNHH_II)),
                           response_stageI =  c(rep(1, trtHH_II), rep(0, trtNHH_II)),
                           treatment_stageII = rep(3, trtHH_II + trtNHH_II),
                           response_stageII = c(rep(1, respHH_II), rep(0, trtHH_II - respHH_II), rep(1, respNHH_II), rep(0, trtNHH_II - respNHH_II)))

  if (trtH_I != sum(trtHL_II, trtNHL_II, trtHH_II, trtNHH_II)){
    data_H <- data.frame(treatment_stageI = rep(3, trtH_I - trtHL_II - trtNHL_II - trtHH_II - trtNHH_II),
                         response_stageI = c(rep(1, respH_I - trtHL_II - trtHH_II),
                                             rep(0, trtH_I - respH_I - trtNHL_II - trtNHH_II)),
                         treatment_stageII = NA,
                         response_stageII = NA)
  }

  data_stageI.II <- rbind(data_P.L.Y,data_P.H.Y,
                          data_L.L.Y,data_L.H.Y,
                          data_H.L.Y,data_H.H.Y)
  return(data_stageI.II)
}



#' BJSM binary (dose level snSMART design)
#'
#' BJSM (Bayesian Joint Stage Modeling) method that borrows information across both stages to estimate the individual response rate of each dose level (placebo, low, high does; binary outcome) based on trial dataset
#' generated by function \code{\link{trial_dataset_dose}} or \code{\link{data_simulation_dose}}
#'
#' @param data trial data generated through function \code{\link{trial_dataset_dose}} or \code{\link{data_simulation_dose}}
#' @param pi_prior vector of two values (a, b). `a` is the parameter `a` of the prior distribution for \code{pi} (response rate) of placebo. `b` is the parameter `b` of the prior distribution for \code{pi} of placebo. Please check the `Details` section for more explanation
#' @param beta_prior vector of two values (a, b). `a` is the parameter `a` of the prior distribution for linkage parameter \code{beta}. `b` is the parameter b of the prior distribution  for linkage parameter \code{beta}
#' @param n_MCMC_chain number of MCMC chains, default to 1
#' @param normal.par vector of two values (normal.mean, normal.var). our function assumes that the logarithm of treatment effect ratio follows a Gaussian prior distribution \eqn{N(\mu, \sigma^2)}, that is \eqn{log(\pi_L/\pi_P)~N(normal.mean, normal.var)}, and \eqn{log(\pi_H/\pi_P)~N(normal.mean, normal.var)}. \code{normal.mean} is the mean of this Gaussian prior. `normal.var` is the variance of this Gaussian prior distribution
#' @param BURN.IN number of burn-in iterations for MCMC
#' @param MCMC_SAMPLE number of iterations for MCMC
#' @param ci coverage probability for credible intervals, default = 0.95
#' @param prior_dist vector of two values ("prior distribution for \code{pi}", "prior distribution for \code{beta}"). User can choose from "gamma", "beta", "pareto". e.g. prior_dist = c("beta", "gamma")
#' @param digits the number of significant digits to use when printing
#'
#' @details
#' For gamma distribution, \code{prior.a} is the shape parameter \code{r}, \code{prior.b} is the rate parameter \code{lambda}. For beta distribution, \code{prior.a} is the shape parameter \code{a}, \code{prior.b} is the shape parameter \code{b}.
#' For pareto distribution, \code{prior.a} is the scale parameter \code{alpha}, \code{prior.b} is the shape parameter \code{c} (see page 29 of the jags user manual version 3.4.0). link: \url{http://www.stats.ox.ac.uk/~nicholls/MScMCMC14/jags_user_manual.pdf}
#'
#' The individual response rate is regarded as a permanent feature of the treatment. The second stage outcome is modeled conditionally on the first stage results linking the first and
#' second stage response probabilities through linkage parameters. The stage 1 response rate for treatment \eqn{k} is denoted as \eqn{\pi_k}. The stage 2 response rate for stage 1 responders to
#' treatment \eqn{k} who receive treatment \eqn{k'} in stage 2 is equal to \eqn{\beta_{1k}\pi_k'}. For nonresponders to treatment \eqn{k} in stage 1 who receive treatment \eqn{k*} in stage 2, the stage
#' 2 response rate is equal to \eqn{\beta_{0k}\pi_k*}
#'
#' Please refer to the paper listed under `reference` section for dose level snSMART design and detailed definition of parameters.
#'
#' @return
#' \describe{
#'   \item{posterior_sample}{posterior samples of the link parameters and response rates generated through the MCMC process}
#'   \item{pi_hat_bjsm}{estimate of response rate/treatment effect}
#'   \item{se_hat_bjsm}{standard error of the estimated response rate}
#'   \item{ci_pi_P, ci_pi_L, ci_pi_H}{x% credible intervals for treatment P, L, H}
#'   \item{diff_PL, diff_LH, diff_PH}{estimate of differences between treatment P and L, L and H, P and H}
#'   \item{se_PL, se_LH, se_PH}{standard error of the estimated differences between treatments}
#'   \item{ci_diff_PL, ci_diff_LH, ci_diff_PH}{x% credible intervals for the differences between treatment P and L, L and H, P and H}
#'   \item{beta_hat}{linkage parameter beta estimates}
#'   \item{se_beta}{standard error of the estimated value of beta}
#'   \item{ci_beta_hat}{linkage parameter beta estimates}}
#'
#'
#' @examples
#' mydata = trial_dataset_dose(trt = c(30, 30, 30), resp = c(5, 10, 15), trt_rep_P = c(3, 2), trt_rep_L = c(5, 5),
#'      trt_rep_H = c(8, 7), resp_rep_P = c(1, 2), resp_rep_L = c(2, 3), resp_rep_H = c(4, 6), trt_nrep_P = c(10, 15),
#'      trt_nrep_L = c(10, 10), trt_nrep_H = 15, resp_nrep_P = c(7, 8), resp_nrep_L = c(7, 6), resp_nrep_H = 10)
#'
#' BJSM_dose_result = BJSM_binary_dose(data = mydata, prior_dist = c("beta", "gamma"),
#'     pi_prior = c(3, 17), normal.par = c(0.2, 100), beta_prior = c(2, 2),
#'     n_MCMC_chain = 2, BURN.IN = 10000, MCMC_SAMPLE = 60000, ci = 0.95)
#'
#' summary(BJSM_dose_result)
#'
#' @references
#' Fang, F., Hochstedler, K.A., Tamura, R.N., Braun, T.M. and Kidwell, K.M., 2021. Bayesian methods to compare dose levels with placebo in a small n,
#' sequential, multiple assignment, randomized trial. Statistics in Medicine, 40(4), pp.963-977.
#'
#' @seealso
#' \code{\link{data_simulation_dose}} \cr
#' \code{\link{trial_dataset_dose}} \cr
#' \code{\link{JSRM_binary_dose}}
#'
#' @rdname BJSM_binary_dose
#' @export

BJSM_binary_dose = function(data, prior_dist, pi_prior, normal.par, beta_prior, n_MCMC_chain, BURN.IN,
                            MCMC_SAMPLE, ci = 0.95){

  pi_prior_dist = prior_dist[1]
  beta_prior_dist = prior_dist[2]

  NUM_ARMS = length(unique(data$treatment_stageI[!is.na(data$treatment_stageI)]))
  pi_prior.a = pi_prior[1]
  pi_prior.b = pi_prior[2]
  beta_prior.a = beta_prior[1]
  beta_prior.b = beta_prior[2]
  normal.mean = normal.par[1]
  normal.var = normal.par[2]
  beta_prior_dist = ifelse(beta_prior_dist == "gamma", "dgamma", ifelse(beta_prior_dist == "beta", "dbeta", "dpar"))
  pi_prior_dist = ifelse(pi_prior_dist == "gamma", "dgamma", ifelse(pi_prior_dist == "beta", "dbeta", "dpar"))

  mydata = data
  mydata$response_status_stageI = mydata$response_stageI + 1

  bugfile <- readLines("inst/BJSM_dose.bug")
  bugfile <- gsub(pattern = "pi_prior_dist", replacement = pi_prior_dist, x = bugfile)
  bugfile2 <- gsub(pattern = "beta_prior_dist", replacement = beta_prior_dist, x = bugfile)

  writeLines(bugfile2, con = "inst/BJSM_dose_new.bug")
  jag.model.name <- "BJSM_dose_new.bug"
  tryCatch({
    jag <- rjags::jags.model(file.path("inst", jag.model.name),
                             data = list(overall_sample_size = nrow(mydata),
                                       num_arms = NUM_ARMS,
                                       response_stageI = mydata$response_stageI,
                                       response_stageII = mydata$response_stageII,
                                       treatment_stageI = mydata$treatment_stageI,
                                       treatment_stageII = mydata$treatment_stageII,
                                       response_discount_status_stageI = mydata$response_status_stageI,

                                       #prior
                                       a_pi = pi_prior.a,
                                       b_pi = pi_prior.b,
                                       a_beta = beta_prior.a,
                                       b_beta = beta_prior.b,
                                       normal.mean = normal.mean,
                                       normal.var = normal.var),
                                       n.chains = n_MCMC_chain)

    posterior_sample <- rjags::coda.samples(jag,
                                            c('pi','beta'),
                                            MCMC_SAMPLE)
  },
  error = function(c) {rbind(error_round,i)
    posterior_sample_burn = window(posterior_sample,start = BURN.IN, end = MCMC_SAMPLE)
    posterior_sample_cmb = do.call(rbind, posterior_sample_burn)
    error_round = rbind(error_round,i)
    error_count = error_count+1
  },
  warning = function(c) {warn_round = rbind(warn_round,i)
  warn_count = warn_count+1},
  finally = {
    posterior_sample_burn = window(posterior_sample,start = BURN.IN, end = MCMC_SAMPLE)
    posterior_sample_cmb = do.call(rbind, posterior_sample_burn)
  }
  )



  out_post = as.data.frame(posterior_sample_cmb)
  colnames(out_post)[c(7:9)] = c("pi_P", "pi_L", "pi_H")



  result = list("posterior_sample" = out_post, # posterior samples of the link parameters and response rates generated through the MCMC process
                "pi_hat_bjsm" = apply(out_post[, 7:9],2,mean),   # estimate of response rate/treatment effect
                "se_hat_bjsm" = apply(out_post[, 7:9],2,sd),     # standard error of the response rate
                "ci_pi_P" = bayestestR::ci(out_post[, 7], ci = ci, method = "HDI"), # x% credible intervals for A
                "ci_pi_L" = bayestestR::ci(out_post[, 8], ci = ci, method = "HDI"), # x% credible intervals for B
                "ci_pi_H" = bayestestR::ci(out_post[, 9], ci = ci, method = "HDI"), # x% credible intervals for C
                "diff_PL" = mean(out_post[, 7] - out_post[, 8]), # estimate of differences between treatments A and B
                "diff_LH" = mean(out_post[, 8] - out_post[, 9]), # estimate of differences between treatments B and C
                "diff_PH" = mean(out_post[, 7] - out_post[, 9]), # estimate of differences between treatments A and C
                "se_PL" = stats::sd(out_post[, 7] - out_post[, 8]),
                "se_LH" = stats::sd(out_post[, 8] - out_post[, 9]),
                "se_PH" = stats::sd(out_post[, 7] - out_post[, 9]),
                "ci_diff_PL" = bayestestR::ci(out_post[, 7] - out_post[, 8], ci = ci, method = "HDI"), # x% credible intervals for the differences between A and B
                "ci_diff_LH" = bayestestR::ci(out_post[, 8] - out_post[, 9], ci = ci, method = "HDI"), # x% credible intervals for the differences between B and C
                "ci_diff_PH" = bayestestR::ci(out_post[, 7] - out_post[, 9], ci = ci, method = "HDI"), # x% credible intervals for the differences between A and C
                "beta_hat" = apply(out_post[, c(1:6)], 2, mean), # linkage parameter beta estimates
                "se_beta" = apply(out_post[, c(1:6)], 2, stats::sd), # linkage parameter beta estimates
                "ci_beta_hat" = HDInterval::hdi(out_post[, 1:6], ci)) # linkage parameter beta estimates

  class(result) = "BJSM_dose_binary"

  return(result)
}



#' Summarizing BJSM fits
#'
#' `summary` method for class "`BJSM_dose_binary`"
#'
#' @param object an object of class "`BJSM_dose_binary`", usually, a result of a call to \code{\link{BJSM_dose_binary}}
#' @param digits the number of significant digits to use when printing
#'
#' @returns
#' \describe{
#'    \item{Treatment Effects Estimate}{a 3 x 5 matrix with columns for the estimated treatment effects, its standard error, coverage probability of its credible interval, lower bound for its credible interval and higher bound for its credible interval}
#'    \item{Differences between Treatments}{a 3 x 5 matrix with columns for the estimated differences in treatment effects between two treatments, its standard error, coverage probability of its credible interval, lower bound and higher bound of the credible interval}
#'    \item{Linkage Parameter Estimate}{a 6 x 5 matrix with columns for the estimated linkage parameters}
#' }
#'
#'
#' @export
summary.BJSM_dose_binary = function(object, digits = 5, ...){
  cat("\nTreatment Effects Estimate:\n")
  trteff = cbind(object$pi_hat_bjsm, object$se_hat_bjsm, rbind(object$ci_pi_P, object$ci_pi_L, object$ci_pi_H))
  rownames(trteff) = c("trtP", "trtL", "trtH")
  colnames(trteff) = c("Estimate", "Std. Error", "C.I.", "CI low", "CI high")
  print(trteff, digits = digits)
  cat("\nDifferences between Treatments:\n")
  trtdiff = cbind(rbind(object$diff_PL, object$diff_LH, object$diff_PH), rbind(object$se_PL, object$se_LH, object$se_PH), rbind(object$ci_diff_PL, object$ci_diff_LH, object$ci_diff_PH))
  rownames(trtdiff) = c("diffPL", "diffLH", "diffPH")
  colnames(trtdiff) = c("Estimate", "Std.Error", "C.I.", "CI low", "CI high")
  print(trtdiff, digits = digits)
  cat("\nLinkage Parameter Estimate:\n")
  betaest = t(rbind(object$beta_hat, object$se_beta, c(rep(trteff[, 3][1], length(object$beta_hat))), object$ci_beta_hat))
  colnames(betaest) = c("Estimate", "Std. Error", "C.I.", "CI low", "CI high")
  print(betaest, digits = digits)
  cat("\n")
}


#' @rdname BJSM_binary_dose
#' @export
print.BJSM_dose_binary = function(object, digits = 5, ...){
  cat("\nTreatment Effects Estimate:\n")
  print(object$pi_hat_bjsm)
  cat("\nDifferences between Treatments:\n")
  trtdiff = rbind(object$diff_PL, object$diff_LH, object$diff_PH)
  colnames(trtdiff) = c("estimate")
  rownames(trtdiff) = c("diffPL", "diffLH", "diffPH")
  print(trtdiff)
  cat("\nLinkage Parameter Estimate:\n")
  print(object$beta_hat)
  cat("\n")
}

