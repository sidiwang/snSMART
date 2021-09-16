#' Trial dataset (dose level snSMART design)
#'
#' Generate trial dataset for dose level snSMART design (placebo, low, high dose; binary outcome) without missing value
#'
#' @param trtP_I number of people who receive placebo in stage 1
#' @param trtL_I number of people who receive low dose treatment in stage 1
#' @param trtH_I number of people who receive high dose treatment in stage 1
#' @param respP_I number of people who respond to placebo in stage 1
#' @param respL_I number of people who respond to low dose treatment in stage 1
#' @param respH_I number of people who respond to high dose treatment in stage 1
#' @param trtPL_II number of 1st stage responders to placebo who receive low dose treatment in stage 2
#' @param trtPH_II number of 1st stage responders to placebo who receive high dose treatment in stage 2
#' @param trtLL_II number of 1st stage responders to low dose treatment who receive low dose treatment in stage 2
#' @param trtLH_II number of 1st stage responders to low dose treatment who receive high dose treatment in stage 2
#' @param trtHL_II number of 1st stage responders to high dose treatment who receive low dose treatment in stage 2
#' @param trtHH_II number of 1st stage responders to high dose treatment who receive high dose treatment in stage 2
#' @param respPL_II number of 1st stage responders to placebo who also respond to low dose treatment in stage 2
#' @param respPH_II number of 1st stage responders to placebo who also respond to high dose treatment in stage 2
#' @param respLL_II number of 1st stage responders to low dose treatment who also respond to low dose treatment  in stage 2
#' @param respLH_II number of 1st stage responders to low dose treatment who also respond to high dose treatment in stage 2
#' @param respHL_II number of 1st stage responders to high dose treatment who also respond to low dose treatment in stage 2
#' @param respHH_II number of 1st stage responders to high dose treatment who also respond to high dose treatment in stage 2
#' @param trtNPL_II number of 1st stage non-responders to placebo who receive low dose treatment in stage 2
#' @param trtNPH_II number of 1st stage non-responders to placebo who receive high dose treatment in stage 2
#' @param trtNLH_II number of 1st stage non-responders to low dose treatment who receive high dose treatment in stage 2
#' @param trtNLL_II number of 1st stage non-responders to low dose treatment who receive low dose treatment again in stage 2
#' @param trtNHH_II number of 1st stage non-responders to high dose treatment who receive high dose treatment again in stage 2
#' @param respNPL_II number of 1st stage non-responders to placebo who respond to low dose treatment in stage 2
#' @param respNPH_II number of 1st stage non-responders to placebo who respond to high dose treatment in stage 2
#' @param respNLH_II number of 1st stage non-responders to low dose treatment who respond to high dose treatment in stage 2
#' @param respNLL_II number of 1st stage non-responders to low dose treatment who respond to low dose treatment in stage 2
#' @param respNHH_II number of 1st stage non-responders to high dose treatment who respond to high dose treatment in stage 2
#'
#' @return a `matrix` of the trial dataset with 4 columns: treatment_stageI, response_stageI, treatment_stageII, response_stageII

#'
#' @examples
#' mydata = trial_dataset_dose(trtP_I = 30, trtL_I = 30, trtH_I = 30, respP_I = 5,
#'      respL_I = 10, respH_I = 15,trtPL_II = 3, trtPH_II = 2, trtLL_II = 5,
#'      trtLH_II = 5, trtHL_II = 8, trtHH_II = 7, respPL_II = 1, respPH_II = 2,
#'      respLL_II = 2, respLH_II = 3, respHL_II = 4, respHH_II = 6, trtNPL_II = 10,
#'      trtNPH_II = 15, trtNLH_II = 10, trtNLL_II = 10, trtNHH_II = 15,
#'      respNPL_II = 7, respNPH_II = 8, respNLH_II = 6, respNLL_II = 7,
#'      respNHH_II = 10)
#'
#' @references
#' Fang, F., Hochstedler, K.A., Tamura, R.N., Braun, T.M. and Kidwell, K.M., 2021. Bayesian methods to compare dose levels with placebo in a small n,
#' sequential, multiple assignment, randomized trial. Statistics in Medicine, 40(4), pp.963-977.
#'
#' @export
#'

trial_dataset_dose <- function(trtP_I, trtL_I, trtH_I, respP_I, respL_I, respH_I,
                          trtPL_II, trtPH_II, trtLL_II, trtLH_II, trtHL_II, trtHH_II,
                          respPL_II, respPH_II, respLL_II, respLH_II, respHL_II, respHH_II,
                          trtNPL_II, trtNPH_II, trtNLH_II, trtNLL_II, trtNHH_II,
                          respNPL_II, respNPH_II, respNLH_II, respNLL_II, respNHH_II){


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
#' generated by function  `trial_dataset_dose`
#'
#' @param data trial data generated through function `trial_dataset_dose`
#' @param NUM_ARMS number of treatment arms
#' @param pi_prior.a parameter a of the prior distribution for \code{pi} (response rate) of placebo. Please check the `Details` section for more explanation
#' @param pi_prior.b parameter b of the prior distribution for \code{pi} (response rate) of placebo. Please check the `Details` section for more explanation
#' @param beta_prior.a parameter a of the prior distribution for linkage parameter \code{beta}
#' @param beta_prior.b parameter b of the prior distribution  for linkage parameter \code{beta}
#' @param n_MCMC_chain number of MCMC chains, default to 1
#' @param normal.mean our function assumes that the logarithm of treatment effect ratio follows a Gaussian prior distribution \eqn{N(\mu, \sigma^2)}, that is \eqn{log(\pi_L/\pi_P)~N(normal.mean, normal.var)}, and \eqn{log(\pi_H/\pi_P)~N(normal.mean, normal.var)}. \code{normal.mean} is the mean of this Gaussian prior
#' @param normal.var similar to above, `normal.var` is the variance of this Gaussian prior distribution
#' @param BURN.IN number of burn-in iterations for MCMC
#' @param MCMC_SAMPLE number of iterations for MCMC
#' @param ci coverage probability for credible intervals, default = 0.95
#' @param pi_prior_dist prior distribution for pi (response rate) of placebo, user can choose from "gamma", "beta", "pareto"
#' @param beta_prior_dist prior distribution for beta, user can choose from gamma, beta, pareto
#'
#' @details
#' For gamma distribution, \code{prior.a} is the shape parameter \code{r}, \code{prior.b} is the rate parameter \code{lambda}. For beta distribution, \code{prior.a} is the shape parameter \code{a}, \code{prior.b} is the shape parameter \code{b}.
#' For pareto distribution, \code{prior.a} is the scale parameter \code{alpha}, \code{prior.b} is the shape parameter \code{c} (see page 29 of the jags user manual version 3.4.0). link: \url{http://www.stats.ox.ac.uk/~nicholls/MScMCMC14/jags_user_manual.pdf}
#'
#' The individual response rate is regarded as a permanent feature of the treatment. The second stage outcome is modeled conditionally on the first stage results linking the first and
#' second stage response probabilities through linkage parameters.
#'
#' @return
#' \strong{`posterior_sample`}: posterior samples of the link parameters and response rates generated through the MCMC process
#' \strong{`pi_hat_bjsm`}: estimate of response rate/treatment effect
#' \strong{`se_hat_bjsm`}: standard error of the response rate
#' \strong{`ci_pi_P`}: x% credible intervals for A
#' \strong{`ci_pi_L`}: x% credible intervals for B
#' \strong{`ci_pi_H`}: x% credible intervals for C
#' \strong{`diff_PL`}: estimate of differences between treatments A and B
#' \strong{`diff_LH`}: estimate of differences between treatments B and C
#' \strong{`diff_PH`}: estimate of differences between treatments A and C
#' \strong{`ci_diff_PL`}: x% credible intervals for the differences between A and B
#' \strong{`ci_diff_LH`}: x% credible intervals for the differences between B and C
#' \strong{`ci_diff_PH`}: x% credible intervals for the differences between A and C
#' \strong{`beta_hat`}: linkage parameter beta estimates
#' \strong{`ci_beta_hat`}: linkage parameter beta estimates
#'
#' @examples
#' mydata = trial_dataset_dose(trtP_I = 30, trtL_I = 30, trtH_I = 30, respP_I = 5,
#'     respL_I = 10, respH_I = 15, trtPL_II = 3, trtPH_II = 2, trtLL_II = 5,
#'     trtLH_II = 5, trtHL_II = 8, trtHH_II = 7, respPL_II = 1, respPH_II = 2,
#'     respLL_II = 2, respLH_II = 3, respHL_II = 4, respHH_II = 6,
#'     trtNPL_II = 10, trtNPH_II = 15, trtNLH_II = 10, trtNLL_II = 10,
#'     trtNHH_II = 15, respNPL_II = 7, respNPH_II = 8, respNLH_II = 6,
#'     respNLL_II = 7, respNHH_II = 10)
#'
#' BJSM_dose_result = BJSM_binary_dose(data = mydata, NUM_ARMS = 3, pi_prior_dist = "beta",
#'     pi_prior.a = 3, pi_prior.b = 17, beta_prior_dist = "gamma",
#'     normal.mean = 0.2, normal.var = 100, beta_prior.a = 2, beta_prior.b = 2,
#'     n_MCMC_chain = 2, BURN.IN = 10000, MCMC_SAMPLE = 60000, ci = 0.95)
#'
#' @references
#' Fang, F., Hochstedler, K.A., Tamura, R.N., Braun, T.M. and Kidwell, K.M., 2021. Bayesian methods to compare dose levels with placebo in a small n,
#' sequential, multiple assignment, randomized trial. Statistics in Medicine, 40(4), pp.963-977.
#'
#' @export

BJSM_binary_dose = function(data, NUM_ARMS, pi_prior_dist, pi_prior.a, pi_prior.b, normal.mean, normal.var, beta_prior_dist, beta_prior.a, beta_prior.b, n_MCMC_chain, BURN.IN,
                            MCMC_SAMPLE, ci = 0.95){

  beta_prior_dist = ifelse(beta_prior_dist == "gamma", "dgamma", ifelse(beta_prior_dist == "beta", "dbeta", "dpar"))
  pi_prior_dist = ifelse(pi_prior_dist == "gamma", "dgamma", ifelse(pi_prior_dist == "beta", "dbeta", "dpar"))

  mydata = data
  mydata$response_status_stageI = mydata$response_stageI + 1

  bugfile  <- readLines("inst/BJSM_dose.bug")
  bugfile  <- gsub(pattern = "pi_prior_dist", replacement = pi_prior_dist, x = bugfile)
  bugfile2  <- gsub(pattern = "beta_prior_dist", replacement = beta_prior_dist, x = bugfile)

  writeLines(bugfile2, con="inst/BJSM_dose_new.bug")
  jag.model.name <- "BJSM_dose_new.bug"
  tryCatch({
    jag <- rjags::jags.model(file.path("inst",jag.model.name),
                             data=list(overall_sample_size = nrow(mydata),
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
    posterior_sample_burn=window(posterior_sample,start=BURN.IN, end=MCMC_SAMPLE)
    posterior_sample_cmb=do.call(rbind, posterior_sample_burn)
    error_round=rbind(error_round,i)
    error_count=error_count+1
  },
  warning = function(c) {warn_round=rbind(warn_round,i)
  warn_count=warn_count+1},
  finally = {
    posterior_sample_burn=window(posterior_sample,start=BURN.IN, end=MCMC_SAMPLE)
    posterior_sample_cmb=do.call(rbind, posterior_sample_burn)
  }
  )



  out_post = posterior_sample_cmb



  result = list("posterior_sample" = out_post, # posterior samples of the link parameters and response rates generated through the MCMC process
                "pi_hat_bjsm" = apply(out_post[,7:9],2,mean),   # estimate of response rate/treatment effect
                "se_hat_bjsm" = apply(out_post[,7:9],2,sd),     # standard error of the response rate
                "ci_pi_P" = bayestestR::ci(out_post[,7], ci = ci, method = "HDI"), # x% credible intervals for A
                "ci_pi_L" = bayestestR::ci(out_post[,8], ci = ci, method = "HDI"), # x% credible intervals for B
                "ci_pi_H" = bayestestR::ci(out_post[,9], ci = ci, method = "HDI"), # x% credible intervals for C
                "diff_PL" = mean(out_post[,7] - out_post[,8]), # estimate of differences between treatments A and B
                "diff_LH" = mean(out_post[,8] - out_post[,9]), # estimate of differences between treatments B and C
                "diff_PH" = mean(out_post[,7] - out_post[,9]), # estimate of differences between treatments A and C
                "ci_diff_PL" = bayestestR::ci(out_post[,7] - out_post[,8], ci = ci, method = "HDI"), # x% credible intervals for the differences between A and B
                "ci_diff_LH" = bayestestR::ci(out_post[,8] - out_post[,9], ci = ci, method = "HDI"), # x% credible intervals for the differences between B and C
                "ci_diff_PH" = bayestestR::ci(out_post[,7] - out_post[,9], ci = ci, method = "HDI"), # x% credible intervals for the differences between A and C
                "beta_hat" = apply(out_post[,c(1:6)],2,mean), # linkage parameter beta estimates
                "ci_beta_hat" = HDInterval::hdi(out_post[,1:6], ci)) # linkage parameter beta estimates



  return(result)
}
