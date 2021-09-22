#' Trial dataset
#'
#' Generate trial dataset for standard snSMART (3 active treatments, non-responders re-randomized; binary outcome) with or without missing value
#'
#' @param trtA_I number of people who receive A in stage 1
#' @param trtB_I number of people who receive B in stage 1
#' @param trtC_I number of people who receive C in stage 1
#' @param respA_I number of people who respond to A in stage 1, who receive same trt in stage 2
#' @param respB_I number of people who respond to B in stage 1, who receive same trt in stage 2
#' @param respC_I number of people who respond to C in stage 1, who receive same trt in stage 2
#' @param trtAA_II number of 1st stage responders to A who receive A in stage 2
#' @param trtBB_II number of 1st stage responders to B who receive B in stage 2
#' @param trtCC_II number of 1st stage responders to C who receive C in stage 2
#' @param respA_II number of 1st stage responders who respond to A again in stage 2
#' @param respB_II number of 1st stage responders who respond to B again in stage 2
#' @param respC_II number of 1st stage responders who respond to C again in stage 2
#' @param trtAB_II number of 1st stage non-responders to A who receive B in stage 2
#' @param trtAC_II number of 1st stage non-responders to A who receive C in stage 2
#' @param trtBA_II number of 1st stage non-responders to B who receive A in stage 2
#' @param trtBC_II number of 1st stage non-responders to B who receive C in stage 2
#' @param trtCA_II number of 1st stage non-responders to C who receive A in stage 2
#' @param trtCB_II number of 1st stage non-responders to C who receive B in stage 2
#' @param respAB_II number of 1st stage non-responders to A who respond to B in stage 2
#' @param respAC_II number of 1st stage non-responders to A who respond to C in stage 2
#' @param respBA_II number of 1st stage non-responders to B who respond to A in stage 2
#' @param respBC_II number of 1st stage non-responders to B who respond to C in stage 2
#' @param respCA_II number of 1st stage non-responders to C who respond to A in stage 2
#' @param respCB_II number of 1st stage non-responders to C who respond to B in stage 2
#'
#' @return a `matrix` of the trial dataset with 4 columns: treatment_stageI, response_stageI, treatment_stageII, response_stageII

#'
#' @examples
#' trial_data = trial_dataset(trtA_I = 9, trtB_I = 12, trtC_I = 9, respA_I = 3,
#'     respB_I = 3, respC_I = 5, trtAA_II = 2, trtBB_II = 2, trtCC_II = 4,
#'     respA_II = 2, respB_II = 2, respC_II = 4, trtAB_II = 3, trtAC_II = 2,
#'     trtBA_II = 3, trtBC_II = 2, trtCA_II = 1, trtCB_II = 2, respAB_II = 3,
#'     respAC_II = 1, respBA_II = 2, respBC_II = 2, respCA_II = 0, respCB_II = 0)
#'
#' @references
#' Wei, B., Braun, T.M., Tamura, R.N. and Kidwell, K.M., 2018. A Bayesian analysis of small n sequential multiple assignment randomized trials (snSMARTs).
#' Statistics in medicine, 37(26), pp.3723-3732.
#'
#' @export
#'

trial_dataset <- function(trtA_I, trtB_I, trtC_I, respA_I, respB_I, respC_I,
                          trtAA_II, trtBB_II, trtCC_II, respA_II, respB_II, respC_II,
                          trtAB_II, trtAC_II, trtBA_II, trtBC_II, trtCA_II, trtCB_II,
                          respAB_II, respAC_II, respBA_II, respBC_II, respCA_II, respCB_II){

  data_A.A.Y <- data.frame(treatment_stageI = rep(1, trtAA_II),
                           response_stageI = rep(1, trtAA_II),
                           treatment_stageII = rep(1, trtAA_II),
                           response_stageII = c(rep(1, respA_II),rep(0, trtAA_II - respA_II)))
  data_A.B.Y <- data.frame(treatment_stageI = rep(1, trtAB_II),
                           response_stageI = rep(0, trtAB_II),
                           treatment_stageII = rep(2, trtAB_II),
                           response_stageII = c(rep(1, respAB_II),rep(0, trtAB_II - respAB_II)))
  data_A.C.Y <- data.frame(treatment_stageI = rep(1, trtAC_II),
                           response_stageI = rep(0, trtAC_II),
                           treatment_stageII = rep(3, trtAC_II),
                           response_stageII = c(rep(1, respAC_II),rep(0, trtAC_II - respAC_II)))

  if (trtA_I != sum(trtAA_II, trtAB_II, trtAC_II)){
    data_A <- data.frame(treatment_stageI = rep(1, trtA_I - trtAA_II - trtAB_II - trtAC_II),
                         response_stageI = c(rep(1, respA_I - trtAA_II),
                                             rep(0, trtA_I - respA_I - trtAB_II - trtAC_II)),
                         treatment_stageII = NA,
                         response_stageII = NA)
  }

  data_B.B.Y <- data.frame(treatment_stageI = rep(2, trtBB_II),
                           response_stageI = rep(1, trtBB_II),
                           treatment_stageII = rep(2, trtBB_II),
                           response_stageII = c(rep(1, trtBB_II),rep(0, trtBB_II - respB_II)))
  data_B.A.Y <- data.frame(treatment_stageI = rep(2, trtBA_II),
                           response_stageI = rep(0, trtBA_II),
                           treatment_stageII = rep(1, trtBA_II),
                           response_stageII = c(rep(1, respBA_II),rep(0, trtBA_II - respBA_II)))
  data_B.C.Y <- data.frame(treatment_stageI = rep(2, trtBC_II),
                           response_stageI = rep(0, trtBC_II),
                           treatment_stageII = rep(3, trtBC_II),
                           response_stageII = c(rep(1, respBC_II),rep(0, trtBC_II - respBC_II)))

  if (trtB_I != sum(trtBB_II, trtBA_II, trtBC_II)){
    data_B <- data.frame(treatment_stageI = rep(2, trtB_I - trtBB_II - trtBA_II - trtBC_II),
                         response_stageI = c(rep(1, respB_I - trtBB_II),
                                             rep(0, trtB_I - respB_I - trtBA_II - trtBC_II)),
                         treatment_stageII = NA,
                         response_stageII = NA)
  }

  data_C.C.Y <- data.frame(treatment_stageI = rep(3, trtCC_II),
                           response_stageI = rep(1, trtCC_II),
                           treatment_stageII = rep(3, trtCC_II),
                           response_stageII = c(rep(1, trtCC_II),rep(0, trtCC_II - respC_II)))
  data_C.A.Y <- data.frame(treatment_stageI = rep(3, trtCA_II),
                           response_stageI = rep(0, trtCA_II),
                           treatment_stageII = rep(1, trtCA_II),
                           response_stageII = c(rep(1, respCA_II),rep(0, trtCA_II - respCA_II)))
  data_C.B.Y <- data.frame(treatment_stageI = rep(3, trtCB_II),
                           response_stageI = rep(0, trtCB_II),
                           treatment_stageII = rep(2, trtCB_II),
                           response_stageII = c(rep(1, respCB_II),rep(0, trtCB_II - respCB_II)))

  if (trtC_I != sum(trtCC_II, trtCA_II, trtCB_II)){
    data_C <- data.frame(treatment_stageI = rep(3, trtC_I - trtCC_II - trtCA_II - trtCB_II),
                         response_stageI = c(rep(1, respC_I - trtCC_II),
                                             rep(0, trtC_I - respC_I - trtCA_II - trtCB_II)),
                         treatment_stageII = NA,
                         response_stageII = NA)
  }

  data_stageI.II <- rbind(data_A.A.Y,data_A.B.Y,data_A.C.Y,data_A,
                          data_B.B.Y,data_B.A.Y,data_B.C.Y,data_B,
                          data_C.C.Y,data_C.A.Y,data_C.B.Y,data_C)
  return(data_stageI.II)
}



#' BJSM binary
#'
#' BJSM (Bayesian Joint Stage Modeling) method that borrows information across both stages to estimate the individual response rate of each treatment (different active treatments, non-responders re-randomized; binary outcome) based on trial dataset
#' generated by function  `trial_data`
#'
#' @param data trial data generated through function `trial_dataset`
#' @param NUM_ARMS number of treatment arms
#' @param pi_prior.a parameter \code{a} of the prior distribution for \code{pi_1K}, a vector with three values, one for each treatment. Please check the `Details` section for more explaination
#' @param pi_prior.b parameter \code{b} of the prior distribution  for \code{pi_1K}, a vector with three values, one for each treatment. Please check the `Details` section for more explaination
#' @param beta0_prior.a parameter \code{a} of the prior distribution for linkage parameter \code{beta0}. Please check the `Details` section for more explaination
#' @param beta0_prior.b parameter \code{b} of the prior distribution  for linkage parameter \code{beta0}. Please check the `Details` section for more explaination
#' @param beta1_prior.a parameter \code{a} of the prior distribution for linkage parameter \code{beta1}. Please check the `Details` section for more explaination
#' @param beta1_prior.c parameter \code{b} of the prior distribution for linkage parameter \code{beta1}. Please check the `Details` section for more explaination
#' @param n_MCMC_chain number of MCMC chains, default to 1.
#' @param BURN.IN number of burn-in iterations for MCMC
#' @param MCMC_SAMPLE number of iterations for MCMC
#' @param ci coverage probability for credible intervals, default = 0.95
#' @param pi_prior_dist prior distribution for \code{pi}, user can choose from "gamma", "beta", "pareto"
#' @param beta0_prior_dist prior distribution for \code{beta0}, user can choose from "gamma", "beta", "pareto"
#' @param beta1_prior_dist prior distribution for \code{beta1}, user can choose from "gamma", "beta", "pareto"
#' @param six TRUE or FALSE. If TRUE, will run the six beta model, if FALSE will run the two beta model. default = TRUE
#' @param DTR TRUE or FALSE. If TRUE, will also return the expected response rate of dynamic treatment regimens. default = TRUE
#'
#' @details
#' For \code{gamma} distribution, \code{prior.a} is the shape parameter \code{r}, \code{prior.b} is the rate parameter \code{lambda}. For \code{beta} distribution, \code{prior.a} is the shape parameter \code{a}, \code{prior.b} is the shape parameter \code{b}.
#' For \code{pareto} distribution, \code{prior.a} is the scale parameter \code{alpha}, \code{prior.b} is the shape parameter \code{c} (see page 29 of the jags user manual version 3.4.0). link: \url{http://www.stats.ox.ac.uk/~nicholls/MScMCMC14/jags_user_manual.pdf}
#'
#' The individual response rate is regarded as a permanent feature of the treatment. The second stage outcome is modeled conditionally on the first stage results linking the first and
#' second stage response probabilities through linkage parameters.
#'
#' @return
#' @param posterior_sample posterior samples of the link parameters and response rates generated through the MCMC process
#' \itemize{
#'
#'
#' \strong{`pi_hat_bjsm`}: estimate of response rate/treatment effect \cr
#'
#' \strong{`se_hat_bjsm`}: standard error of the response rate \cr
#'
#' \strong{`ci_pi_A`}: x% credible intervals for A \cr
#'
#' \strong{`ci_pi_B`}: x% credible intervals for B \cr
#'
#' \strong{`ci_pi_C`}: x% credible intervals for C \cr
#'
#' \strong{`diff_AB`}: estimate of differences between treatments A and B \cr
#'
#' \strong{`diff_BC`}: estimate of differences between treatments B and C \cr
#'
#' \strong{`diff_AC`}: estimate of differences between treatments A and C \cr
#'
#' \strong{`ci_diff_AB`}: x% credible intervals for the differences between A and B \cr
#'
#' \strong{`ci_diff_BC`}: x% credible intervals for the differences between B and C \cr
#'
#' \strong{`ci_diff_AC`}: x% credible intervals for the differences between A and C \cr
#'
#' \strong{`beta0_hat`}: linkage parameter \code{beta0} estimates \cr
#'
#' \strong{`beta1_hat`}: linkage parameter \code{beta1} estimates \cr
#'
#' \strong{`ci_beta0_hat`}: linkage parameter \code{beta0} credible interval \cr
#'
#' \strong{`ci_beta1_hat`}: linkage parameter \code{beta1} credible interval \cr
#'
#' \strong{`pi_DTR_est`}: expected response rate of dynamic treatment regimens (DTRs) \cr
#' }
#'
#' @examples
#' mydata = trial_dataset(trtA_I = 9, trtB_I = 12, trtC_I = 9, respA_I = 3,
#'     respB_I = 3, respC_I = 5, trtAA_II = 2, trtBB_II = 2, trtCC_II = 4,
#'     respA_II = 2, respB_II = 2, respC_II = 4, trtAB_II = 3, trtAC_II = 2,
#'     trtBA_II = 3, trtBC_II = 2, trtCA_II = 1, trtCB_II = 2, respAB_II = 3,
#'     respAC_II = 1, respBA_II = 2, respBC_II = 2, respCA_II = 0, respCB_II = 0)
#'
#' BJSM_result = BJSM_binary(data = mydata, NUM_ARMS = 3, pi_prior_dist = "beta",
#'     pi_prior.a = c(0.4, 0.4, 0.4), pi_prior.b = c(1.6, 1.6, 1.6),
#'     beta0_prior_dist = "beta", beta0_prior.a = 1.6, beta0_prior.b = 0.4,
#'     beta1_prior_dist = "pareto", beta1_prior.a = 3, beta1_prior.c = 1,
#'     n_MCMC_chain = 1, BURN.IN = 10000, MCMC_SAMPLE = 60000, ci = 0.95,
#'     six = TRUE, DTR = TRUE)
#'
#' @references
#' Wei, B., Braun, T.M., Tamura, R.N. and Kidwell, K.M., 2018. A Bayesian analysis of small n sequential multiple assignment randomized trials (snSMARTs).
#' Statistics in medicine, 37(26), pp.3723-3732.
#'
#' Chao, Y.C., Trachtman, H., Gipson, D.S., Spino, C., Braun, T.M. and Kidwell, K.M., 2020. Dynamic treatment regimens in small n, sequential, multiple assignment,
#' randomized trials: An application in focal segmental glomerulosclerosis. Contemporary clinical trials, 92, p.105989.
#'
#' @export

BJSM_binary = function(data, NUM_ARMS, pi_prior_dist, pi_prior.a, pi_prior.b, beta0_prior_dist, beta0_prior.a, beta0_prior.b, beta1_prior_dist, beta1_prior.a, beta1_prior.c, n_MCMC_chain, BURN.IN,
                       MCMC_SAMPLE, ci = 0.95, six = TRUE, DTR = TRUE){
  # NUM_ARMS number of treatment arms
  # pi_prior.a  alpha parameter of the prior beta distribution for pi_1K, a vector with three values, one for each treatment
  # pi_prior.b  beta parameter of the prior beta distribution  for pi_1K, a vector with three values, one for each treatment
  # beta0_prior.a alpha parameter of the prior beta distribution for linkage parameter beta0
  # beta0_prior.b beta parameter of the prior beta distribution  for linkage parameter beta0
  # beta1_prior.a scale parameter of the prior pareto distribution for linkage parameter beta1
  # beta1_prior.c shape parameter of the prior pareto distribution for linkage parameter beta1
  # n_MCMC_chain number of MCMC chains, default to 1. If this is set to a number more than 1
  # BURN.IN number of burn-in iterations for MCMC
  # MCMC_SAMPLE number of iterations for MCMC
  # ci coverage probability for credible intervals
  # pi_prior_dist, beta0_prior_dist, beta1_prior_dist are prior distributions for pi, beta0, and beta1, user can choose from gamma, beta, pareto
  # for gamma, prior.a is r, prior.b is lambda, for beta, prior.a is alpha, prior.b is beta, for pareto, prior.a is alpha, prior.b is c (page 29 of the jags user manual version 3.4.0)
  # six, if TRUE, will run the six beta model, if FALSE will run the two beta model

  pi_prior_dist = ifelse(pi_prior_dist == "gamma", "dgamma", ifelse(pi_prior_dist == "beta", "dbeta", "dpar"))
  beta0_prior_dist = ifelse(beta0_prior_dist == "gamma", "dgamma", ifelse(beta0_prior_dist == "beta", "dbeta", "dpar"))
  beta1_prior_dist = ifelse(beta1_prior_dist == "gamma", "dgamma", ifelse(beta1_prior_dist == "beta", "dbeta", "dpar"))

  mydata = data
  mydata$disc <- 2 * mydata$treatment_stageI - (mydata$response_stageI == 0)

  if (six == TRUE){
    # If using 6-betas model
    bugfile  <- readLines("inst/BJSM_6betas_missing.bug")
    bugfile  <- gsub(pattern = "pi_prior_dist", replacement = pi_prior_dist, x = bugfile)
    bugfile  <- gsub(pattern = "beta0_prior_dist", replacement = beta0_prior_dist, x = bugfile)
    bugfile2  <- gsub(pattern = "beta1_prior_dist", replacement = beta1_prior_dist, x = bugfile)

    writeLines(bugfile2, con="inst/BJSM_6betas_missing_new.bug")
    jag.model.name <- "BJSM_6betas_missing_new.bug"
    tryCatch({
      jag <- rjags::jags.model(file.path("inst",jag.model.name),
                        data=list(n1 = nrow(mydata),
                                  n2 = nrow(mydata[!is.na(mydata$response_stageII),]),
                                  num_arms = NUM_ARMS,
                                  Y1 = mydata$response_stageI,
                                  Y2 = mydata$response_stageII[!is.na(mydata$response_stageII)],
                                  treatment_stageI = mydata$treatment_stageI,
                                  treatment_stageII = mydata$treatment_stageII[!is.na(mydata$response_stageII)],
                                  response_stageI_disc = mydata$disc[!is.na(mydata$response_stageII)],
                                  #prior
                                  pi_prior.a = pi_prior.a,
                                  pi_prior.b = pi_prior.b,
                                  beta0_prior.a = beta0_prior.a,
                                  beta0_prior.b = beta0_prior.b,
                                  beta1_prior.a = beta1_prior.a,    # pareto
                                  beta1_prior.c = beta1_prior.c     # pareto
                        ),
                        n.chains=n_MCMC_chain, n.adapt = BURN.IN)
      posterior_sample <- rjags::coda.samples(jag,
                                       c('pi','beta'),
                                       MCMC_SAMPLE)
    },
    warning = function(war){},
    error = function(err){},
    finally = {}
    )

  } else {
    bugfile  <- readLines("inst/BJSM_2beta_missing.bug")
    bugfile  <- gsub(pattern = "pi_prior_dist", replacement = pi_prior_dist, x = bugfile)
    bugfile  <- gsub(pattern = "beta0_prior_dist", replacement = beta0_prior_dist, x = bugfile)
    bugfile2  <- gsub(pattern = "beta1_prior_dist", replacement = beta1_prior_dist, x = bugfile)

    writeLines(bugfile2, con="inst/BJSM_2beta_missing_new.bug")

    # If using 2-betas model
    jag.model.name <- "BJSM_2beta_missing_new.bug"  # beta1 ~ pareto
    tryCatch({
      jag <- rjags::jags.model(file.path("inst",jag.model.name),
                        data=list(n1 = nrow(mydata),
                                  n2 = nrow(mydata[!is.na(mydata$response_stageII),]),
                                  num_arms = NUM_ARMS,
                                  Y1 = mydata$response_stageI,
                                  Y2 = mydata$response_stageII[!is.na(mydata$response_stageII)],
                                  treatment_stageI = mydata$treatment_stageI,
                                  treatment_stageII = mydata$treatment_stageII[!is.na(mydata$response_stageII)],
                                  response1 = mydata$response_stageI[!is.na(mydata$response_stageII)] + 1,
                                  #prior
                                  pi_prior.a = pi_prior.a,
                                  pi_prior.b = pi_prior.b,
                                  beta0_prior.a = beta0_prior.a,
                                  beta0_prior.b = beta0_prior.b,
                                  beta1_prior.a = beta1_prior.a,
                                  beta1_prior.c = beta1_prior.c),
                        n.chains=n_MCMC_chain,n.adapt = BURN.IN)
      posterior_sample <- rjags::coda.samples(jag,
                                       c('pi','beta'),
                                       MCMC_SAMPLE)
    },
    warning = function(war){},
    error = function(err){},
    finally = {}
    )

  }

  out_post = posterior_sample[[1]]


  if (six == TRUE){

    pi_DTR_est = c()
    pi_AB_tt <- out_post[,7]^2*out_post[,2]+(1-out_post[,7])*out_post[,8]*out_post[,1]
    pi_AC_tt <- out_post[,7]^2*out_post[,2]+(1-out_post[,7])*out_post[,9]*out_post[,1]
    pi_BA_tt <- out_post[,8]^2*out_post[,4]+(1-out_post[,8])*out_post[,7]*out_post[,3]
    pi_BC_tt <- out_post[,8]^2*out_post[,4]+(1-out_post[,8])*out_post[,9]*out_post[,3]
    pi_CA_tt <- out_post[,9]^2*out_post[,6]+(1-out_post[,9])*out_post[,7]*out_post[,5]
    pi_CB_tt <- out_post[,9]^2*out_post[,6]+(1-out_post[,9])*out_post[,8]*out_post[,5]
    pi_DTR_est <- rbind(pi_DTR_est,c(mean(pi_AB_tt),mean(pi_AC_tt),mean(pi_BA_tt),mean(pi_BC_tt),mean(pi_CA_tt),mean(pi_CB_tt)))

    if (DTR == TRUE){

      result = list("posterior_sample" = out_post, # posterior samples of the link parameters and response rates generated through the MCMC process
                    "pi_hat_bjsm" = apply(out_post[,7:9],2,mean),   # estimate of response rate/treatment effect
                    "se_hat_bjsm" = apply(out_post[,7:9],2,stats::sd),     # standard error of the response rate
                    "ci_pi_A" = bayestestR::ci(out_post[,7], ci = ci, method = "HDI"), # x% credible intervals for A
                    "ci_pi_B" = bayestestR::ci(out_post[,8], ci = ci, method = "HDI"), # x% credible intervals for B
                    "ci_pi_C" = bayestestR::ci(out_post[,9], ci = ci, method = "HDI"), # x% credible intervals for C
                    "diff_AB" = mean(out_post[,7] - out_post[,8]), # estimate of differences between treatments A and B
                    "diff_BC" = mean(out_post[,8] - out_post[,9]), # estimate of differences between treatments B and C
                    "diff_AC" = mean(out_post[,7] - out_post[,9]), # estimate of differences between treatments A and C
                    "ci_diff_AB" = bayestestR::ci(out_post[,7] - out_post[,8], ci = ci, method = "HDI"), # x% credible intervals for the differences between A and B
                    "ci_diff_BC" = bayestestR::ci(out_post[,8] - out_post[,9], ci = ci, method = "HDI"), # x% credible intervals for the differences between B and C
                    "ci_diff_AC" = bayestestR::ci(out_post[,7] - out_post[,9], ci = ci, method = "HDI"), # x% credible intervals for the differences between A and C
                    "beta0_hat" = apply(out_post[,c(1,3,5)],2,mean), # linkage parameter beta0 estimates
                    "beta1_hat" = apply(out_post[,c(2,4,6)],2,mean),  # linkage parameter beta1 estimates
                    "ci_beta0_hat" = HDInterval::hdi(out_post[,c(1,3,5)], ci), # linkage parameter beta0 credible interval
                    "ci_beta1_hat" = HDInterval::hdi(out_post[,c(2,4,6)], ci), # linkage parameter beta1 credible interval
                    "pi_DTR_est" = pi_DTR_est) # expected response rate of dynamic treatment regimens (DTRs)
    }else{
      result = list("posterior_sample" = out_post, # posterior samples of the link parameters and response rates generated through the MCMC process
                    "pi_hat_bjsm" = apply(out_post[,7:9],2,mean),   # estimate of response rate/treatment effect
                    "se_hat_bjsm" = apply(out_post[,7:9],2,stats::sd),     # standard error of the response rate
                    "ci_pi_A" = bayestestR::ci(out_post[,7], ci = ci, method = "HDI"), # x% credible intervals for A
                    "ci_pi_B" = bayestestR::ci(out_post[,8], ci = ci, method = "HDI"), # x% credible intervals for B
                    "ci_pi_C" = bayestestR::ci(out_post[,9], ci = ci, method = "HDI"), # x% credible intervals for C
                    "diff_AB" = mean(out_post[,7] - out_post[,8]), # estimate of differences between treatments A and B
                    "diff_BC" = mean(out_post[,8] - out_post[,9]), # estimate of differences between treatments B and C
                    "diff_AC" = mean(out_post[,7] - out_post[,9]), # estimate of differences between treatments A and C
                    "ci_diff_AB" = bayestestR::ci(out_post[,7] - out_post[,8], ci = ci, method = "HDI"), # x% credible intervals for the differences between A and B
                    "ci_diff_BC" = bayestestR::ci(out_post[,8] - out_post[,9], ci = ci, method = "HDI"), # x% credible intervals for the differences between B and C
                    "ci_diff_AC" = bayestestR::ci(out_post[,7] - out_post[,9], ci = ci, method = "HDI"), # x% credible intervals for the differences between A and C
                    "beta0_hat" = apply(out_post[,c(1,3,5)],2,mean), # linkage parameter beta0 estimates
                    "beta1_hat" = apply(out_post[,c(2,4,6)],2,mean),  # linkage parameter beta1 estimates
                    "ci_beta0_hat" = HDInterval::hdi(out_post[,c(1,3,5)], ci), # linkage parameter beta0 credible interval
                    "ci_beta1_hat" = HDInterval::hdi(out_post[,c(2,4,6)], ci)) # linkage parameter beta1 credible interval
                    #    "pi_DTR_est" = pi_DTR_est) # expected response rate of dynamic treatment regimens (DTRs)
    }
  } else {

    pi_DTR_est = c()
    pi_AB_tt <- out_post[,3]^2*out_post[,2]+(1-out_post[,3])*out_post[,4]*out_post[,1]
    pi_AC_tt <- out_post[,3]^2*out_post[,2]+(1-out_post[,3])*out_post[,5]*out_post[,1]
    pi_BA_tt <- out_post[,4]^2*out_post[,2]+(1-out_post[,4])*out_post[,3]*out_post[,1]
    pi_BC_tt <- out_post[,4]^2*out_post[,2]+(1-out_post[,4])*out_post[,5]*out_post[,1]
    pi_CA_tt <- out_post[,5]^2*out_post[,2]+(1-out_post[,5])*out_post[,3]*out_post[,1]
    pi_CB_tt <- out_post[,5]^2*out_post[,2]+(1-out_post[,5])*out_post[,4]*out_post[,1]
    pi_DTR_est <- rbind(pi_DTR_est,c(mean(pi_AB_tt),mean(pi_AC_tt),mean(pi_BA_tt),mean(pi_BC_tt),mean(pi_CA_tt),mean(pi_CB_tt)))

    if (DTR == TRUE){
      result = list("posterior_sample" = out_post, # posterior samples of the link parameters and response rates generated through the MCMC process
                    "pi_hat_bjsm" = apply(out_post[,3:5],2,mean), # estimate of treatment response rate
                    "se_hat_bjsm" = apply(out_post[,3:5],2,stats::sd), # estimated standard error of the response rate
                    "ci_pi_A" = bayestestR::ci(out_post[,3], ci = ci, method = "HDI"), # x% credible intervals for the response rate of treatment A
                    "ci_pi_B" = bayestestR::ci(out_post[,4], ci = ci, method = "HDI"), # x% credible intervals for the response rate of treatment B
                    "ci_pi_C" = bayestestR::ci(out_post[,5], ci = ci, method = "HDI"), # x% credible intervals for the response rate of treatment C
                    "diff_AB" = mean(out_post[,3] - out_post[,4]), # estimate of differences between treatments A and B
                    "diff_BC" = mean(out_post[,4] - out_post[,5]), # estimate of differences between treatments B and C
                    "diff_AC" = mean(out_post[,3] - out_post[,5]), # estimate of differences between treatments A and C
                    "ci_diff_AB" = bayestestR::ci(out_post[,3] - out_post[,4], ci = ci, method = "HDI"), # x% credible intervals for the differences between A and B
                    "ci_diff_BC" = bayestestR::ci(out_post[,4] - out_post[,5], ci = ci, method = "HDI"), # x% credible intervals for the differences between B and C
                    "ci_diff_AC" = bayestestR::ci(out_post[,3] - out_post[,5], ci = ci, method = "HDI"), # x% credible intervals for the differences between A and C
                    "beta0_hat" = mean(out_post[,1]), # linkage parameter beta0 estimates
                    "beta1_hat" = mean(out_post[,2]), # linkage parameter beta1 estimates
                    "ci_beta0_hat" = bayestestR::ci(out_post[,1], ci = ci, method = "HDI"), # linkage parameter beta0 credible interval
                    "ci_beta1_hat"  = bayestestR::ci(out_post[,2], ci = ci, method = "HDI"),
                    "pi_DTR_est" = pi_DTR_est) # expected response rate of dynamic treatment regimens (DTRs)
    }else{
      result = list("posterior_sample" = out_post, # posterior samples of the link parameters and response rates generated through the MCMC process
                    "pi_hat_bjsm" = apply(out_post[,3:5],2,mean), # estimate of treatment response rate
                    "se_hat_bjsm" = apply(out_post[,3:5],2,stats::sd), # estimated standard error of the response rate
                    "ci_pi_A" = bayestestR::ci(out_post[,3], ci = ci, method = "HDI"), # x% credible intervals for the response rate of treatment A
                    "ci_pi_B" = bayestestR::ci(out_post[,4], ci = ci, method = "HDI"), # x% credible intervals for the response rate of treatment B
                    "ci_pi_C" = bayestestR::ci(out_post[,5], ci = ci, method = "HDI"), # x% credible intervals for the response rate of treatment C
                    "diff_AB" = apply(out_post[,3] - out_post[,4],2,mean), # estimate of differences between treatments A and B
                    "diff_BC" = apply(out_post[,4] - out_post[,5],2,mean), # estimate of differences between treatments B and C
                    "diff_AC" = apply(out_post[,3] - out_post[,5],2,mean), # estimate of differences between treatments A and C
                    "ci_diff_AB" = bayestestR::ci(out_post[,3] - out_post[,4], ci = ci, method = "HDI"), # x% credible intervals for the differences between A and B
                    "ci_diff_BC" = bayestestR::ci(out_post[,4] - out_post[,5], ci = ci, method = "HDI"), # x% credible intervals for the differences between B and C
                    "ci_diff_AC" = bayestestR::ci(out_post[,3] - out_post[,5], ci = ci, method = "HDI"), # x% credible intervals for the differences between A and C
                    "beta0_hat" = mean(out_post[,1]), # linkage parameter beta0 estimates
                    "beta1_hat" = mean(out_post[,2]), # linkage parameter beta1 estimates
                    "ci_beta0_hat" = bayestestR::ci(out_post[,1], ci = ci, method = "HDI"), # linkage parameter beta0 credible interval
                    "ci_beta1_hat"  = bayestestR::ci(out_post[,2], ci = ci, method = "HDI"))

    }
  }


  return(result)
}
