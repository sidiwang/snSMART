#' Trial dataset
#'
#' Generate trial dataset based on the exact number of responders/non-responders for standard snSMART (3 active treatments, non-responders re-randomized; binary outcome) with or without missing value.
#' Useful for recording real snSMART trial result or generating specific simulation scenario. Not to be confused with \code{\link{data_simulation}}. Missing data in stage 2 is allowed.
#'
#' @param trt vector of 3 values c(number of people who receive A in stage 1, number of people who receive B in stage 1, number of people who receive C in stage 1)
#' @param resp vector of 3 values c(number of people who respond to A in stage 1, number of people who respond to B in stage 1, number of people who respond to C in stage 1)
#' @param trt_same_II vector of 3 values c(number of 1st stage responders to A who receive A in stage 2, number of 1st stage responders to B who receive B in stage 2, number of 1st stage responders to C who receive C in stage 2)
#' @param resp_same_II vector of 3 values c(number of 1st stage responders who respond to A again in stage 2, number of 1st stage responders who respond to B again in stage 2, number of 1st stage responders who respond to C again in stage 2)
#' @param trt_negA vector of 2 values c(number of 1st stage non-responders to A who receive B in stage 2, number of 1st stage non-responders to A who receive C in stage 2)
#' @param trt_negB vector of 2 values c(number of 1st stage non-responders to B who receive A in stage 2, number of 1st stage non-responders to B who receive C in stage 2)
#' @param trt_negc vector of 2 values c(number of 1st stage non-responders to C who receive A in stage 2, number of 1st stage non-responders to C who receive B in stage 2)
#' @param resp_negA vector of 2 values c(number of 1st stage non-responders to A who respond to B in stage 2, number of 1st stage non-responders to A who respond to C in stage 2)
#' @param resp_negB vector of 2 values c(number of 1st stage non-responders to B who respond to A in stage 2, number of 1st stage non-responders to B who respond to C in stage 2)
#' @param resp_negC vector of 2 values c(number of 1st stage non-responders to C who respond to A in stage 2, number of 1st stage non-responders to C who respond to B in stage 2)
#'
#' @return a `matrix` of the trial dataset with 4 columns: treatment_stageI, response_stageI, treatment_stageII, response_stageII

#' @examples
#' trial_data = trial_dataset(trt = c(9, 12, 9), resp = c(3, 3, 5), trt_same_II = c(2, 2, 4),
#'     resp_same_II = c(2, 2, 4), trt_negA = c(3, 2), trt_negB = c(3, 2), trt_negc = c(1, 2),
#'     resp_negA = c(3, 1), resp_negB = c(2, 2), resp_negC = c(0, 0))
#'
#' @references
#' Wei, B., Braun, T.M., Tamura, R.N. and Kidwell, K.M., 2018. A Bayesian analysis of small n sequential multiple assignment randomized trials (snSMARTs).
#' Statistics in medicine, 37(26), pp.3723-3732.
#'
#' @seealso
#' \code{\link{data_simulation}} \cr
#' \code{\link{BJSM_binary}} \cr
#' \code{\link{JSRM_binary}} \cr
#' \code{\link{sample_size}}
#'
#' @export
#'

trial_dataset <- function(trt, resp, trt_same_II, resp_same_II,
                          trt_negA, trt_negB, trt_negc, resp_negA, resp_negB, resp_negC){

  trtA_I = trt[1]
  trtB_I = trt[2]
  trtC_I = trt[3]

  respA_I = resp[1]
  respB_I = resp[2]
  respC_I = resp[3]

  trtAA_II = trt_same_II[1]
  trtBB_II = trt_same_II[2]
  trtCC_II = trt_same_II[3]

  respA_II = resp_same_II[1]
  respB_II = resp_same_II[2]
  respC_II = resp_same_II[3]

  trtAB_II = trt_negA[1]
  trtAC_II = trt_negA[2]

  trtBA_II = trt_negB[1]
  trtBC_II = trt_negB[2]

  trtCA_II = trt_negc[1]
  trtCB_II = trt_negc[2]

  respAB_II = resp_negA[1]
  respAC_II = resp_negA[2]

  respBA_II = resp_negB[1]
  respBC_II = resp_negB[2]

  respCA_II = resp_negC[1]
  respCB_II = resp_negC[2]

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

  data_A = NULL

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
                           response_stageII = c(rep(1, respB_II),rep(0, trtBB_II - respB_II)))
  data_B.A.Y <- data.frame(treatment_stageI = rep(2, trtBA_II),
                           response_stageI = rep(0, trtBA_II),
                           treatment_stageII = rep(1, trtBA_II),
                           response_stageII = c(rep(1, respBA_II),rep(0, trtBA_II - respBA_II)))
  data_B.C.Y <- data.frame(treatment_stageI = rep(2, trtBC_II),
                           response_stageI = rep(0, trtBC_II),
                           treatment_stageII = rep(3, trtBC_II),
                           response_stageII = c(rep(1, respBC_II),rep(0, trtBC_II - respBC_II)))
  data_B = NULL

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
                           response_stageII = c(rep(1, respC_II),rep(0, trtCC_II - respC_II)))
  data_C.A.Y <- data.frame(treatment_stageI = rep(3, trtCA_II),
                           response_stageI = rep(0, trtCA_II),
                           treatment_stageII = rep(1, trtCA_II),
                           response_stageII = c(rep(1, respCA_II),rep(0, trtCA_II - respCA_II)))
  data_C.B.Y <- data.frame(treatment_stageI = rep(3, trtCB_II),
                           response_stageI = rep(0, trtCB_II),
                           treatment_stageII = rep(2, trtCB_II),
                           response_stageII = c(rep(1, respCB_II),rep(0, trtCB_II - respCB_II)))

  data_C = NULL

  if (trtC_I != sum(trtCC_II, trtCA_II, trtCB_II)){
    data_C <- data.frame(treatment_stageI = rep(3, trtC_I - trtCC_II - trtCA_II - trtCB_II),
                         response_stageI = c(rep(1, respC_I - trtCC_II),
                                             rep(0, trtC_I - respC_I - trtCA_II - trtCB_II)),
                         treatment_stageII = NA,
                         response_stageII = NA)
  }

  data_stageI.II <- rbind(data_A.A.Y, data_A.B.Y, data_A.C.Y, data_A,
                          data_B.B.Y, data_B.A.Y, data_B.C.Y, data_B,
                          data_C.C.Y, data_C.A.Y, data_C.B.Y, data_C)
  return(data_stageI.II)
}



#' BJSM binary (standard snSMART design)
#'
#' This function implements the BJSM (Bayesian Joint Stage Modeling) method which borrows information across both stages to estimate the individual response rate of each treatment in a standard snSMART design (different active treatments, non-responders re-randomized; binary outcome) based on trial dataset
#' generated by function \code{\link{trial_dataset}} or \code{\link{data_simulation}}
#'
#' @param data trial data generated through function \code{\link{trial_dataset}} or \code{\link{data_simulation}}
#' @param pi_prior.a parameter \code{a} of the prior distribution for \code{pi_1K}, where K = treatment A, treatment B, or treatment C. a vector with three values, one for each treatment. Please check the `Details` section for more explanation
#' @param pi_prior.b parameter \code{b} of the prior distribution  for \code{pi_1K}, where K = treatment A, treatment B, or treatment C. a vector with three values, one for each treatment. Please check the `Details` section for more explanation
#' @param beta0_prior vector of two values (a, b). `a` is the value of parameter \code{a} of the prior distribution for linkage parameter \code{beta0}, `b` is the value of parameter \code{b} of the prior distribution for linkage parameter \code{beta0}. Please check the `Details` section for more explanation
#' @param beta1_prior vector of two values (a, b). `a` is the value of parameter \code{a} of the prior distribution for linkage parameter \code{beta1}. `b` is the value of parameter \code{b} of the prior distribution for linkage parameter \code{beta1}. Please check the `Details` section for more explanation
#' @param n_MCMC_chain number of MCMC chains, default to 1.
#' @param BURN.IN number of burn-in iterations for MCMC
#' @param MCMC_SAMPLE number of iterations for MCMC
#' @param ci coverage probability for credible intervals, default = 0.95
#' @param prior_dist vector of three values ("prior distribution for \code{pi}", "prior distribution for \code{beta0}", "prior distribution for \code{beta1}"). User can choose from "gamma", "beta", "pareto". e.g. prior_dist = c("beta", "beta", "pareto")
#' @param six TRUE or FALSE. If TRUE, will run the six beta model (allow for estimating `beta_0k` and `beta_1k` values that differ among different treatments k), if FALSE will run the two beta model. default = TRUE.
#' @param DTR TRUE or FALSE. If TRUE, will also return the expected response rate of dynamic treatment regimens. default = TRUE
#' @param digits the number of significant digits to use when printing
#'
#' @details
#' For \code{gamma} distribution, \code{prior.a} is the shape parameter \code{r}, \code{prior.b} is the rate parameter \code{lambda}. For \code{beta} distribution, \code{prior.a} is the shape parameter \code{a}, \code{prior.b} is the shape parameter \code{b}.
#' For \code{pareto} distribution, \code{prior.a} is the scale parameter \code{alpha}, \code{prior.b} is the shape parameter \code{c} (see page 29 of the jags user manual version 3.4.0). link: \url{http://www.stats.ox.ac.uk/~nicholls/MScMCMC14/jags_user_manual.pdf}
#'
#' The individual response rate is regarded as a permanent feature of the treatment. The second stage outcome is modeled conditionally on the first stage results linking the first and
#' second stage response probabilities through linkage parameters. The first stage response rate is denoted as \eqn{\pi_k} for treatment \eqn{k}. In the two \eqn{\beta} model, the second stage response rate for first stage responders is equal to \eqn{\beta_1\pi_k}. For nonresponders to treatment \eqn{k} in the first stage who
#' receive treatment \eqn{k'} in the second the stage, the second stage response rate in the second stage is equal to \eqn{\beta_0\pi_{k'}}. In the six \eqn{\beta} model, the second stage response rate of the first stage responders to treatment k is denoted by \eqn{\beta_{1k}\pi_k}, and the second stage response rate of the non-responders
#' to first stage treatment $k$ who receive treatment \eqn{k'} in the second stage is denoted by \eqn{\beta_{0k}\pi_{k'}}. All the \eqn{\beta}s are linkage parameters.
#'
#' Please refer to the paper listed under `reference` section for standard snSMART trial design and detailed definition of parameters.
#'
#' @return
#' \describe{
#'    \item{posterior_sample}{posterior samples of the link parameters and response rates generated through the MCMC process}
#'    \item{pi_hat_bjsm}{estimate of response rate/treatment effect}
#'
#' \item{se_hat_bjsm}{standard error of the response rate}
#'
#' \item{ci_pi_A, ci_pi_B, ci_pi_C}{x% credible intervals for treatment A, B, C}
#'
#' \item{diff_AB, diff_BC. diff_AC}{estimate of differences between treatments A and B, B and C, A and C}
#'
#' \item{ci_diff_AB, ci_diff_BC, ci_diff_AC}{x% credible intervals for the estimated differences between treatments A and B, B and C, A and C}
#'
#' \item{se_AB, se_BC, se_AC}{standard error for the estimated differences between treatments A and B, B and C, A and C}
#'
#' \item{beta0_hat, beta1_hat}{linkage parameter \code{beta0} and \code{beta1} estimates}
#'
#' \item{se_beta0_hat, se_beta1_hat}{standard error of the estimated value of linkage parameter \code{beta0} and \code{beta1}}
#'
#' \item{ci_beta0_hat, ci_beta1_hat}{linkage parameter \code{beta0} and \code{beta1} credible interval}
#'
#' \item{pi_DTR_est}{expected response rate of dynamic treatment regimens (DTRs)}
#' }
#'
#' @examples
#' mydata = trial_dataset(trt = c(9, 12, 9), resp = c(3, 3, 5), trt_same_II = c(2, 2, 4),
#'     resp_same_II = c(2, 2, 4), trt_negA = c(3, 2), trt_negB = c(3, 2), trt_negc = c(1, 2),
#'     resp_negA = c(3, 1), resp_negB = c(2, 2), resp_negC = c(0, 0))
#'
#' BJSM_result = BJSM_binary(data = mydata, prior_dist = c("beta", "beta", "pareto"),
#'     pi_prior.a = c(0.4, 0.4, 0.4), pi_prior.b = c(1.6, 1.6, 1.6),
#'     beta0_prior = c(1.6, 0.4), beta1_prior = c(3, 1), n_MCMC_chain = 1, BURN.IN = 10000,
#'     MCMC_SAMPLE = 60000, ci = 0.95, six = TRUE, DTR = TRUE)
#'
#' BJSM_result2 = BJSM_binary(data = mydata, prior_dist = c("beta", "beta", "pareto"),
#'     pi_prior.a = c(0.4, 0.4, 0.4), pi_prior.b = c(1.6, 1.6, 1.6),
#'     beta0_prior = c(1.6, 0.4), beta1_prior = c(3, 1), n_MCMC_chain = 1, BURN.IN = 10000,
#'     MCMC_SAMPLE = 60000, ci = 0.95, six = FALSE, DTR = FALSE)
#'
#' summary(BJSM_result)
#' summary(BJSM_result2)
#'
#' @references
#' Wei, B., Braun, T.M., Tamura, R.N. and Kidwell, K.M., 2018. A Bayesian analysis of small n sequential multiple assignment randomized trials (snSMARTs).
#' Statistics in medicine, 37(26), pp.3723-3732.
#'
#' Chao, Y.C., Trachtman, H., Gipson, D.S., Spino, C., Braun, T.M. and Kidwell, K.M., 2020. Dynamic treatment regimens in small n, sequential, multiple assignment, randomized trials: An application in focal segmental glomerulosclerosis. Contemporary clinical trials, 92, p.105989.
#'
#' @seealso
#' \code{\link{data_simulation}} \cr
#' \code{\link{trial_dataset}} \cr
#' \code{\link{JSRM_binary}} \cr
#' \code{\link{sample_size}}
#'
#' @rdname BJSM_binary
#' @export

BJSM_binary = function(data, prior_dist, pi_prior.a, pi_prior.b, beta0_prior, beta1_prior, n_MCMC_chain, BURN.IN,
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

  NUM_ARMS = length(unique(data$treatment_stageI[!is.na(data$treatment_stageI)]))
  beta0_prior.a = beta0_prior[1]
  beta0_prior.b = beta0_prior[2]
  beta1_prior.a = beta1_prior[1]
  beta1_prior.c = beta1_prior[2]
  pi_prior_dist = prior_dist[1]
  beta0_prior_dist = prior_dist[2]
  beta1_prior_dist = prior_dist[3]

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
                        n.chains = n_MCMC_chain, n.adapt = BURN.IN)
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
                        data = list(n1 = nrow(mydata),
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
                        n.chains = n_MCMC_chain, n.adapt = BURN.IN)
      posterior_sample <- rjags::coda.samples(jag,
                                       c('pi', 'beta'),
                                       MCMC_SAMPLE)
    },
    warning = function(war){},
    error = function(err){},
    finally = {}
    )

  }

  out_post = as.data.frame(posterior_sample[[1]])


  if (six == TRUE){

    colnames(out_post) = c("beta0A", "beta1A", "beta0B", "beta1B", "beta0C", "beta1C", "pi_A", "pi_B", "pi_C")

    pi_DTR_est = c()
    pi_AB_tt <- out_post[, 7]^2 * out_post[, 2] + (1 - out_post[, 7]) * out_post[, 8] * out_post[, 1]
    pi_AC_tt <- out_post[, 7]^2 * out_post[, 2] + (1 - out_post[, 7]) * out_post[, 9] * out_post[, 1]
    pi_BA_tt <- out_post[, 8]^2 * out_post[, 4] + (1 - out_post[, 8]) * out_post[, 7] * out_post[, 3]
    pi_BC_tt <- out_post[, 8]^2 * out_post[, 4] + (1 - out_post[, 8]) * out_post[, 9] * out_post[, 3]
    pi_CA_tt <- out_post[, 9]^2 * out_post[, 6] + (1 - out_post[, 9]) * out_post[, 7] * out_post[, 5]
    pi_CB_tt <- out_post[, 9]^2 * out_post[, 6] + (1 - out_post[, 9]) * out_post[, 8] * out_post[, 5]
    pi_DTR_est <- rbind(pi_DTR_est, c(mean(pi_AB_tt), mean(pi_AC_tt), mean(pi_BA_tt), mean(pi_BC_tt), mean(pi_CA_tt), mean(pi_CB_tt)))
    colnames(pi_DTR_est) = c("rep_AB", "rep_AC", "rep_BA", "rep_BC", "rep_CA", "rep_CB")
    rownames(pi_DTR_est) = c("result")

    if (DTR == TRUE){

      result = list("posterior_sample" = out_post, # posterior samples of the link parameters and response rates generated through the MCMC process
                    "pi_hat_bjsm" = apply(out_post[, 7:9], 2, mean),   # estimate of response rate/treatment effect
                    "se_hat_bjsm" = apply(out_post[, 7:9], 2, stats::sd),     # standard error of the response rate
                    "ci_pi_A" = bayestestR::ci(out_post[, 7], ci = ci, method = "HDI"), # x% credible intervals for A
                    "ci_pi_B" = bayestestR::ci(out_post[, 8], ci = ci, method = "HDI"), # x% credible intervals for B
                    "ci_pi_C" = bayestestR::ci(out_post[, 9], ci = ci, method = "HDI"), # x% credible intervals for C
                    "diff_AB" = mean(out_post[, 7] - out_post[, 8]), # estimate of differences between treatments A and B
                    "se_AB" = stats::sd(out_post[, 7] - out_post[, 8]),
                    "diff_BC" = mean(out_post[, 8] - out_post[, 9]), # estimate of differences between treatments B and C
                    "se_BC" = stats::sd(out_post[, 8] - out_post[, 9]),
                    "diff_AC" = mean(out_post[, 7] - out_post[, 9]), # estimate of differences between treatments A and C
                    "se_AC" = stats::sd(out_post[, 7] - out_post[, 9]),
                    "ci_diff_AB" = bayestestR::ci(out_post[, 7] - out_post[, 8], ci = ci, method = "HDI"), # x% credible intervals for the differences between A and B
                    "ci_diff_BC" = bayestestR::ci(out_post[, 8] - out_post[, 9], ci = ci, method = "HDI"), # x% credible intervals for the differences between B and C
                    "ci_diff_AC" = bayestestR::ci(out_post[, 7] - out_post[, 9], ci = ci, method = "HDI"), # x% credible intervals for the differences between A and C
                    "beta0_hat" = apply(out_post[, c(1, 3, 5)], 2, mean), # linkage parameter beta0 estimates
                    "beta1_hat" = apply(out_post[, c(2, 4, 6)], 2, mean),  # linkage parameter beta1 estimates
                    "se_beta0_hat" = apply(out_post[, c(1, 3, 5)], 2, stats::sd),
                    "se_beta1_hat" = apply(out_post[, c(2, 4, 6)], 2, stats::sd),
                    "ci_beta0_hat" = HDInterval::hdi(out_post[, c(1, 3, 5)], ci), # linkage parameter beta0 credible interval
                    "ci_beta1_hat" = HDInterval::hdi(out_post[, c(2, 4, 6)], ci), # linkage parameter beta1 credible interval
                    "pi_DTR_est" = t(pi_DTR_est)) # expected response rate of dynamic treatment regimens (DTRs)
    }else{
      result = list("posterior_sample" = out_post, # posterior samples of the link parameters and response rates generated through the MCMC process
                    "pi_hat_bjsm" = apply(out_post[, 7:9], 2, mean),   # estimate of response rate/treatment effect
                    "se_hat_bjsm" = apply(out_post[, 7:9], 2, stats::sd),     # standard error of the response rate
                    "ci_pi_A" = bayestestR::ci(out_post[, 7], ci = ci, method = "HDI"), # x% credible intervals for A
                    "ci_pi_B" = bayestestR::ci(out_post[, 8], ci = ci, method = "HDI"), # x% credible intervals for B
                    "ci_pi_C" = bayestestR::ci(out_post[, 9], ci = ci, method = "HDI"), # x% credible intervals for C
                    "diff_AB" = mean(out_post[, 7] - out_post[, 8]), # estimate of differences between treatments A and B
                    "se_AB" = stats::sd(out_post[, 7] - out_post[, 8]),
                    "diff_BC" = mean(out_post[, 8] - out_post[, 9]), # estimate of differences between treatments B and C
                    "se_BC" = stats::sd(out_post[, 8] - out_post[, 9]),
                    "diff_AC" = mean(out_post[, 7] - out_post[, 9]), # estimate of differences between treatments A and C
                    "se_AC" = stats::sd(out_post[, 7] - out_post[, 9]),
                    "ci_diff_AB" = bayestestR::ci(out_post[, 7] - out_post[, 8], ci = ci, method = "HDI"), # x% credible intervals for the differences between A and B
                    "ci_diff_BC" = bayestestR::ci(out_post[, 8] - out_post[, 9], ci = ci, method = "HDI"), # x% credible intervals for the differences between B and C
                    "ci_diff_AC" = bayestestR::ci(out_post[, 7] - out_post[, 9], ci = ci, method = "HDI"), # x% credible intervals for the differences between A and C
                    "beta0_hat" = apply(out_post[, c(1, 3, 5)], 2, mean), # linkage parameter beta0 estimates
                    "beta1_hat" = apply(out_post[, c(2, 4, 6)], 2, mean),  # linkage parameter beta1 estimates
                    "se_beta0_hat" = apply(out_post[, c(1, 3, 5)], 2, stats::sd),
                    "se_beta1_hat" = apply(out_post[, c(2, 4, 6)], 2, stats::sd),
                    "ci_beta0_hat" = HDInterval::hdi(out_post[, c(1, 3, 5)], ci), # linkage parameter beta0 credible interval
                    "ci_beta1_hat" = HDInterval::hdi(out_post[, c(2, 4, 6)], ci)) # linkage parameter beta1 credible interval
                    #    "pi_DTR_est" = pi_DTR_est) # expected response rate of dynamic treatment regimens (DTRs)
    }
  } else {

    colnames(out_post) = c("beta0", "beta1", "pi_A", "pi_B", "pi_C")

    pi_DTR_est = c()
    pi_AB_tt <- out_post[, 3]^2 * out_post[, 2] + (1 - out_post[, 3]) * out_post[, 4] * out_post[, 1]
    pi_AC_tt <- out_post[, 3]^2 * out_post[, 2] + (1 - out_post[, 3]) * out_post[, 5] * out_post[, 1]
    pi_BA_tt <- out_post[, 4]^2 * out_post[, 2] + (1 - out_post[, 4]) * out_post[, 3] * out_post[, 1]
    pi_BC_tt <- out_post[, 4]^2 * out_post[, 2] + (1 - out_post[, 4]) * out_post[, 5] * out_post[, 1]
    pi_CA_tt <- out_post[, 5]^2 * out_post[, 2] + (1 - out_post[, 5]) * out_post[, 3] * out_post[, 1]
    pi_CB_tt <- out_post[, 5]^2 * out_post[, 2] + (1 - out_post[, 5]) * out_post[, 4] * out_post[, 1]
    pi_DTR_est <- rbind(pi_DTR_est, c(mean(pi_AB_tt), mean(pi_AC_tt), mean(pi_BA_tt), mean(pi_BC_tt), mean(pi_CA_tt), mean(pi_CB_tt)))
    colnames(pi_DTR_est) = c("rep_AB", "rep_AC", "rep_BA", "rep_BC", "rep_CA", "rep_CB")
    rownames(pi_DTR_est) = c("result")

    if (DTR == TRUE){
      result = list("posterior_sample" = out_post, # posterior samples of the link parameters and response rates generated through the MCMC process
                    "pi_hat_bjsm" = apply(out_post[, 3:5], 2, mean), # estimate of treatment response rate
                    "se_hat_bjsm" = apply(out_post[, 3:5], 2, stats::sd), # estimated standard error of the response rate
                    "ci_pi_A" = bayestestR::ci(out_post[, 3], ci = ci, method = "HDI"), # x% credible intervals for the response rate of treatment A
                    "ci_pi_B" = bayestestR::ci(out_post[, 4], ci = ci, method = "HDI"), # x% credible intervals for the response rate of treatment B
                    "ci_pi_C" = bayestestR::ci(out_post[, 5], ci = ci, method = "HDI"), # x% credible intervals for the response rate of treatment C
                    "diff_AB" = mean(out_post[, 3] - out_post[, 4]), # estimate of differences between treatments A and B
                    "se_AB" = stats::sd(out_post[, 3] - out_post[, 3]),
                    "diff_BC" = mean(out_post[, 4] - out_post[, 5]), # estimate of differences between treatments B and C
                    "se_BC" = stats::sd(out_post[, 4] - out_post[, 5]),
                    "diff_AC" = mean(out_post[, 3] - out_post[, 5]), # estimate of differences between treatments A and C
                    "se_AC" = stats::sd(out_post[, 3] - out_post[, 5]),
                    "ci_diff_AB" = bayestestR::ci(out_post[, 3] - out_post[, 4], ci = ci, method = "HDI"), # x% credible intervals for the differences between A and B
                    "ci_diff_BC" = bayestestR::ci(out_post[, 4] - out_post[, 5], ci = ci, method = "HDI"), # x% credible intervals for the differences between B and C
                    "ci_diff_AC" = bayestestR::ci(out_post[, 3] - out_post[, 5], ci = ci, method = "HDI"), # x% credible intervals for the differences between A and C
                    "beta0_hat" = mean(out_post[, 1]), # linkage parameter beta0 estimates
                    "beta1_hat" = mean(out_post[, 2]), # linkage parameter beta1 estimates
                    "se_beta0_hat" = stats::sd(out_post[, 1]),
                    "se_beta1_hat" = stats::sd(out_post[, 2]),
                    "ci_beta0_hat" = bayestestR::ci(out_post[, 1], ci = ci, method = "HDI"), # linkage parameter beta0 credible interval
                    "ci_beta1_hat"  = bayestestR::ci(out_post[, 2], ci = ci, method = "HDI"),
                    "pi_DTR_est" = t(pi_DTR_est)) # expected response rate of dynamic treatment regimens (DTRs)
    }else{
      result = list("posterior_sample" = out_post, # posterior samples of the link parameters and response rates generated through the MCMC process
                    "pi_hat_bjsm" = apply(out_post[, 3:5], 2, mean), # estimate of treatment response rate
                    "se_hat_bjsm" = apply(out_post[, 3:5], 2, stats::sd), # estimated standard error of the response rate
                    "ci_pi_A" = bayestestR::ci(out_post[, 3], ci = ci, method = "HDI"), # x% credible intervals for the response rate of treatment A
                    "ci_pi_B" = bayestestR::ci(out_post[, 4], ci = ci, method = "HDI"), # x% credible intervals for the response rate of treatment B
                    "ci_pi_C" = bayestestR::ci(out_post[, 5], ci = ci, method = "HDI"), # x% credible intervals for the response rate of treatment C
                    "diff_AB" = mean(out_post[, 3] - out_post[, 4]), # estimate of differences between treatments A and B
                    "se_AB" = stats::sd(out_post[, 3] - out_post[, 3]),
                    "diff_BC" = mean(out_post[, 4] - out_post[, 5]), # estimate of differences between treatments B and C
                    "se_BC" = stats::sd(out_post[, 4] - out_post[, 5]),
                    "diff_AC" = mean(out_post[, 3] - out_post[, 5]), # estimate of differences between treatments A and C
                    "se_AC" = stats::sd(out_post[, 3] - out_post[, 5]),
                    "ci_diff_AB" = bayestestR::ci(out_post[, 3] - out_post[, 4], ci = ci, method = "HDI"), # x% credible intervals for the differences between A and B
                    "ci_diff_BC" = bayestestR::ci(out_post[, 4] - out_post[, 5], ci = ci, method = "HDI"), # x% credible intervals for the differences between B and C
                    "ci_diff_AC" = bayestestR::ci(out_post[, 3] - out_post[, 5], ci = ci, method = "HDI"), # x% credible intervals for the differences between A and C
                    "beta0_hat" = mean(out_post[, 1]), # linkage parameter beta0 estimates
                    "beta1_hat" = mean(out_post[, 2]), # linkage parameter beta1 estimates
                    "se_beta0_hat" = stats::sd(out_post[, 1]),
                    "se_beta1_hat" = stats::sd(out_post[, 2]),
                    "ci_beta0_hat" = bayestestR::ci(out_post[, 1], ci = ci, method = "HDI"), # linkage parameter beta0 credible interval
                    "ci_beta1_hat"  = bayestestR::ci(out_post[, 2], ci = ci, method = "HDI"))

    }
  }

  class(result) = "BJSM_binary"

  return(result)
}


#' Summarizing BJSM fits
#'
#' `summary` method for class "`BJSM_binary`"
#'
#' @param object an object of class "`BJSM_binary`", usually, a result of a call to \code{\link{BJSM_binary}}
#' @param digits the number of significant digits to use when printing
#'
#' @returns
#' \describe{
#'    \item{Treatment Effects Estimate}{a 3 x 5 matrix with columns for the estimated treatment effects, its standard error, coverage probability of its credible interval, lower bound for its credible interval and higher bound for its credible interval}
#'    \item{Differences between Treatments}{a 3 x 5 matrix with columns for the estimated differences in treatment effects between two treatments, its standard error, coverage probability of its credible interval, lower bound and higher bound of the credible interval}
#'    \item{Linkage Parameter Estimate}{a 2 x 5 matrix, if the two beta model is fitted, or a 6 x 5 matrix, if the six beta model is fitted, with columns for the estimated linkage parameters}
#'    \item{Expected Response Rate of Dynamic Treatment Regimens (DTR)}{}
#' }
#'
#'
#' @export
summary.BJSM_binary = function(object, digits = 5, ...){
  cat("\nTreatment Effects Estimate:\n")
  trteff = cbind(object$pi_hat_bjsm, object$se_hat_bjsm, rbind(object$ci_pi_A, object$ci_pi_B, object$ci_pi_C))
  rownames(trteff) = c("trtA", "trtB", "trtC")
  colnames(trteff) = c("Estimate", "Std. Error", "C.I.", "CI low", "CI high")
  print(trteff, digits = digits)
  cat("\nDifferences between Treatments:\n")
  trtdiff = cbind(rbind(object$diff_AB, object$diff_BC, object$diff_AC), rbind(object$se_AB, object$se_BC, object$se_AC), rbind(object$ci_diff_AB, object$ci_diff_BC, object$ci_diff_AC))
  rownames(trtdiff) = c("diffAB", "diffBC", "diffAC")
  colnames(trtdiff) = c("Estimate", "Std.Error", "C.I.", "CI low", "CI high")
  print(trtdiff, digits = digits)
  cat("\nLinkage Parameter Estimate:\n")
  if (length(object$beta0_hat) == 1){
    betaest = rbind(as.matrix(cbind(object$beta0_hat, object$se_beta0_hat, object$ci_beta0_hat)), as.matrix(cbind(object$beta1_hat, object$se_beta1_hat, object$ci_beta1_hat)))
    colnames(betaest) = c("Estimate", "Std. Error", "C.I.", "CI low", "CI high")
    rownames(betaest) = c("beta0", "beta1")
    } else {
    betaest = rbind(cbind(object$beta0_hat, object$se_beta0_hat, c(rep(trteff$C.I.[1], length(object$beta0_hat))), t(object$ci_beta0_hat)), cbind(object$beta1_hat, object$se_beta1_hat, c(rep(trteff$C.I.[1], length(object$beta1_hat))), t(object$ci_beta1_hat)))
    colnames(betaest) = c("Estimate", "Std. Error", "C.I.", "CI low", "CI high")
  }

  print(betaest, digits = digits)
  if (!is.null(object$pi_DTR_est) == TRUE){
    cat("\nExpected Response Rate of Dynamic Treatment Regimens (DTR):\n")
    print(object$pi_DTR_est, digits = digits)
  }
  cat("\n")
}


#' @rdname BJSM_binary
#' @export
print.BJSM_binary = function(object, digits = 5, ...){
  cat("\nTreatment Effects Estimate:\n")
  print(object$pi_hat_bjsm)
  cat("\nDifferences between Treatments:\n")
  trtdiff = rbind(object$diff_AB, object$diff_BC, object$diff_AC)
  colnames(trtdiff) = c("estimate")
  rownames(trtdiff) = c("diffAB", "diffBC", "diffAC")
  print(trtdiff)
  cat("\nLinkage Parameter Estimate:\n")
  ret = cbind(object$beta0_hat, object$beta1_hat)
  colnames(ret) = c("Estimate", "Std. Error")
  print(ret)
  cat("\n")
}

