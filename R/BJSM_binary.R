#'
#'
#'
#'
#'
#'
#'
library(rjags)
library(bayestestR)

# generate trial dataset (with or without missing value)
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




BJSM_binary = function(data, NUM_ARMS, pi_prior.a, pi_prior.b, beta0_prior.a, beta0_prior.b, beta1_prior.a, beta1_prior.c, n_MCMC_chain, BURN.IN,
                       MCMC_SAMPLE, six = TRUE){
  # NUM_ARMS number of treatment arms
  # pi_prior.a  alpha parameter of the prior beta distribution for pi_1K
  # pi_prior.b  beta parameter of the prior beta distribution  for pi_1K
  # beta0_prior.a alpha parameter of the prior beta distribution for linkage parameter beta0
  # beta0_prior.b beta parameter of the prior beta distribution  for linkage parameter beta0
  # beta1_prior.a scale parameter of the prior pareto distribution for linkage parameter beta1
  # beta1_prior.c shape parameter of the prior pareto distribution for linkage parameter beta1
  # n_MCMC_chain number of MCMC chains, default to 1. If this is set to a number more than 1
  # BURN.IN number of burn-in iterations for MCMC
  # MCMC_SAMPLE number of iterations for MCMC
  mydata = data
  mydata$disc <- 2 * mydata$treatment_stageI - (mydata$response_stageI == 0)

  if (six == TRUE){
    # If using 6-betas model
    jag.model.name <- "BJSM_6betas_missing.bug"  # beta1 ~ pareto
    tryCatch({
      jag <- rjags::jags.model(file.path("data",jag.model.name),
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

    # If using 2-betas model
    jag.model.name <- "BJSM_2beta_missing.bug"  # beta1 ~ pareto
    tryCatch({
      jag <- rjags::jags.model(file.path("data",jag.model.name),
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
    result = list("posterior_sample" = out_post,
                  "pi_hat_bjsm" = apply(out_post[,7:9],2,mean),   # estimate of response rate
                  "se_hat_bjsm" = apply(out_post[,7:9],2,sd),     # standard error of the response rate
                  "ci_pi_A" = bayestestR::ci(out_post[,7], ci = 0.95, method = "HDI"),
                  "ci_pi_B" = bayestestR::ci(out_post[,8], ci = 0.95, method = "HDI"),
                  "ci_pi_C" = bayestestR::ci(out_post[,9], ci = 0.95, method = "HDI"),
                  "beta0_hat" = apply(out_post[,c(1,3,5)],2,mean),
                  "beta1_hat" = apply(out_post[,c(2,4,6)],2,mean))

  } else {
    result = list("posterior_sample" = out_post,
                  "pi_hat_bjsm" = apply(out_post[,3:5],2,mean),
                  "se_hat_bjsm" = apply(out_post[,3:5],2,sd),
                  "ci_pi_A" = bayestestR::ci(out_post[,3], ci = 0.95, method = "HDI"),
                  "ci_pi_B" = bayestestR::ci(out_post[,4], ci = 0.95, method = "HDI"),
                  "ci_pi_C" = bayestestR::ci(out_post[,5], ci = 0.95, method = "HDI"),
                  "beta0_hat" = mean(out_post[,1]),
                  "beta1_hat" = mean(out_post[,2]))
  }


  return(result)
}
