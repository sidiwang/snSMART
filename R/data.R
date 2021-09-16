
#' Patient Entry
#'
#' sample dataset of group sequential trial design snSMART, can be used for interim analysis
#'
#' @name patient_entry
#' @docType data
#' @keywords data
#' @examples
#' mydata = patient_entry
#'
#' result1 = group_seq(data = mydata, rule.type = 1, interim = TRUE, drop_threshold=0.5, NUM_ARMS = 3, pi_prior_dist = "beta", pi_prior.a =  c(0.4,0.4,0.4), pi_prior.b = c(1.6, 1.6, 1.6), beta0_prior_dist = "beta",
#' beta0_prior.a = 1.6, beta0_prior.b = 0.4, beta1_prior_dist = "pareto", beta1_prior.a = 3, beta1_prior.c = 1, MCMC_SAMPLE = 60000, BURN.IN = 10000,
#' n_MCMC_chain = 1)
#'
#' result2 = group_seq(data = mydata, rule.type = 2, interim = TRUE, drop_threshold_large=0.5, drop_threshold_small = 0.4, NUM_ARMS = 3, pi_prior_dist = "beta", pi_prior.a =  c(0.4,0.4,0.4), pi_prior.b = c(1.6, 1.6, 1.6), beta0_prior_dist = "beta",
#' beta0_prior.a = 1.6, beta0_prior.b = 0.4, beta1_prior_dist = "pareto", beta1_prior.a = 3, beta1_prior.c = 1, MCMC_SAMPLE = 60000, BURN.IN = 10000,
#' n_MCMC_chain = 1)
#'

NULL

#' Patient Entry (full)
#'
#' sample dataset of group sequential trial design snSMART, can be used for final analysis
#' @name patient_entry_full
#' @docType data
#' @keywords data
#' @examples
#' mydata = patient_entry_full
#' result3 = group_seq(data = mydata, rule.type = 2, interim = FALSE, NUM_ARMS = 3, pi_prior_dist = "beta", pi_prior.a =  c(0.4,0.4,0.4), pi_prior.b = c(1.6, 1.6, 1.6), beta0_prior_dist = "beta",
#' beta0_prior.a = 1.6, beta0_prior.b = 0.4, beta1_prior_dist = "pareto", beta1_prior.a = 3, beta1_prior.c = 1, MCMC_SAMPLE = 60000, BURN.IN = 10000,
#' n_MCMC_chain = 1, ci = 0.95, DTR = TRUE)
#'
NULL
