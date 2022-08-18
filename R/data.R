
#' group sequential data look 1
#'
#' sample dataset of group sequential trial design snSMART, can be used for interim analysis
#'
#' @name groupseqDATA_look1
#' @docType data
#' @keywords data
#' @examples
#' mydata = groupseqDATA_look1
#'
#' result1 = group_seq(data = mydata, interim = TRUE, drop_threshold_pair = c(0.5, 0.4),
#'   prior_dist = c("beta", "beta", "pareto"), pi_prior =  c(0.4, 1.6, 0.4, 1.6, 0.4, 1.6),
#'   beta_prior = c(1.6, 0.4, 3, 1), MCMC_SAMPLE = 6000, BURN.IN = 1000, n_MCMC_chain = 1)
#'
#'

NULL

#' group sequential full data
#'
#' sample dataset of group sequential trial design snSMART, can be used for final analysis
#' @name groupseqDATA_full
#' @docType data
#' @keywords data
#' @examples
#' mydata = groupseqDATA_full
#' result2 = group_seq(data = mydata, interim = FALSE, prior_dist = c("beta",
#'    "beta", "pareto"), pi_prior =  c(0.4, 1.6, 0.4, 1.6, 0.4, 1.6),
#'    beta_prior = c(1.6, 0.4, 3, 1), MCMC_SAMPLE = 6000, BURN.IN = 1000,
#'    n_MCMC_chain = 1, ci = 0.95, DTR = TRUE)
#'
NULL

#' Data Binary
#'
#' sample dataset of snSMART (3 active treatment) with binary outcomes
#' @name data_binary
#' @docType data
#' @keywords data
#' @examples
#' mydata = data_binary
#' LPJSM_result = LPJSM_binary(data = mydata, six = TRUE, DTR = TRUE)
#'
NULL

#' Data Dose Level
#'
#' sample dataset of snSMART (dose level treatment) with binary outcomes
#' @name data_dose
#' @docType data
#' @keywords data
#' @examples
#' mydata = data_dose
#' BJSM_dose_result = BJSM_binary(data = data_dose, prior_dist = c("beta", "gamma"),
#'     pi_prior = c(3, 17), normal.par = c(0.2, 100), beta_prior = c(2, 2),
#'     n_MCMC_chain = 2, BURN.IN = 1000, MCMC_SAMPLE = 6000, ci = 0.95)
#'
NULL

#' Data Mapping Function
#'
#' sample dataset of snSMART (mapping function) with continuous outcomes
#' @name trialDataMF
#' @docType data
#' @keywords data
#' @examples
#' trialData = trialDataMF
#'
#' BJSM_result = BJSM_c(data = trialData, xi_prior.mean = c(50, 50, 50),
#'     xi_prior.sd = c(50, 50, 50), phi3_prior.sd = 20, n_MCMC_chain = 1,
#'     n.adapt = 1000, MCMC_SAMPLE = 5000, ci = 0.95, n.digits = 5)
#'
#' summary(BJSM_result)
#' print(BJSM_result)
#'
NULL
