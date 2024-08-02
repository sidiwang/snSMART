
#' Group sequential data look 1
#'
#' sample synthetic dataset of group sequential trial design snSMART, can be used for interim analysis
#'
#' @name groupseqDATA_look1
#' @usage groupseqDATA_look1
#' @format This data frame contains the following columns:
#' \describe{
#'  \item{time.1st.trt}{first treatment time}
#'  \item{time.1st.resp}{first response time}
#'  \item{time.2nd.trt}{second treatment time}
#'  \item{time.2nd.resp}{second response time}
#'  \item{trt.1st}{treatment arm for first treatment}
#'  \item{resp.1st}{response for first treatment}
#'  \item{trt.2nd}{treatment arm for second treatment}
#'  \item{resp.2nd}{response for second treatment}
#' }
#' @docType data
#' @keywords data
#' @examples
#' mydata <- groupseqDATA_look1
#'
#' result1 <- group_seq(
#'   data = mydata, interim = TRUE, drop_threshold_pair = c(0.5, 0.4),
#'   prior_dist = c("beta", "beta", "pareto"), pi_prior = c(0.4, 1.6, 0.4, 1.6, 0.4, 1.6),
#'   beta_prior = c(1.6, 0.4, 3, 1), MCMC_SAMPLE = 6000, n.adapt = 1000, n_MCMC_chain = 1
#' )
#'
NULL

#' Group sequential full data
#'
#' sample synthetic dataset of group sequential trial design snSMART, can be used for final analysis
#' @name groupseqDATA_full
#' @usage groupseqDATA_full
#' @format This data frame contains the following columns:
#' \describe{
#'  \item{time.1st.trt}{first treatment time}
#'  \item{time.1st.resp}{first response time}
#'  \item{time.2nd.trt}{second treatment time}
#'  \item{time.2nd.resp}{second response time}
#'  \item{trt.1st}{treatment arm for first treatment}
#'  \item{resp.1st}{response for first treatment}
#'  \item{trt.2nd}{treatment arm for second treatment}
#'  \item{resp.2nd}{response for second treatment}
#' }
#' @docType data
#' @keywords data
#' @examples
#' mydata <- groupseqDATA_full
#' result2 <- group_seq(
#'   data = mydata, interim = FALSE, prior_dist = c(
#'     "beta", "beta", "pareto"
#'   ), pi_prior = c(0.4, 1.6, 0.4, 1.6, 0.4, 1.6),
#'   beta_prior = c(1.6, 0.4, 3, 1), MCMC_SAMPLE = 6000, n.adapt = 1000,
#'   n_MCMC_chain = 1, ci = 0.95, DTR = TRUE
#' )
NULL

#' Dataset with binary outcomes
#'
#' sample synthetic dataset of snSMART (3 active treatment) with binary outcomes
#' @name data_binary
#' @usage data_binary
#' @format This data frame contains the following columns:
#' \describe{
#'  \item{treatment_stageI}{treatment received in stage 1 - possible values: 1 (placebo), 2, 3}
#'  \item{response_stageI}{whether patients respond to stage 1 treatment - possible values: 0 (nonresponder), 1 (responder)}
#'  \item{treatment_stageII}{treatment received in stage 2 - possible values: 2, 3}
#'  \item{response_stageII}{whether patients respond to stage 2 treatment - possible values: 0 (nonresponder), 1 (responder)}
#' }
#' @docType data
#' @keywords data
#' @examples
#' mydata <- data_binary
#' LPJSM_result <- LPJSM_binary(data = mydata, six = TRUE, DTR = TRUE)
#'
NULL

#' Dose Level dataset with binary outcomes
#'
#' sample synthetic dataset of snSMART (dose level treatment) with binary outcomes
#' @name data_dose
#' @usage data_dose
#' @format This data frame contains the following columns:
#' \describe{
#'  \item{treatment_stageI}{treatment received in stage 1 - possible values: 1 (placebo), 2, 3}
#'  \item{response_stageI}{whether patients respond to stage 1 treatment - possible values: 0 (nonresponder), 1 (responder)}
#'  \item{treatment_stageII}{treatment received in stage 2 - possible values: 2, 3}
#'  \item{response_stageII}{whether patients respond to stage 2 treatment - possible values: 0 (nonresponder), 1 (responder)}
#' }
#' @docType data
#' @keywords data
#' @examples
#' mydata <- data_dose
#' BJSM_dose_result <- BJSM_binary(
#'   data = data_dose, prior_dist = c("beta", "gamma"),
#'   pi_prior = c(3, 17), normal.par = c(0.2, 100), beta_prior = c(2, 2),
#'   n_MCMC_chain = 2, n.adapt = 100, MCMC_SAMPLE = 2000, ci = 0.95
#' )
#'
NULL

#' Dataset with continuous outcomes
#'
#' sample synthetic dataset of snSMART (mapping function) with continuous outcomes
#' @name trialDataMF
#' @usage trialDataMF
#' @format This data frame contains the following columns:
#' \describe{
#'  \item{id}{participant ID}
#'  \item{trt1}{treatment received in stage 1 - possible values: 1 (placebo), 2, 3}
#'  \item{stage1outcome}{a number between 0-100 that represents the stage 1 treatment effect}
#'  \item{stay}{indicates whether the participant stayed on the same treatment arm in stage 2 - possible values: 0 (didn't stay), 1 (stayed)}
#'  \item{trt2}{treatment received in stage 2 - possible values: 2, 3}
#'  \item{stage2outcome}{a number between 0-100 that represents the stage 2 treatment effect}
#' }
#' @docType data
#' @keywords data
#' @examples
#' trialData <- trialDataMF
#'
#' BJSM_result <- BJSM_c(
#'   data = trialData, xi_prior.mean = c(50, 50, 50),
#'   xi_prior.sd = c(50, 50, 50), phi3_prior.sd = 20, n_MCMC_chain = 1,
#'   n.adapt = 1000, MCMC_SAMPLE = 5000, ci = 0.95, n.digits = 5
#' )
#'
#' summary(BJSM_result)
#' print(BJSM_result)
#'
NULL
