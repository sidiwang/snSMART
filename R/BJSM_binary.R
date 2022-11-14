#' BJSM for snSMART (3 active treatments/placebo and 2 dose level) with binary outcome
#'
#' This function implements the BJSM (Bayesian Joint Stage Modeling) method which
#' borrows information across both stages to estimate the individual response rate
#' of each treatment/dose level in a snSMART design with binary outcomes.
#'
#' @param data trial data with 4 columns: \code{treatment_stageI, response_stageI,
#'  treatment_stageII} and \code{response_stageII}. Missing data is allowed in stage 2.
#' @param pi_prior for 3 active treatment design: vector of six values (a, b, c, d, e, f),
#'  where a and b are the parameter \code{a} and parameter \code{b} of the prior distribution
#'  for \code{pi_1A}, c and d are the parameter \code{a} and parameter \code{b} of the prior
#'  distribution for \code{pi_1B}, and e and f are the parameter \code{a} and parameter
#'  \code{b} of the prior distribution for \code{pi_1C}. for dose level design: vector
#'  of two values (a, b). `a` is the parameter `a` of the prior distribution for
#'  \code{pi} (response rate) of placebo. `b` is the parameter `b` of the prior
#'  distribution for \code{pi} of placebo. Please check the `Details` section for
#'  more explanation
#' @param beta_prior for 3 active treatment design: vector of four values (a, b, c, d).
#'  `a` is the value of parameter \code{a} of the prior distribution for linkage parameter
#'  \code{beta_0} or \code{beta_0m}, `b` is the value of parameter \code{b} of the
#'  prior distribution for linkage parameter \code{beta_0} or \code{beta_0m}. `c`
#'  is the value of parameter \code{a} of the prior distribution for linkage parameter
#'  \code{beta_1} or \code{beta_1m}. `d` is the value of parameter \code{b} of the
#'  prior distribution for linkage parameter \code{beta_1} or \code{beta_1m}. for
#'  dose level design: vector of two values (a, b). `a` is the parameter `a` of the
#'  prior distribution for linkage parameter \code{beta}. `b` is the parameter b of
#'  the prior distribution  for linkage parameter \code{beta}. Please check the `Details`
#'  section for more explanation
#' @param n_MCMC_chain number of MCMC chains, default to 1.
#' @param normal.par for dose level design: vector of two values (normal.mean, normal.var).
#'  our function assumes that the logarithm of treatment effect ratio follows a Gaussian
#'  prior distribution \eqn{N(\mu, \sigma^2)}, that is \eqn{log(\pi_L/\pi_P)~N(normal.mean, normal.var)},
#'  and \eqn{log(\pi_H/\pi_P)~N(normal.mean, normal.var)}. \code{normal.mean} is the mean of
#'  this Gaussian prior. `normal.var` is the variance of this Gaussian prior distribution
#' @param n.adapt the number of iterations for adaptation
#' @param BURN.IN number of burn-in iterations for MCMC
#' @param thin thinning interval for monitors
#' @param MCMC_SAMPLE number of iterations for MCMC
#' @param ci coverage probability for credible intervals, default = 0.95
#' @param prior_dist for 3 active treatment design: vector of three values
#'  ("prior distribution for \code{pi}", "prior distribution for \code{beta0}",
#'  "prior distribution for \code{beta1}"). User can choose from "gamma", "beta", "pareto".
#'  e.g. prior_dist = c("beta", "beta", "pareto"); for dose level design: vector of two
#'  values ("prior distribution for \code{pi_P}", "prior distribution for \code{beta}")
#' @param six TRUE or FALSE. If TRUE, will run the six beta model (allow for estimating
#'  `beta_0m` and `beta_1m` values that differ among different treatments m), if FALSE
#'  will run the two beta model. default = TRUE. Only need to specify this for 3 active
#'  treatment design.
#' @param DTR TRUE or FALSE. If TRUE, will also return the expected response rate of
#'  dynamic treatment regimens. default = TRUE. Only need to specify this for 3 active
#'  treatment design.
#' @param verbose TRUE or FALSE. If FALSE, no function message and progress bar will be
#'  printed.
#' @param jags.model_options a list of optional arguments that are passed to \code{jags.model()} function.
#' @param coda.samples_options a list of optional arguments that are passed to \code{coda.samples()} function.
#'
#' @details
#' For \code{gamma} distribution, \code{prior.a} is the shape parameter \code{r}, \code{prior.b} is the rate parameter \code{lambda}. For \code{beta} distribution, \code{prior.a} is the shape parameter \code{a}, \code{prior.b} is the shape parameter \code{b}.
#' For \code{pareto} distribution, \code{prior.a} is the scale parameter \code{alpha}, \code{prior.b} is the shape parameter \code{c} (see page 29 of the jags user manual version 3.4.0). link: \url{http://www.stats.ox.ac.uk/~nicholls/MScMCMC14/jags_user_manual.pdf}
#'
#' The individual response rate is regarded as a permanent feature of the treatment. The second stage outcome is modeled conditionally on the first stage results linking the first and
#' second stage response probabilities through linkage parameters. The first stage response rate is denoted as \eqn{\pi_m} for treatment \eqn{m}. In the two \eqn{\beta} model, the second stage response rate for first stage responders is equal to \eqn{\beta_1\pi_m}. For nonresponders to treatment \eqn{m} in the first stage who
#' receive treatment \eqn{m'} in the second the stage, the second stage response rate in the second stage is equal to \eqn{\beta_0\pi_{m'}}. In the six \eqn{\beta} model, the second stage response rate of the first stage responders to treatment m is denoted by \eqn{\beta_{1m}\pi_m}, and the second stage response rate of the non-responders
#' to first stage treatment $m$ who receive treatment \eqn{m'} in the second stage is denoted by \eqn{\beta_{0m}\pi_{m'}}. All the \eqn{\beta}s are linkage parameters.
#'
#' Please refer to the paper listed under `reference` section for standard snSMART trial design and detailed definition of parameters.
#'
#' Note that this package does not include the JAGS library, users need to install JAGS separately. Please check this page for more details: \url{https://sourceforge.net/projects/mcmc-jags/}
#' @return
#' \item{posterior_sample}{an \code{mcmc.list} object generated through the \code{coda.samples()} function,
#'    which includes posterior samples of the link parameters and response rates generated through the MCMC
#'    process}
#' \item{pi_hat_bjsm}{estimate of response rate/treatment effect}
#'
#' \item{se_hat_bjsm}{standard error of the response rate}
#'
#' \item{ci_pi_A(P), ci_pi_B(L), ci_pi_C(H)}{x% credible intervals for treatment A(P), B(L), C(H)}
#'
#' \item{diff_AB(PL), diff_BC(LH). diff_AC(PH)}{estimate of differences between
#' treatments A(P) and B(L), B(L) and C(H), A(P) and C(H)}
#'
#' \item{ci_diff_AB(PL), ci_diff_BC(LH), ci_diff_AC(PH)}{x% credible intervals
#' for the estimated differences between treatments A(P) and B(L), B(L) and C(H), A(P) and C(H)}
#'
#' \item{se_AB(PL), se_BC(LH), se_AC(PH)}{standard error for the estimated differences
#' between treatments A(P) and B(L), B(L) and C(H), A(P) and C(H)}
#'
#' \item{beta0_hat, beta1_hat}{linkage parameter \code{beta0} and \code{beta1} estimates}
#'
#' \item{se_beta0_hat, se_beta1_hat}{standard error of the estimated value of linkage
#' parameter \code{beta0} and \code{beta1}}
#'
#' \item{ci_beta0_hat, ci_beta1_hat}{linkage parameter \code{beta0} and \code{beta1}
#' credible interval}
#'
#' \item{pi_DTR_est}{expected response rate of dynamic treatment regimens (DTRs)}
#'
#' \item{pi_DTR_se}{standard error for the estimated DTR response rate}
#'
#' \item{ci_pi_AB, ci_pi_AC, ci_pi_BA, ci_pi_BC, ci_pi_CA, ci_pi_CB}{x% credible intervals for the estimated DTR response rate}
#'
#' @importFrom stats update
#'
#' @examples
#' mydata <- data_binary
#'
#' BJSM_result <- BJSM_binary(
#'   data = mydata, prior_dist = c("beta", "beta", "pareto"),
#'   pi_prior = c(0.4, 1.6, 0.4, 1.6, 0.4, 1.6), beta_prior = c(1.6, 0.4, 3, 1),
#'   n_MCMC_chain = 1, n.adapt = 1000, MCMC_SAMPLE = 2000, ci = 0.95,
#'   six = TRUE, DTR = TRUE, verbose = FALSE
#' )
#'
#' \donttest{
#' BJSM_result2 <- BJSM_binary(
#'   data = mydata, prior_dist = c("beta", "beta", "pareto"),
#'   pi_prior = c(0.4, 1.6, 0.4, 1.6, 0.4, 1.6), beta_prior = c(1.6, 0.4, 3, 1),
#'   n_MCMC_chain = 1, n.adapt = 10000, MCMC_SAMPLE = 60000, ci = 0.95,
#'   six = FALSE, DTR = FALSE, verbose = FALSE
#' )
#'
#' summary(BJSM_result)
#' summary(BJSM_result2)
#' }
#'
#' data <- data_dose
#' BJSM_dose_result <- BJSM_binary(
#'   data = data_dose, prior_dist = c("beta", "gamma"),
#'   pi_prior = c(3, 17), normal.par = c(0.2, 100), beta_prior = c(2, 2),
#'   n_MCMC_chain = 2, n.adapt = 1000, MCMC_SAMPLE = 6000, ci = 0.95, verbose = FALSE
#' )
#'
#' summary(BJSM_dose_result)
#'
#' @references
#' Wei, B., Braun, T.M., Tamura, R.N. and Kidwell, K.M., 2018. A Bayesian analysis of small n sequential multiple assignment randomized trials (snSMARTs).
#' Statistics in medicine, 37(26), pp.3723-3732.
#'
#' Chao, Y.C., Trachtman, H., Gipson, D.S., Spino, C., Braun, T.M. and Kidwell, K.M., 2020. Dynamic treatment regimens in small n, sequential, multiple assignment, randomized trials: An application in focal segmental glomerulosclerosis. Contemporary clinical trials, 92, p.105989.
#'
#' Fang, F., Hochstedler, K.A., Tamura, R.N., Braun, T.M. and Kidwell, K.M., 2021. Bayesian methods to compare dose levels with placebo in a small n,
#' sequential, multiple assignment, randomized trial. Statistics in Medicine, 40(4), pp.963-977.
#'
#' @seealso
#' \code{\link{LPJSM_binary}} \cr
#' \code{\link{sample_size}}
#'
#' @rdname BJSM_binary
#' @export

BJSM_binary <- function(data, prior_dist, pi_prior, normal.par, beta_prior, n_MCMC_chain, n.adapt, BURN.IN = 100,
                        thin = 1, MCMC_SAMPLE, ci = 0.95, six = TRUE, DTR = TRUE,
                        jags.model_options = NULL, coda.samples_options = NULL, verbose = FALSE, ...) {
  quiet <- FALSE
  progress.bar <- "text"

  if (verbose == FALSE) {
    quiet <- TRUE
    progress.bar <- "none"
  }
  # bug files written to temporary directory on function call to satisfy CRAN
  # requirements of not accessing user's system files

  # "BJSM_6betas_missing.bug"
  BJSM_6betas_missing_file <- tempfile(fileext = ".bug")
  writeLines(BJSM_6betas_missing_text(), con = BJSM_6betas_missing_file)

  # "BJSM_6betas_missing_new.bug"
  BJSM_6betas_missing_new_file <- tempfile(fileext = ".bug")
  writeLines(BJSM_6betas_missing_new_text(), con = BJSM_6betas_missing_new_file)

  # "BJSM_2beta_missing.bug"
  BJSM_2beta_missing_file <- tempfile(fileext = ".bug")
  writeLines(BJSM_2beta_missing_text(), con = BJSM_2beta_missing_file)

  # "BJSM_2beta_missing_new.bug"
  BJSM_2beta_missing_new_file <- tempfile(fileext = ".bug")
  writeLines(BJSM_2beta_missing_new_text(), con = BJSM_2beta_missing_new_file)

  # "BJSM_dose.bug"
  BJSM_dose_file <- tempfile(fileext = ".bug")
  writeLines(BJSM_dose_text(), con = BJSM_dose_file)

  # "BJSM_dose_new.bug"
  BJSM_dose_new_file <- tempfile(fileext = ".bug")
  writeLines(BJSM_dose_new_text(), con = BJSM_dose_new_file)

  if (length(pi_prior) > 2) {
    pi_prior.a <- c(pi_prior[1], pi_prior[3], pi_prior[5])
    pi_prior.b <- c(pi_prior[2], pi_prior[4], pi_prior[6])
    NUM_ARMS <- length(unique(data$treatment_stageI[!is.na(data$treatment_stageI)]))
    beta0_prior.a <- beta_prior[1]
    beta0_prior.b <- beta_prior[2]
    beta1_prior.a <- beta_prior[3]
    beta1_prior.c <- beta_prior[4]
    pi_prior_dist <- prior_dist[1]
    beta0_prior_dist <- prior_dist[2]
    beta1_prior_dist <- prior_dist[3]

    pi_prior_dist <- ifelse(pi_prior_dist == "gamma", "dgamma", ifelse(pi_prior_dist == "beta", "dbeta", "dpar"))
    beta0_prior_dist <- ifelse(beta0_prior_dist == "gamma", "dgamma", ifelse(beta0_prior_dist == "beta", "dbeta", "dpar"))
    beta1_prior_dist <- ifelse(beta1_prior_dist == "gamma", "dgamma", ifelse(beta1_prior_dist == "beta", "dbeta", "dpar"))

    mydata <- data
    mydata$disc <- 2 * mydata$treatment_stageI - (mydata$response_stageI == 0)

    if (six == TRUE) {
      # If using 6-betas model
      # bugfile  <- readLines(system.file("BJSM_6betas_missing.bug", package = "snSMART"))
      bugfile <- readLines(BJSM_6betas_missing_file)
      bugfile <- gsub(pattern = "pi_prior_dist", replacement = pi_prior_dist, x = bugfile)
      bugfile <- gsub(pattern = "beta0_prior_dist", replacement = beta0_prior_dist, x = bugfile)
      bugfile2 <- gsub(pattern = "beta1_prior_dist", replacement = beta1_prior_dist, x = bugfile)

      bugfile2_file <- tempfile(fileext = ".bug")
      writeLines(bugfile2, con = bugfile2_file)
      jag.model.name <- bugfile2_file
      tryCatch(
        {
          jag <- do.call(rjags::jags.model, c(list(
            file = file.path(jag.model.name),
            data = list(
              n1 = nrow(mydata),
              n2 = nrow(mydata[!is.na(mydata$response_stageII), ]),
              num_arms = NUM_ARMS,
              Y1 = mydata$response_stageI,
              Y2 = mydata$response_stageII[!is.na(mydata$response_stageII)],
              treatment_stageI = mydata$treatment_stageI,
              treatment_stageII = mydata$treatment_stageII[!is.na(mydata$response_stageII)],
              response_stageI_disc = mydata$disc[!is.na(mydata$response_stageII)],
              # prior
              pi_prior.a = pi_prior.a,
              pi_prior.b = pi_prior.b,
              beta0_prior.a = beta0_prior.a,
              beta0_prior.b = beta0_prior.b,
              beta1_prior.a = beta1_prior.a, # pareto
              beta1_prior.c = beta1_prior.c # pareto
            ),
            n.chains = n_MCMC_chain, n.adapt = n.adapt, quiet = quiet, jags.model_options
          )))
          update(jag, BURN.IN, progress.bar = progress.bar)
          posterior_sample <- do.call(rjags::coda.samples, c(list(
            model = jag,
            variable.names = c("pi", "beta"),
            n.iter = MCMC_SAMPLE, thin = thin, progress.bar = progress.bar, coda.samples_options
          )))
        },
        warning = function(war) {},
        error = function(err) {},
        finally = {}
      )
    } else {
      # bugfile  <- readLines(system.file("BJSM_2beta_missing.bug", package = "snSMART"))
      bugfile <- readLines(BJSM_2beta_missing_file)
      bugfile <- gsub(pattern = "pi_prior_dist", replacement = pi_prior_dist, x = bugfile)
      bugfile <- gsub(pattern = "beta0_prior_dist", replacement = beta0_prior_dist, x = bugfile)
      bugfile2 <- gsub(pattern = "beta1_prior_dist", replacement = beta1_prior_dist, x = bugfile)

      bugfile2_file <- tempfile(fileext = ".bug")
      writeLines(bugfile2, con = bugfile2_file)
      # writeLines(bugfile2, con=system.file("BJSM_2beta_missing_new.bug", package = "snSMART"))

      # If using 2-betas model
      # jag.model.name <- system.file("BJSM_2beta_missing_new.bug", package = "snSMART")  # beta1 ~ pareto
      jag.model.name <- bugfile2_file
      tryCatch(
        {
          jag <- do.call(rjags::jags.model, c(list(
            file = file.path(jag.model.name),
            data = list(
              n1 = nrow(mydata),
              n2 = nrow(mydata[!is.na(mydata$response_stageII), ]),
              num_arms = NUM_ARMS,
              Y1 = mydata$response_stageI,
              Y2 = mydata$response_stageII[!is.na(mydata$response_stageII)],
              treatment_stageI = mydata$treatment_stageI,
              treatment_stageII = mydata$treatment_stageII[!is.na(mydata$response_stageII)],
              response1 = mydata$response_stageI[!is.na(mydata$response_stageII)] + 1,
              # prior
              pi_prior.a = pi_prior.a,
              pi_prior.b = pi_prior.b,
              beta0_prior.a = beta0_prior.a,
              beta0_prior.b = beta0_prior.b,
              beta1_prior.a = beta1_prior.a,
              beta1_prior.c = beta1_prior.c
            ),
            n.chains = n_MCMC_chain, n.adapt = n.adapt, quiet = quiet, jags.model_options
          )))
          update(jag, BURN.IN, progress.bar = progress.bar)
          posterior_sample <- do.call(rjags::coda.samples, c(list(
            model = jag,
            variable.names = c("pi", "beta"),
            n.iter = MCMC_SAMPLE, thin = thin, progress.bar = progress.bar, coda.samples_options
          )))
        },
        warning = function(war) {},
        error = function(err) {},
        finally = {}
      )
    }

    out_post <- as.data.frame(posterior_sample[[1]])


    if (six == TRUE) {
      colnames(out_post) <- c("beta0A", "beta1A", "beta0B", "beta1B", "beta0C", "beta1C", "pi_A", "pi_B", "pi_C")

      pi_AB_tt <- out_post[, 7]^2 * out_post[, 2] + (1 - out_post[, 7]) * out_post[, 8] * out_post[, 1]
      pi_AC_tt <- out_post[, 7]^2 * out_post[, 2] + (1 - out_post[, 7]) * out_post[, 9] * out_post[, 1]
      pi_BA_tt <- out_post[, 8]^2 * out_post[, 4] + (1 - out_post[, 8]) * out_post[, 7] * out_post[, 3]
      pi_BC_tt <- out_post[, 8]^2 * out_post[, 4] + (1 - out_post[, 8]) * out_post[, 9] * out_post[, 3]
      pi_CA_tt <- out_post[, 9]^2 * out_post[, 6] + (1 - out_post[, 9]) * out_post[, 7] * out_post[, 5]
      pi_CB_tt <- out_post[, 9]^2 * out_post[, 6] + (1 - out_post[, 9]) * out_post[, 8] * out_post[, 5]
      pi_DTR <- cbind(pi_AB_tt, pi_AC_tt, pi_BA_tt, pi_BC_tt, pi_CA_tt, pi_CB_tt)
      colnames(pi_DTR) <- c("rep_AB", "rep_AC", "rep_BA", "rep_BC", "rep_CA", "rep_CB")

      if (DTR == TRUE) {
        result <- list(
          "posterior_sample" = posterior_sample, # posterior samples of the link parameters and response rates generated through the MCMC process
          "pi_hat_bjsm" = apply(out_post[, 7:9], 2, mean), # estimate of response rate/treatment effect
          "se_hat_bjsm" = apply(out_post[, 7:9], 2, stats::sd), # standard error of the response rate
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
          "beta1_hat" = apply(out_post[, c(2, 4, 6)], 2, mean), # linkage parameter beta1 estimates
          "se_beta0_hat" = apply(out_post[, c(1, 3, 5)], 2, stats::sd),
          "se_beta1_hat" = apply(out_post[, c(2, 4, 6)], 2, stats::sd),
          "ci_beta0_hat" = HDInterval::hdi(out_post[, c(1, 3, 5)], ci), # linkage parameter beta0 credible interval
          "ci_beta1_hat" = HDInterval::hdi(out_post[, c(2, 4, 6)], ci), # linkage parameter beta1 credible interval
          "pi_DTR_est" = apply(pi_DTR, 2, mean), # expected response rate of dynamic treatment regimens (DTRs)
          "pi_DTR_se" = apply(pi_DTR, 2, stats::sd),
          "ci_pi_AB" = bayestestR::ci(pi_DTR[, 1], ci = ci, method = "HDI"),
          "ci_pi_AC" = bayestestR::ci(pi_DTR[, 2], ci = ci, method = "HDI"),
          "ci_pi_BA" = bayestestR::ci(pi_DTR[, 3], ci = ci, method = "HDI"),
          "ci_pi_BC" = bayestestR::ci(pi_DTR[, 4], ci = ci, method = "HDI"),
          "ci_pi_CA" = bayestestR::ci(pi_DTR[, 5], ci = ci, method = "HDI"),
          "ci_pi_CB" = bayestestR::ci(pi_DTR[, 6], ci = ci, method = "HDI")
        )
      } else {
        result <- list(
          "posterior_sample" = posterior_sample, # posterior samples of the link parameters and response rates generated through the MCMC process
          "pi_hat_bjsm" = apply(out_post[, 7:9], 2, mean), # estimate of response rate/treatment effect
          "se_hat_bjsm" = apply(out_post[, 7:9], 2, stats::sd), # standard error of the response rate
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
          "beta1_hat" = apply(out_post[, c(2, 4, 6)], 2, mean), # linkage parameter beta1 estimates
          "se_beta0_hat" = apply(out_post[, c(1, 3, 5)], 2, stats::sd),
          "se_beta1_hat" = apply(out_post[, c(2, 4, 6)], 2, stats::sd),
          "ci_beta0_hat" = HDInterval::hdi(out_post[, c(1, 3, 5)], ci), # linkage parameter beta0 credible interval
          "ci_beta1_hat" = HDInterval::hdi(out_post[, c(2, 4, 6)], ci)
        ) # linkage parameter beta1 credible interval
        #    "pi_DTR_est" = pi_DTR_est) # expected response rate of dynamic treatment regimens (DTRs)
      }
    } else {
      colnames(out_post) <- c("beta0", "beta1", "pi_A", "pi_B", "pi_C")

      pi_AB_tt <- out_post[, 3]^2 * out_post[, 2] + (1 - out_post[, 3]) * out_post[, 4] * out_post[, 1]
      pi_AC_tt <- out_post[, 3]^2 * out_post[, 2] + (1 - out_post[, 3]) * out_post[, 5] * out_post[, 1]
      pi_BA_tt <- out_post[, 4]^2 * out_post[, 2] + (1 - out_post[, 4]) * out_post[, 3] * out_post[, 1]
      pi_BC_tt <- out_post[, 4]^2 * out_post[, 2] + (1 - out_post[, 4]) * out_post[, 5] * out_post[, 1]
      pi_CA_tt <- out_post[, 5]^2 * out_post[, 2] + (1 - out_post[, 5]) * out_post[, 3] * out_post[, 1]
      pi_CB_tt <- out_post[, 5]^2 * out_post[, 2] + (1 - out_post[, 5]) * out_post[, 4] * out_post[, 1]
      pi_DTR <- cbind(pi_AB_tt, pi_AC_tt, pi_BA_tt, pi_BC_tt, pi_CA_tt, pi_CB_tt)
      colnames(pi_DTR) <- c("rep_AB", "rep_AC", "rep_BA", "rep_BC", "rep_CA", "rep_CB")

      if (DTR == TRUE) {
        result <- list(
          "posterior_sample" = posterior_sample, # posterior samples of the link parameters and response rates generated through the MCMC process
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
          "ci_beta1_hat" = bayestestR::ci(out_post[, 2], ci = ci, method = "HDI"),
          "pi_DTR_est" = apply(pi_DTR, 2, mean), # expected response rate of dynamic treatment regimens (DTRs)
          "pi_DTR_se" = apply(pi_DTR, 2, stats::sd),
          "ci_pi_AB" = bayestestR::ci(pi_DTR[, 1], ci = ci, method = "HDI"),
          "ci_pi_AC" = bayestestR::ci(pi_DTR[, 2], ci = ci, method = "HDI"),
          "ci_pi_BA" = bayestestR::ci(pi_DTR[, 3], ci = ci, method = "HDI"),
          "ci_pi_BC" = bayestestR::ci(pi_DTR[, 4], ci = ci, method = "HDI"),
          "ci_pi_CA" = bayestestR::ci(pi_DTR[, 5], ci = ci, method = "HDI"),
          "ci_pi_CB" = bayestestR::ci(pi_DTR[, 6], ci = ci, method = "HDI")
        )
      } else {
        result <- list(
          "posterior_sample" = posterior_sample, # posterior samples of the link parameters and response rates generated through the MCMC process
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
          "ci_beta1_hat" = bayestestR::ci(out_post[, 2], ci = ci, method = "HDI")
        )
      }
    }

    class(result) <- "BJSM_binary"
  } else {
    pi_prior_dist <- prior_dist[1]
    beta_prior_dist <- prior_dist[2]

    NUM_ARMS <- length(unique(data$treatment_stageI[!is.na(data$treatment_stageI)]))
    pi_prior.a <- pi_prior[1]
    pi_prior.b <- pi_prior[2]
    beta_prior.a <- beta_prior[1]
    beta_prior.b <- beta_prior[2]
    normal.mean <- normal.par[1]
    normal.var <- normal.par[2]
    beta_prior_dist <- ifelse(beta_prior_dist == "gamma", "dgamma", ifelse(beta_prior_dist == "beta", "dbeta", "dpar"))
    pi_prior_dist <- ifelse(pi_prior_dist == "gamma", "dgamma", ifelse(pi_prior_dist == "beta", "dbeta", "dpar"))

    mydata <- data
    mydata$response_status_stageI <- mydata$response_stageI + 1

    bugfile <- readLines(BJSM_dose_file)
    bugfile <- gsub(pattern = "pi_prior_dist", replacement = pi_prior_dist, x = bugfile)
    bugfile2 <- gsub(pattern = "beta_prior_dist", replacement = beta_prior_dist, x = bugfile)

    ##
    bugfile2_file <- tempfile(fileext = ".bug")
    writeLines(bugfile2, con = bugfile2_file)

    # writeLines(bugfile2, con = system.file("BJSM_dose_new.bug", package = "snSMART"))
    # jag.model.name <- system.file("BJSM_dose_new.bug", package = "snSMART")
    jag.model.name <- bugfile2_file

    tryCatch({
      jag <- do.call(rjags::jags.model, c(list(
        file = file.path(jag.model.name),
        data = list(
          overall_sample_size = nrow(mydata),
          num_arms = NUM_ARMS,
          response_stageI = mydata$response_stageI,
          response_stageII = mydata$response_stageII,
          treatment_stageI = mydata$treatment_stageI,
          treatment_stageII = mydata$treatment_stageII,
          response_discount_status_stageI = mydata$response_status_stageI,

          # prior
          a_pi = pi_prior.a,
          b_pi = pi_prior.b,
          a_beta = beta_prior.a,
          b_beta = beta_prior.b,
          normal.mean = normal.mean,
          normal.var = normal.var
        ),
        n.chains = n_MCMC_chain, n.adapt = n.adapt, quiet = quiet, jags.model_options
      )))
      update(jag, BURN.IN, progress.bar = progress.bar)
      posterior_sample <- do.call(rjags::coda.samples, c(list(
        model = jag,
        variable.names = c("pi", "beta"),
        n.iter = MCMC_SAMPLE, thin = thin, progress.bar = progress.bar, coda.samples_options
      )))
    })

    out_post <- as.data.frame(posterior_sample[[1]])
    colnames(out_post)[c(7:9)] <- c("pi_P", "pi_L", "pi_H")



    result <- list(
      "posterior_sample" = posterior_sample, # posterior samples of the link parameters and response rates generated through the MCMC process
      "pi_hat_bjsm" = apply(out_post[, 7:9], 2, mean), # estimate of response rate/treatment effect
      "se_hat_bjsm" = apply(out_post[, 7:9], 2, sd), # standard error of the response rate
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
      "ci_beta_hat" = HDInterval::hdi(out_post[, 1:6], ci)
    ) # linkage parameter beta estimates

    class(result) <- "BJSM_dose_binary"
  }

  return(result)
}


#' Summarizing BJSM fits
#'
#' `summary` method for class "`BJSM_binary`"
#'
#' @param object an object of class "`BJSM_binary`", usually, a result of a call to \code{\link{BJSM_binary}}
#' @param ... further arguments. Not currently used.
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
summary.BJSM_binary <- function(object, ...) {
  trteff <- cbind(object$pi_hat_bjsm, object$se_hat_bjsm, rbind(object$ci_pi_A, object$ci_pi_B, object$ci_pi_C))
  rownames(trteff) <- c("trtA", "trtB", "trtC")
  colnames(trteff) <- c("Estimate", "Std. Error", "C.I.", "CI low", "CI high")

  trtdiff <- cbind(rbind(object$diff_AB, object$diff_BC, object$diff_AC), rbind(object$se_AB, object$se_BC, object$se_AC), rbind(object$ci_diff_AB, object$ci_diff_BC, object$ci_diff_AC))
  rownames(trtdiff) <- c("diffAB", "diffBC", "diffAC")
  colnames(trtdiff) <- c("Estimate", "Std.Error", "C.I.", "CI low", "CI high")

  if (length(object$beta0_hat) == 1) {
    betaest <- rbind(as.matrix(cbind(object$beta0_hat, object$se_beta0_hat, object$ci_beta0_hat)), as.matrix(cbind(object$beta1_hat, object$se_beta1_hat, object$ci_beta1_hat)))
    colnames(betaest) <- c("Estimate", "Std. Error", "C.I.", "CI low", "CI high")
    rownames(betaest) <- c("beta0", "beta1")
  } else {
    betaest <- rbind(cbind(object$beta0_hat, object$se_beta0_hat, c(rep(trteff$C.I.[1], length(object$beta0_hat))), t(object$ci_beta0_hat)), cbind(object$beta1_hat, object$se_beta1_hat, c(rep(trteff$C.I.[1], length(object$beta1_hat))), t(object$ci_beta1_hat)))
    colnames(betaest) <- c("Estimate", "Std. Error", "C.I.", "CI low", "CI high")
  }

  if (!is.null(object$pi_DTR_est) == TRUE) {
    dtreff <- cbind(object$pi_DTR_est, object$pi_DTR_se, rbind(object$ci_pi_AB, object$ci_pi_AC, object$ci_pi_BA, object$ci_pi_BC, object$ci_pi_CA, object$ci_pi_CB))
    colnames(dtreff) <- c("Estimate", "Std. Error", "C.I.", "CI low", "CI high")
    obj <- list(trteff = trteff, trtdiff = trtdiff, betaest = betaest, dtreff = dtreff)
  } else {
    obj <- list(trteff = trteff, trtdiff = trtdiff, betaest = betaest)
  }


  class(obj) <- "summary.BJSM_binary"
  obj
}

#' @rdname BJSM_binary
#' @param x object to summarize.
#' @param ... further arguments. Not currently used.
#' @export
#' @export print.summary.BJSM_binary
print.summary.BJSM_binary <- function(x, ...) {
  cat("\nTreatment Effects Estimate:\n")
  print(x$trteff)
  cat("\nDifferences between Treatments:\n")
  print(x$trtdiff)
  cat("\nLinkage Parameter Estimate:\n")
  print(x$betaest)
  if (length(x) == 4) {
    cat("\nExpected Response Rate of Dynamic Treatment Regimens (DTR):\n")
    print(x$dtreff)
  }
  cat("\n")
}


#' @rdname BJSM_binary
#' @param x object to summarize.
#' @param ... further arguments. Not currently used.
#' @export
#' @export print.BJSM_binary

print.BJSM_binary <- function(x, ...) {
  cat("\nTreatment Effects Estimate:\n")
  trteff <- cbind(x$pi_hat_bjsm, x$se_hat_bjsm, rbind(x$ci_pi_A, x$ci_pi_B, x$ci_pi_C))
  rownames(trteff) <- c("trtA", "trtB", "trtC")
  colnames(trteff) <- c("Estimate", "Std. Error", "C.I.", "CI low", "CI high")
  print(trteff)
  cat("\nDifferences between Treatments:\n")
  trtdiff <- cbind(rbind(x$diff_AB, x$diff_BC, x$diff_AC), rbind(x$se_AB, x$se_BC, x$se_AC), rbind(x$ci_diff_AB, x$ci_diff_BC, x$ci_diff_AC))
  rownames(trtdiff) <- c("diffAB", "diffBC", "diffAC")
  colnames(trtdiff) <- c("Estimate", "Std.Error", "C.I.", "CI low", "CI high")
  print(trtdiff)
}




#' Summarizing BJSM fits
#'
#' `summary` method for class `BJSM_dose_binary`
#'
#' @param object an object of class `BJSM_dose_binary`, usually, a result of a call to \code{\link{BJSM_binary}}
#' @param ... further arguments. Not currently used.
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
summary.BJSM_dose_binary <- function(object, ...) {
  trteff <- cbind(object$pi_hat_bjsm, object$se_hat_bjsm, rbind(object$ci_pi_P, object$ci_pi_L, object$ci_pi_H))
  rownames(trteff) <- c("trtP", "trtL", "trtH")
  colnames(trteff) <- c("Estimate", "Std. Error", "C.I.", "CI low", "CI high")

  trtdiff <- cbind(rbind(object$diff_PL, object$diff_LH, object$diff_PH), rbind(object$se_PL, object$se_LH, object$se_PH), rbind(object$ci_diff_PL, object$ci_diff_LH, object$ci_diff_PH))
  rownames(trtdiff) <- c("diffPL", "diffLH", "diffPH")
  colnames(trtdiff) <- c("Estimate", "Std.Error", "C.I.", "CI low", "CI high")

  betaest <- t(rbind(object$beta_hat, object$se_beta, c(rep(trteff[, 3][1], length(object$beta_hat))), object$ci_beta_hat))
  colnames(betaest) <- c("Estimate", "Std. Error", "C.I.", "CI low", "CI high")

  obj <- list(trteff = trteff, trtdiff = trtdiff, betaest = betaest)
  class(obj) <- "summary.BJSM_dose_binary"
  obj
}

#' @rdname BJSM_binary
#' @param x object to summarize.
#' @param ... further arguments. Not currently used.
#' @export
#' @export print.summary.BJSM_dose_binary
print.summary.BJSM_dose_binary <- function(x, ...) {
  cat("\nTreatment Effects Estimate:\n")
  print(x$trteff)
  cat("\nDifferences between Treatments:\n")
  print(x$trtdiff)
  cat("\nLinkage Parameter Estimate:\n")
  print(x$betaest)
  cat("\n")
}


#' @rdname BJSM_binary
#' @param x object to summarize.
#' @param ... further arguments. Not currently used.
#' @export
#' @export print.BJSM_dose_binary
print.BJSM_dose_binary <- function(x, ...) {
  cat("\nTreatment Effects Estimate:\n")
  trteff <- cbind(x$pi_hat_bjsm, x$se_hat_bjsm, rbind(x$ci_pi_P, x$ci_pi_L, x$ci_pi_H))
  rownames(trteff) <- c("trtP", "trtL", "trtH")
  colnames(trteff) <- c("Estimate", "Std. Error", "C.I.", "CI low", "CI high")
  print(trteff)
  cat("\nDifferences between Treatments:\n")
  trtdiff <- cbind(rbind(x$diff_PL, x$diff_LH, x$diff_PH), rbind(x$se_PL, x$se_LH, x$se_PH), rbind(x$ci_diff_PL, x$ci_diff_LH, x$ci_diff_PH))
  rownames(trtdiff) <- c("diffPL", "diffLH", "diffPH")
  colnames(trtdiff) <- c("Estimate", "Std.Error", "C.I.", "CI low", "CI high")
  print(trtdiff)
}
