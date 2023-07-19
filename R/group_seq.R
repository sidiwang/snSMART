#' BJSM method for interim analysis and final analysis of group sequential trial design
#'
#' After obtain real trial data, this function can be used to decide which arm to
#' drop in an interim analysis or provide a full final analysis.
#'
#'
#' @param data dataset should include 8 columns: `time.1st.trt` (first treatment
#'  starts time), `time.1st.resp` (first response time), `time.2nd.trt` (second
#'  treatment starts time), `time.2nd.resp` (second response time),
#'  `trt.1st` (treatment arm for first treatment), `resp.1st` (response for first
#'  treatment), `trt.2nd` (treatment arm for second treatment), `resp.2nd` (response
#'  for second
#'  treatment) data yet to be observed should be marked as "`NA`"
#' @param interim indicates whether user is conducting an interim analysis via BJSM (`interim` = TRUE) or an final analysis via BJSM (`interim` = FALSE)
#' @param drop_threshold_pair a vector of 2 values (`drop_threshold_tau_l`, `drop_threshold_psi_l`). Both `drop_threshold_tau_l` and `drop_threshold_psi_l` should be between 0 and 1. only assign value to this parameter when `interim = TRUE`. See the details section for more explanation
#' @param pi_prior vector of six values (a, b, c, d, e, f), where a and b are the parameter \code{a} and parameter \code{b} of the prior distribution for \code{pi_1A}, c and d are the parameter \code{a} and parameter \code{b} of the prior distribution for \code{pi_1B}, and e and f are the parameter \code{a} and parameter \code{b} of the prior distribution for \code{pi_1C}. Please check the `Details` section for more explanation
#' @param beta_prior vector of four values (`beta0_prior.a`, `beta0_prior.b`, `beta1_prior.a`, `beta1_prior.c`).  `beta0_prior.a` is the parameter a of the prior distribution for linkage parameter `beta0`. `beta0_prior.b` is the parameter b of the prior distribution  for linkage parameter `beta0`. `beta1_prior.a` is the parameter a of the prior distribution for linkage parameter `beta1`. `beta1_prior.c` is the parameter b of the prior distribution for linkage parameter `beta1`. Please check the `Details` section for more explanation
#' @param n_MCMC_chain number of MCMC chains, default to 1
#' @param n.adapt the number of iterations for adaptation
#' @param BURN.IN number of burn-in iterations for MCMC
#' @param thin thinning interval for monitors
#' @param MCMC_SAMPLE number of iterations for MCMC
#' @param prior_dist vector of three values ("prior distribution for \code{pi}",
#'  "prior distribution for \code{beta0}", "prior distribution for \code{beta1}"),
#'  user can choose from "gamma", "beta", "pareto". e.g. prior_dist = c("beta",
#'  "beta", "pareto")
#' @param ci coverage probability for credible intervals, default = 0.95. only
#'  assign value to this parameter when `interim = FALSE`.
#' @param DTR, if TRUE, will also return the expected response rate of dynamic
#'  treatment regimens. default = TRUE. only assign value to this parameter when
#'  `interim = FALSE`.
#' @param verbose TRUE or FALSE. If FALSE, no function message and progress bar will be
#'  printed.
#' @param jags.model_options a list of optional arguments that are passed to \code{jags.model()} function.
#' @param coda.samples_options a list of optional arguments that are passed to \code{coda.samples()} function.
#'
#' @details
#' For \code{gamma} distribution, \code{prior.a} is the shape parameter \code{r},
#' \code{prior.b} is the rate parameter \code{lambda}. For \code{beta} distribution,
#' \code{prior.a} is the shape parameter \code{a}, \code{prior.b} is the shape parameter
#' \code{b}.
#' For \code{pareto} distribution, \code{prior.a} is the scale parameter \code{alpha},
#' \code{prior.b} is the shape parameter \code{c} (see jags user manual).
#' The individual response rate is regarded as a permanent feature of the treatment.
#' The second stage outcome is modeled conditionally on the first stage results
#' linking the first and
#' second stage response probabilities through linkage parameters.
#'
#' (paper provided in the reference section, section 2.2.2 Bayesian decision rules. drop_threshold_tau_l and drop_threshold_psi_l correspond to \eqn{tau_l} and \eqn{psi_l} respectively)
#'
#' Please refer to the paper listed under `reference` section for detailed definition of parameters.
#' Note that this package does not include the JAGS library, users need to install JAGS separately. Please check this page for more details: \url{https://sourceforge.net/projects/mcmc-jags/}

#' @examples
#' mydata <- groupseqDATA_look1
#'
#' result1 <- group_seq(
#'   data = mydata, interim = TRUE, drop_threshold_pair = c(0.5, 0.4),
#'   prior_dist = c("beta", "beta", "pareto"), pi_prior = c(0.4, 1.6, 0.4, 1.6, 0.4, 1.6),
#'   beta_prior = c(1.6, 0.4, 3, 1), MCMC_SAMPLE = 6000, n.adapt = 1000, n_MCMC_chain = 1
#' )
#'
#' summary(result1)
#'
#'
#' mydata <- groupseqDATA_full
#' result2 <- group_seq(
#'   data = mydata, interim = FALSE, prior_dist = c("beta", "beta", "pareto"),
#'   pi_prior = c(0.4, 1.6, 0.4, 1.6, 0.4, 1.6),
#'   beta_prior = c(1.6, 0.4, 3, 1), MCMC_SAMPLE = 6000, n.adapt = 1000,
#'   n_MCMC_chain = 1, ci = 0.95, DTR = TRUE
#' )
#'
#' summary(result2)
#'
#' @return
#' if `interim = TRUE`, this function returns either 0 - no arm is dropped,
#' or A/B/C - arm A/B/C is dropped \cr
#'
#' if `interim = FALSE`, this function returns:
#'
#' \item{posterior_sample}{an \code{mcmc.list} object generated through the \code{coda.samples()} function,
#'    which includes posterior samples of the link parameters and response rates generated through the MCMC
#'    process}
#' \item{pi_hat_bjsm}{estimate of response rate/treatment effect}
#'
#' \item{se_hat_bjsm}{standard error of the response rate}
#'
#' \item{ci_pi_A, ci_pi_B, ci_pi_C}{x% credible intervals for treatment A, B, C}
#'
#' \item{diff_AB, diff_BC. diff_AC}{estimate of differences between treatments A
#' and B, B and C, A and C}
#'
#' \item{ci_diff_AB, ci_diff_BC, ci_diff_AC}{x% credible intervals for the differences
#' between treatments A and B, B and C, A and C}
#'
#' \item{se_AB, se_BC, se_AC}{standard error for the differences between treatments
#' A and B, B and C, A and C}
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
#'
#' @references
#' Chao, Y.C., Braun, T.M., Tamura, R.N. and Kidwell, K.M., 2020. A Bayesian group
#' sequential small n sequential multiple‐assignment randomized trial. Journal of
#' the Royal Statistical Society: Series C (Applied Statistics), 69(3), pp.663-680.
#'
# #' @seealso
# #' \code{\link{sim_group_seq}}
#'
#' @export
#'
#' @rdname group_seq
group_seq <- function(data, interim = TRUE, drop_threshold_pair = NULL, prior_dist, pi_prior,
                      beta_prior, MCMC_SAMPLE, n.adapt, thin = 1, BURN.IN = 100, n_MCMC_chain,
                      ci = 0.95, DTR = TRUE,
                      jags.model_options = NULL, coda.samples_options = NULL, verbose = FALSE, ...) {
  quiet <- FALSE
  progress.bar <- "text"

  if (verbose == FALSE) {
    quiet <- TRUE
    progress.bar <- "none"
  }

  # bug files written to temporary directory on function call to satisfy CRAN
  # requirements of not accessing user's system files

  # "Bayes_AR.bug"
  Bayes_AR_file <- tempfile(fileext = ".bug")
  writeLines(Bayes_AR_text(), con = Bayes_AR_file)

  # "Bayes_AR_new.bug"
  Bayes_AR_new_file <- tempfile(fileext = ".bug")
  writeLines(Bayes_AR_new_text(), con = Bayes_AR_new_file)

  # "Bayes.bug"
  Bayes_file <- tempfile(fileext = ".bug")
  writeLines(Bayes_text(), con = Bayes_file)

  # "Bayes_new.bug"
  Bayes_new_file <- tempfile(fileext = ".bug")
  writeLines(Bayes_new_text(), con = Bayes_new_file)


  if (!is.null(drop_threshold_pair)) {
    drop_threshold_large <- drop_threshold_pair[1]
    drop_threshold_small <- drop_threshold_pair[2]
  }

  NUM_ARMS <- length(unique(data$trt.1st[!is.na(data$trt.1st)]))

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


  assn.stage2 <- function(i, trt, y, rand.prob) # Function that assigns the second stage treatment
  {
    alltrt <- 1:3
    if (y[i] == 1) newtrt <- trt[i]
    if (y[i] == 0) newtrt <- sample(alltrt[alltrt != trt[i]], 1, prob = rand.prob[alltrt != trt[i]])
    return(newtrt)
  }

  if (interim == TRUE) {

    # pi_hat <- matrix(NA,nrow=n.update,ncol=3)
    patient_entry <- data

    patient_entry$disc <- 2 * patient_entry$trt.1st - (patient_entry$resp.1st == 0)

    bugfile <- readLines(Bayes_AR_file)
    bugfile <- gsub(pattern = "pi_prior_dist", replacement = pi_prior_dist, x = bugfile)
    bugfile <- gsub(pattern = "beta0_prior_dist", replacement = beta0_prior_dist, x = bugfile)
    bugfile2 <- gsub(pattern = "beta1_prior_dist", replacement = beta1_prior_dist, x = bugfile)

    bugfile2_file <- tempfile(fileext = ".bug")
    writeLines(bugfile2, con = bugfile2_file)

    error_ind <- 0
    error_count <- 0
    tryCatch(
      {
        jags <- do.call(rjags::jags.model, c(list(
          file = bugfile2_file,
          data = list(
            n1 = nrow(patient_entry[!is.na(patient_entry$resp.1st), ]),
            n2 = nrow(patient_entry[!is.na(patient_entry$resp.2nd), ]),
            num_arms = NUM_ARMS,
            Y1 = patient_entry$resp.1st,
            Y2 = patient_entry$resp.2nd[!is.na(patient_entry$resp.2nd)],
            treatment_stageI = patient_entry$trt.1st,
            treatment_stageII = patient_entry$trt.2nd[!is.na(patient_entry$resp.2nd)],
            response_stageI_disc = patient_entry$disc[!is.na(patient_entry$resp.2nd)],
            # prior
            pi_prior.a = pi_prior[c(1, 3, 5)],
            pi_prior.b = pi_prior[c(2, 4, 6)],
            beta0_prior.a = beta0_prior.a,
            beta0_prior.b = beta0_prior.b,
            beta1_prior.a = beta1_prior.a,
            beta1_prior.c = beta1_prior.c
          ),
          n.chains = n_MCMC_chain, n.adapt = n.adapt, quiet = quiet, jags.model_options
        )))
        update(jags, BURN.IN, progress.bar = progress.bar)
        posterior_sample <- do.call(rjags::coda.samples, c(list(
          model = jags,
          variable.names = c("pi", "beta"),
          n.iter = MCMC_SAMPLE,
          thin = thin, progress.bar = progress.bar, coda.samples_options
        )))
      },
      warning = function(war) {
        warning_count <- warning_count + 1
        err_war_message <- rbind(paste("The warning ", warning_count, " is: ", war))
      },
      error = function(err) {
        error_count <- error_count + 1
        err_war_message <- rbind(paste("The error ", error_count, " is: ", err))
        error_ind <- 1
      },
      finally = {
      }
    )
    out_post <- posterior_sample[[1]]

    min_A <- mean(apply(out_post[, 7:9], 1, function(x) {
      x[1] == min(x)
    })) # posterior probability that A has smallest response rate
    min_B <- mean(apply(out_post[, 7:9], 1, function(x) {
      x[2] == min(x)
    })) # posterior probability that B has smallest response rate
    min_C <- mean(apply(out_post[, 7:9], 1, function(x) {
      x[3] == min(x)
    })) # posterior probability that C has smallest response rate
    max_A <- mean(apply(out_post[, 7:9], 1, function(x) {
      x[1] == max(x)
    })) # posterior probability that A has largest response rate
    max_B <- mean(apply(out_post[, 7:9], 1, function(x) {
      x[2] == max(x)
    })) # posterior probability that B has largest response rate
    max_C <- mean(apply(out_post[, 7:9], 1, function(x) {
      x[3] == max(x)
    })) # posterior probability that C has largest response rate
    keep_A <- (max_A > drop_threshold_large)
    message("\nInterim Analysis Outcome:\n")
    message("Threshold tau_l is set to: ")
    message(drop_threshold_large)
    message("\nThreshold psi_l is set to: ")
    message(drop_threshold_small)
    if (keep_A == 1) {
      message("\nStep 1: Arm A's interim posterior probability of having the greatest response is bigger than threshold ")
      message(drop_threshold_large)
      message("\n")
    }
    keep_B <- (max_B > drop_threshold_large)
    if (keep_B == 1) {
      message("\nStep 1: Arm B's interim posterior probability of having the greatest response is bigger than threshold ")
      message(drop_threshold_large)
      message("\n")
    }
    keep_C <- (max_C > drop_threshold_large)
    if (keep_C == 1) {
      message("\nStep 1: Arm C's interim posterior probability of having the greatest response is bigger than threshold ")
      message(drop_threshold_large)
      message("\n")
    }
    if (any(c(keep_A, keep_B, keep_C) > 0)) {
      drop_A <- ((keep_A == 0) * (min_A == max(min_A, min_B, min_C)) == 1)
      drop_B <- ((keep_B == 0) * (min_B == max(min_A, min_B, min_C)) == 1)
      drop_C <- ((keep_C == 0) * (min_C == max(min_A, min_B, min_C)) == 1)
      if (sum(drop_A, drop_B, drop_C) > 1) { # if more than one arm is dropped
        randompick <- sample(1:3, 1, prob = c(drop_A, drop_B, drop_C) / sum(drop_A, drop_B, drop_C))
        drop_A <- (randompick == 1)
        drop_B <- (randompick == 2)
        drop_C <- (randompick == 3)
      }
      if (drop_A == 1) {
        message("Step 2: Arm A's interim posterior probability of having the lowest response is higher\n")
        message("Arm A is dropped\n")
      } else if (drop_B == 1) {
        message("Step 2: Arm B's interim posterior probability of having the lowest response is higher\n")
        message("Arm B is dropped\n")
      } else if (drop_C == 1) {
        message("Step 2: Arm C's interim posterior probability of having the lowest response is higher\n")
        message("Arm B is dropped\n")
      }
    } else {
      message("Step 1: No treatment has P_{m,l} bigger than threshold ")
      message(drop_threshold_large)
      message("\n")
      drop_A <- (min_A > drop_threshold_small)
      if (drop_A == 1) {
        message("Step 2: Arm A's interim posterior probability of having the lowest response is higher than threshold ")
        message(drop_threshold_small)
        message("\n")
        message("Arm A is dropped\n")
      }
      drop_B <- (min_B > drop_threshold_small)
      if (drop_B == 1) {
        message("Step 2: Arm B's interim posterior probability of having the lowest response is higher than threshold ")
        message(drop_threshold_small)
        message("\n")
        message("Arm B is dropped\n")
      }
      drop_C <- (min_C > drop_threshold_small)
      if (drop_A == 1) {
        message("Step 2: Arm B's interim posterior probability of having the lowest response is higher than threshold ")
        message(drop_threshold_small)
        message("\n")
        message("Arm B is dropped\n")
      }
    }

    dropped_arm <- (drop_A == 1) * 1 + (drop_B == 1) * 2 + (drop_C == 1) * 3
    if (all(c(drop_A, drop_B, drop_C) == 0)) { # if none of the arm is dropped, move on to next update
      message("none of the arm is removed, move on to next update\n")
    }

    message("\n")
    result <- list("dropped_arm" = dropped_arm, "posterior_sample" = posterior_sample)
    class(result) <- "group_seq"
    return(result)
  } else {
    mydata <- data
    mydata$disc <- 2 * mydata$trt.1st - (mydata$resp.1st == 0)

    bugfile <- readLines(Bayes_file)
    bugfile <- gsub(pattern = "pi_prior_dist", replacement = pi_prior_dist, x = bugfile)
    bugfile <- gsub(pattern = "beta0_prior_dist", replacement = beta0_prior_dist, x = bugfile)
    bugfile2 <- gsub(pattern = "beta1_prior_dist", replacement = beta1_prior_dist, x = bugfile)

    bugfile2_file <- tempfile(fileext = ".bug")
    writeLines(bugfile2, con = bugfile2_file)
    error_ind <- 0
    tryCatch(
      {
        jags <- do.call(rjags::jags.model, c(list(
          file = bugfile2_file,
          data = list(
            n = nrow(mydata),
            num_arms = NUM_ARMS,
            Y1 = mydata$resp.1st,
            Y2 = mydata$resp.2nd,
            treatment_stageI = mydata$trt.1st,
            treatment_stageII = mydata$trt.2nd,
            response_stageI_disc = mydata$disc,
            # prior
            pi_prior.a = pi_prior[c(1, 3, 5)],
            pi_prior.b = pi_prior[c(2, 4, 6)],
            beta0_prior.a = beta0_prior.a,
            beta0_prior.b = beta0_prior.b,
            beta1_prior.a = beta1_prior.a,
            beta1_prior.c = beta1_prior.c
          ),
          n.chains = n_MCMC_chain, n.adapt = BURN.IN, quiet = quiet, jags.model_options
        )))
        update(jags, BURN.IN, progress.bar = progress.bar)
        posterior_sample <- do.call(
          rjags::coda.samples,
          c(list(
            model = jags,
            variable.names = c("pi", "beta"),
            thin = thin,
            n.iter = MCMC_SAMPLE,
            progress.bar = progress.bar,
            coda.samples_options
          ))
        )
      },
      warning = function(war) {
        warning_count <- warning_count + 1
        err_war_message <- rbind(paste("The warning ", warning_count, " is: ", war))
      }
    )
    out_post <- as.data.frame(posterior_sample[[1]])
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
    class(result) <- "group_seq"
    return(result)
  }
}



#' Summarizing BJSM fits
#'
#' `summary` method for class "`group_seq`"
#'
#' @param object an object of class "`group_seq`", usually, a result of a call to \code{\link{group_seq}}
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
summary.group_seq <- function(object, ...) {
  if (length(object) != 2) {
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

    if (!is.null(object$pi_DTR_est)) {
      dtreff <- cbind(object$pi_DTR_est, object$pi_DTR_se, rbind(object$ci_pi_AB, object$ci_pi_AC, object$ci_pi_BA, object$ci_pi_BC, object$ci_pi_CA, object$ci_pi_CB))
      colnames(dtreff) <- c("Estimate", "Std. Error", "C.I.", "CI low", "CI high")
      obj <- list(trteff = trteff, trtdiff = trtdiff, betaest = betaest, dtreff = dtreff)
    }

    obj <- list(trteff = trteff, trtdiff = trtdiff, betaest = betaest)
  } else {
    obj <- list(dropped_arm = object$dropped_arm, MCMC_result = summary(object$posterior_sample))
  }

  class(obj) <- "summary.group_seq"
  obj
}


#' @rdname group_seq
#' @param x object to print
#' @param ... further arguments. Not currently used.
#' @export
#' @export print.summary.group_seq
print.summary.group_seq <- function(x, ...) {
  if (length(x) != 2) {
    cat("\nTreatment Effects Estimate:\n")
    print(x$trteff)
    cat("\nDifferences between Treatments:\n")
    print(x$trtdiff)
    cat("\nLinkage Parameter Estimate:\n")
    print(x$betaest)
    if (!is.null(x$dtreff)) {
      cat("\nExpected Response Rate of Dynamic Treatment Regimens (DTR):\n")
      print(x$dtreff)
    }
    cat("\n")
  } else {
    if (x$dropped_arm == 0) { # if none of the arm is dropped, move on to next update
      cat("none of the arm is removed, move on to next update\n")
    } else if (x$dropped_arm == 1) {
      cat("Arm A is dropped\n")
    } else if (x$dropped_arm == 2) {
      cat("Arm B is dropped\n")
    } else {
      cat("Arm C is dropped\n")
    }

    print(x$MCMC_result)
  }
}


#' @rdname group_seq
#' @param x object to summarize.
#' @param ... further arguments. Not currently used.
#' @export
#' @export print.group_seq
print.group_seq <- function(x, ...) {
  if (length(x) != 2) {
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
  } else {
    if (x$dropped_arm == 0) { # if none of the arm is dropped, move on to next update
      cat("none of the arm is removed, move on to next update\n")
    } else if (x$dropped_arm == 1) {
      cat("Arm A is dropped\n")
    } else if (x$dropped_arm == 2) {
      cat("Arm B is dropped\n")
    } else {
      cat("Arm C is dropped\n")
    }
  }
}
