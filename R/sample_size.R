#' Sample size calculation for snSMART with 3 active treatments and a binary outcome
#'
#' conduct Bayesian sample size calculation for a snSMART design with 3 active
#' treatments and a binary outcome to distinguish the best treatment from the second-best
#' treatment using the Bayesian joint stage model.
#'
#' @param pi a vector with 3 values (`piA`, `piB`, `piC`). `piA` is the the response
#'  rate (ranges from 0.01 to 0.99) for treatment A, `piB` is the response rate
#'  (ranges from 0.01 to 0.99) for treatment B, `piC` is the response rate (ranges
#'  from 0.01 to 0.99) for treatment C
#' @param beta1 the linkage parameter (ranges from 1.00 to 1/largest response rate)
#'  for first stage responders. (A smaller value leads to more conservative sample
#'  size calculation because two stages are less correlated)
#' @param beta0 the linkage parameter (ranges from 0.01 to 0.99) for first stage
#'  non-responders. A larger value leads to a more conservative sample size calculation
#'  because two stages are less correlated
#' @param coverage the coverage rate (ranges from 0.01 to 0.99) for the posterior
#'  difference of top two treatments
#' @param power the probability (ranges from 0.01 to 0.99) for identify the best treatment
#' @param mu a vector with 3 values (`muA`, `muB`, `muC`). `muA` is the prior mean
#'  (ranges from 0.01 to 0.99) for treatment A, `muB` is the prior mean (ranges from
#'  0.01 to 0.99) for treatment B, `muC` is the prior mean (ranges from 0.01 to 0.99)
#'  for treatment C
#' @param n a vector with 3 values (`nA`, `nB`, `nC`). `nA` is the prior sample size
#'  (larger than 0) for treatment A. `nB` is the prior sample size (larger than 0)
#'  for treatment B. `nC` is the prior sample size (larger than 0) for treatment C
#' @param test for testing purposes only. Defaults to `FALSE`.
#' @param verbose TRUE or FALSE. If FALSE, no function message and progress bar will be
#'  printed.
#'
#' @details
#' Note that this package does not include the JAGS library, users need to install JAGS separately. Please check this page for more details: \url{https://sourceforge.net/projects/mcmc-jags/}
#' Please load the \code{EnvStats} package before calculating sample size.
#' This function may take a few minutes to run
#'
#' @return
#' \describe{
#'     \item{final_N}{the estimated sample size per arm for this snSMART}
#'     \item{critical_value}{ critical value based on the provided coverage value}
#'     \item{grid_result}{for each iteration we calculate \code{l}, where \code{l} belongs to \code{{2 * (pi_(1) - pi_(2)), ..., 0.02, 0.01}}; \code{E(D)}: the mean of the posterior distribution of \code{D},
#'     , where \code{D = pi_(1) = pi_(2)}; \code{Var(D)}: the variance of the posterior distribution of \code{D}; \code{N}: the corresponding
#'     sample size; and \code{power}: the resulting power of this iteration}
#'
#' }
#' @examples
#' require(EnvStats)
#'
#' sampleSize <- sample_size(
#'   pi = c(0.7, 0.5, 0.25), beta1 = 1.4, beta0 = 0.5, coverage = 0.9,
#'   power = 0.8, mu = c(0.65, 0.55, 0.25), n = c(4, 2, 3), test = TRUE
#' )
#'
#' # sampleSize = sample_size(pi = c(0.7, 0.5, 0.25), beta1 = 1.4, beta0 = 0.5, coverage = 0.9,
#' #    power = 0.8, mu = c(0.65, 0.55, 0.25), n = c(4, 2, 3), test = FALSE)
#'
#' @references
#' Wei, B., Braun, T.M., Tamura, R.N. and Kidwell, K.M., 2018. A Bayesian analysis of
#' small n sequential multiple assignment randomized trials (snSMARTs).
#' Statistics in medicine, 37(26), pp.3723-3732.
#'
#' Wei, B., Braun, T.M., Tamura, R.N. and Kidwell, K., 2020. Sample size determination
#' for Bayesian analysis of small n sequential, multiple assignment, randomized trials
#' (snSMARTs) with three agents. Journal of Biopharmaceutical Statistics, 30(6), pp.1109-1120.
#'
#' @seealso
# #' \code{\link{JSRM_binary}} \cr
#' \code{\link{BJSM_binary}}
#'
#' @importFrom EnvStats qpareto ppareto
#'
#' @export
#'


sample_size <- function(pi, beta1, beta0, coverage, power, mu, n, test = FALSE, verbose = FALSE) {
  if (test == TRUE) {
    return()
  }

  beta_prior_generator <- function(info_level, prior_mean) {
    alpha0 <- prior_mean * info_level
    beta0 <- info_level * (1 - prior_mean)
    alpha0_beta0 <- list(
      alpha0 = alpha0,
      beta0 = beta0
    )
  }

  # cFunc[a_, b_, beta1_] = (beta1 a^2)/(a + b)
  #
  cFunc <- function(a, b, beta1) {
    c <- (beta1 * a^2) / (a + b)
    return(c)
  }
  # cFunc(1,2,3)
  #
  # dFunc[a_, b_, beta0_] = (a b - beta1 a^2)/(a + b)
  #
  dFunc <- function(a, b, beta1) {
    d <- (a^2 + a * b - beta1 * a^2) / (a + b)
    return(d)
  }
  # dFunc(1,2,3)
  #
  # eFunc[a_, b_, beta1_, K_] = (beta0 a b)/((a + b) (K - 1))
  #
  eFunc <- function(a, b, beta0, K) {
    e <- (beta0 * a * b) / ((a + b) * (K - 1))
    return(e)
  }
  # eFunc(1,2,3,4)
  #
  # fFunc[a_, b_, beta1_, K_] = (beta0 a b + b^2)/((a + b) (K - 1))
  fFunc <- function(a, b, beta0, K) {
    f <- (beta0 * a * b + b^2) / ((a + b) * (K - 1))
    return(f)
  }

  sigmaSqABC <- function(K,
                         piA, piB, piC,
                         beta1, beta0,
                         aA, bA, cA, dA, eA, fA,
                         n) {
    sigmaSqABC <-
      1 / ((1 / 1 / beta0^2 * ((eA + ((n - piB * n) / (K - 1) * beta0 * piA + (n - piC * n) / (K - 1) * beta0 * piA)) * ((2 * n - piB * n - piC * n) / 2 - ((n - piB * n) / (K - 1) * beta0 * piA + (n - piC * n) / (K - 1) * beta0 * piA) +
        fA)) / ((eA + fA + (2 * n - piB * n - piC * n) / 2)^2 * (eA + fA + (2 * n - piB * n - piC * n) / 2 +
        1))) +
        1 / (((aA + piA * n) * (bA + n - piA * n)) / ((aA + bA + n)^2 * (aA + bA + n + 1))) +
        1 / (1 / beta1^2 * ((cA + piA * n * piA * beta1) * (dA + piA * n - piA * n * piA * beta1)) / ((cA + dA + piA * n)^2 * (cA + dA + piA * n + 1))))

    return(sigmaSqABC)
  }

  normal_cdf <- function(x, mu, sigmaSq) {
    normal_cdf <- 1 / 2 * (1 + pracma::erf((x - mu) / (sqrt(sigmaSq * 2))))
    return(normal_cdf)
  }

  # normal_cdf(0,1,Inf)

  order1_pdf <- function(mu1, sigmaSq1, mu2, sigmaSq2, mu3, sigmaSq3, x) {
    order1_pdf <- dnorm(x, mean = mu1, sd = sqrt(sigmaSq1)) * normal_cdf(x, mu2, sigmaSq2) * normal_cdf(x, mu3, sigmaSq3) +
      dnorm(x, mean = mu2, sd = sqrt(sigmaSq2)) * normal_cdf(x, mu1, sigmaSq1) * normal_cdf(x, mu3, sigmaSq3) +
      dnorm(x, mean = mu3, sd = sqrt(sigmaSq3)) * normal_cdf(x, mu1, sigmaSq1) * normal_cdf(x, mu2, sigmaSq2)
    return(order1_pdf)
  }

  order2_pdf <- function(mu1, sigmaSq1, mu2, sigmaSq2, mu3, sigmaSq3, x) {
    order2_pdf <- (1 - normal_cdf(x, mu3, sigmaSq3)) * (normal_cdf(x, mu1, sigmaSq1) * dnorm(x, mean = mu2, sd = sqrt(sigmaSq2)) + normal_cdf(x, mu2, sigmaSq2) * dnorm(x, mean = mu1, sd = sqrt(sigmaSq1))) +
      (1 - normal_cdf(x, mu2, sigmaSq2)) * (normal_cdf(x, mu1, sigmaSq1) * dnorm(x, mean = mu3, sd = sqrt(sigmaSq3)) + normal_cdf(x, mu3, sigmaSq3) * dnorm(x, mean = mu1, sd = sqrt(sigmaSq1))) +
      (1 - normal_cdf(x, mu1, sigmaSq1)) * (normal_cdf(x, mu2, sigmaSq2) * dnorm(x, mean = mu3, sd = sqrt(sigmaSq3)) + normal_cdf(x, mu3, sigmaSq3) * dnorm(x, mean = mu2, sd = sqrt(sigmaSq2)))
    return(order2_pdf)
  }

  prod_o1d <- function(mu1, sigmaSq1, mu2, sigmaSq2, mu3, sigmaSq3, x, d) {
    prod_o1d <- order1_pdf(mu1, sigmaSq1, mu2, sigmaSq2, mu3, sigmaSq3, x + d) * d
    return(prod_o1d)
  }

  prod_o1dd <- function(mu1, sigmaSq1, mu2, sigmaSq2, mu3, sigmaSq3, x, d) {
    prod_o1dd <- order1_pdf(mu1, sigmaSq1, mu2, sigmaSq2, mu3, sigmaSq3, x + d) * d * d
    return(prod_o1dd)
  }

  o1d_integral <- function(mu1, sigmaSq1, mu2, sigmaSq2, mu3, sigmaSq3, x) {
    o1d_integral <- cubature::hcubature(prod_o1d,
      mu1 = mu1, sigmaSq1 = sigmaSq1,
      mu2 = mu2, sigmaSq2 = sigmaSq2, mu3 = mu3, sigmaSq3 = sigmaSq3, x = x, lowerLimit = 0, upperLimit = Inf
    )$integral
    return(o1d_integral)
  }

  o1dd_integral <- function(mu1, sigmaSq1, mu2, sigmaSq2, mu3, sigmaSq3, x) {
    o1dd_integral <- cubature::hcubature(prod_o1dd,
      mu1 = mu1, sigmaSq1 = sigmaSq1,
      mu2 = mu2, sigmaSq2 = sigmaSq2, mu3 = mu3, sigmaSq3 = sigmaSq3, x = x, lowerLimit = 0, upperLimit = Inf
    )$integral
    return(o1dd_integral)
  }

  diff_o1o2_mean_integrand <- function(mu1, sigmaSq1, mu2, sigmaSq2, mu3, sigmaSq3, x) {
    diff_o1o2_mean_integrand <- order2_pdf(
      mu1 = mu1, sigmaSq1 = sigmaSq1,
      mu2 = mu2, sigmaSq2 = sigmaSq2,
      mu3 = mu3, sigmaSq3 = sigmaSq3, x = x
    ) * o1d_integral(
      mu1 = mu1,
      sigmaSq1 = sigmaSq1,
      mu2 = mu2, sigmaSq2 = sigmaSq2,
      mu3 = mu3, sigmaSq3 = sigmaSq3, x = x
    )
    return(diff_o1o2_mean_integrand)
  }

  diff_o1o2_variance1_integrand <- function(mu1, sigmaSq1, mu2, sigmaSq2, mu3, sigmaSq3, x) {
    diff_o1o2_variance1_integrand <- order2_pdf(
      mu1 = mu1, sigmaSq1 = sigmaSq1,
      mu2 = mu2, sigmaSq2 = sigmaSq2,
      mu3 = mu3, sigmaSq3 = sigmaSq3, x = x
    ) * o1dd_integral(
      mu1 = mu1, sigmaSq1 = sigmaSq1,
      mu2 = mu2, sigmaSq2 = sigmaSq2,
      mu3 = mu3, sigmaSq3 = sigmaSq3, x = x
    )
    return(diff_o1o2_variance1_integrand)
  }

  mean_o1o2_diff <- function(K,
                             piA, piB, piC,
                             beta1, beta0,
                             aA, bA, cA, dA, eA, fA,
                             aB, bB, cB, dB, eB, fB,
                             aC, bC, cC, dC, eC, fC,
                             n) {
    sigmaSq1 <- sigmaSqABC(
      K,
      piA, piB, piC,
      beta1, beta0,
      aA, bA, cA, dA, eA, fA,
      n
    )
    sigmaSq2 <- sigmaSqABC(
      K,
      piB, piA, piC,
      beta1, beta0,
      aB, bB, cB, dB, eB, fB,
      n
    )
    sigmaSq3 <- sigmaSqABC(
      K,
      piC, piA, piB,
      beta1, beta0,
      aC, bC, cC, dC, eC, fC,
      n
    )
    mean_o1o2_diff <- cubature::hcubature(diff_o1o2_mean_integrand,
      mu1 = piA, sigmaSq1 = sigmaSq1,
      mu2 = piB, sigmaSq2 = sigmaSq2, mu3 = piC, sigmaSq3 = sigmaSq3, lowerLimit = -Inf, upperLimit = Inf
    )$integral

    return(mean_o1o2_diff)
  }

  var1_o1o2_diff <- function(K,
                             piA, piB, piC,
                             beta1, beta0,
                             aA, bA, cA, dA, eA, fA,
                             aB, bB, cB, dB, eB, fB,
                             aC, bC, cC, dC, eC, fC,
                             n) {
    sigmaSq1 <- sigmaSqABC(
      K,
      piA, piB, piC,
      beta1, beta0,
      aA, bA, cA, dA, eA, fA,
      n
    )
    sigmaSq2 <- sigmaSqABC(
      K,
      piB, piA, piC,
      beta1, beta0,
      aB, bB, cB, dB, eB, fB,
      n
    )
    sigmaSq3 <- sigmaSqABC(
      K,
      piC, piA, piB,
      beta1, beta0,
      aC, bC, cC, dC, eC, fC,
      n
    )
    # o1_o2_mean=mean_o1o2_diff(K,
    #                           piA, piB, piC,
    #                           beta1, beta0,
    #                           aA, bA, cA, dA, eA, fA,
    #                           aB, bB, cB, dB, eB, fB,
    #                           aC, bC, cC, dC, eC, fC,
    #                           n)

    var1_o1o2_diff <- cubature::hcubature(diff_o1o2_variance1_integrand,
      mu1 = piA, sigmaSq1 = sigmaSq1,
      mu2 = piB, sigmaSq2 = sigmaSq2, mu3 = piC, sigmaSq3 = sigmaSq3, lowerLimit = -Inf, upperLimit = Inf
    )$integral


    return(var1_o1o2_diff)
  }

  # sample_size_equation <- function (K,
  #                                   piA, piB, piC,
  #                                   beta1, beta0,
  #                                   aA, bA, cA, dA, eA, fA,
  #                                   aB, bB, cB, dB, eB, fB,
  #                                   aC, bC, cC, dC, eC, fC,
  #                                   ciL,
  #                                   n) {
  #   2*ciZ*sqrt(var1_o1o2_diff(K,
  #                              piA, piB, piC,
  #                              beta1, beta0,
  #                              aA, bA, cA, dA, eA, fA,
  #                              aB, bB, cB, dB, eB, fB,
  #                              aC, bC, cC, dC, eC, fC,
  #                              n))-ciL
  # }





  piA <- pi[1]
  piB <- pi[2]
  piC <- pi[3]

  sortT <- sort(c(piA, piB, piC))
  if (sortT[-1][1] == sortT[-1][2]) {
    stop("Top 2 treatments have the same expected response rate, only a unique best treatment is allowed")
  }

  muA <- mu[1]
  muB <- mu[2]
  muC <- mu[3]

  nA <- n[1]
  nB <- n[2]
  nC <- n[3]


  max_beta1 <- 1 / min(c(max(piA, piB, piC), 0.99))
  str_beta1 <- paste0("The value of ", "\u03B21", " should be within 1~", max_beta1)


  if (!piA < 1 & piA > 0) {
    message("The value of \u03C0A should be within 0.01~0.99")
  }
  if (!piB < 1 & piB > 0) {
    message("The value of \u03C0B should be within 0.01~0.99")
  }
  if (!piC < 1 & piC > 0) {
    message("The value of \u03C0C should be within 0.01~0.99")
  }
  if (!(beta1 <= max_beta1 & beta1 > 1) |
    is.na(max_beta1) |
    piA > 1 | piA < 0 |
    piB > 1 | piB < 0 |
    piC > 1 | piC < 0) {
    message(str_beta1)
  }
  if (!beta0 < 1 & beta0 > 0) {
    message("The value of \u03B20 should be within 0.01~0.99")
  }
  if (!coverage < 1 & coverage > 0) {
    message("The value of 1-\u03B1 should be within 0.01~0.99")
  }
  if (!power < 1 & power > 0) {
    message("The value of 1-\u03BE should be within 0.01~0.99")
  }

  if (!muA < 1 & muA > 0) {
    message("The value of \u03bcA should be within 0.01~0.99")
  }
  if (!muB < 1 & muB > 0) {
    message("The value of \u03bcB should be within 0.01~0.99")
  }
  if (!muC < 1 & muC > 0) {
    message("The value of \u03bcC should be within 0.01~0.99")
  }
  if (!nA > 0) {
    message("The value of nA should be greater than 0")
  }
  if (!nB > 0) {
    message("The value of nB should be greater than 0")
  }
  if (!nC > 0) {
    message("The value of nC should be greater than 0")
  }


  K <- 3
  magic_Z <- 1.5
  pi_A <- as.numeric(piA)
  pi_B <- as.numeric(piB)
  pi_C <- as.numeric(piC)

  muA_input <- as.numeric(muA)
  muB_input <- as.numeric(muB)
  muC_input <- as.numeric(muC)

  nA_input <- as.numeric(nA)
  nB_input <- as.numeric(nB)
  nC_input <- as.numeric(nC)

  LIST_OF_PIS <- c(pi_A, pi_B, pi_C)
  LIST_OF_PIS <- LIST_OF_PIS[order(LIST_OF_PIS, decreasing = T)]

  COVRAGE <- as.numeric(coverage)
  # COVRAGE_BONFF=1-(1-COVRAGE)/(K-1)
  # COVRAGE=COVRAGE_BONFF
  POW <- as.numeric(power) # when POW is small, we need to decrease sample size lower limit, need to estimate the lower limit or there will not be a opposite sign

  # SAMPLE_SIZE_LLIMIT=1
  # SAMPLE_SIZE_ULIMIT=1000

  # CIL_MIN=0.1# can not be too small
  # CIL_MAX=2
  # CIL_STEP=0.01

  SS_LOW <- 1
  SS_HIGH <- 300
  CONVERGE_TOL <- 0.001
  CIL_MIN <- 0.01 # can not be too small
  CIL_STEP_I <- 0.01

  ###################
  # check if pis are all the same
  # see if the first two treatment are the same-if not the same then do computation, if the same, do computation with the first and last
  # need all treatment arm for computation
  if (LIST_OF_PIS[1] == LIST_OF_PIS[2]) {
    piA <- LIST_OF_PIS[1]
    piB <- LIST_OF_PIS[K]
    piC <- LIST_OF_PIS[setdiff(1:K, c(1, K))]
  } else {
    piA <- LIST_OF_PIS[1]
    piB <- LIST_OF_PIS[2]
    piC <- LIST_OF_PIS[setdiff(1:K, c(1, 2))]
  }

  # generate truncated beta1 and beta0
  # set.seed(199)
  # generate pareto truncated at 1.5 with location 1 and scale 3, mean 1.5
  beta1_mean <- as.numeric(beta1)
  pareto_beta <- 1 / (1 - 1 / beta1_mean)
  # beta1_sample=rtrunc(LINKAGE_SAMPLE, spec="pareto",a=1,b=1/max(c(piA,piB,piC)),1,pareto_beta)
  # beta1_sample=rep(1,LINKAGE_SAMPLE)
  beta0_mean <- as.numeric(beta0)
  # beta0_sample=rbeta(LINKAGE_SAMPLE,beta0_mean*2,2-beta0_mean*2)
  # beta0_sample=rep(1,LINKAGE_SAMPLE)
  # plot(hist(beta0))
  # load functions
  # require(rmutil)
  CIL_MAX <- (piA - piB) * magic_Z
  # require(EnvStats)
  beta1_sample <- round(mean(truncdist::rtrunc(99999, spec = "pareto", a = 1, b = 1 / max(c(piA, piB, piC)), 1, pareto_beta)), 3)
  beta0_sample <- rep(beta0_mean, 1)

  sample_size_list_pair1 <- NULL

  error_count_pair1 <- 0
  warn_count_pair1 <- 0

  sim_count <- 0

  error_round_pair1 <- NULL
  warn_round_pair1 <- NULL
  error_mesg_pair1 <- NULL
  warn_mesg_pair1 <- NULL

  error_round_pair1 <- NULL
  warn_round_pair1 <- NULL
  error_mesg_pair1 <- NULL
  warn_mesg_pair1 <- NULL

  # calculate priors

  beta1 <- beta1_sample
  beta0 <- beta0_sample

  aA <- muA_input * nA_input
  bA <- nA_input - aA
  # aA = 0.25*2
  # bA = 2-aA
  # source(system.file("functions.R", package = "snSMART"))

  cA <- cFunc(aA, bA, beta1)
  dA <- dFunc(aA, bA, beta1)
  eA <- eFunc(aA, bA, beta0, K)
  fA <- fFunc(aA, bA, beta0, K)

  aB <- muB_input * nB_input
  bB <- nB_input - aB
  # aB = 0.25 * 2
  # bB = 2-aB
  cB <- cFunc(aB, bB, beta1)
  dB <- dFunc(aB, bB, beta1)
  eB <- eFunc(aB, bB, beta0, K)
  fB <- fFunc(aB, bB, beta0, K)

  aC <- muC_input * nC_input
  bC <- nC_input - aC
  # aC = 0.25 * 2
  # bC = 2-aC
  cC <- cFunc(aC, bC, beta1)
  dC <- dFunc(aC, bC, beta1)
  eC <- eFunc(aC, bC, beta0, K)
  fC <- fFunc(aC, bC, beta0, K)

  ciZ <- qnorm(1 - (1 - COVRAGE) / 2) # critical value of coverage

  # i=0
  if (verbose == TRUE) {
    pb <- txtProgressBar(min = 0, max = 1, initial = 0, style = 3)
    setTxtProgressBar(pb, 0.3)
  }

  grid_result <- NULL

  for (CIL_I in seq(CIL_MAX, CIL_MIN, by = -CIL_STEP_I)) {
    # i=i+1
    # message(CIL_I)
    # get sample size solved
    # message(CIL_I)
    # CIL_I=0.26

    if (CIL_I / 2 == (max(c(piA, piB, piC)) - min(c(piA, piB, piC)))) {
      next
    }
    tryCatch(
      {
        fun <- function(x) {
          2 * ciZ * sqrt(var1_o1o2_diff(
            K = K,
            piA = piA, piB = piB, piC = piC,
            beta1 = beta1, beta0 = beta0,
            aA = aA, bA = bA, cA = cA, dA = dA, eA = eA, fA = fA,
            aB = aB, bB = bB, cB = cB, dB = dB, eB = eB, fB = fB,
            aC = aC, bC = bC, cC = cC, dC = dC, eC = eC, fC = fC,
            n = x
          ) - mean_o1o2_diff(
            K = K,
            piA = piA, piB = piB, piC = piC,
            beta1 = beta1, beta0 = beta0,
            aA = aA, bA = bA, cA = cA, dA = dA, eA = eA, fA = fA,
            aB = aB, bB = bB, cB = cB, dB = dB, eB = eB, fB = fB,
            aC = aC, bC = bC, cC = cC, dC = dC, eC = eC, fC = fC,
            n = x
          )^2) - CIL_I
        }

        # for(i in 1:300){
        #   message(fun(i))
        # }

        # ciL=0.28

        # sample_size_tmp_pair1=ceiling(uniroot(fun, c(SS_LOW, SS_HIGH),tol = CONVERGE_TOL)$root)
        sample_size_tmp_pair1 <- ceiling(uniroot(fun, c(SS_LOW, SS_HIGH), tol = CONVERGE_TOL)$root)
        # All <- uniroot.all(sample_size_equation, c(0, 8))
        # calculate powerpower
        # 1-pnorm((ciL/2-(max(c(piA,piB)-min(piA,piB))))/sqrt(sigmaSqABDiffFunc(sample_size_tmp_pair1)))+pnorm((-ciL/2-(max(c(piA,piB)-min(piA,piB))))/sqrt(sigmaSqABDiffFunc(sample_size_tmp_pair1)))
        mu_o1o2_diff <- mean_o1o2_diff(
          K = K,
          piA = piA, piB = piB, piC = piC,
          beta1 = beta1, beta0 = beta0,
          aA = aA, bA = bA, cA = cA, dA = dA, eA = eA, fA = fA,
          aB = aB, bB = bB, cB = cB, dB = dB, eB = eB, fB = fB,
          aC = aC, bC = bC, cC = cC, dC = dC, eC = eC, fC = fC,
          n = sample_size_tmp_pair1
        )
        mu_o1o2_sq_diff <- var1_o1o2_diff(
          K = K,
          piA = piA, piB = piB, piC = piC,
          beta1 = beta1, beta0 = beta0,
          aA = aA, bA = bA, cA = cA, dA = dA, eA = eA, fA = fA,
          aB = aB, bB = bB, cB = cB, dB = dB, eB = eB, fB = fB,
          aC = aC, bC = bC, cC = cC, dC = dC, eC = eC, fC = fC,
          n = sample_size_tmp_pair1
        )
        var_o1o2_diff <- mu_o1o2_sq_diff - mu_o1o2_diff^2
        # pow_pair1 = (1-pnorm((ciL/2-mu_o1o2_diff)/sqrt(var_o1o2_diff)))
        pow_pair1 <- (1 - pnorm((CIL_I / 2 - mu_o1o2_diff) / sqrt(var_o1o2_diff)))

        # total_process=length(seq(CIL_MAX,CIL_MIN,by=-CIL_STEP_I))

        if (verbose == TRUE) {
          setTxtProgressBar(pb, pow_pair1 / POW)
        }

        grid_result <- rbind(grid_result, c(CIL_I, mu_o1o2_diff, mu_o1o2_sq_diff, pow_pair1, sample_size_tmp_pair1))


        if (pow_pair1 > POW) {
          if (verbose == TRUE) {
            setTxtProgressBar(pb, 1)
          }
          break
        }
      },
      error = function(c) {
        # error_round_tmp_pair1=cbind(i,beta1,beta0,CIL_I)
        # error_round_pair1=rbind(error_round_pair1,error_round_tmp_pair1)
        # next
      },
      warning = function(c) {
        # warn_round_tmp_pair1=cbind(i,beta1,beta0,CIL_I)
        # warn_round_pair1=rbind(warn_round_pair1,warn_round_tmp_pair1)
        # message(i)
        # message(CIL_I)
        # next
      },
      finally = { # posterior_sample_burn=window(posterior_sample,start=BURNING, end=MCMC_SAMPLE)
        # posterior_sample_cmb=do.call(rbind, posterior_sample_burn)
      }
    )
    # message(CIL_I)
    # message(pow_pair1)
    # message(mu_o1o2_diff)
    # message(mu_o1o2_sq_diff)
    # message(sample_size_tmp_pair1)
  }

  if (verbose == TRUE){
    close(pb)
  }

  colnames(grid_result) <- c("l", "E(D)", "Var(D)", "power", "N")
  result <- list(critical_value = ciZ, grid_result = as.data.frame(grid_result), final_N = sample_size_tmp_pair1)

  message(paste0(
    "With given settings, the estimated sample size per arm for an snSMART is: ", sample_size_tmp_pair1, "\n",
    "This implies that for an snSMART with sample size of ", sample_size_tmp_pair1, " per arm (", 3 * sample_size_tmp_pair1, " in total for three agents):", "\n",
    "The probability of successfully identifying the best treatment is ", power, " when the difference of response rates between the best and second best treatment is at least ", LIST_OF_PIS[1] - LIST_OF_PIS[2], ", and the response rate of the best treatment is ", LIST_OF_PIS[1], "\n"
  ))

  class(result) <- "sample_size"
  return(result)
}


#' @rdname sample_size
#' @param object object to summarize.
#' @param ... further arguments. Not currently used.
#' @export
#'
summary.sample_size <- function(object, ...) {
  obj <- list(final_N = object$final_N, grid_result = object$grid_result)
  class(obj) <- "summary.sample_size"
  obj
}

#' @rdname sample_size
#' @param x object to print
#' @param ... further arguments. Not currently used.
#' @export
#' @export print.summary.sample_size
#'
print.summary.sample_size <- function(x, ...) {
  cat("With given settings, the estimated sample size per arm for an snSMART is: ")
  cat(as.numeric(x$final_N))
  cat(paste0("In total ", nrow(x$grid_result), " iterations were taken:\n"))
  message(x$grid_result)
}

#' @rdname sample_size
#' @param x object to print
#' @param ... further arguments. Not currently used.
#' @export
#' @export print.sample_size
#'
print.sample_size <- function(x, ...) {
  cat("With given settings, the estimated sample size per arm for an snSMART is:\n")
  cat(as.numeric(x$final_N))
}
