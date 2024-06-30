# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html

library(testthat)
library(snSMART)

test_that("BJSM_binary 1", {
  mydata <- data_binary

  BJSM_result <- BJSM_binary(
    data = mydata, prior_dist = c("beta", "beta", "pareto"),
    pi_prior = c(0.4, 1.6, 0.4, 1.6, 0.4, 1.6), beta_prior = c(1.6, 0.4, 3, 1),
    n_MCMC_chain = 1, n.adapt = 1000, MCMC_SAMPLE = 2000, ci = 0.95,
    six = TRUE, DTR = TRUE, verbose = FALSE
  )

  summary(BJSM_result)
  print(summary(BJSM_result))
  print(BJSM_result)

  result = c(0.3889, 0.4371, 0.5793)
  names(result) = c("pi_A", "pi_B", "pi_C")

  expect_equal(BJSM_result$pi_hat_bjsm, result, tolerance = 1e-1)
})

test_that("BJSM_binary 2", {
  mydata <- data_binary

  BJSM_result <- BJSM_binary(
    data = mydata, prior_dist = c("beta", "beta", "pareto"),
    pi_prior = c(0.4, 1.6, 0.4, 1.6, 0.4, 1.6), beta_prior = c(1.6, 0.4, 3, 1),
    n_MCMC_chain = 1, n.adapt = 1000, MCMC_SAMPLE = 2000, ci = 0.95,
    six = TRUE, DTR = FALSE, verbose = FALSE
  )

  summary(BJSM_result)
  print(summary(BJSM_result))
  print(BJSM_result)

  result = c(0.3889, 0.4371, 0.5793)
  names(result) = c("pi_A", "pi_B", "pi_C")

  expect_equal(BJSM_result$pi_hat_bjsm, result, tolerance = 1e-1)
})

test_that("BJSM_binary 3", {
  mydata <- data_binary

  BJSM_result <- BJSM_binary(
    data = mydata, prior_dist = c("beta", "beta", "pareto"),
    pi_prior = c(0.4, 1.6, 0.4, 1.6, 0.4, 1.6), beta_prior = c(1.6, 0.4, 3, 1),
    n_MCMC_chain = 1, n.adapt = 1000, MCMC_SAMPLE = 2000, ci = 0.95,
    six = TRUE, DTR = TRUE, verbose = TRUE
  )

  summary(BJSM_result)
  print(summary(BJSM_result))
  print(BJSM_result)

  result = c(0.3889, 0.4371, 0.5793)
  names(result) = c("pi_A", "pi_B", "pi_C")

  expect_equal(BJSM_result$pi_hat_bjsm, result, tolerance = 1e-1)
})

test_that("BJSM_binary 4", {
  mydata <- data_binary

  BJSM_result <- BJSM_binary(
    data = mydata, prior_dist = c("beta", "beta", "pareto"),
    pi_prior = c(0.4, 1.6, 0.4, 1.6, 0.4, 1.6), beta_prior = c(1.6, 0.4, 3, 1),
    n_MCMC_chain = 1, n.adapt = 1000, MCMC_SAMPLE = 2000, ci = 0.95,
    six = TRUE, DTR = FALSE, verbose = FALSE
  )

  summary(BJSM_result)
  print(summary(BJSM_result))
  print(BJSM_result)

  result = c(0.3889, 0.4371, 0.5793)
  names(result) = c("pi_A", "pi_B", "pi_C")

  expect_equal(BJSM_result$pi_hat_bjsm, result, tolerance = 1e-1)
})

test_that("BJSM_binary 5", {
  mydata <- data_binary

  BJSM_result2 <- BJSM_binary(
    data = mydata, prior_dist = c("beta", "beta", "pareto"),
    pi_prior = c(0.4, 1.6, 0.4, 1.6, 0.4, 1.6), beta_prior = c(1.6, 0.4, 3, 1),
    n_MCMC_chain = 1, n.adapt = 10000, MCMC_SAMPLE = 60000, ci = 0.95,
    six = FALSE, DTR = FALSE, verbose = FALSE
  )

  summary(BJSM_result2)
  print(BJSM_result2)
  print(summary(BJSM_result2))

  result = c(0.3993, 0.4250, 0.5411)
  names(result) = c("pi_A", "pi_B", "pi_C")

  expect_equal(BJSM_result2$pi_hat_bjsm, result, tolerance = 1e-1)
})

test_that("BJSM_binary 6", {
  mydata <- data_binary

  BJSM_result2 <- BJSM_binary(
    data = mydata, prior_dist = c("beta", "beta", "pareto"),
    pi_prior = c(0.4, 1.6, 0.4, 1.6, 0.4, 1.6), beta_prior = c(1.6, 0.4, 3, 1),
    n_MCMC_chain = 1, n.adapt = 10000, MCMC_SAMPLE = 60000, ci = 0.95,
    six = FALSE, DTR = TRUE, verbose = FALSE
  )

  summary(BJSM_result2)
  print(BJSM_result2)
  print(summary(BJSM_result2))

  result = c(0.3993, 0.4250, 0.5411)
  names(result) = c("pi_A", "pi_B", "pi_C")

  expect_equal(BJSM_result2$pi_hat_bjsm, result, tolerance = 1e-1)
})

test_that("BJSM_binary 7", {
  mydata <- data_binary

  BJSM_result2 <- BJSM_binary(
    data = mydata, prior_dist = c("beta", "beta", "pareto"),
    pi_prior = c(0.4, 1.6, 0.4, 1.6, 0.4, 1.6), beta_prior = c(1.6, 0.4, 3, 1),
    n_MCMC_chain = 1, n.adapt = 10000, MCMC_SAMPLE = 60000, ci = 0.95,
    six = FALSE, DTR = TRUE, verbose = TRUE
  )

  summary(BJSM_result2)
  print(BJSM_result2)
  print(summary(BJSM_result2))

  result = c(0.3993, 0.4250, 0.5411)
  names(result) = c("pi_A", "pi_B", "pi_C")

  expect_equal(BJSM_result2$pi_hat_bjsm, result, tolerance = 1e-1)
})


test_that("BJSM_binary 7", {
  mydata <- data_binary

  BJSM_result2 <- BJSM_binary(
    data = mydata, prior_dist = c("beta", "beta", "pareto"),
    pi_prior = c(0.4, 1.6, 0.4, 1.6, 0.4, 1.6), beta_prior = c(1.6, 0.4, 3, 1),
    n_MCMC_chain = 1, n.adapt = 10000, MCMC_SAMPLE = 60000, ci = 0.95,
    six = FALSE, DTR = FALSE, verbose = TRUE
  )

  summary(BJSM_result2)
  print(BJSM_result2)
  print(summary(BJSM_result2))

  result = c(0.3993, 0.4250, 0.5411)
  names(result) = c("pi_A", "pi_B", "pi_C")

  expect_equal(BJSM_result2$pi_hat_bjsm, result, tolerance = 1e-1)
})


test_that("BJSM_binary 8", {

  data <- data_dose
  BJSM_dose_result <- BJSM_binary(
    data = data_dose, prior_dist = c("beta", "gamma"),
    pi_prior = c(3, 17), normal.par = c(0.2, 100), beta_prior = c(2, 2),
    n_MCMC_chain = 2, n.adapt = 1000, MCMC_SAMPLE = 6000, ci = 0.95, verbose = FALSE
  )

  summary(BJSM_dose_result)
  print(BJSM_dose_result)
  print(summary(BJSM_dose_result))

  result = c(0.06971, 0.40131, 0.73859)
  names(result) = c("pi_P", "pi_L", "pi_H")

  expect_equal(BJSM_dose_result$pi_hat_bjsm, result, tolerance = 1e-1)
})


test_that("BJSM_binary 9", {

  data <- data_dose
  BJSM_dose_result <- BJSM_binary(
    data = data_dose, prior_dist = c("beta", "gamma"),
    pi_prior = c(3, 17), normal.par = c(0.2, 100), beta_prior = c(2, 2),
    n_MCMC_chain = 2, n.adapt = 1000, MCMC_SAMPLE = 6000, ci = 0.95, verbose = TRUE
  )

  summary(BJSM_dose_result)
  print(BJSM_dose_result)
  print(summary(BJSM_dose_result))

  result = c(0.06971, 0.40131, 0.73859)
  names(result) = c("pi_P", "pi_L", "pi_H")

  expect_equal(BJSM_dose_result$pi_hat_bjsm, result, tolerance = 1e-1)
})


test_that("BJSM_c 1", {

  trialData <- trialDataMF

  BJSM_result <- BJSM_c(
    data = trialData, xi_prior.mean = c(50, 50, 50),
    xi_prior.sd = c(50, 50, 50), phi3_prior.sd = 20, n_MCMC_chain = 1,
    n.adapt = 1000, MCMC_SAMPLE = 5000, BURIN.IN = 1000, ci = 0.95, n.digits = 5, verbose = FALSE
  )

  summary(BJSM_result)
  print(BJSM_result)
  print(summary(BJSM_result))

  result = 51.12
  names(result) = c("xi_[1]")

  expect_equal(BJSM_result$mean_estimate[c("xi_[1]")], result, tolerance = 1e-1)
})

test_that("BJSM_c 2", {

  trialData <- trialDataMF

  BJSM_result <- BJSM_c(
    data = trialData, xi_prior.mean = c(50, 50, 50),
    xi_prior.sd = c(50, 50, 50), phi3_prior.sd = 20, n_MCMC_chain = 1,
    n.adapt = 1000, MCMC_SAMPLE = 5000, BURIN.IN = 1000, ci = 0.95, n.digits = 5, verbose = TRUE
  )

  summary(BJSM_result)
  print(BJSM_result)
  print(summary(BJSM_result))

  result = 51.12
  names(result) = c("xi_[1]")

  expect_equal(BJSM_result$mean_estimate[c("xi_[1]")], result, tolerance = 1e-1)
})

test_that("group_seq 1", {

  mydata <- groupseqDATA_look1

  result1 <- group_seq(
    data = mydata, interim = TRUE, drop_threshold_pair = c(0.5, 0.4),
    prior_dist = c("beta", "beta", "pareto"), pi_prior = c(0.4, 1.6, 0.4, 1.6, 0.4, 1.6),
    beta_prior = c(1.6, 0.4, 3, 1), MCMC_SAMPLE = 6000, n.adapt = 1000, n_MCMC_chain = 1
  )

  summary(result1)
  print(result1)
  print(summary(result1))

  expect_equal(result1$dropped_arm, 1, tolerance = 1e-1)
})

test_that("group_seq 2", {

  mydata <- groupseqDATA_look1

  result1 <- group_seq(
    data = mydata, interim = TRUE, drop_threshold_pair = c(0.5, 0.4),
    prior_dist = c("beta", "beta", "pareto"), pi_prior = c(0.4, 1.6, 0.4, 1.6, 0.4, 1.6),
    beta_prior = c(1.6, 0.4, 3, 1), MCMC_SAMPLE = 6000, n.adapt = 1000, n_MCMC_chain = 1,
    DTR = TRUE, verbose = TRUE
  )

  summary(result1)
  print(result1)
  print(summary(result1))

  expect_equal(result1$dropped_arm, 1, tolerance = 1e-1)
})

test_that("group_seq 3", {

  mydata <- groupseqDATA_look1

  result1 <- group_seq(
    data = mydata, interim = TRUE, drop_threshold_pair = c(0.5, 0.4),
    prior_dist = c("beta", "beta", "pareto"), pi_prior = c(0.4, 1.6, 0.4, 1.6, 0.4, 1.6),
    beta_prior = c(1.6, 0.4, 3, 1), MCMC_SAMPLE = 6000, n.adapt = 1000, n_MCMC_chain = 1,
    DTR = FALSE, verbose = TRUE
  )

  summary(result1)
  print(result1)
  print(summary(result1))

  expect_equal(result1$dropped_arm, 1, tolerance = 1e-1)
})


test_that("group_seq 4", {

  mydata <- groupseqDATA_look1

  result1 <- group_seq(
    data = mydata, interim = TRUE, drop_threshold_pair = c(0.5, 0.4),
    prior_dist = c("beta", "beta", "pareto"), pi_prior = c(0.4, 1.6, 0.4, 1.6, 0.4, 1.6),
    beta_prior = c(1.6, 0.4, 3, 1), MCMC_SAMPLE = 6000, n.adapt = 1000, n_MCMC_chain = 1,
    DTR = FALSE, verbose = FALSE
  )

  summary(result1)
  print(result1)
  print(summary(result1))

  expect_equal(result1$dropped_arm, 1, tolerance = 1e-1)
})

test_that("group_seq 5", {

  mydata <- groupseqDATA_full
  result2 <- group_seq(
    data = mydata, interim = FALSE, prior_dist = c("beta", "beta", "pareto"),
    pi_prior = c(0.4, 1.6, 0.4, 1.6, 0.4, 1.6),
    beta_prior = c(1.6, 0.4, 3, 1), MCMC_SAMPLE = 6000, n.adapt = 1000,
    n_MCMC_chain = 1, ci = 0.95, DTR = TRUE
  )

  summary(result2)
  print(result2)
  print(summary(result2))

  result = c(0.3014, 0.4729, 0.6761)
  names(result) = c("pi_A", "pi_B", "pi_C")

  expect_equal(result2$pi_hat_bjsm, result, tolerance = 1e-1)
})

test_that("group_seq 6", {

  mydata <- groupseqDATA_full
  result2 <- group_seq(
    data = mydata, interim = FALSE, prior_dist = c("beta", "beta", "pareto"),
    pi_prior = c(0.4, 1.6, 0.4, 1.6, 0.4, 1.6),
    beta_prior = c(1.6, 0.4, 3, 1), MCMC_SAMPLE = 6000, n.adapt = 1000,
    n_MCMC_chain = 1, ci = 0.95, DTR = TRUE, verbose = TRUE
  )

  summary(result2)
  print(result2)
  print(summary(result2))

  result = c(0.3014, 0.4729, 0.6761)
  names(result) = c("pi_A", "pi_B", "pi_C")

  expect_equal(result2$pi_hat_bjsm, result, tolerance = 1e-1)
})

test_that("group_seq 7", {

  mydata <- groupseqDATA_full
  result2 <- group_seq(
    data = mydata, interim = FALSE, prior_dist = c("beta", "beta", "pareto"),
    pi_prior = c(0.4, 1.6, 0.4, 1.6, 0.4, 1.6),
    beta_prior = c(1.6, 0.4, 3, 1), MCMC_SAMPLE = 6000, n.adapt = 1000,
    n_MCMC_chain = 1, ci = 0.95, DTR = FALSE, verbose = TRUE
  )

  summary(result2)
  print(result2)
  print(summary(result2))

  result = c(0.3014, 0.4729, 0.6761)
  names(result) = c("pi_A", "pi_B", "pi_C")

  expect_equal(result2$pi_hat_bjsm, result, tolerance = 1e-1)
})


test_that("group_seq 8", {

  mydata <- groupseqDATA_full
  result2 <- group_seq(
    data = mydata, interim = FALSE, prior_dist = c("beta", "beta", "pareto"),
    pi_prior = c(0.4, 1.6, 0.4, 1.6, 0.4, 1.6),
    beta_prior = c(1.6, 0.4, 3, 1), MCMC_SAMPLE = 6000, n.adapt = 1000,
    n_MCMC_chain = 1, ci = 0.95, DTR = FALSE, verbose = FALSE
  )

  summary(result2)
  print(result2)
  print(summary(result2))

  result = c(0.3014, 0.4729, 0.6761)
  names(result) = c("pi_A", "pi_B", "pi_C")

  expect_equal(result2$pi_hat_bjsm, result, tolerance = 1e-1)
})

test_that("group_seq 9", {

  mydata <- groupseqDATA_look1
  mydata$trt.1st = ifelse(mydata$trt.1st == 1, 4, mydata$trt.1st)
  mydata$trt.1st = ifelse(mydata$trt.1st == 2, 1, mydata$trt.1st)
  mydata$trt.1st = ifelse(mydata$trt.1st == 4, 2, mydata$trt.1st)

  mydata$trt.2nd = ifelse(mydata$trt.2nd == 1, 4, mydata$trt.2nd)
  mydata$trt.2nd = ifelse(mydata$trt.2nd == 2, 1, mydata$trt.2nd)
  mydata$trt.2nd = ifelse(mydata$trt.2nd == 4, 2, mydata$trt.2nd)

  result1 <- group_seq(
    data = mydata, interim = TRUE, drop_threshold_pair = c(0.5, 0.4),
    prior_dist = c("beta", "beta", "pareto"), pi_prior = c(0.4, 1.6, 0.4, 1.6, 0.4, 1.6),
    beta_prior = c(1.6, 0.4, 3, 1), MCMC_SAMPLE = 6000, n.adapt = 1000, n_MCMC_chain = 1
  )

  summary(result1)
  print(result1)
  print(summary(result1))

  expect_equal(result1$dropped_arm, 2, tolerance = 1e-1)
})

test_that("group_seq 10", {

  mydata <- groupseqDATA_look1
  mydata$trt.1st = ifelse(mydata$trt.1st == 1, 4, mydata$trt.1st)
  mydata$trt.1st = ifelse(mydata$trt.1st == 3, 1, mydata$trt.1st)
  mydata$trt.1st = ifelse(mydata$trt.1st == 4, 3, mydata$trt.1st)

  result1 <- group_seq(
    data = mydata, interim = TRUE, drop_threshold_pair = c(0.5, 0.4),
    prior_dist = c("beta", "beta", "pareto"), pi_prior = c(0.4, 1.6, 0.4, 1.6, 0.4, 1.6),
    beta_prior = c(1.6, 0.4, 3, 1), MCMC_SAMPLE = 6000, n.adapt = 1000, n_MCMC_chain = 1
  )

  summary(result1)
  print(result1)
  print(summary(result1))

  expect_equal(result1$dropped_arm, 3, tolerance = 1e-1)
})


test_that("group_seq 11", {

  mydata <- groupseqDATA_full
  mydata$trt.1st = ifelse(mydata$trt.1st == 1, 4, mydata$trt.1st)
  mydata$trt.1st = ifelse(mydata$trt.1st == 3, 1, mydata$trt.1st)
  mydata$trt.1st = ifelse(mydata$trt.1st == 4, 3, mydata$trt.1st)

  mydata$trt.2nd = ifelse(mydata$trt.2nd == 1, 4, mydata$trt.2nd)
  mydata$trt.2nd = ifelse(mydata$trt.2nd == 3, 1, mydata$trt.2nd)
  mydata$trt.2nd = ifelse(mydata$trt.2nd == 4, 3, mydata$trt.2nd)

  result2 <- group_seq(
    data = mydata, interim = FALSE, prior_dist = c("beta", "beta", "pareto"),
    pi_prior = c(0.4, 1.6, 0.4, 1.6, 0.4, 1.6),
    beta_prior = c(1.6, 0.4, 3, 1), MCMC_SAMPLE = 6000, n.adapt = 1000,
    n_MCMC_chain = 1, ci = 0.95, DTR = FALSE, verbose = FALSE
  )

  summary(result2)
  print(result2)
  print(summary(result2))

  result = c(0.6878912, 0.4730337, 0.3022561)
  names(result) = c("pi_A", "pi_B", "pi_C")

  expect_equal(result2$pi_hat_bjsm, result, tolerance = 1e-1)
})


test_that("group_seq 12", {

  mydata <- groupseqDATA_full
  mydata$trt.1st = ifelse(mydata$trt.1st == 1, 4, mydata$trt.1st)
  mydata$trt.1st = ifelse(mydata$trt.1st == 2, 1, mydata$trt.1st)
  mydata$trt.1st = ifelse(mydata$trt.1st == 4, 2, mydata$trt.1st)

  mydata$trt.2nd = ifelse(mydata$trt.2nd == 1, 4, mydata$trt.2nd)
  mydata$trt.2nd = ifelse(mydata$trt.2nd == 2, 1, mydata$trt.2nd)
  mydata$trt.2nd = ifelse(mydata$trt.2nd == 4, 2, mydata$trt.2nd)

  result2 <- group_seq(
    data = mydata, interim = FALSE, prior_dist = c("beta", "beta", "pareto"),
    pi_prior = c(0.4, 1.6, 0.4, 1.6, 0.4, 1.6),
    beta_prior = c(1.6, 0.4, 3, 1), MCMC_SAMPLE = 6000, n.adapt = 1000,
    n_MCMC_chain = 1, ci = 0.95, DTR = FALSE, verbose = FALSE
  )

  summary(result2)
  print(result2)
  print(summary(result2))

  result = c(0.4727878, 0.3003522, 0.6828518)
  names(result) = c("pi_A", "pi_B", "pi_C")

  expect_equal(result2$pi_hat_bjsm, result, tolerance = 1e-1)
})

test_that("LPJSM_binary 3", {

  data <- data_binary

  LPJSM_result <- LPJSM_binary(data = data, six = TRUE, DTR = TRUE)

  summary(LPJSM_result)
  print(LPJSM_result)
  print(summary(LPJSM_result))

  result = c(0.2966, 0.3736, 0.4298)
  names(result) = c("alphaA", "alphaB", "alphaC")

  expect_equal(LPJSM_result$pi_hat, result, tolerance = 1e-1)
})

test_that("LPJSM_binary 4", {

  data <- data_binary

  LPJSM_result <- LPJSM_binary(data = data, six = FALSE, DTR = TRUE)

  summary(LPJSM_result)
  print(LPJSM_result)
  print(summary(LPJSM_result))

  result = c(0.2966, 0.3736, 0.4298)
  names(result) = c("alphaA", "alphaB", "alphaC")

  expect_equal(LPJSM_result$pi_hat, result, tolerance = 1e-1)
})

test_that("LPJSM_binary 5", {

  data <- data_binary

  LPJSM_result <- LPJSM_binary(data = data, six = FALSE, DTR = FALSE)

  summary(LPJSM_result)
  print(LPJSM_result)
  print(summary(LPJSM_result))

  result = c(0.2966, 0.3736, 0.4298)
  names(result) = c("alphaA", "alphaB", "alphaC")

  expect_equal(LPJSM_result$pi_hat, result, tolerance = 1e-1)
})

test_that("LPJSM_binary 5", {

  data <- data_binary

  LPJSM_result <- LPJSM_binary(data = data, six = TRUE, DTR = FALSE)

  result = c(0.2966, 0.3736, 0.4298)
  names(result) = c("alphaA", "alphaB", "alphaC")

  summary(LPJSM_result)
  print(LPJSM_result)
  print(summary(LPJSM_result))

  expect_equal(LPJSM_result$pi_hat, result, tolerance = 1e-1)
})

test_that("sampleSize 1", {
  sampleSize <- sample_size(
    pi = c(0.7, 0.5, 0.25), beta1 = 1.4, beta0 = 0.5, coverage = 0.9,
    power = 0.3, mu = c(0.65, 0.55, 0.25), n = c(10, 10, 10)
  )
  result = 17

  summary(sampleSize)
  print(sampleSize)
  print(summary(sampleSize))

  expect_equal(sampleSize$final_N, result, tolerance = 1e-1)
})

test_that("sampleSize 2", {
  sampleSize <- sample_size(
    pi = c(0.7, 0.5, 0.25), beta1 = 1.4, beta0 = 0.5, coverage = 0.9,
    power = 0.3, mu = c(0.65, 0.55, 0.25), n = c(10, 10, 10), verbose = TRUE
  )
  result = 17

  summary(sampleSize)
  print(sampleSize)
  print(summary(sampleSize))

  expect_equal(sampleSize$final_N, result, tolerance = 1e-1)
})


test_that("sampleSize 2", {
  try({sampleSize <- sample_size(
    pi = c(2, 2, 2), beta1 = -2, beta0 = -2, coverage = 2,
    power = -2, mu = c(2, 2, 2), n = c(-2, -2, -2), verbose = TRUE
  )}, silent = TRUE)
})
