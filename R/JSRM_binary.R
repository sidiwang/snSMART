#' JSRM for standard snSMART
#'
#' A joint-stage regression model (JSRM) is a frequentist modeling approach that incorporates the responses of both stages as repeated measurements for each subject.
#' Generalized estimating equations (GEE) are used to estimate the response rates of each treatment. The marginal response rates for each DTR can also be obtained based on the GEE results
#'
#' @param data data format produced by \code{\link{trial_dataset}} or \code{\link{data_simulation}}
#' @param six if TRUE, will run the six beta model, if FALSE will run the two
#' beta model. Default is `six = TRUE`
#' @param DTR if TRUE, will also return the expected response rate and its standard error of dynamic treatment regimens
#' @param digits the number of significant digits to use when printing
#'
#' @return a `list` containing
#' \itemize{
#'   \item{`GEE_output`}{ - original output of the GEE (geeglm) model}
#'   \item{`pi_hat`}{ - estimate of response rate/treatment effect}
#'   \item{`sd_pi_hat`}{ - standard error of the response rate}
#'   \item{`pi_DTR_hat`}{ - expected response rate of dynamic treatment regimens (DTRs)}
#'   \item{`pi_DTR_se`}{ - standard deviation of DTR estimates}
#' }
#'
#' @examples
#' data = trial_dataset(trt = c(9, 12, 9), resp = c(3, 3, 5), trt_same_II = c(2, 2, 4),
#'     resp_same_II = c(2, 2, 4), trt_negA = c(3, 2), trt_negB = c(3, 2), trt_negc = c(1, 2),
#'     resp_negA = c(3, 1), resp_negB = c(2, 2), resp_negC = c(0, 0))
#'
#' JSRM_result = JSRM_binary(data = data, six = TRUE, DTR = TRUE)
#'
#' summary(JSRM_result)
#'
#' @references
#' Wei, B., Braun, T.M., Tamura, R.N. and Kidwell, K.M., 2018. A Bayesian analysis of small n sequential multiple assignment randomized trials (snSMARTs).
#' Statistics in medicine, 37(26), pp.3723-3732.
#'
#' Chao, Y.C., Trachtman, H., Gipson, D.S., Spino, C., Braun, T.M. and Kidwell, K.M., 2020. Dynamic treatment regimens in small n, sequential, multiple assignment,
#' randomized trials: An application in focal segmental glomerulosclerosis. Contemporary clinical trials, 92, p.105989.
#'
#' @seealso
#' \code{\link{data_simulation}} \cr
#' \code{\link{trial_dataset}} \cr
#' \code{\link{BJSM_binary}} \cr
#' \code{\link{sample_size}}
#'
#' @rdname JSRM_binary
#' @export
#'
JSRM_binary = function(data, six = TRUE, DTR = TRUE){

  # data, same format as the bjsm_binary.R file trial dataset format
  # six, if TRUE, will run the six beta model, if FALSE will run the two beta model


  mydata = data
  mydata$disc <- 2 * mydata$treatment_stageI - (mydata$response_stageI == 0)

  Y <- c(mydata$response_stageI, mydata$response_stageII)
  XA1 <- as.numeric(mydata$treatment_stageI == 1)
  XB1 <- as.numeric(mydata$treatment_stageI == 2)
  XC1 <- as.numeric(mydata$treatment_stageI == 3)
  XA2 <- as.numeric(mydata$treatment_stageII == 1)
  XB2 <- as.numeric(mydata$treatment_stageII == 2)
  XC2 <- as.numeric(mydata$treatment_stageII == 3)
  XA <- c(XA1, XA2)
  XB <- c(XB1, XB2)
  XC <- c(XC1, XC2)
  Z1A <- c(rep(0, nrow(mydata)), ifelse(mydata$treatment_stageI == 1, mydata$response_stageI ,0))
  Z2A <- c(rep(0, nrow(mydata)), ifelse(mydata$treatment_stageI == 1, 1 - mydata$response_stageI ,0))
  Z1B <- c(rep(0, nrow(mydata)), ifelse(mydata$treatment_stageI == 2, mydata$response_stageI ,0))
  Z2B <- c(rep(0, nrow(mydata)), ifelse(mydata$treatment_stageI == 2, 1 - mydata$response_stageI ,0))
  Z1C <- c(rep(0, nrow(mydata)), ifelse(mydata$treatment_stageI == 3, mydata$response_stageI ,0))
  Z2C <- c(rep(0, nrow(mydata)), ifelse(mydata$treatment_stageI == 3, 1 - mydata$response_stageI ,0))
  ptid <- rep(1:nrow(mydata), 2)

  if (six == TRUE){
    geedata <- data.frame(ptid, XA, XB, XC, Z1A, Z2A, Z1B, Z2B, Z1C, Z2C, Y)
    geedata <- geedata[order(geedata$ptid),]
    rm(ptid, Y, Z1A, Z2A, Z1B, Z2B, Z1C, Z2C, XA, XB, XC)
    try({
      mod1 <- geepack::geeglm(Y ~ XA + XB + XC + Z1A + Z2A + Z1B + Z2B + Z1C + Z2C - 1, family = poisson(link = "log"), data = geedata, id = ptid,corstr = "independence")
      beta_hat <- mod1$coefficients[1:3];
      sd_beta_hat <- summary(mod1)$coef[1:3, 2];
      pi_hat <- exp(beta_hat);
      sd_pi_hat <- exp(beta_hat) * sd_beta_hat;
      b_hat <- coef(mod1);
      grad <- exp(rbind(c(2 * b_hat[1] + b_hat[4], b_hat[2] + b_hat[5], b_hat[1] + b_hat[2] + b_hat[5]),
                        c(2 * b_hat[1] + b_hat[4], b_hat[3] + b_hat[5], b_hat[1] + b_hat[3] + b_hat[5]),
                        c(2 * b_hat[2] + b_hat[6], b_hat[1] + b_hat[7], b_hat[2] + b_hat[1] + b_hat[7]),
                        c(2 * b_hat[2] + b_hat[6], b_hat[3] + b_hat[7], b_hat[2] + b_hat[3] + b_hat[7]),
                        c(2 * b_hat[3] + b_hat[8], b_hat[1] + b_hat[9], b_hat[3] + b_hat[1] + b_hat[9]),
                        c(2 * b_hat[3] + b_hat[8], b_hat[2] + b_hat[9], b_hat[3] + b_hat[2] + b_hat[9])));
      pi_DTR_hat <- c(1, 1, -1) %*% t(grad)
      grad[, 3] <- -grad[, 3]
      sigma_b <- mod1$geese$vbeta
      L1 <- matrix(c(2, 0, 0, 1, 0, 0, 0, 0, 0,
                     0, 1, 0, 0, 1, 0, 0, 0, 0,
                     1, 1, 0, 0, 1, 0, 0, 0, 0), nrow = 3, ncol = 9, byrow = T)
      sigma_g <- L1 %*% sigma_b %*% t(L1)
      seAB <- sqrt(grad[1,] %*% sigma_g %*% grad[1,])

      L2 <- matrix(c(2, 0, 0, 1, 0, 0, 0, 0, 0,
                     0, 0, 1, 0, 1, 0, 0, 0, 0,
                     1, 0, 1, 0, 1, 0, 0, 0, 0), nrow = 3, ncol = 9, byrow = T)
      sigma_g <- L2 %*% sigma_b %*% t(L2)
      seAC <- sqrt(grad[2,] %*% sigma_g %*% grad[2,])

      L3 <- matrix(c(0, 2, 0, 0, 0, 1, 0, 0, 0,
                     1, 0, 0, 0, 0, 0, 1, 0, 0,
                     1, 1, 0, 0, 0, 0, 1, 0, 0), nrow = 3, ncol = 9, byrow = T)
      sigma_g <- L3 %*% sigma_b %*% t(L3)
      seBA <- sqrt(grad[3,] %*% sigma_g %*% grad[3,])

      L4 <- matrix(c(0, 2, 0, 0, 0, 1, 0, 0, 0,
                     0, 0, 1, 0, 0, 0, 1, 0, 0,
                     0, 1, 1, 0, 0, 0, 1, 0, 0), nrow = 3, ncol = 9, byrow = T)
      sigma_g <- L4 %*% sigma_b %*% t(L4)
      seBC <- sqrt(grad[4,] %*% sigma_g %*% grad[4,])

      L5 <- matrix(c(0, 0, 2, 0, 0, 0, 0, 1, 0,
                     1, 0, 0, 0, 0, 0, 0, 0, 1,
                     1, 0, 1, 0, 0, 0, 0, 0, 1), nrow = 3, ncol = 9, byrow = T)
      sigma_g <- L5 %*% sigma_b %*% t(L5)
      seCA <- sqrt(grad[5,] %*% sigma_g %*% grad[5,])

      L6 <- matrix(c(0, 0, 2, 0, 0, 0, 0, 1, 0,
                     0, 1, 0, 0, 0, 0, 0, 0, 1,
                     0, 1, 1, 0, 0, 0, 0, 0, 1), nrow = 3, ncol = 9, byrow = T)
      sigma_g <- L6 %*% sigma_b %*% t(L6)
      seCB <- sqrt(grad[6,] %*% sigma_g %*% grad[6,])
      pi_DTR_se <- c(seAB, seAC, seBA, seBC, seCA, seCB)
    })

  } else {
    Z1 = Z1A + Z1B + Z1C
    Z2 = Z2A + Z2B + Z2C
    geedata <- data.frame(ptid, XA, XB, XC, Z1, Z2, Y)
    geedata <- geedata[order(geedata$ptid),]
    rm(ptid, Y, Z1A, Z2A, Z1B, Z2B, Z1C, Z2C, XA, XB, XC, Z1, Z2)
    try({
      mod1 <- geepack::geeglm(Y ~ XA + XB + XC + Z1 + Z2 - 1, family = poisson(link = "log"), data = geedata, id = ptid,corstr = "independence")
      beta_hat <- mod1$coefficients[1:3];
      sd_beta_hat <- summary(mod1)$coef[1:3,2];
      pi_hat <- exp(beta_hat);
      sd_pi_hat <- exp(beta_hat) * sd_beta_hat;
      b_hat <- coef(mod1);
      grad <- exp(rbind(c(2 * b_hat[1] + b_hat[4], b_hat[2] + b_hat[5], b_hat[1] + b_hat[2] + b_hat[5]),
                        c(2 * b_hat[1] + b_hat[4], b_hat[3] + b_hat[5], b_hat[1] + b_hat[3] + b_hat[5]),
                        c(2 * b_hat[2] + b_hat[4], b_hat[1] + b_hat[5], b_hat[2] + b_hat[1] + b_hat[5]),
                        c(2 * b_hat[2] + b_hat[4], b_hat[3] + b_hat[5], b_hat[2] + b_hat[3] + b_hat[5]),
                        c(2 * b_hat[3] + b_hat[4], b_hat[1] + b_hat[5], b_hat[3] + b_hat[1] + b_hat[5]),
                        c(2 * b_hat[3] + b_hat[4], b_hat[2] + b_hat[5], b_hat[3] + b_hat[2] + b_hat[5])));
      pi_DTR_hat <- c(1, 1, -1) %*% t(grad)
      grad[,3] <- -grad[,3]
      sigma_b <- mod1$geese$vbeta
      L1 <- matrix(c(2, 0, 0, 1, 0,
                     0, 1, 0, 0, 1,
                     1, 1, 0, 0, 1), nrow = 3, ncol = 5, byrow = T)
      sigma_g <- L1 %*% sigma_b %*% t(L1)
      seAB <- sqrt(grad[1,] %*% sigma_g %*% grad[1,])

      L2 <- matrix(c(2, 0, 0, 1, 0,
                     0, 0, 1, 0, 1,
                     1, 0, 1, 0, 1), nrow = 3, ncol = 5, byrow = T)
      sigma_g <- L2 %*% sigma_b %*% t(L2)
      seAC <- sqrt(grad[2,] %*% sigma_g %*% grad[2,])

      L3 <- matrix(c(0, 2, 0, 1, 0,
                     1, 0, 0, 0, 1,
                     1, 1, 0, 0, 1), nrow = 3, ncol = 5, byrow = T)
      sigma_g <- L3 %*% sigma_b %*% t(L3)
      seBA <- sqrt(grad[3,] %*% sigma_g %*% grad[3,])

      L4 <- matrix(c(0, 2, 0, 1, 0,
                     0, 0, 1, 0, 1,
                     0, 1, 1, 0, 1), nrow = 3, ncol = 5, byrow = T)
      sigma_g <- L4 %*% sigma_b %*% t(L4)
      seBC <- sqrt(grad[4,] %*% sigma_g %*% grad[4,])

      L5 <- matrix(c(0, 0, 2, 1, 0,
                     1, 0, 0, 0, 1,
                     1, 0, 1, 0, 1), nrow = 3, ncol = 5, byrow = T)
      sigma_g <- L5 %*% sigma_b %*% t(L5)
      seCA <- sqrt(grad[5,] %*% sigma_g %*% grad[5,])

      L6 <- matrix(c(0, 0, 2, 1, 0,
                     0, 1, 0, 0, 1,
                     0, 1, 1, 0, 1), nrow = 3, ncol = 5, byrow = T)
      sigma_g <- L6 %*% sigma_b %*% t(L6)
      seCB <- sqrt(grad[6,] %*% sigma_g %*% grad[6,])
      pi_DTR_se <- c(seAB, seAC, seBA, seBC, seCA, seCB)
    })
  }

  names(pi_DTR_hat) = c("rep_AB", "rep_AC", "rep_BA", "rep_BC", "rep_CA", "rep_CB")
  names(pi_DTR_se) = c("se_AB", "se_AC", "se_BA", "se_BC", "se_CA", "se_CB")

  if (DTR == TRUE){
    result = list("GEE_output" = mod1, # original output of the GEE (geeglm) model
                  "pi_hat" = pi_hat, # estimate of response rate/treatment effect
                  "sd_pi_hat" = sd_pi_hat, # standard error of the response rate
                  "pi_DTR_hat" = pi_DTR_hat, # expected response rate of dynamic treatment regimens (DTRs)
                  "pi_DTR_se" = pi_DTR_se) # standard deviation of DTR estimates
  }else{

    result = list("GEE_output" = mod1, # original output of the GEE (geeglm) model
                  "pi_hat" = pi_hat, # estimate of response rate/treatment effect
                  "sd_pi_hat" = sd_pi_hat) # standard error of the response rate


  }

  class(result) = "JSRM_binary"
  return(result)
}



#' @rdname JSRM_binary
#' @export
summary.JSRM_binary = function(object, digits = 5, ...){
  cat("\nGEE output:\n")
  print(summary(object$GEE_output))
  cat("\nTreatment Effect Estimate:\n")
  trteff = cbind(object$pi_hat, object$sd_pi_hat)
  rownames(trteff) = c("trtA", "trtB", "trtC")
  colnames(trteff) = c("Estimate", "Std. Error")
  print(trteff, digits = digits)

  if (!is.null(object$pi_DTR_hat) == TRUE){
    cat("\nExpected Response Rate of Dynamic Treatment Regimens (DTR):\n")
    DTRest = t(rbind(object$pi_DTR_hat, object$pi_DTR_se))
    colnames(DTRest) = c("Estimate", "Std. Error")
    rownames(DTRest) = c("rep_AB", "rep_AC", "rep_BA", "rep_BC", "rep_CA", "rep_CB")
    print(DTRest, digits = digits)
  }
  cat("\n")
}



#' @rdname JSRM_binary
#' @export
print.JSRM_binary = function(object, digits = 5, ...){
  cat("\nTreatment Effect Estimate\n")
  trteff = cbind(object$pi_hat, object$sd_pi_hat)
  rownames(trteff) = c("trtA", "trtB", "trtC")
  colnames(trteff) = c("Estimate", "Std. Error")
  print(trteff, digits = digits)
  cat("\n")
}


