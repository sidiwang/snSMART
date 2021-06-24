#' JSRM
#'
#' A joint-stage regression model (JSRM) is a frequentist modeling approach that incorporates the responses of both stages as repeated measurements for each subject.
#' Generalized estimating equations (GEE) are used to estimate the response rates of each treatment. The marginal response rates for each DTR can also be obtained based on the GEE results/
#'
#' @param data data format produced by the `trial_dataset_dose` function
#'
#' @return a `list` containing
#' \itemize{
#'   \item{GEE_output}{ - original output of the GEE (geeglm) model}
#'   \item{pi_hat}{ - estimate of response rate/treatment effect}
#'   \item{sd_pi_hat}{ - standard error of the response rate}
#' }
#'
#' @examples
#' #data
#' mydata = trial_dataset_dose(trtP_I = 30, trtL_I = 30, trtH_I = 30, respP_I = 5, respL_I = 10, respH_I = 15,
#' trtPL_II = 3, trtPH_II = 2, trtLL_II = 5, trtLH_II = 5, trtHL_II = 8, trtHH_II = 7,
#' respPL_II = 1, respPH_II = 2, respLL_II = 2, respLH_II = 3, respHL_II = 4, respHH_II = 6,
#' trtNPL_II = 10, trtNPH_II = 15, trtNLH_II = 10, trtNLL_II = 10, trtNHH_II = 15,
#' respNPL_II = 7, respNPH_II = 8, respNLH_II = 6, respNLL_II = 7, respNHH_II = 10)
#'
#' #JSRM_result = JSRM_dose_binary(data = mydata)
#'
#' @references
#' Fang, F., Hochstedler, K.A., Tamura, R.N., Braun, T.M. and Kidwell, K.M., 2021. Bayesian methods to compare dose levels with placebo in a small n,
#' sequential, multiple assignment, randomized trial. Statistics in Medicine, 40(4), pp.963-977.
#'
JSRM_dose_binary = function(data){

  # data, same format as the bjsm_binary.R file trial dataset format
  # six, if TRUE, will run the six beta model, if FALSE will run the two beta model


  mydata = data
  mydata$disc <- 2 * mydata$treatment_stageI - (mydata$response_stageI == 0)

  Y <- c(mydata$response_stageI, mydata$response_stageII)
  XP1 <- as.numeric(mydata$treatment_stageI==1)
  XL1 <- as.numeric(mydata$treatment_stageI==2)
  XH1 <- as.numeric(mydata$treatment_stageI==3)
  XP2 <- as.numeric(mydata$treatment_stageII==1)
  XL2 <- as.numeric(mydata$treatment_stageII==2)
  XH2 <- as.numeric(mydata$treatment_stageII==3)
  XP <- c(XP1, XP2)
  XL <- c(XL1, XL2)
  XH <- c(XH1, XH2)
  X0P <- c(rep(0, nrow(mydata)), ifelse(mydata$treatment_stageI==1,1-mydata$response_stageI,0))
  X1P <- c(rep(0, nrow(mydata)), ifelse(mydata$treatment_stageI==1,mydata$response_stageI,0))
  X0L <- c(rep(0, nrow(mydata)), ifelse(mydata$treatment_stageI==2,1-mydata$response_stageI,0))
  X1L <- c(rep(0, nrow(mydata)), ifelse(mydata$treatment_stageI==2,mydata$response_stageI,0))
  X0H <- c(rep(0, nrow(mydata)), ifelse(mydata$treatment_stageI==3,1-mydata$response_stageI,0))
  X1H <- c(rep(0, nrow(mydata)), ifelse(mydata$treatment_stageI==3,mydata$response_stageI,0))
  ptid <- rep(1:nrow(mydata), 2)

  geedata <- data.frame(ptid, XP, XL, XH, X0P, X1P, X0L, X1L, X0H, X1H,  Y)
  geedata <- geedata[order(geedata$ptid),]
  rm(ptid, XP, XL, XH, X0P, X1P, X0L, X1L, X0H, X1H,  Y)
  try({
    mod1=NULL
    mod1 <- geepack::geeglm(Y~XP+XL+XH+X0P+X1P+X0L+X1L+X0H+X1H-1, family=poisson(link="log"), data=geedata, id=ptid,corstr = "independence")
    beta_hat <- mod1$coefficients[1:3];
    sd_beta_hat <- summary(mod1)$coef[1:3,2];
    pi_hat <- exp(beta_hat);
    sd_pi_hat <- exp(beta_hat)*sd_beta_hat;
    b_hat <- coef(mod1);

  })



  result = list("GEE_output" = mod1, # original output of the GEE (geeglm) model
                "pi_hat" = pi_hat, # estimate of response rate/treatment effect
                "sd_pi_hat" = sd_pi_hat) # standard error of the response rate



  return(result)
}
