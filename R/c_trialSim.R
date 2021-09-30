#' Trial Characteristics Assessment based on Simulation (continuous snSMART design)
#'
#' assess the characteristics of trial design based on datasets simulated by \code{\link{c_trialDataGen}}
#'
#' @param stage1effects a vector of 3 values (`mean effect of treatment A`, `mean effect of treatment B`, `main effect of treatment C`). e.g. `stage1effects = c(40, 50, 70)`
#' @param stage2weights a vector of 3 values (`weight of stage I`, `weight of stage II`, `weight of staying on the same treatment`). The second stage mean outcome is a weighted average of the treatment effects from stage one and stage two with an additional effect if the patient stays on the same treatment. e.g. `stage2weights = c(0.2, 0.7, 0.3)`
#' @param stagecorrs a vector of 2 values (`correlations between two stages if patient stays on the same treatment`, `correlations between two stages if patent doesn't stay on the same treatment`). This allows those who stay on the same treatmetn to have a different correlation between stage one and stage two outcomes than those who switch treatments.
#' @param variance the variance of the treatment effects
#' @param n vector of number of patients on each treatment, e.g. \code{n <- c(100, 100, 100)}
#' @param nsim number of simulations. Use the \code{\link{c_trialDataGen}} function to generate \code{nsim} simulated datasets
#' @param stay.ethical numerical value, if stage 1 outcome is bigger than `stay.ethical` value, patient has probability of 1 of staying on the same treatment in stage 2
#' @param switch.safety numerical value, if stage 1 outcome is smaller than `switch.safety` value, patient has probability of 0 of staying on the same treatment in stage 2

#' @return
#' \describe{
#' \item{switch.safety}{input by user or estimated based on data if `switch.safety` = NULL}
#' \item{stay.ethical}{input by user or estimated based on data if `stay.ethical` = NULL}
#' \item{stage1Eff}{average simulated stage 1 treatment effect }
#' \item{stage2EFF}{average simulated stage 2 treatment effect }
#' \item{stage2best}{percentage of patients assigned to the best treatment in stage 2; overallbest: percentage of patients assigned to the best treatment at least once in the trial}
#' \item{stage2worst}{percentage of patients assigned to the worst treatment in stage 2; overallworst: percentage of patients assigned to the worst treatment at least once in the trial}
#' \item{improve}{percentage of patients received an improved treatment outcome in stage 2}
#' \item{switchBetterOrStayedBest}{percentage of patients received more effective treatments in stage 2 or stayed on the best treatment throughout the trial}
#' \item{stayedSame}{percentage of patients received the same treatment in stage 2 after randomization}
#' \item{gotBetter}{percentage of patients that received more effective treatments in stage 2}
#' \item{gotWorse}{percentage of patients that received less effective treatments in stage 2}
#' }
#'
#' @examples
#' treat.a<-c(70, 15, c(0,0,0))
#' treat.b<-c(50, 10, c(0,0,0))
#' treat.c<-c(60, 12, c(0,0,0))
#' treatInfo = rbind(treat.a, treat.b, treat.c)
#'
#' treatCors<-diag(.9, 3)
#'
#' switch.safety <- NULL
#' stay.ethical<- NULL

#' na<-100
#' nb<-100
#' nc<-100
#' n<-c(na,nb,nc)

#' k = c_trialSim(stage1effects, stage2weights, stagecorrs, variance, n, 100)
#'
#' @references
#' Hartman, H., Tamura, R.N., Schipper, M.J. and Kidwell, K.M., 2021. Design and analysis considerations for utilizing a mapping function in a small sample,
#' sequential, multiple assignment, randomized trials with continuous outcomes. Statistics in Medicine, 40(2), pp.312-326.
#'
#'
#' @seealso
#' \code{\link{BJSM_c}} \cr
#' \code{\link{c_trialDataGen}}
#'
#' @export




c_trialSim = function(stage1effects, stage2weights, stagecorrs, variance, stay.ethical = NULL, switch.safety = NULL, n, nsim){
  x<-c()

  treatInfo = rbind(c(stage1effects[1], rep(0, 4)), c(stage1effects[2], rep(0, 4)), c(stage1effects[3], rep(0, 4)))

    for(i in 1:nsim){
    trialData = c_trialDataGen(stage1effects = stage1effects, stage2weights = stage2weights, stagecorrs = stagecorrs, variance = variance, n = n, stay.ethical = stay.ethical, switch.safety = switch.safety, wideForm = T)
    stage1EFF = mean(trialData$stage1outcome)
    stage2EFF = mean(trialData$stage2outcome)

    stage2best = sum(trialData$trt2 == which.max(treatInfo[,1]))/sum(n)

    overallbest = sum(trialData$trt2 == which.max(treatInfo[,1]) | trialData$trt1 == which.max(treatInfo[,1]))/sum(n)

    stage2worst = sum(trialData$trt2 == which.min(treatInfo[,1]))/sum(n)

    overallworst = sum(trialData$trt2 == which.min(treatInfo[,1]) | trialData$trt1 == which.min(treatInfo[,1]))/sum(n)

    improve = sum(trialData$stage1outcome < trialData$stage2outcome)/sum(n)

    switchBetterOrStayedBest = sum((trialData$trt1==which.min(treatInfo[,1]) & trialData$trt2!=which.min(treatInfo[,1])) | trialData$trt2 == which.max (treatInfo[,1]))/sum(n)
    stayedSame = sum(trialData$trt1!=which.max(treatInfo[,1]) & trialData$trt1==trialData$trt2)/sum(n)

    gotBetter =  sum((trialData$trt1==which.min(treatInfo[,1]) & trialData$trt2!=which.min(treatInfo[,1])) | (trialData$trt2 != which.max(treatInfo[,1]) & trialData$trt2 == which.max(treatInfo[,1])))/sum(n)
    gotWorse = sum((trialData$trt1==which.max(treatInfo[,1]) & trialData$trt2!=which.max(treatInfo[,1])) | (trialData$trt2 != which.min(treatInfo[,1]) & trialData$trt2 == which.min(treatInfo[,1])))/sum(n)

    x = rbind(x, cbind(stage1EFF, stage2EFF, stage2best, overallbest, stage2worst, overallworst, improve, switchBetterOrStayedBest, stayedSame, gotBetter, gotWorse))
  }
  z = x
  summary = round(c(switch.safety, stay.ethical, apply(data.frame(x),2,mean)), 2)
  names(summary)[1:2] = c("switch.safety", "stay.ethical")
  return(summary)
}
