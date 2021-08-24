#' Trial Characteristics Assessment based on Simulation (continuous snSMART design)
#'
#' assess the characteristics of trial design based on simulated data
#'
#' @param treatInfo Treatment distribution information. Each row represents one treatment, format: (mean of treatment effects, sd of treatment effects, c(priming effect on trt A,  trt B, trt C)), see example for detail
#' @param treatCors Treatment correlations matrix. Format example:
#' \[\[ (A,A) (A,B) (A,C)
#'    (B,A) (B,B) (B,C)
#'    (C,A) (C,B) (C,C) \]\], (A,B) denotes the treatment correlation between first stage treatment A and second stage treatment B
#' @param n vector of number of patients on each treatment e.g. n <- c(100, 100, 100)
#' @param nsim number of simulations
#' @param stay.ethical numerical value, if stage 1 outcome is bigger than `stay.ethical` value, patient has probability of 1 of staying on the same treatment in stage 2
#' @param switch.safety numerical value, if stage 1 outcome is smaller than `switch.safety` value, patient has probability of 0 of staying on the same treatment in stage 2

#' @return switch.safety: input by user; stay.ethical: input by user; stage1Eff: average simulated stage 1 treatment effect; stage2EFF: average simulated stage 2 treatment effect
#' stage2best: percentage of patients assigned to the best treatment in stage 2; overallbest: percentage of patients assigned to the best treatment at least once in the trial;
#' stage2worst: percentage of patients assigned to the worst treatment in stage 2; overallworst: percentage of patients assigned to the worst treatment at least once in the trial;
#' improve: percentage of patients received an improved treatment outcome in stage 2;
#' switchBetterOrStayedBest: percentage of patients received more effective treatments in stage 2 or stayed on the best treatment throughout the trial;
#' stayedSame: percentage of patients received the same treatment in stage 2 after randomization;
#' gotBetter: percentage of patients that received more effective treatments in stage 2;
#' gotWorse: percentage of patients that received less effective treatments in stage 2;
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

#' k = c_trialSim(treatInfo, treatCors, stay.ethical, switch.safety, n, 100)
#'
#' @references
#' Hartman, H., Tamura, R.N., Schipper, M.J. and Kidwell, K.M., 2021. Design and analysis considerations for utilizing a mapping function in a small sample,
#' sequential, multiple assignment, randomized trials with continuous outcomes. Statistics in Medicine, 40(2), pp.312-326.
#'




c_trialSim = function(treatInfo, treatCors, stay.ethical = NULL, switch.safety=NULL, n, nsim){
  x<-c()
    for(i in 1:nsim){
    trialData = c_trialDataGen(treatInfo, treatCors, stay.ethical = 100, switch.safety=0, n, wideForm = T)
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
