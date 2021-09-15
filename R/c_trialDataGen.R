#' Data Simulation (continuous snSMART design)
#'
#' simulate data for snSMART with continuous outcome (non-responders re-randomized)
#'
#' @param treatInfo Treatment distribution information. Each row represents one treatment, format: (mean of treatment effects, sd of treatment effects, c(priming effect on trt A,  trt B, trt C)), see example for detail
#' @param treatCors Treatment correlations matrix. Format example:
#' \[\[ (A,A) (A,B) (A,C)
#'    (B,A) (B,B) (B,C)
#'    (C,A) (C,B) (C,C) \]\], (A,B) denotes the treatment correlation between first stage treatment A and second stage treatment B
#' @param n vector of number of patients on each treatment e.g. n <- c(100, 100, 100)
#' @param stay.ethical numerical value, if stage 1 outcome is bigger than `stay.ethical` value, patient has probability of 1 of staying on the same treatment in stage 2
#' @param switch.safety numerical value, if stage 1 outcome is smaller than `switch.safety` value, patient has probability of 0 of staying on the same treatment in stage 2
#' @param wideForm whether to output the simulated dataset in wide form

#' @return The simulated dataset with 7 variables: patient id, treatment 1, stage 1 outcome, probability of stay on the same treatment, realization of staying flag (1 = stay on the same treatment, 0 = rerandomization in stage 2), treatment 2, stage 2 outcome
#'
#' @examples
#' treat.a<-c(70, 15, c(0,0,0))
#' treat.b<-c(50, 10, c(0,0,0))
#' treat.c<-c(60, 12, c(0,0,0))
#' treatInfo = rbind(treat.a, treat.b, treat.c)
#'
#' treatCors<-diag(.9, 3)
#'
#' switch.safety<-NULL
#' stay.ethical<- NULL

#' na<-100
#' nb<-100
#' nc<-100
#' n<-c(na,nb,nc)

#' trialData = c_trialDataGen(treatInfo, treatCors, n, wideForm = FALSE)
#'
#' @references
#' Hartman, H., Tamura, R.N., Schipper, M.J. and Kidwell, K.M., 2021. Design and analysis considerations for utilizing a mapping function in a small sample,
#' sequential, multiple assignment, randomized trials with continuous outcomes. Statistics in Medicine, 40(2), pp.312-326.
#'
#' @export
#'

library(tidyr)

c_trialDataGen = function(treatInfo, treatCors, n,
                        stay.ethical = NULL, switch.safety = NULL,
                        wideForm = TRUE){

  #number of treatments based on input values
  n.trt<-dim(treatInfo)[1]

  stage1outcome = rnorm(n[1], mean=treatInfo[1,1], sd=treatInfo[1,2])
  for(i in 2:n.trt){
    stage1outcome = c(stage1outcome, rnorm(n[i], mean=treatInfo[i,1], sd=treatInfo[i,2]))
  }

  #probability of staying on same trt for each patient
  pstay<-ifelse(stage1outcome/100 > 1, 1,
                ifelse(stage1outcome/100 < 0, 0 ,
                       stage1outcome/100))
  if(!is.null(switch.safety)){
    pstay<-ifelse(stage1outcome < switch.safety, 0, pstay)
  }
  if(!is.null(stay.ethical)){
    pstay<-ifelse(stage1outcome > stay.ethical, 1, pstay)
  }

  #realization of staying for each patient
  stay<-rbinom(sum(n),1,pstay)
  #vector of rerandomization treatments
  rerand = sample(2:n.trt,n[1],replace=T)
  trt.vect=1:n.trt
  for(i in 2:n.trt){
    rerand = c(rerand, sample(trt.vect[!trt.vect == i],n[i],replace=T))
  }
  trt1<-rep(1, n[1])
  for(i in 2:n.trt){
    trt1 <- c(trt1, rep(i, n[i]))
  }
  #realization of trt2 for each patient
  trt2<-ifelse(stay==1, trt1, rerand)

  stage2outcomeDistn = condMVNorm::condMVN(mean=c(treatInfo[trt1[1],1],treatInfo[trt2[1],1]+treatInfo[trt1[1],trt2[1]+2]),
                               sigma = matrix(c(treatInfo[trt1[1],2],
                                                treatCors[trt1[1], trt2[1]],
                                                treatCors[trt1[1], trt2[1]],
                                                treatInfo[trt2[1],2]),
                                              nrow=2),
                               dependent.ind = 2,
                               given.ind = 1,
                               X.given = stage1outcome[1])
  stage2outcome = rnorm(1, stage2outcomeDistn$condMean,stage2outcomeDistn$condVar)
  for(i in 2:sum(n)){
    stage2outcomeDistn = condMVNorm::condMVN(mean=c(treatInfo[trt1[i],1],treatInfo[trt2[i],1]+treatInfo[trt1[i],trt2[i]+2]),
                                 sigma = matrix(c(treatInfo[trt1[i],2],
                                                  treatCors[trt1[i], trt2[i]],
                                                  treatCors[trt1[i], trt2[i]],
                                                  treatInfo[trt2[i],2]),
                                                nrow=2),
                                 dependent.ind = 2,
                                 given.ind = 1,
                                 X.given = stage1outcome[i])
    stage2outcome = c(stage2outcome, rnorm(1, stage2outcomeDistn$condMean,stage2outcomeDistn$condVar))
  }

  id= seq(1:sum(n))
  trial.data=data.frame(id, trt1, stage1outcome, pstay, stay, trt2, stage2outcome, stage1outcome)

  if(wideForm == FALSE){
    trial.data.long = trial.data %>%
      gather(stage, outcome, c(stage1outcome, stage2outcome))
    trial.data.long$stage = ifelse(trial.data.long$stage == "stage1outcome", 1, 2)
    trial.data.long$trt = ifelse(trial.data.long$stage == 1, trial.data.long$trt1, trial.data.long$trt2)
    return(trial.data.long)
  }else{

  return(trial.data)
  }
}


