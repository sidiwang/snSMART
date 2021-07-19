#' Data Simulation (continuous snSMART design)
#'
#' simulate data for snSMART with continuous outcome (non-responders re-randomized)
#'
#' @param stage1effects vector of mean stage 1 effect of each treatment e.g. stage1effects<-c(50,60,70)
#' @param stage2weights first stage response rate of B
#' @param stagecorrs vector of correlations of stage effects. The first number is the correlation of effects between the same treatment in stage 1 and 2, and the second number is the correlation of effects between different treatments in stage 1 and 2. e.g. stagecorrs <- c(0.8, 0.3)
#' @param variance variance of stage 1 effect among patients
#' @param n vector of number of patients on each treatment e.g. n <- c(100, 100, 100)
#' @param stay.ethical numerical value, if stage 1 outcome is bigger than `stay.ethical` value, patient has probability of 1 of staying on the same treatment in stage 2
#' @param switch.safety numerical value, if stage 1 outcome is smaller than `switch.safety` value, patient has probability of 0 of staying on the same treatment in stage 2
#' @param wideForm whether to output the simulated dataset in wide form

#' @return The simulated dataset with 7 variables: patient id, treatment 1, stage 1 outcome, probability of stay on the same treatment, realization of staying flag (1 = stay on the same treatment, 0 = rerandomization in stage 2), treatment 2, stage 2 outcome
#'
#' @examples
#' stage1effects<-c(50,60,70)
#' stage2weights<-c(.2, .8, 5)
#' stagecorrs<-c(.8,.3)
#' variance<-10
#' switch.safety<-NULL
#' stay.ethical<- NULL

#' na<-100
#' nb<-100
#' nc<-100
#' n<-c(na,nb,nc)

#' trialData = trialDataGen(stage1effects, stage2weights, stagecorrs, variance, n)
#'
#' @references
#' Hartman, H., Tamura, R.N., Schipper, M.J. and Kidwell, K.M., 2021. Design and analysis considerations for utilizing a mapping function in a small sample,
#' sequential, multiple assignment, randomized trials with continuous outcomes. Statistics in Medicine, 40(2), pp.312-326.
#'
#' @export
#'

library(tidyr)

trialDataGen = function(stage1effects, stage2weights, stagecorrs, variance, n,
                        stay.ethical = NULL, switch.safety = NULL,
                        wideForm = TRUE){

  #number of treatments based on input values
  n.trt<-length(stage1effects)

  stage1outcome = c()
  for(i in 1:n.trt){
    stage1outcome = c(stage1outcome, rnorm(n[i], mean=stage1effects[i], sd=sqrt(variance)))
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
  rerand = c()
  trt.vect=1:n.trt
  for(i in 1:n.trt){
    rerand = c(rerand, sample(trt.vect[!trt.vect == i],n[i],replace=T))
  }
  trt1<-c()
  for(i in 1:n.trt){
    trt1 <- c(trt1, rep(i, n[i]))
  }
  #realization of trt2 for each patient
  trt2<-ifelse(stay==1, trt1, rerand)

  stage2outcome<-c()
  for(i in 1:sum(n)){
    stage2outcomeDistn = condMVNorm::condMVN(mean=c(stage1effects[trt1[i]],stage2weights[1]*stage1effects[trt1[i]]+stage2weights[2]*stage1effects[trt2[i]] + stage2weights[3]*stay[i]),
                                 sigma = variance * matrix(c(1,
                                                             stagecorrs[1]*stay[i]+ stagecorrs[2]*(1-stay[i]),
                                                             stagecorrs[1]*stay[i]+ stagecorrs[2]*(1-stay[i]),
                                                             1),
                                                           nrow=2),
                                 dependent.ind = 2,
                                 given.ind = 1,
                                 X.given = stage1outcome[i])
    stage2outcome = c(stage2outcome, rnorm(1, stage2outcomeDistn$condMean,stage2outcomeDistn$condVar))
  }

  id= seq(1:sum(n))
  trial.data=data.frame(id, trt1, stage1outcome, pstay, stay, trt2, stage2outcome)

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


