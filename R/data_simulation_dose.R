#' Data Simulation (dose level design)
#'
#' simulate data for the dose level design of snSMART (placebo, low, high dose; binary outcome)
#'
#' @param p_IA first stage response rate of placebo
#' @param p_IB first stage response rate of low dose treatment
#' @param p_IC first stage response rate of high dose treatment
#' @param RESPONSE_RATE_DISCOUNT_P_Y_0 linkage parameters, for nonresponders to placebo in the first stage who received who received treatment `k'` in the second stage, the second stage response rate in the second stage is equal to `RESPONSE_RATE_DISCOUNT_P_Y_0 * pi_IK'`
#' @param RESPONSE_RATE_DISCOUNT_P_Y_1 linkage parameters, for responders to placebo in the first stage who received who received treatment `k'` in the second stage, the second stage response rate in the second stage is equal to `RESPONSE_RATE_DISCOUNT_P_Y_1 * pi_IK'`
#' @param RESPONSE_RATE_DISCOUNT_L_Y_0 linkage parameters, for nonresponders to low dose treatment in the first stage who received who received treatment `k'` in the second stage, the second stage response rate in the second stage is equal to `RESPONSE_RATE_DISCOUNT_L_Y_0 * pi_IK'`
#' @param RESPONSE_RATE_DISCOUNT_L_Y_1 linkage parameters, for responders to low dose treatment in the first stage who received who received treatment `k'` in the second stage, the second stage response rate in the second stage is equal to `RESPONSE_RATE_DISCOUNT_L_Y_1 * pi_IK'`
#' @param RESPONSE_RATE_DISCOUNT_H_Y_0 similar to `RESPONSE_RATE_DISCOUNT_L_Y_0`
#' @param RESPONSE_RATE_DISCOUNT_H_Y_1 similar to `RESPONSE_RATE_DISCOUNT_L_Y_1`
#' @param RAND_PROB_POS_P_TRT_H probability for first stage responders to placebo being randomized to high dose treatment in second stage
#' @param RAND_PROB_POS_P_TRT_L probability for first stage responders to placebo being randomized to low dose treatment in second stage
#' @param RAND_PROB_NEG_P_TRT_L probability for first stage nonresponders to placebo being randomized to low dose treatment in second stage
#' @param RAND_PROB_NEG_P_TRT_H probability for first stage nonresponders to placebo being randomized to high dose treatment in second stage
#' @param RAND_PROB_POS_L_TRT_L similar to above
#' @param RAND_PROB_POS_L_TRT_H similar to above
#' @param RAND_PROB_NEG_L_TRT_L similar to above
#' @param RAND_PROB_NEG_L_TRT_H similar to above
#' @param RAND_PROB_POS_H_TRT_L similar to above
#' @param RAND_PROB_POS_H_TRT_H similar to above
#' @param n_A first stage sample size of A
#' @param n_B first stage sample size of B
#' @param n_C first stage sample size of C
#'
#' @return The simulated dataset with five columns: response_stageI, treatment_stageI, response_stageII, treatment_stageII, response_status_stageI (responders will be two and non-responders will be one)
#'
#' @examples
#' data = data_simulation_dose(p_IP = 0.15, p_IL = 0.15, p_IH = 0.15, RESPONSE_RATE_DISCOUNT_P_Y_0 = 0.9, RESPONSE_RATE_DISCOUNT_P_Y_1 = 1.3, RESPONSE_RATE_DISCOUNT_L_Y_0 = 0.8, RESPONSE_RATE_DISCOUNT_L_Y_1 = 1.2,
#' RESPONSE_RATE_DISCOUNT_H_Y_0 = 0.8, RESPONSE_RATE_DISCOUNT_H_Y_1 = 1.2, RAND_PROB_POS_P_TRT_H = 0.5, RAND_PROB_POS_P_TRT_L = 0.5, RAND_PROB_NEG_P_TRT_L = 0.5,
#' RAND_PROB_NEG_P_TRT_H = 0.5, RAND_PROB_POS_L_TRT_L = 0.5, RAND_PROB_POS_L_TRT_H = 0.5, RAND_PROB_NEG_L_TRT_L = 0.5, RAND_PROB_NEG_L_TRT_H = 0.5,
#' RAND_PROB_POS_H_TRT_L = 0.5, RAND_PROB_POS_H_TRT_H = 0.5, n_P = 30, n_L = 30, n_H = 30)
#'
#' @references
#' Fang, F., Hochstedler, K.A., Tamura, R.N., Braun, T.M. and Kidwell, K.M., 2021. Bayesian methods to compare dose levels with placebo in a small n, sequential, multiple assignment, randomized trial.
#' Statistics in Medicine, 40(4), pp.963-977.

#'
#' @export
#'

data_simulation_dose <- function(p_IP, p_IL, p_IH, RESPONSE_RATE_DISCOUNT_P_Y_0, RESPONSE_RATE_DISCOUNT_P_Y_1, RESPONSE_RATE_DISCOUNT_L_Y_0, RESPONSE_RATE_DISCOUNT_L_Y_1,
                                 RESPONSE_RATE_DISCOUNT_H_Y_0, RESPONSE_RATE_DISCOUNT_H_Y_1, RAND_PROB_POS_P_TRT_H, RAND_PROB_POS_P_TRT_L, RAND_PROB_NEG_P_TRT_L,
                                 RAND_PROB_NEG_P_TRT_H, RAND_PROB_POS_L_TRT_L, RAND_PROB_POS_L_TRT_H, RAND_PROB_NEG_L_TRT_L, RAND_PROB_NEG_L_TRT_H,
                                 RAND_PROB_POS_H_TRT_L, RAND_PROB_POS_H_TRT_H, n_P, n_L, n_H){
  #stage I
  # Placebo
  num_responseP_stageI=rbinom(n=1, size=n_P, prob=p_IP)
  num_negativeP_stageI=n_P-num_responseP_stageI

  sample_size_armL_negativeP_stageII=rbinom(n=1, size=num_negativeP_stageI, prob=RAND_PROB_NEG_P_TRT_L)
  sample_size_armH_negativeP_stageII=num_negativeP_stageI-sample_size_armL_negativeP_stageII


  sample_size_armH_responseP_stageII=rbinom(n=1, size=num_responseP_stageI, prob=RAND_PROB_POS_P_TRT_H)
  sample_size_armL_responseP_stageII=num_responseP_stageI-sample_size_armH_responseP_stageII

  # Low dose
  num_responseL_stageI=rbinom(n=1, size=n_L, prob=p_IL)
  num_negativeL_stageI=n_L-num_responseL_stageI

  sample_size_armL_negativeL_stageII=rbinom(n=1, size=num_negativeL_stageI, prob=RAND_PROB_NEG_L_TRT_L)
  sample_size_armH_negativeL_stageII=num_negativeL_stageI-sample_size_armL_negativeL_stageII

  sample_size_armL_responseL_stageII=rbinom(n=1, size=num_responseL_stageI, prob=RAND_PROB_POS_L_TRT_L)
  sample_size_armH_responseL_stageII=num_responseL_stageI-sample_size_armL_responseL_stageII

  # High dose
  num_responseH_stageI=rbinom(n=1, size=n_H, prob=p_IH)
  num_negativeH_stageI=n_H-num_responseH_stageI

  sample_size_armL_responseH_stageII=rbinom(n=1, size=num_responseH_stageI, prob=RAND_PROB_POS_H_TRT_L)
  sample_size_armH_responseH_stageII=num_responseH_stageI-sample_size_armL_responseH_stageII

  sample_size_armH_negativeH_stageII=num_negativeH_stageI


  # Placebo: 0P, 1P

  RESPONSE_RATE_STAGE_II_PL_Y_0=p_IL*RESPONSE_RATE_DISCOUNT_P_Y_0
  RESPONSE_RATE_STAGE_II_PH_Y_0=p_IH*RESPONSE_RATE_DISCOUNT_P_Y_0

  RESPONSE_RATE_STAGE_II_PL_Y_1=p_IL*RESPONSE_RATE_DISCOUNT_P_Y_1
  RESPONSE_RATE_STAGE_II_PH_Y_1=p_IH*RESPONSE_RATE_DISCOUNT_P_Y_1

  # Placebo: 0L, 1L

  RESPONSE_RATE_STAGE_II_LL_Y_0=p_IL*RESPONSE_RATE_DISCOUNT_L_Y_0
  RESPONSE_RATE_STAGE_II_LH_Y_0=p_IH*RESPONSE_RATE_DISCOUNT_L_Y_0

  RESPONSE_RATE_STAGE_II_LL_Y_1=p_IL*RESPONSE_RATE_DISCOUNT_L_Y_1
  RESPONSE_RATE_STAGE_II_LH_Y_1=p_IH*RESPONSE_RATE_DISCOUNT_L_Y_1

  # Placebo: 0H, 1H

  RESPONSE_RATE_STAGE_II_HH_Y_0=p_IH*RESPONSE_RATE_DISCOUNT_H_Y_0

  RESPONSE_RATE_STAGE_II_HL_Y_1=p_IL*RESPONSE_RATE_DISCOUNT_H_Y_1
  RESPONSE_RATE_STAGE_II_HH_Y_1=p_IH*RESPONSE_RATE_DISCOUNT_H_Y_1

  #stage II

  #0PL
  num_responseL_negativeP_stageII=rbinom(n=1, size=sample_size_armL_negativeP_stageII, prob=RESPONSE_RATE_STAGE_II_PL_Y_0)
  data_negativeP_stageI_treatmentL_stageII=data.frame(response_stageI=c(rep(0,sample_size_armL_negativeP_stageII)),
                                                      treatment_stageI=c(rep(1,sample_size_armL_negativeP_stageII)),
                                                      response_stageII=c(rep(1,num_responseL_negativeP_stageII),rep(0,sample_size_armL_negativeP_stageII-num_responseL_negativeP_stageII)),
                                                      treatment_stageII=c(rep(2,sample_size_armL_negativeP_stageII)))
  #0PH
  num_responseH_negativeP_stageII=rbinom(n=1, size=sample_size_armH_negativeP_stageII, prob=RESPONSE_RATE_STAGE_II_PH_Y_0)
  data_nagativeP_stageI_treatmentH_stageII=data.frame(response_stageI=c(rep(0,sample_size_armH_negativeP_stageII)),
                                                      treatment_stageI=c(rep(1,sample_size_armH_negativeP_stageII)),
                                                      response_stageII=c(rep(1,num_responseH_negativeP_stageII),rep(0,sample_size_armH_negativeP_stageII-num_responseH_negativeP_stageII)),
                                                      treatment_stageII=c(rep(3,sample_size_armH_negativeP_stageII)))

  # 1PL
  num_responseL_responseP_stageII=rbinom(n=1, size=sample_size_armL_responseP_stageII, prob=RESPONSE_RATE_STAGE_II_PL_Y_1)
  data_responseP_stageI_treatmentL_stageII=data.frame(response_stageI=c(rep(1,sample_size_armL_responseP_stageII)),
                                                      treatment_stageI=c(rep(1,sample_size_armL_responseP_stageII)),
                                                      response_stageII=c(rep(1,num_responseL_responseP_stageII),rep(0,sample_size_armL_responseP_stageII-num_responseL_responseP_stageII)),
                                                      treatment_stageII=c(rep(2,sample_size_armL_responseP_stageII)))
  # 1PH
  num_responseH_responseP_stageII=rbinom(n=1, size=sample_size_armH_responseP_stageII, prob=RESPONSE_RATE_STAGE_II_PH_Y_1)
  data_responseP_stageI_treatmentH_stageII=data.frame(response_stageI=c(rep(1,sample_size_armH_responseP_stageII)),
                                                      treatment_stageI=c(rep(1,sample_size_armH_responseP_stageII)),
                                                      response_stageII=c(rep(1,num_responseH_responseP_stageII),rep(0,sample_size_armH_responseP_stageII-num_responseH_responseP_stageII)),
                                                      treatment_stageII=c(rep(3,sample_size_armH_responseP_stageII)))
  #0LL
  num_responseL_negativeL_stageII=rbinom(n=1, size=sample_size_armL_negativeL_stageII, prob=RESPONSE_RATE_STAGE_II_LL_Y_0)
  data_negativeL_stageI_treatmentL_stageII=data.frame(response_stageI=c(rep(0,sample_size_armL_negativeL_stageII)),
                                                      treatment_stageI=c(rep(2,sample_size_armL_negativeL_stageII)),
                                                      response_stageII=c(rep(1,num_responseL_negativeL_stageII),rep(0,sample_size_armL_negativeL_stageII-num_responseL_negativeL_stageII)),
                                                      treatment_stageII=c(rep(2,sample_size_armL_negativeL_stageII)))
  #0LH
  num_responseH_negativeL_stageII=rbinom(n=1, size=sample_size_armH_negativeL_stageII, prob=RESPONSE_RATE_STAGE_II_LH_Y_0)
  data_nagativeL_stageI_treatmentH_stageII=data.frame(response_stageI=c(rep(0,sample_size_armH_negativeL_stageII)),
                                                      treatment_stageI=c(rep(2,sample_size_armH_negativeL_stageII)),
                                                      response_stageII=c(rep(1,num_responseH_negativeL_stageII),rep(0,sample_size_armH_negativeL_stageII-num_responseH_negativeL_stageII)),
                                                      treatment_stageII=c(rep(3,sample_size_armH_negativeL_stageII)))
  #1LL
  num_responseL_responseL_stageII=rbinom(n=1, size=sample_size_armL_responseL_stageII, prob=RESPONSE_RATE_STAGE_II_LL_Y_1)
  data_responseL_stageI_treatmentL_stageII=data.frame(response_stageI=c(rep(1,sample_size_armL_responseL_stageII)),
                                                      treatment_stageI=c(rep(2,sample_size_armL_responseL_stageII)),
                                                      response_stageII=c(rep(1,num_responseL_responseL_stageII),rep(0,sample_size_armL_responseL_stageII-num_responseL_responseL_stageII)),
                                                      treatment_stageII=c(rep(2,sample_size_armL_responseL_stageII)))
  #1LH
  num_responseH_responseL_stageII=rbinom(n=1, size=sample_size_armH_responseL_stageII, prob=RESPONSE_RATE_STAGE_II_LH_Y_1)
  data_responseL_stageI_treatmentH_stageII=data.frame(response_stageI=c(rep(1,sample_size_armH_responseL_stageII)),
                                                      treatment_stageI=c(rep(2,sample_size_armH_responseL_stageII)),
                                                      response_stageII=c(rep(1,num_responseH_responseL_stageII),rep(0,sample_size_armH_responseL_stageII-num_responseH_responseL_stageII)),
                                                      treatment_stageII=c(rep(3,sample_size_armH_responseL_stageII)))
  #0HH

  num_responseH_negativeH_stageII=rbinom(n=1, size=sample_size_armH_negativeH_stageII, prob=RESPONSE_RATE_STAGE_II_HH_Y_0)
  data_negativeH_stageI_treatmentH_stageII=data.frame(response_stageI=c(rep(0,sample_size_armH_negativeH_stageII)),
                                                      treatment_stageI=c(rep(3,sample_size_armH_negativeH_stageII)),
                                                      response_stageII=c(rep(1,num_responseH_negativeH_stageII),rep(0,sample_size_armH_negativeH_stageII-num_responseH_negativeH_stageII)),
                                                      treatment_stageII=c(rep(3,sample_size_armH_negativeH_stageII)))
  #1HL
  num_responseL_responseH_stageII=rbinom(n=1, size=sample_size_armL_responseH_stageII, prob=RESPONSE_RATE_STAGE_II_HL_Y_1)
  data_responseH_stageI_treatmentL_stageII=data.frame(response_stageI=c(rep(1,sample_size_armL_responseH_stageII)),
                                                      treatment_stageI=c(rep(3,sample_size_armL_responseH_stageII)),
                                                      response_stageII=c(rep(1,num_responseL_responseH_stageII),rep(0,sample_size_armL_responseH_stageII-num_responseL_responseH_stageII)),
                                                      treatment_stageII=c(rep(2,sample_size_armL_responseH_stageII)))
  #1HH

  num_responseH_responseH_stageII=rbinom(n=1, size=sample_size_armH_responseH_stageII, prob=RESPONSE_RATE_STAGE_II_HH_Y_1)
  data_responseH_stageI_treatmentH_stageII=data.frame(response_stageI=c(rep(1,sample_size_armH_responseH_stageII)),
                                                      treatment_stageI=c(rep(3,sample_size_armH_responseH_stageII)),
                                                      response_stageII=c(rep(1,num_responseH_responseH_stageII),rep(0,sample_size_armH_responseH_stageII-num_responseH_responseH_stageII)),
                                                      treatment_stageII=c(rep(3,sample_size_armH_responseH_stageII)))

  data_stageI_stageII_tmp=rbind(data_responseP_stageI_treatmentH_stageII,
                                data_negativeP_stageI_treatmentL_stageII,
                                data_responseP_stageI_treatmentL_stageII,
                                data_nagativeP_stageI_treatmentH_stageII,
                                data_negativeL_stageI_treatmentL_stageII,
                                data_responseL_stageI_treatmentL_stageII,
                                data_nagativeL_stageI_treatmentH_stageII,
                                data_responseL_stageI_treatmentH_stageII,
                                data_negativeH_stageI_treatmentH_stageII,
                                data_responseH_stageI_treatmentH_stageII,
                                data_responseH_stageI_treatmentL_stageII
  )
  data_stageI_stageII_tmp$response_status_stageI<-data_stageI_stageII_tmp$response_stageI+1

  return(data_stageI_stageII_tmp)
}
