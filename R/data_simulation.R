#' Data Simulation ((snSMART with 3 active treatments and a binary outcome))
#'
#' simulate data for the standard design of snSMART (3 active treatments, non-responders re-randomized; binary outcome) based on response rate of each treatment. Useful for generating large number of simulations Not to be confused with \code{\link{trial_dataset}}
#'
#' @param p_trt vector of 3 values (first stage response rate of trt A, first stage response rate of trt B, first stage response rate of trt C). e.g. `p_trt = c(0.2, 0.3, 0.5)`
#' @param discount_y0 linkage parameters, for nonresponders to treatment k in the first stage who received treatment `k'` in the second stage, the second stage response rate is equal to `discount_y0 * pi_IK'`
#' @param discount_y1 linkage parameters, the second stage response rate for first stage responders is equal to `discount_y1 * p_IK`
#' @param p_1nA_2B probability for first stage non-responders to A being randomized to arm B in second stage
#' @param p_1nB_2A probability for first stage non-responders to B being randomized to arm A in second stage
#' @param p_1nC_2A probability for first stage non-responders to C being randomized to arm A in second stage
#' @param n vector of 3 values (first stage sample size of A, first stage sample size of B, first stage sample size of C). e.g. `n = c(30, 30, 30)`
#'
#' @return The simulated dataset with four columns: `response_stageI`, `treatment_stageI`, `response_stageII`, `treatment_stageII`
#'
#' @examples
#' data = data_simulation(p_trt = c(0.2, 0.3, 0.4),
#'     discount_y0 = c(0.6, 0.6, 0.6), discount_y1 = c(1.5, 1.5, 1.5),
#'     p_1nA_2B = 0.5, p_1nB_2A = 0.5, p_1nC_2A = 0.5, n = c(30, 30, 30))
#'
#' @references
#' Wei, B., Braun, T.M., Tamura, R.N. and Kidwell, K.M., 2018. A Bayesian analysis of small n sequential multiple assignment randomized trials (snSMARTs).
#' Statistics in medicine, 37(26), pp.3723-3732.
#'
#' Chao, Y.C., Trachtman, H., Gipson, D.S., Spino, C., Braun, T.M. and Kidwell, K.M., 2020. Dynamic treatment regimens in small n, sequential, multiple assignment, randomized trials: An application in focal segmental glomerulosclerosis. Contemporary clinical trials, 92, p.105989.

#' @seealso
#' \code{\link{trial_dataset}} \cr
#' \code{\link{BJSM_binary}} \cr
#' \code{\link{JSRM_binary}} \cr
#' \code{\link{sample_size}}
#'
#' @export
#'

data_simulation <- function(p_trt, discount_y0, discount_y1, p_1nA_2B, p_1nB_2A, p_1nC_2A, n){

  p_IA = p_trt[1]
  p_IB = p_trt[2]
  p_IC = p_trt[3]

  n_A = n[1]
  n_B = n[2]
  n_C = n[3]

  RESPONSE_RATE_STAGE_II_A_Y_0=p_IA*discount_y0[2]
  RESPONSE_RATE_STAGE_II_A_Y_1=p_IA*discount_y1[2]
  RESPONSE_RATE_STAGE_II_B_Y_0=p_IB*discount_y0[2]
  RESPONSE_RATE_STAGE_II_B_Y_1=p_IB*discount_y1[2]
  RESPONSE_RATE_STAGE_II_C_Y_0=p_IC*discount_y0[2]
  RESPONSE_RATE_STAGE_II_C_Y_1=p_IC*discount_y1[2]

  num_responseA_stageI=rbinom(n=1, size=n_A, prob=p_IA)
  num_negativeA_stageI=n_A-num_responseA_stageI

  sample_size_armB_negativeA_stageII=rbinom(n=1, size=num_negativeA_stageI, prob=p_1nA_2B)
  sample_size_armC_negativeA_stageII=num_negativeA_stageI-sample_size_armB_negativeA_stageII

  #   response_trtA_stageI=c(rep(1,num_responseA_stageI),rep(0,num_negativeA_stageI))
  #   data_response_trtA_stageI=data.frame(response=response_trtA_stageI,trt=rep(1,n_A))

  # trt B
  num_responseB_stageI=rbinom(n=1, size=n_B, prob=p_IB)
  num_negativeB_stageI=n_B-num_responseB_stageI

  sample_size_armA_negativeB_stageII=rbinom(n=1, size=num_negativeB_stageI, prob=p_1nB_2A)
  sample_size_armC_negativeB_stageII=num_negativeB_stageI-sample_size_armA_negativeB_stageII

  #   response_trtB_stageI=c(rep(1,num_responseB_stageI),rep(0,num_negativeB_stageI))
  #   data_response_trtB_stageI=data.frame(response=response_trtB_stageI,trt=rep(2,n_B))
  #

  # trt C
  num_responseC_stageI=rbinom(n=1, size=n_C, prob=p_IC)
  num_negativeC_stageI=n_C-num_responseC_stageI

  sample_size_armA_negativeC_stageII=rbinom(n=1, size=num_negativeC_stageI, prob=p_1nC_2A)
  sample_size_armB_negativeC_stageII=num_negativeC_stageI-sample_size_armA_negativeC_stageII

  #   response_trtC_stageI=c(rep(1,num_responseC_stageI),rep(0,num_negativeC_stageI))
  #   data_response_trtC_stageI=data.frame(response=response_trtC_stageI,trt=rep(3,n_C))
  #

  #stage II
  # trt A
  num_responseA_negativeB_stageII=rbinom(n=1, size=sample_size_armA_negativeB_stageII, prob=RESPONSE_RATE_STAGE_II_A_Y_0)
  data_negativeB_stageI_treatmentA_stageII=data.frame(response_stageI=c(rep(0,sample_size_armA_negativeB_stageII)),
                                                      treatment_stageI=c(rep(2,sample_size_armA_negativeB_stageII)),
                                                      response_stageII=c(rep(1,num_responseA_negativeB_stageII),rep(0,sample_size_armA_negativeB_stageII-num_responseA_negativeB_stageII)),
                                                      treatment_stageII=c(rep(1,sample_size_armA_negativeB_stageII)))

  num_responseA_negativeC_stageII=rbinom(n=1, size=sample_size_armA_negativeC_stageII, prob=RESPONSE_RATE_STAGE_II_A_Y_0)
  data_negativeC_stageI_treatmentA_stageII=data.frame(response_stageI=c(rep(0,sample_size_armA_negativeC_stageII)),
                                                      treatment_stageI=c(rep(3,sample_size_armA_negativeC_stageII)),
                                                      response_stageII=c(rep(1,num_responseA_negativeC_stageII),rep(0,sample_size_armA_negativeC_stageII-num_responseA_negativeC_stageII)),
                                                      treatment_stageII=c(rep(1,sample_size_armA_negativeC_stageII)))

  num_responseA_positiveA_stageII=rbinom(n=1, size=num_responseA_stageI, prob=RESPONSE_RATE_STAGE_II_A_Y_1)
  data_positiveA_stageI_treatmentA_stageII=data.frame(response_stageI=c(rep(1,num_responseA_stageI)),
                                                      treatment_stageI=c(rep(1,num_responseA_stageI)),
                                                      response_stageII=c(rep(1,num_responseA_positiveA_stageII),rep(0,num_responseA_stageI-num_responseA_positiveA_stageII)),
                                                      treatment_stageII=c(rep(1,num_responseA_stageI)))

  # trt B
  num_responseB_negativeA_stageII=rbinom(n=1, size=sample_size_armB_negativeA_stageII, prob=RESPONSE_RATE_STAGE_II_B_Y_0)
  data_negativeA_stageI_treatmentB_stageII=data.frame(response_stageI=c(rep(0,sample_size_armB_negativeA_stageII)),
                                                      treatment_stageI=c(rep(1,sample_size_armB_negativeA_stageII)),
                                                      response_stageII=c(rep(1,num_responseB_negativeA_stageII),rep(0,sample_size_armB_negativeA_stageII-num_responseB_negativeA_stageII)),
                                                      treatment_stageII=c(rep(2,sample_size_armB_negativeA_stageII)))

  num_responseB_negativeC_stageII=rbinom(n=1, size=sample_size_armB_negativeC_stageII, prob=RESPONSE_RATE_STAGE_II_B_Y_0)
  data_negativeC_stageI_treatmentB_stageII=data.frame(response_stageI=c(rep(0,sample_size_armB_negativeC_stageII)),
                                                      treatment_stageI=c(rep(3,sample_size_armB_negativeC_stageII)),
                                                      response_stageII=c(rep(1,num_responseB_negativeC_stageII),rep(0,sample_size_armB_negativeC_stageII-num_responseB_negativeC_stageII)),
                                                      treatment_stageII=c(rep(2,sample_size_armB_negativeC_stageII)))
  num_responseB_positiveB_stageII=rbinom(n=1, size=num_responseB_stageI, prob=RESPONSE_RATE_STAGE_II_B_Y_1)
  data_positiveB_stageI_treatmentB_stageII=data.frame(response_stageI=c(rep(1,num_responseB_stageI)),
                                                      treatment_stageI=c(rep(2,num_responseB_stageI)),
                                                      response_stageII=c(rep(1,num_responseB_positiveB_stageII),rep(0,num_responseB_stageI-num_responseB_positiveB_stageII)),
                                                      treatment_stageII=c(rep(2,num_responseB_stageI)))

  # trt C
  num_responseC_negativeA_stageII=rbinom(n=1, size=sample_size_armC_negativeA_stageII, prob=RESPONSE_RATE_STAGE_II_C_Y_0)
  data_negativeA_stageI_treatmentC_stageII=data.frame(response_stageI=c(rep(0,sample_size_armC_negativeA_stageII)),
                                                      treatment_stageI=c(rep(1,sample_size_armC_negativeA_stageII)),
                                                      response_stageII=c(rep(1,num_responseC_negativeA_stageII),rep(0,sample_size_armC_negativeA_stageII-num_responseC_negativeA_stageII)),
                                                      treatment_stageII=c(rep(3,sample_size_armC_negativeA_stageII)))
  num_responseC_negativeB_stageII=rbinom(n=1, size=sample_size_armC_negativeB_stageII, prob=RESPONSE_RATE_STAGE_II_C_Y_0)
  data_negativeB_stageI_treatmentC_stageII=data.frame(response_stageI=c(rep(0,sample_size_armC_negativeB_stageII)),
                                                      treatment_stageI=c(rep(2,sample_size_armC_negativeB_stageII)),
                                                      response_stageII=c(rep(1,num_responseC_negativeB_stageII),rep(0,sample_size_armC_negativeB_stageII-num_responseC_negativeB_stageII)),
                                                      treatment_stageII=c(rep(3,sample_size_armC_negativeB_stageII)))
  num_responseC_positiveC_stageII=rbinom(n=1, size=num_responseC_stageI, prob=RESPONSE_RATE_STAGE_II_C_Y_1)
  data_positiveC_stageI_treatmentC_stageII=data.frame(response_stageI=c(rep(1,num_responseC_stageI)),
                                                      treatment_stageI=c(rep(3,num_responseC_stageI)),
                                                      response_stageII=c(rep(1,num_responseC_positiveC_stageII),rep(0,num_responseC_stageI-num_responseC_positiveC_stageII)),
                                                      treatment_stageII=c(rep(3,num_responseC_stageI)))

  # formulate data
  # when formulate the data, if we don't want to use event driven simulation,
  # we need to formulate the data based on the last stage bz we only know all the stage outcome
  # at the last stage, which has already been done

  data_stageI_stageII_tmp=rbind(data_negativeB_stageI_treatmentA_stageII,
                                data_negativeC_stageI_treatmentA_stageII,
                                data_positiveA_stageI_treatmentA_stageII,
                                data_negativeA_stageI_treatmentB_stageII,
                                data_negativeC_stageI_treatmentB_stageII,
                                data_positiveB_stageI_treatmentB_stageII,
                                data_negativeA_stageI_treatmentC_stageII,
                                data_negativeB_stageI_treatmentC_stageII,
                                data_positiveC_stageI_treatmentC_stageII
  )
  # responders will be two and non-responders will be one
  data_stageI_stageII=cbind(data_stageI_stageII_tmp)

  return(data_stageI_stageII)
}

