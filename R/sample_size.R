#' snSMART Sample Size Calculation
#'
#' conduct Bayesian sample size calculation for a standard snSMART design to distinguish the best treatment from the second-best treatment using the Bayesian joint stage model
#'
#' @param piA the response rate (ranges from 0.01 to 0.99) for treatment A
#' @param piB the response rate (ranges from 0.01 to 0.99) for treatment B
#' @param piC the response rate (ranges from 0.01 to 0.99) for treatment C
#' @param beta1 the linkage parameter (ranges from 1.00 to 1/largest response rate) for first stage responders. (A smaller value leads to more conservative sample size calculation because two stages are less correlated)
#' @param beta0 the linkage parameter (ranges from 0.01 to 0.99) for first stage non-responders. A larger value leads to a more conservative sample size calculation because two stages are less correlated
#' @param coverage the coverage rate (ranges from 0.01 to 0.99) for the posterior difference of top two treatments
#' @param power the probability (ranges from 0.01 to 0.99) for identify the best treatment
#' @param muA the prior mean (ranges from 0.01 to 0.99) for treatment A
#' @param muB the prior mean (ranges from 0.01 to 0.99) for treatment B
#' @param muC the prior mean (ranges from 0.01 to 0.99) for treatment C
#' @param nA the prior sample size (larger than 0) for treatment A
#' @param nB the prior sample size (larger than 0) for treatment B
#' @param nC the prior sample size (larger than 0) for treatment C
#'
#' @return the estimated sample size per arm for an snSMART
#'
#' @examples
#' sampleSize = sample_size(piA = 0.7, piB = 0.5, piC = 0.25, beta1 = 1.4, beta0 = 0.5, coverage = 0.9, power = 0.8, muA = 0.65, muB = 0.55, muC = 0.25, nA = 4, nB = 2, nC = 3)
#'
#' @references
#' Wei, B., Braun, T.M., Tamura, R.N. and Kidwell, K.M., 2018. A Bayesian analysis of small n sequential multiple assignment randomized trials (snSMARTs).
#' Statistics in medicine, 37(26), pp.3723-3732.
#'
#' Wei, B., Braun, T.M., Tamura, R.N. and Kidwell, K., 2020. Sample size determination for Bayesian analysis of small n sequential, multiple assignment, randomized trials (snSMARTs) with three agents. Journal of Biopharmaceutical Statistics, 30(6), pp.1109-1120.
#'
#' @export
#'


sample_size <- function(piA, piB, piC, beta1, beta0, coverage, power, muA, muB, muC, nA, nB, nC){

  max_beta1=1/min(c(max(piA,piB,piC),0.99))
  str_beta1=paste0("The value of ","\u03B21"," should be within 1~",max_beta1)


  if(!piA<1 & piA>0){print("The value of \u03C0A should be within 0.01~0.99")}
  if(!piB<1 & piB>0){print("The value of \u03C0B should be within 0.01~0.99")}
  if(!piC<1 & piC>0){print("The value of \u03C0C should be within 0.01~0.99")}
  if(!(beta1<=max_beta1 & beta1>1)|
     is.na(max_beta1)|
     piA>1 | piA<0|
     piB>1 | piB<0|
     piC>1 | piC<0){
    print(str_beta1)}
  if(!beta0<1 & beta0>0){print("The value of \u03B20 should be within 0.01~0.99")}
  if(!coverage<1 & coverage>0){print("The value of 1-\u03B1 should be within 0.01~0.99")}
  if(!power<1 & power>0){print("The value of 1-\u03BE should be within 0.01~0.99")}

  if(!muA<1 & muA>0){print("The value of \u03bcA should be within 0.01~0.99")}
  if(!muB<1 & muB>0){print("The value of \u03bcB should be within 0.01~0.99")}
  if(!muC<1 & muC>0){print("The value of \u03bcC should be within 0.01~0.99")}
  if(!nA>0){print("The value of nA should be greater than 0")}
  if(!nB>0){print("The value of nB should be greater than 0")}
  if(!nC>0){print("The value of nC should be greater than 0")}


  K=3
  magic_Z=1.5
  pi_A=as.numeric(piA)
  pi_B=as.numeric(piB)
  pi_C=as.numeric(piC)

  muA_input=as.numeric(muA)
  muB_input=as.numeric(muB)
  muC_input=as.numeric(muC)

  nA_input=as.numeric(nA)
  nB_input=as.numeric(nB)
  nC_input=as.numeric(nC)

  LIST_OF_PIS=c(pi_A,pi_B,pi_C)
  LIST_OF_PIS=LIST_OF_PIS[order(LIST_OF_PIS,decreasing = T)]

  COVRAGE=as.numeric(coverage)
  # COVRAGE_BONFF=1-(1-COVRAGE)/(K-1)
  # COVRAGE=COVRAGE_BONFF
  POW = as.numeric(power) #when POW is small, we need to decrease sample size lower limit, need to estimate the lower limit or there will not be a opposite sign

  # SAMPLE_SIZE_LLIMIT=1
  # SAMPLE_SIZE_ULIMIT=1000

  # CIL_MIN=0.1# can not be too small
  # CIL_MAX=2
  # CIL_STEP=0.01

  SS_LOW=1
  SS_HIGH=300
  CONVERGE_TOL=0.001
  CIL_MIN=0.01# can not be too small
  CIL_STEP_I=0.01

  ###################
  # check if pis are all the same
  # see if the first two treatment are the same-if not the same then do computation, if the same, do computation with the first and last
  # need all treatment arm for computation
  if(LIST_OF_PIS[1]==LIST_OF_PIS[2]){
    piA=LIST_OF_PIS[1]
    piB=LIST_OF_PIS[K]
    piC=LIST_OF_PIS[setdiff(1:K,c(1,K))]
  }else{
    piA=LIST_OF_PIS[1]
    piB=LIST_OF_PIS[2]
    piC=LIST_OF_PIS[setdiff(1:K,c(1,2))]
  }

  # generate truncated beta1 and beta0
  set.seed(199)
  # generate pareto truncated at 1.5 with location 1 and scale 3, mean 1.5
  beta1_mean=as.numeric(beta1)
  pareto_beta=1/(1-1/beta1_mean)
  # beta1_sample=rtrunc(LINKAGE_SAMPLE, spec="pareto",a=1,b=1/max(c(piA,piB,piC)),1,pareto_beta)
  # beta1_sample=rep(1,LINKAGE_SAMPLE)
  beta0_mean=as.numeric(beta0)
  # beta0_sample=rbeta(LINKAGE_SAMPLE,beta0_mean*2,2-beta0_mean*2)
  # beta0_sample=rep(1,LINKAGE_SAMPLE)
  # plot(hist(beta0))
  # load functions
  #require(rmutil)
  CIL_MAX=(piA-piB)*magic_Z
  beta1_sample=round(mean(truncdist::rtrunc(99999, spec="pareto",a=1,b=1/max(c(piA,piB,piC)),1,pareto_beta)),3)
  beta0_sample=rep(beta0_mean,1)

  sample_size_list_pair1=NULL

  error_count_pair1=0
  warn_count_pair1=0

  sim_count=0

  error_round_pair1=NULL
  warn_round_pair1=NULL
  error_mesg_pair1=NULL
  warn_mesg_pair1=NULL

  error_round_pair1=NULL
  warn_round_pair1=NULL
  error_mesg_pair1=NULL
  warn_mesg_pair1=NULL

  # calculate priors

  beta1=beta1_sample
  beta0=beta0_sample

  aA = muA_input*nA_input
  bA = nA_input-aA
  # aA = 0.25*2
  # bA = 2-aA
  source(file = "inst/functions.R")

  cA = cFunc(aA, bA, beta1)
  dA = dFunc(aA, bA, beta1)
  eA = eFunc(aA, bA, beta0, K)
  fA = fFunc(aA, bA, beta0, K)

  aB = muB_input*nB_input
  bB = nB_input-aB
  # aB = 0.25*2
  # bB = 2-aB
  cB = cFunc(aB, bB, beta1)
  dB = dFunc(aB, bB, beta1)
  eB = eFunc(aB, bB, beta0, K)
  fB = fFunc(aB, bB, beta0, K)

  aC = muC_input*nC_input
  bC = nC_input-aC
  # aC = 0.25*2
  # bC = 2-aC
  cC = cFunc(aC, bC, beta1)
  dC = dFunc(aC, bC, beta1)
  eC = eFunc(aC, bC, beta0, K)
  fC = fFunc(aC, bC, beta0, K)

  ciZ=qnorm(1-(1-COVRAGE)/2) # critical value of coverage

  # i=0

  for(CIL_I in seq(CIL_MAX,CIL_MIN,by=-CIL_STEP_I)){
    # i=i+1
    # print(CIL_I)
    # print(CIL_I)
    # get sample size solved
    # print(CIL_I)
    # CIL_I=0.26
    if(CIL_I/2==(max(c(piA,piB,piC))-min(c(piA,piB,piC)))){
      next
    }
    tryCatch({

      fun <- function (x) {
        2*ciZ*sqrt(var1_o1o2_diff(K=K,
                                  piA=piA, piB=piB, piC=piC,
                                  beta1=beta1, beta0=beta0,
                                  aA=aA, bA=bA, cA=cA, dA=dA, eA=eA, fA=fA,
                                  aB=aB, bB=bB, cB=cB, dB=dB, eB=eB, fB=fB,
                                  aC=aC, bC=bC, cC=cC, dC=dC, eC=eC, fC=fC,
                                  n=x)-mean_o1o2_diff(K=K,
                                                      piA=piA, piB=piB, piC=piC,
                                                      beta1=beta1, beta0=beta0,
                                                      aA=aA, bA=bA, cA=cA, dA=dA, eA=eA, fA=fA,
                                                      aB=aB, bB=bB, cB=cB, dB=dB, eB=eB, fB=fB,
                                                      aC=aC, bC=bC, cC=cC, dC=dC, eC=eC, fC=fC,
                                                      n=x)^2)-CIL_I
      }

      # for(i in 1:300){
      #   print(fun(i))
      # }

      # ciL=0.28

      # sample_size_tmp_pair1=ceiling(uniroot(fun, c(SS_LOW, SS_HIGH),tol = CONVERGE_TOL)$root)
      sample_size_tmp_pair1=ceiling(uniroot(fun, c(SS_LOW, SS_HIGH),tol = CONVERGE_TOL)$root)
      # All <- uniroot.all(sample_size_equation, c(0, 8))
      # calculate powerpower
      # 1-pnorm((ciL/2-(max(c(piA,piB)-min(piA,piB))))/sqrt(sigmaSqABDiffFunc(sample_size_tmp_pair1)))+pnorm((-ciL/2-(max(c(piA,piB)-min(piA,piB))))/sqrt(sigmaSqABDiffFunc(sample_size_tmp_pair1)))
      mu_o1o2_diff=mean_o1o2_diff(K=K,
                                  piA=piA, piB=piB, piC=piC,
                                  beta1=beta1, beta0=beta0,
                                  aA=aA, bA=bA, cA=cA, dA=dA, eA=eA, fA=fA,
                                  aB=aB, bB=bB, cB=cB, dB=dB, eB=eB, fB=fB,
                                  aC=aC, bC=bC, cC=cC, dC=dC, eC=eC, fC=fC,
                                  n=sample_size_tmp_pair1)
      mu_o1o2_sq_diff=var1_o1o2_diff(K=K,
                                     piA=piA, piB=piB, piC=piC,
                                     beta1=beta1, beta0=beta0,
                                     aA=aA, bA=bA, cA=cA, dA=dA, eA=eA, fA=fA,
                                     aB=aB, bB=bB, cB=cB, dB=dB, eB=eB, fB=fB,
                                     aC=aC, bC=bC, cC=cC, dC=dC, eC=eC, fC=fC,
                                     n=sample_size_tmp_pair1)
      var_o1o2_diff=mu_o1o2_sq_diff-mu_o1o2_diff^2
      # pow_pair1=(1-pnorm((ciL/2-mu_o1o2_diff)/sqrt(var_o1o2_diff)))
      pow_pair1=(1-pnorm((CIL_I/2-mu_o1o2_diff)/sqrt(var_o1o2_diff)))

      # total_process=length(seq(CIL_MAX,CIL_MIN,by=-CIL_STEP_I))


      if(pow_pair1>POW) break
    },
    error = function(c) {
      # error_round_tmp_pair1=cbind(i,beta1,beta0,CIL_I)
      # error_round_pair1=rbind(error_round_pair1,error_round_tmp_pair1)
      # next
    },
    warning = function(c) {
      # warn_round_tmp_pair1=cbind(i,beta1,beta0,CIL_I)
      # warn_round_pair1=rbind(warn_round_pair1,warn_round_tmp_pair1)
      # print(i)
      # print(CIL_I)
      # next
    },
    finally = {# posterior_sample_burn=window(posterior_sample,start=BURNING, end=MCMC_SAMPLE)
      # posterior_sample_cmb=do.call(rbind, posterior_sample_burn)
    }
    )
    # print(CIL_I)
    # print(pow_pair1)
    # print(mu_o1o2_diff)
    # print(mu_o1o2_sq_diff)
    # print(sample_size_tmp_pair1)


  }

  print(cat(paste0("With given settings, the estimated sample size per arm for an snSMART is: ", sample_size_tmp_pair1, "\n",
               "This implies that for an snSMART with sample size of ", sample_size_tmp_pair1, " per arm (", 3*sample_size_tmp_pair1, " in total for three agents):", "\n",
               "The probability of successfully identifying the best treatment is ", power, " when the difference of response rates between the best and second best treatment is at least ", LIST_OF_PIS[1] - LIST_OF_PIS[2], ", and the response rate of the best treatment is ", LIST_OF_PIS[1], "\n")))

  return(sample_size_tmp_pair1)

}
