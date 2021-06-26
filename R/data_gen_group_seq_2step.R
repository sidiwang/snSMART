#' Data generation of a group sequential snSMART with a two-step rule
#'
#' generate data for the Group Sequential snSMART, which is designed based on the standard design of snSMART (3 treatments, non-responders re-randomized; binary outcome)
#'
#' @param pi_1A first stage response rate of A
#' @param pi_1B first stage response rate of B
#' @param pi_1C first stage response rate of C
#' @param discount_y linkage parameters, for responders to treatment k in the first stage who received treatment k again in the second stage, the second stage response rate is equal to `discount_y * pi_IK`
#' @param discount_n1 linkage parameters, vector of three values: 1) nonresponders to A receive treatment B in second stage, 2) nonresponders to B receive treatment A in second stage, 3) nonresponders to C receive A in the second stage
#' @param discount_n2 linkage parameters, vector of three values: 1) nonresponders to A receive treatment C in second stage, 2) nonresponders to B receive treatment C in second stage, 3) nonresponders to C receive B in the second stage
#' @param rate accrual rate. Rate equals to 5 if participants are enrolled in the study at the rate of five people per month
#' @param n.month number of months to enroll subjects
#' @param n.update number of updates during the study. Default is 1 (number of interim analysis? need to check with Kelley)
#' @param drop_threshold_large a number between 0 and 1. Default value is 0.5, the number of dropping threshold should be the same as the number update in the study. See the details section for more explaination
#' @param drop_threshold_small a number between 0 and 1. Default value is 0.5, the number of dropping threshold should be the same as the number update in the study. See the details section for more explaination
#' @param NUM_ARMS number of treatment arms
#' @param pi_prior.a parameter a of the prior distribution for pi_1K, a vector with three values, one for each treatment. Please check the `Details` section for more explaination
#' @param pi_prior.b parameter b of the prior distribution  for pi_1K, a vector with three values, one for each treatment. Please check the `Details` section for more explaination
#' @param beta0_prior.a parameter a of the prior distribution for linkage parameter beta0
#' @param beta0_prior.b parameter b of the prior distribution  for linkage parameter beta0
#' @param beta1_prior.a parameter a of the prior distribution for linkage parameter beta1
#' @param beta1_prior.c parameter b of the prior distribution for linkage parameter beta1
#' @param n_MCMC_chain number of MCMC chains, default to 1. If this is set to a number more than 1
#' @param BURN.IN number of burn-in iterations for MCMC
#' @param MCMC_SAMPLE number of iterations for MCMC
#' @param pi_prior_dist prior distribution for pi, user can choose from gamma, beta, pareto
#' @param beta0_prior_dist prior distribution for beta0, user can choose from gamma, beta, pareto
#' @param beta1_prior_dist prior distribution for beta1, user can choose from gamma, beta, pareto
#'
#' @details
#' (paper provided in the reference section, section 2.2.2 Bayesian decision rules. drop_threshold_large and drop_threshold_small are corresponding to `\tau_l` and `\phi_l` respectively {need to check with Kelley})
#'
#' @return
#' return the simulated dataset, randomization probabilities before each update (the number of rows equals to the number of updates plus 1), removed arm and the round of the update that removal of treatment arm occurs
#'
#'
#' @example
#' simu1 = data_gen_group_seq_2step(pi_1A = 0.25, pi_1B = 0.45, pi_1C = 0.65, discount_y = c(1.5,1.5,1.5), discount_n1 = c(0.8, 0.8, 0.8), discount_n2 = c(0.8, 0.8, 0.8),
#' rate = 5, n.month = 60, n.update = 1, drop_threshold_large = 0.5, drop_threshold_small = 0.5, NUM_ARMS  = 3, pi_prior_dist = "beta", pi_prior.a =  c(0.4,0.4,0.4), pi_prior.b = c(1.6, 1.6, 1.6), beta0_prior_dist = "beta",
#' beta0_prior.a = 1.6, beta0_prior.b = 0.4, beta1_prior_dist = "pareto", beta1_prior.a = 3, beta1_prior.c = 1, MCMC_SAMPLE = 60000, BURN.IN = 10000, n_MCMC_chain = 1)
#'
#' simu2 = data_gen_group_seq_2step(pi_1A = 0.1, pi_1B = 0.1, pi_1C = 0.11, discount_y = c(1.5,1.5,1.5), discount_n1 = c(0.8, 0.8, 0.8), discount_n2 = c(0.8, 0.8, 0.8),
#' rate = 5, n.month = 60, n.update = 2, drop_threshold_large = c(0.5, 0.5), drop_threshold_small = c(0.5, 0.5), NUM_ARMS  = 3, pi_prior_dist = "beta", pi_prior.a =  c(0.4,0.4,0.4), pi_prior.b = c(1.6, 1.6, 1.6), beta0_prior_dist = "beta",
#' beta0_prior.a = 1.6, beta0_prior.b = 0.4, beta1_prior_dist = "pareto", beta1_prior.a = 3, beta1_prior.c = 1, MCMC_SAMPLE = 60000, BURN.IN = 10000, n_MCMC_chain = 1)

#' @references
#' Chao, Y.C., Braun, T.M., Tamura, R.N. and Kidwell, K.M., 2020. A Bayesian group sequential small n sequential multiple‐assignment randomized trial. Journal of the Royal Statistical Society: Series C (Applied Statistics), 69(3), pp.663-680.
#'
#' @export



## Data generation of a group sequential snSMART with a two-step rule
data_gen_group_seq_2step <- function(pi_1A,pi_1B,pi_1C,discount_y,discount_n1,discount_n2, rate, n.month,n.update=1,drop_threshold_large=0.5,drop_threshold_small=0.5,
                                     NUM_ARMS,pi_prior_dist, pi_prior.a, pi_prior.b, beta0_prior_dist, beta0_prior.a, beta0_prior.b, beta1_prior_dist,
                                     beta1_prior.a, beta1_prior.c, MCMC_SAMPLE, BURN.IN, n_MCMC_chain){

  pi_prior_dist = ifelse(pi_prior_dist == "gamma", "dgamma", ifelse(pi_prior_dist == "beta", "dbeta", "dpar"))
  beta0_prior_dist = ifelse(beta0_prior_dist == "gamma", "dgamma", ifelse(beta0_prior_dist == "beta", "dbeta", "dpar"))
  beta1_prior_dist = ifelse(beta1_prior_dist == "gamma", "dgamma", ifelse(beta1_prior_dist == "beta", "dbeta", "dpar"))

  assn.stage2 <- function(i, trt, y, rand.prob)    # Function that assigns the second stage treatment
  {
    alltrt <- 1:3
    if (y[i]==1) newtrt <- trt[i]
    if (y[i]==0) newtrt <- sample(alltrt[alltrt!=trt[i]], 1, prob=rand.prob[alltrt!=trt[i]])
    return(newtrt)
  }


  if (length(drop_threshold_large)!=length(drop_threshold_small))
    stop("Please assign same number of dropping thresholds for pi_largest and pi_smallest")
  if (length(drop_threshold_large)!=n.update)
    stop("Please assign same number of dropping threshold as the number of update")
  ## Bayes randomization setup
  file_path <- 'inst'
  jags.model.name.update <- 'Bayes_AR_new.bug'

  pi_2A_y <<- pi_1A * discount_y[1]        # Second stage response rate of responders to A
  pi_2B.A_n <<- pi_1A * discount_n1[2]     # Second stage response rate of non-responders to B who receive A in the second stage
  pi_2C.A_n <<- pi_1A * discount_n1[3]     # Second stage response rate of non-responders to C who receive A in the second stage
  pi_2B_y <<- pi_1B * discount_y[2]        # Second stage response rate of responders to B
  pi_2A.B_n <<- pi_1B * discount_n1[1]     # Second stage response rate of non-responders to A who receive B in the second stage
  pi_2C.B_n <<- pi_1B * discount_n2[3]     # Second stage response rate of non-responders to C who receive B in the second stage
  pi_2C_y <<- pi_1C * discount_y[3]        # Second stage response rate of responders to C
  pi_2A.C_n <<- pi_1C * discount_n2[1]     # Second stage response rate of non-responders to A who receive C in the second stage
  pi_2B.C_n <<- pi_1C * discount_n2[2]     # Second stage response rate of non-responders to B who receive C in the second stage

  N <- n.month * rate  # total number of subjects

  t <- sapply(1:n.month,function(x) runif(rate,min=(x-1)*30,max=x*30))
  time.1st.trt <- round(t[order(t)],0)
  time.1st.resp <- time.1st.trt + 180
  time.2nd.trt <- time.1st.trt + 181
  time.2nd.resp <- time.2nd.trt + 180
  patient_entry <- data.frame(time.1st.trt = time.1st.trt,
                              time.1st.resp = time.1st.resp,
                              time.2nd.trt = time.2nd.trt,
                              time.2nd.resp = time.2nd.resp,
                              trt.1st = rep(NA,N),
                              resp.1st = rep(NA,N),
                              trt.2nd = rep(NA,N),
                              resp.2nd = rep(NA,N))
  n.group <- n.update + 1    # n.update is the number of times of update, n.group is the number of groups that we divide the subjects into
  subject.update <- rate*n.month*seq(1,n.group,1)/n.group
  time.update <- c(-1,time.1st.resp[subject.update])
  rand.prob <- matrix(NA,nrow=n.update+1,ncol=3)
  rand.prob[1,] <- c(1/3,1/3,1/3)
  pi_hat <- matrix(NA,nrow=n.update,ncol=3)

  for (k in 1:n.update){
    # assignment
    temp.ind <- which(time.update[k]<time.1st.trt&time.1st.trt<=time.update[k+1])
    patient_entry$trt.1st[temp.ind] <- sample(1:3,length(temp.ind),replace=T,prob=rand.prob[k,])
    patient_entry$resp.1st[temp.ind] <- (patient_entry$trt.1st[temp.ind]==1)*rbinom(length(temp.ind),1,prob=pi_1A) +
      (patient_entry$trt.1st[temp.ind]==2)*rbinom(length(temp.ind),1,prob=pi_1B) + (patient_entry$trt.1st[temp.ind]==3)*rbinom(length(temp.ind),1,prob=pi_1C)

    temp.ind <- which(time.update[k]<time.2nd.trt&time.2nd.trt<=time.update[k+1])
    patient_entry$trt.2nd[temp.ind] <- sapply(temp.ind,assn.stage2,trt=patient_entry$trt.1st,
                                              y=patient_entry$resp.1st,rand.prob=rand.prob[k,])

    # table(patient_entry[temp.ind,]$trt.2nd[patient_entry[temp.ind,]$trt.1st==1&patient_entry[temp.ind,]$resp.1st==0])
    # prop.table(table(patient_entry[temp.ind,]$trt.2nd[patient_entry[temp.ind,]$trt.1st==1&patient_entry[temp.ind,]$resp.1st==0]))

    patient_entry$resp.2nd[temp.ind] <- (patient_entry$trt.1st[temp.ind]==1&patient_entry$trt.2nd[temp.ind]==1)*rbinom(length(temp.ind),1,prob=pi_2A_y) +
      (patient_entry$trt.1st[temp.ind]==1&patient_entry$trt.2nd[temp.ind]==2)*rbinom(length(temp.ind),1,prob=pi_2A.B_n) +
      (patient_entry$trt.1st[temp.ind]==1&patient_entry$trt.2nd[temp.ind]==3)*rbinom(length(temp.ind),1,prob=pi_2A.C_n) +
      (patient_entry$trt.1st[temp.ind]==2&patient_entry$trt.2nd[temp.ind]==2)*rbinom(length(temp.ind),1,prob=pi_2B_y) +
      (patient_entry$trt.1st[temp.ind]==2&patient_entry$trt.2nd[temp.ind]==1)*rbinom(length(temp.ind),1,prob=pi_2B.A_n) +
      (patient_entry$trt.1st[temp.ind]==2&patient_entry$trt.2nd[temp.ind]==3)*rbinom(length(temp.ind),1,prob=pi_2B.C_n) +
      (patient_entry$trt.1st[temp.ind]==3&patient_entry$trt.2nd[temp.ind]==3)*rbinom(length(temp.ind),1,prob=pi_2C_y) +
      (patient_entry$trt.1st[temp.ind]==3&patient_entry$trt.2nd[temp.ind]==1)*rbinom(length(temp.ind),1,prob=pi_2C.A_n) +
      (patient_entry$trt.1st[temp.ind]==3&patient_entry$trt.2nd[temp.ind]==2)*rbinom(length(temp.ind),1,prob=pi_2C.B_n)

    # randomization update

    patient_entry$disc <- 2 * patient_entry$trt.1st - (patient_entry$resp.1st == 0)
    error_ind <- 0

    bugfile  <- readLines("inst/Bayes_AR.bug")
    bugfile  <- gsub(pattern = "pi_prior_dist", replace = pi_prior_dist, x = bugfile)
    bugfile  <- gsub(pattern = "beta0_prior_dist", replace = beta0_prior_dist, x = bugfile)
    bugfile2  <- gsub(pattern = "beta1_prior_dist", replace = beta1_prior_dist, x = bugfile)

    writeLines(bugfile2, con="inst/Bayes_AR_new.bug")

    tryCatch({
      jags <- jags.model(file.path(file_path,jags.model.name.update),
                         data=list(n1 = nrow(patient_entry[patient_entry$time.1st.resp<=time.update[k+1],]),
                                   n2 = nrow(patient_entry[patient_entry$time.2nd.resp<=time.update[k+1],]),
                                   num_arms = NUM_ARMS,
                                   Y1 = patient_entry$resp.1st[patient_entry$time.1st.resp<=time.update[k+1]],
                                   Y2 = patient_entry$resp.2nd[patient_entry$time.2nd.resp<=time.update[k+1]],
                                   treatment_stageI = patient_entry$trt.1st[patient_entry$time.1st.resp<=time.update[k+1]],
                                   treatment_stageII = patient_entry$trt.2nd[patient_entry$time.2nd.resp<=time.update[k+1]],
                                   response_stageI_disc = patient_entry$disc[patient_entry$time.2nd.resp<=time.update[k+1]],
                                   #prior
                                   pi_prior.a = pi_prior.a,
                                   pi_prior.b = pi_prior.b,
                                   beta0_prior.a = beta0_prior.a,
                                   beta0_prior.b = beta0_prior.b,
                                   beta1_prior.a = beta1_prior.a,
                                   beta1_prior.c = beta1_prior.c),
                         n.chains=n_MCMC_chain,n.adapt = BURN.IN)
      posterior_sample <- coda.samples(jags,
                                       c('pi','beta'),
                                       MCMC_SAMPLE)
    },
    warning = function(war){
      warning_count <- warning_count + 1
      err_war_message <- rbind(paste("The warning ", warning_count, " is: ", war))
    },
    error = function(err){
      error_count <- error_count + 1
      err_war_message <- rbind(paste("The error ", error_count, " is: ", err))
      error_ind <- 1
    },
    finally = {
      print(k)     # show the number of iterations run
    }
    )
    out_post <- posterior_sample[[1]]
    min_A <- mean(apply(out_post[,7:9],1,function(x) {x[1]==min(x)}))  # posterior probability that A has smallest response rate
    min_B <- mean(apply(out_post[,7:9],1,function(x) {x[2]==min(x)}))  # posterior probability that B has smallest response rate
    min_C <- mean(apply(out_post[,7:9],1,function(x) {x[3]==min(x)}))  # posterior probability that C has smallest response rate
    max_A <- mean(apply(out_post[,7:9],1,function(x) {x[1]==max(x)}))  # posterior probability that A has largest response rate
    max_B <- mean(apply(out_post[,7:9],1,function(x) {x[2]==max(x)}))  # posterior probability that B has largest response rate
    max_C <- mean(apply(out_post[,7:9],1,function(x) {x[3]==max(x)}))  # posterior probability that C has largest response rate
    keep_A <- (max_A>drop_threshold_large[k])
    keep_B <- (max_B>drop_threshold_large[k])
    keep_C <- (max_C>drop_threshold_large[k])
    if(any(c(keep_A,keep_B,keep_C)>0)){
      drop_A <- ((keep_A==0)*(min_A==max(min_A,min_B,min_C))==1)
      drop_B <- ((keep_B==0)*(min_B==max(min_A,min_B,min_C))==1)
      drop_C <- ((keep_C==0)*(min_C==max(min_A,min_B,min_C))==1)
      if (sum(drop_A,drop_B,drop_C)>1){    # if more than one arm is dropped
        randompick <- sample(1:3,1,prob=c(drop_A,drop_B,drop_C)/sum(drop_A,drop_B,drop_C))
        drop_A <- (randompick==1)
        drop_B <- (randompick==2)
        drop_C <- (randompick==3)
      }
    } else {
      drop_A <- (min_A>drop_threshold_small[k])
      drop_B <- (min_B>drop_threshold_small[k])
      drop_C <- (min_C>drop_threshold_small[k])
    }

    if(all(c(drop_A,drop_B,drop_C)==0)){     # if none of the arm is dropped, move on to next update
      rand.prob[k+1,] <- rand.prob[k,]
      n.round <- k
    } else if(sum(drop_A,drop_B,drop_C)==1){  # if one of the arms is dropped, move on to last assignment
      rand.prob[k+1,c(drop_A,drop_B,drop_C)] <- 0
      rand.prob[k+1,!c(drop_A,drop_B,drop_C)] <- c(1/2,1/2)
      n.round <- k
      break
    }
  }

  # last assignment
  if (sum(drop_A,drop_B,drop_C)<=1){
    temp.ind <- which(time.update[n.round+1]<time.1st.trt)
    patient_entry$trt.1st[temp.ind] <- sample(1:3,length(temp.ind),replace=T,prob=rand.prob[n.round+1,])
    patient_entry$resp.1st[temp.ind] <- (patient_entry$trt.1st[temp.ind]==1)*rbinom(length(temp.ind),1,prob=pi_1A) +
      (patient_entry$trt.1st[temp.ind]==2)*rbinom(length(temp.ind),1,prob=pi_1B) + (patient_entry$trt.1st[temp.ind]==3)*rbinom(length(temp.ind),1,prob=pi_1C)

    temp.ind <- which(time.update[n.round+1]<time.2nd.trt)
    patient_entry$trt.2nd[temp.ind] <- sapply(temp.ind,assn.stage2,trt=patient_entry$trt.1st,
                                              y=patient_entry$resp.1st,rand.prob=rand.prob[n.round+1,])
    patient_entry$resp.2nd[temp.ind] <- (patient_entry$trt.1st[temp.ind]==1&patient_entry$trt.2nd[temp.ind]==1)*rbinom(length(temp.ind),1,prob=pi_2A_y) +
      (patient_entry$trt.1st[temp.ind]==1&patient_entry$trt.2nd[temp.ind]==2)*rbinom(length(temp.ind),1,prob=pi_2A.B_n) +
      (patient_entry$trt.1st[temp.ind]==1&patient_entry$trt.2nd[temp.ind]==3)*rbinom(length(temp.ind),1,prob=pi_2A.C_n) +
      (patient_entry$trt.1st[temp.ind]==2&patient_entry$trt.2nd[temp.ind]==2)*rbinom(length(temp.ind),1,prob=pi_2B_y) +
      (patient_entry$trt.1st[temp.ind]==2&patient_entry$trt.2nd[temp.ind]==1)*rbinom(length(temp.ind),1,prob=pi_2B.A_n) +
      (patient_entry$trt.1st[temp.ind]==2&patient_entry$trt.2nd[temp.ind]==3)*rbinom(length(temp.ind),1,prob=pi_2B.C_n) +
      (patient_entry$trt.1st[temp.ind]==3&patient_entry$trt.2nd[temp.ind]==3)*rbinom(length(temp.ind),1,prob=pi_2C_y) +
      (patient_entry$trt.1st[temp.ind]==3&patient_entry$trt.2nd[temp.ind]==1)*rbinom(length(temp.ind),1,prob=pi_2C.A_n) +
      (patient_entry$trt.1st[temp.ind]==3&patient_entry$trt.2nd[temp.ind]==2)*rbinom(length(temp.ind),1,prob=pi_2C.B_n)
    dropped_arm <- (drop_A==1) * 1 + (drop_B==1) * 2 + (drop_C==1) * 3
  }

  if (sum(drop_A,drop_B,drop_C)==0){
    dropped_round <- 0
  } else if (sum(drop_A,drop_B,drop_C)==1){
    dropped_round <- n.round
  }

  colnames(rand.prob) = c("trtA", "trtB", "trtC")

  return(list("simulated.data" = patient_entry[,-9], "randomization.probabilities" = rand.prob, "dropped.arm" = dropped_arm, "dropped.round" = dropped_round))

}
