
BJSM_6betas_missing_text <- function() {
  return(
    "model{
  for (i in 1:n1){   # n1 is number of first-stage outcomes
  # likelihood
    Y1[i]~dbern(pi_1[i])
  # explaining
    pi_1[i] <- pi[treatment_stageI[i]]
  }
  for (k in 1:n2){   # n2 is number of second-stage outcomes
  # likelihood
    Y2[k]~dbern(pi_2[k])
  # explaining
    pi_2[k] <- pi[treatment_stageII[k]] * beta[response_stageI_disc[k]]
  }
  for (j in 1:num_arms){
    pi[j]~pi_prior_dist(pi_prior.a[j],pi_prior.b[j])
  }
  beta[1]~beta0_prior_dist(beta0_prior.a,beta0_prior.b)
  beta[2]~beta1_prior_dist(beta1_prior.a,beta1_prior.c)
  beta[3]~beta0_prior_dist(beta0_prior.a,beta0_prior.b)
  beta[4]~beta1_prior_dist(beta1_prior.a,beta1_prior.c)
  beta[5]~beta0_prior_dist(beta0_prior.a,beta0_prior.b)
  beta[6]~beta1_prior_dist(beta1_prior.a,beta1_prior.c)
}"
  )
}


BJSM_6betas_missing_new_text <- function() {
  return(
    "model{
  for (i in 1:n1){   # n1 is number of first-stage outcomes
  # likelihood
    Y1[i]~dbern(pi_1[i])
  # explaining
    pi_1[i] <- pi[treatment_stageI[i]]
  }
  for (k in 1:n2){   # n2 is number of second-stage outcomes
  # likelihood
    Y2[k]~dbern(pi_2[k])
  # explaining
    pi_2[k] <- pi[treatment_stageII[k]] * beta[response_stageI_disc[k]]
  }
  for (j in 1:num_arms){
    pi[j]~dbeta(pi_prior.a[j],pi_prior.b[j])
  }
  beta[1]~dbeta(beta0_prior.a,beta0_prior.b)
  beta[2]~dpar(beta1_prior.a,beta1_prior.c)
  beta[3]~dbeta(beta0_prior.a,beta0_prior.b)
  beta[4]~dpar(beta1_prior.a,beta1_prior.c)
  beta[5]~dbeta(beta0_prior.a,beta0_prior.b)
  beta[6]~dpar(beta1_prior.a,beta1_prior.c)
}"
  )
}

BJSM_2beta_missing_text <- function() {
  return(
    "model{
  for (i in 1:n1){   # n1 is number of first-stage outcomes
  # likelihood
    Y1[i]~dbern(pi_1[i])
  # explaining
    pi_1[i] <- pi[treatment_stageI[i]]
  }
  for (k in 1:n2){   # n2 is number of second-stage outcomes
  # likelihood
    Y2[k]~dbern(pi_2[k])
  # explaining
    pi_2[k] <- pi[treatment_stageII[k]] * beta[response1[k]]
  }

  for (j in 1:num_arms){
    pi[j]~pi_prior_dist(pi_prior.a[j],pi_prior.b[j])
  }
  beta[1]~beta0_prior_dist(beta0_prior.a,beta0_prior.b)
  beta[2]~beta1_prior_dist(beta1_prior.a,beta1_prior.c)
}"
  )
}

BJSM_2beta_missing_new_text <- function() {
  return(
    "model{
  for (i in 1:n1){   # n1 is number of first-stage outcomes
  # likelihood
    Y1[i]~dbern(pi_1[i])
  # explaining
    pi_1[i] <- pi[treatment_stageI[i]]
  }
  for (k in 1:n2){   # n2 is number of second-stage outcomes
  # likelihood
    Y2[k]~dbern(pi_2[k])
  # explaining
    pi_2[k] <- pi[treatment_stageII[k]] * beta[response1[k]]
  }

  for (j in 1:num_arms){
    pi[j]~dbeta(pi_prior.a[j],pi_prior.b[j])
  }
  beta[1]~dbeta(beta0_prior.a,beta0_prior.b)
  beta[2]~dpar(beta1_prior.a,beta1_prior.c)
}"
  )
}

BJSM_dose_text <- function() {
  return(
    "model{
	for(i in 1:overall_sample_size) {
  #likelihood
    response_stageI[i]~dbern(response_rate_stageI[i])
    response_stageII[i]~dbern(response_rate_stageII[i])
  #explaining
    response_rate_stageI[i]<-pi[treatment_stageI[i]]
    response_rate_stageII[i]<-pi[treatment_stageII[i]]*beta[response_discount_status_stageI[i],treatment_stageI[i]]


}
  #give prior for first stage

  pi[1]~pi_prior_dist(a_pi[1],b_pi[1])

  #Consider the relationship between pi_L, pi_H and pi_P

  pi[2] <- pi[1]*exp(tau1)
  tau1 ~dnorm(normal.mean,1/normal.var)
  pi[3] <- pi[1]*exp(tau2)
  tau2 ~dnorm(normal.mean,1/normal.var)


  for(i in 1:2){
    for(j in 1:num_arms){
      beta[i,j]~beta_prior_dist(a_beta,b_beta)
      }
  }

}"
  )
}

BJSM_dose_new_text <- function() {
  return(
    "model{
	for(i in 1:overall_sample_size) {
  #likelihood
    response_stageI[i]~dbern(response_rate_stageI[i])
    response_stageII[i]~dbern(response_rate_stageII[i])
  #explaining
    response_rate_stageI[i]<-pi[treatment_stageI[i]]
    response_rate_stageII[i]<-pi[treatment_stageII[i]]*beta[response_discount_status_stageI[i],treatment_stageI[i]]


}
  #give prior for first stage

  pi[1]~dbeta(a_pi[1],b_pi[1])

  #Consider the relationship between pi_L, pi_H and pi_P

  pi[2] <- pi[1]*exp(tau1)
  tau1 ~dnorm(normal.mean,1/normal.var)
  pi[3] <- pi[1]*exp(tau2)
  tau2 ~dnorm(normal.mean,1/normal.var)


  for(i in 1:2){
    for(j in 1:num_arms){
      beta[i,j]~dgamma(a_beta,b_beta)
      }
  }

}"
  )
}

csnSMART_text <- function() {
  return(
    "#continuous snSMART bayesian analysis

#Holly Hartman
#March 18, 2019

#Inputs:
##Y (matrix), trt1, trt2, stay,ntrt,
##priors: xi, phi, tau, sigma

#Estimates:
##xi, phi, tau, sigma


model{
  for (i in 1:n){   # n is total sample size
    # likelihood
    Y[i, 1:2]~dmnorm(mu[i,1:2], rho[stay1[i], , ]) #multivariate normal
      #rho is variance-covariance matrix

    mu[i,1] <- xi_[trt1[i]]
    mu[i,2] <-phi1*xi_[trt1[i]] + (1-phi1)*xi_[trt2[i]]+ phi3*stay[i]
	 }

  #Priors
  for (j in 1:ntrt){
    xi_[j]~dnorm(xi_prior.mean[j], xi_prior.sd[j])
  }
  for(k in 1:2){
	  rho[k,1:2,1:2] ~ dwish(Omega1[,], 2)
  }
  phi1~dunif(0, 0.5)
  phi3~dnorm(0, phi3_prior.sd) T(0,)

 Omega1[1,1] <- 1
 Omega1[2,2] <- 1
 Omega1[1,2] <- 0
 Omega1[2,1] <- 0
}"
  )
}

Bayes_AR_text <- function() {
  return(
    "model{
  for (i in 1:n1){   # n1 is number of first-stage outcomes
  # likelihood
    Y1[i]~dbern(pi_1[i])
  # explaining
    pi_1[i] <- pi[treatment_stageI[i]]
  }
  for (k in 1:n2){   # n2 is number of second-stage outcomes
  # likelihood
    Y2[k]~dbern(pi_2[k])
  # explaining
    pi_2[k] <- pi[treatment_stageII[k]] * beta[response_stageI_disc[k]]
  }
  for (j in 1:num_arms){
    pi[j]~pi_prior_dist(pi_prior.a[j],pi_prior.b[j])
  }
  beta[1]~beta0_prior_dist(beta0_prior.a,beta0_prior.b)
  beta[2]~beta1_prior_dist(beta1_prior.a,beta1_prior.c)
  beta[3]~beta0_prior_dist(beta0_prior.a,beta0_prior.b)
  beta[4]~beta1_prior_dist(beta1_prior.a,beta1_prior.c)
  beta[5]~beta0_prior_dist(beta0_prior.a,beta0_prior.b)
  beta[6]~beta1_prior_dist(beta1_prior.a,beta1_prior.c)
}"
  )
}

Bayes_AR_new_text <- function() {
  return(
    "model{
  for (i in 1:n1){   # n1 is number of first-stage outcomes
  # likelihood
    Y1[i]~dbern(pi_1[i])
  # explaining
    pi_1[i] <- pi[treatment_stageI[i]]
  }
  for (k in 1:n2){   # n2 is number of second-stage outcomes
  # likelihood
    Y2[k]~dbern(pi_2[k])
  # explaining
    pi_2[k] <- pi[treatment_stageII[k]] * beta[response_stageI_disc[k]]
  }
  for (j in 1:num_arms){
    pi[j]~dbeta(pi_prior.a[j],pi_prior.b[j])
  }
  beta[1]~dbeta(beta0_prior.a,beta0_prior.b)
  beta[2]~dpar(beta1_prior.a,beta1_prior.c)
  beta[3]~dbeta(beta0_prior.a,beta0_prior.b)
  beta[4]~dpar(beta1_prior.a,beta1_prior.c)
  beta[5]~dbeta(beta0_prior.a,beta0_prior.b)
  beta[6]~dpar(beta1_prior.a,beta1_prior.c)
}"
  )
}

Bayes_text <- function() {
  return(
    "model{
  for (i in 1:n){   # n is total sample size
  # likelihood
    Y1[i]~dbern(pi_1[i])
    Y2[i]~dbern(pi_2[i])
  # explaining
    pi_1[i] <- pi[treatment_stageI[i]]
    pi_2[i] <- pi[treatment_stageII[i]] * beta[response_stageI_disc[i]]
  }

  for (j in 1:num_arms){
    pi[j]~pi_prior_dist(pi_prior.a[j],pi_prior.b[j])
  }
  beta[1]~beta0_prior_dist(beta0_prior.a,beta0_prior.b)
  beta[2]~beta1_prior_dist(beta1_prior.a,beta1_prior.c)
  beta[3]~beta0_prior_dist(beta0_prior.a,beta0_prior.b)
  beta[4]~beta1_prior_dist(beta1_prior.a,beta1_prior.c)
  beta[5]~beta0_prior_dist(beta0_prior.a,beta0_prior.b)
  beta[6]~beta1_prior_dist(beta1_prior.a,beta1_prior.c)
}"
  )
}

Bayes_new_text <- function() {
  return(
    "model{
  for (i in 1:n){   # n is total sample size
  # likelihood
    Y1[i]~dbern(pi_1[i])
    Y2[i]~dbern(pi_2[i])
  # explaining
    pi_1[i] <- pi[treatment_stageI[i]]
    pi_2[i] <- pi[treatment_stageII[i]] * beta[response_stageI_disc[i]]
  }

  for (j in 1:num_arms){
    pi[j]~dbeta(pi_prior.a[j],pi_prior.b[j])
  }
  beta[1]~dbeta(beta0_prior.a,beta0_prior.b)
  beta[2]~dpar(beta1_prior.a,beta1_prior.c)
  beta[3]~dbeta(beta0_prior.a,beta0_prior.b)
  beta[4]~dpar(beta1_prior.a,beta1_prior.c)
  beta[5]~dbeta(beta0_prior.a,beta0_prior.b)
  beta[6]~dpar(beta1_prior.a,beta1_prior.c)
}"
  )
}
