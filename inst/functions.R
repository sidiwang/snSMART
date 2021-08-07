beta_prior_generator=function(info_level, prior_mean){
  alpha0=prior_mean*info_level
  beta0=info_level*(1-prior_mean)
  alpha0_beta0=list(alpha0=alpha0,
                    beta0=beta0)
}

# cFunc[a_, b_, beta1_] = (beta1 a^2)/(a + b)
#
cFunc=function(a, b, beta1){
  c= (beta1*a^2)/(a + b)
  return(c)
}
# cFunc(1,2,3)
#
# dFunc[a_, b_, beta0_] = (a b - beta1 a^2)/(a + b)
#
dFunc=function(a, b, beta1){
  d= (a^2+ a*b - beta1*a^2)/(a + b)
  return(d)
}
# dFunc(1,2,3)
#
# eFunc[a_, b_, beta1_, K_] = (beta0 a b)/((a + b) (K - 1))
#
eFunc=function(a, b, beta0, K){
  e=(beta0*a*b)/((a + b)*(K - 1))
  return(e)
}
# eFunc(1,2,3,4)
#
# fFunc[a_, b_, beta1_, K_] = (beta0 a b + b^2)/((a + b) (K - 1))
fFunc=function(a, b, beta0, K){
  f= (beta0*a*b + b^2)/((a + b)*(K - 1))
  return(f)
}

sigmaSqABC=function(K,
                    piA, piB, piC,
                    beta1, beta0,
                    aA, bA, cA, dA, eA, fA,
                    n){
  sigmaSqABC=

    1/((1/1/beta0^2*((eA + ((n - piB*n)/(K - 1)*beta0*piA+(n - piC*n)/(K - 1)*beta0*piA))*((2*n - piB*n - piC*n)/2 - ((n - piB*n)/(K - 1)*beta0*piA+(n - piC*n)/(K - 1)*beta0*piA) +
                                                                                             fA))/((eA + fA + (2*n - piB*n - piC*n)/2)^2*(eA + fA + (2*n - piB*n - piC*n)/2 +
                                                                                                                                            1))) +
         1/(((aA + piA*n)*(bA + n - piA*n))/((aA + bA + n)^2*(aA + bA + n + 1)) )+
         1/(1/beta1^2*((cA + piA*n*piA*beta1)*(dA + piA*n - piA*n*piA*beta1))/((cA + dA + piA*n)^2*(cA + dA + piA*n + 1))))

  return(sigmaSqABC)
}

normal_cdf=function(x,mu,sigmaSq){
  normal_cdf=1/2*(1+pracma::erf((x-mu)/(sqrt(sigmaSq*2))))
  return(normal_cdf)
}

# normal_cdf(0,1,Inf)

order1_pdf=function(mu1,sigmaSq1,mu2,sigmaSq2,mu3,sigmaSq3,x){
  order1_pdf=dnorm(x,mean=mu1,sd=sqrt(sigmaSq1))*normal_cdf(x,mu2,sigmaSq2)*normal_cdf(x,mu3,sigmaSq3)+
    dnorm(x,mean=mu2,sd=sqrt(sigmaSq2))*normal_cdf(x,mu1,sigmaSq1)*normal_cdf(x,mu3,sigmaSq3)+
    dnorm(x,mean=mu3,sd=sqrt(sigmaSq3))*normal_cdf(x,mu1,sigmaSq1)*normal_cdf(x,mu2,sigmaSq2)
  return(order1_pdf)
}

order2_pdf=function(mu1,sigmaSq1,mu2,sigmaSq2,mu3,sigmaSq3,x){
  order2_pdf=(1-normal_cdf(x,mu3,sigmaSq3))*(normal_cdf(x,mu1,sigmaSq1)*dnorm(x,mean=mu2,sd=sqrt(sigmaSq2))+normal_cdf(x,mu2,sigmaSq2)*dnorm(x,mean=mu1,sd=sqrt(sigmaSq1)))+
    (1-normal_cdf(x,mu2,sigmaSq2))*(normal_cdf(x,mu1,sigmaSq1)*dnorm(x,mean=mu3,sd=sqrt(sigmaSq3))+normal_cdf(x,mu3,sigmaSq3)*dnorm(x,mean=mu1,sd=sqrt(sigmaSq1)))+
    (1-normal_cdf(x,mu1,sigmaSq1))*(normal_cdf(x,mu2,sigmaSq2)*dnorm(x,mean=mu3,sd=sqrt(sigmaSq3))+normal_cdf(x,mu3,sigmaSq3)*dnorm(x,mean=mu2,sd=sqrt(sigmaSq2)))
  return(order2_pdf)
}

prod_o1d=function(mu1,sigmaSq1,mu2,sigmaSq2,mu3,sigmaSq3,x,d){
  prod_o1d=order1_pdf(mu1,sigmaSq1,mu2,sigmaSq2,mu3,sigmaSq3,x+d)*d
  return(prod_o1d)
}

prod_o1dd=function(mu1,sigmaSq1,mu2,sigmaSq2,mu3,sigmaSq3,x,d){
  prod_o1dd=order1_pdf(mu1,sigmaSq1,mu2,sigmaSq2,mu3,sigmaSq3,x+d)*d*d
  return(prod_o1dd)
}

o1d_integral=function(mu1,sigmaSq1,mu2,sigmaSq2,mu3,sigmaSq3,x){
  o1d_integral=cubature::hcubature(prod_o1d,mu1=mu1,sigmaSq1=sigmaSq1,
                  mu2=mu2,sigmaSq2=sigmaSq2,mu3=mu3,sigmaSq3=sigmaSq3,x=x,lowerLimit = 0, upperLimit = Inf)$integral
  return(o1d_integral)
}

o1dd_integral=function(mu1,sigmaSq1,mu2,sigmaSq2,mu3,sigmaSq3,x){
  o1dd_integral=cubature::hcubature(prod_o1dd,mu1=mu1,sigmaSq1=sigmaSq1,
                         mu2=mu2,sigmaSq2=sigmaSq2,mu3=mu3,sigmaSq3=sigmaSq3,x=x,lowerLimit = 0, upperLimit = Inf)$integral
  return(o1dd_integral)
}

diff_o1o2_mean_integrand=function(mu1,sigmaSq1,mu2,sigmaSq2,mu3,sigmaSq3,x){
  diff_o1o2_mean_integrand=order2_pdf(mu1=mu1,sigmaSq1=sigmaSq1,
                                      mu2=mu2,sigmaSq2=sigmaSq2,
                                      mu3=mu3,sigmaSq3=sigmaSq3,x=x)*o1d_integral(mu1=mu1,
                                                                              sigmaSq1=sigmaSq1,
                                                                             mu2=mu2,sigmaSq2=sigmaSq2,
                                                                             mu3=mu3,sigmaSq3=sigmaSq3,x=x)
  return(diff_o1o2_mean_integrand)
}

diff_o1o2_variance1_integrand=function(mu1,sigmaSq1,mu2,sigmaSq2,mu3,sigmaSq3,x){
  diff_o1o2_variance1_integrand=order2_pdf(mu1=mu1,sigmaSq1=sigmaSq1,
                                           mu2=mu2,sigmaSq2=sigmaSq2,
                                           mu3=mu3,sigmaSq3=sigmaSq3,x=x)*o1dd_integral(mu1=mu1,sigmaSq1=sigmaSq1,
                                        mu2=mu2,sigmaSq2=sigmaSq2,
                                        mu3=mu3,sigmaSq3=sigmaSq3,x=x)
  return(diff_o1o2_variance1_integrand)
}

mean_o1o2_diff=function(K,
                        piA, piB, piC,
                        beta1, beta0,
                        aA, bA, cA, dA, eA, fA,
                        aB, bB, cB, dB, eB, fB,
                        aC, bC, cC, dC, eC, fC,
                        n){
  sigmaSq1=sigmaSqABC(K,
                      piA, piB, piC,
                      beta1, beta0,
                      aA, bA, cA, dA, eA, fA,
                      n)
  sigmaSq2=sigmaSqABC(K,
                      piB, piA, piC,
                      beta1, beta0,
                      aB, bB, cB, dB, eB, fB,
                      n)
  sigmaSq3=sigmaSqABC(K,
                      piC, piA, piB,
                      beta1, beta0,
                      aC, bC, cC, dC, eC, fC,
                      n)
  mean_o1o2_diff=cubature::hcubature(diff_o1o2_mean_integrand, mu1=piA,sigmaSq1=sigmaSq1,
                            mu2=piB,sigmaSq2=sigmaSq2,mu3=piC,sigmaSq3=sigmaSq3,lowerLimit = -Inf, upperLimit =Inf)$integral

  return(mean_o1o2_diff)

}

var1_o1o2_diff=function(K,
                       piA, piB, piC,
                       beta1, beta0,
                       aA, bA, cA, dA, eA, fA,
                       aB, bB, cB, dB, eB, fB,
                       aC, bC, cC, dC, eC, fC,
                       n){
  sigmaSq1=sigmaSqABC(K,
                      piA, piB, piC,
                      beta1, beta0,
                      aA, bA, cA, dA, eA, fA,
                      n)
  sigmaSq2=sigmaSqABC(K,
                      piB, piA, piC,
                      beta1, beta0,
                      aB, bB, cB, dB, eB, fB,
                      n)
  sigmaSq3=sigmaSqABC(K,
                      piC, piA, piB,
                      beta1, beta0,
                      aC, bC, cC, dC, eC, fC,
                      n)
  # o1_o2_mean=mean_o1o2_diff(K,
  #                           piA, piB, piC,
  #                           beta1, beta0,
  #                           aA, bA, cA, dA, eA, fA,
  #                           aB, bB, cB, dB, eB, fB,
  #                           aC, bC, cC, dC, eC, fC,
  #                           n)

  var1_o1o2_diff=cubature::hcubature(diff_o1o2_variance1_integrand,mu1=piA,sigmaSq1=sigmaSq1,
                               mu2=piB,sigmaSq2=sigmaSq2,mu3=piC,sigmaSq3=sigmaSq3, lowerLimit = -Inf, upperLimit = Inf)$integral


  return(var1_o1o2_diff)

}

# sample_size_equation <- function (K,
#                                   piA, piB, piC,
#                                   beta1, beta0,
#                                   aA, bA, cA, dA, eA, fA,
#                                   aB, bB, cB, dB, eB, fB,
#                                   aC, bC, cC, dC, eC, fC,
#                                   ciL,
#                                   n) {
#   2*ciZ*sqrt(var1_o1o2_diff(K,
#                              piA, piB, piC,
#                              beta1, beta0,
#                              aA, bA, cA, dA, eA, fA,
#                              aB, bB, cB, dB, eB, fB,
#                              aC, bC, cC, dC, eC, fC,
#                              n))-ciL
# }



