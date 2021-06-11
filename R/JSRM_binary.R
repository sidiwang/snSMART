#'
#'
#'
#'
#'

# JSRM

JSRM_binary = function(data, six = TRUE){

  # data, same format as the bjsm_binary.R file trial dataset format
  # six, if TRUE, will run the six beta model, if FALSE will run the two beta model


  mydata = data
  mydata$disc <- 2 * mydata$treatment_stageI - (mydata$response_stageI == 0)

  Y <- c(mydata$response_stageI, mydata$response_stageII)
  XA1 <- as.numeric(mydata$treatment_stageI==1)
  XB1 <- as.numeric(mydata$treatment_stageI==2)
  XC1 <- as.numeric(mydata$treatment_stageI==3)
  XA2 <- as.numeric(mydata$treatment_stageII==1)
  XB2 <- as.numeric(mydata$treatment_stageII==2)
  XC2 <- as.numeric(mydata$treatment_stageII==3)
  XA <- c(XA1, XA2)
  XB <- c(XB1, XB2)
  XC <- c(XC1, XC2)
  Z1A <- c(rep(0, nrow(mydata)), ifelse(mydata$treatment_stageI==1,mydata$response_stageI,0))
  Z2A <- c(rep(0, nrow(mydata)), ifelse(mydata$treatment_stageI==1,1-mydata$response_stageI,0))
  Z1B <- c(rep(0, nrow(mydata)), ifelse(mydata$treatment_stageI==2,mydata$response_stageI,0))
  Z2B <- c(rep(0, nrow(mydata)), ifelse(mydata$treatment_stageI==2,1-mydata$response_stageI,0))
  Z1C <- c(rep(0, nrow(mydata)), ifelse(mydata$treatment_stageI==3,mydata$response_stageI,0))
  Z2C <- c(rep(0, nrow(mydata)), ifelse(mydata$treatment_stageI==3,1-mydata$response_stageI,0))
  ptid <- rep(1:nrow(mydata), 2)

  if (six == TRUE){
    geedata <- data.frame(ptid, XA, XB, XC, Z1A, Z2A, Z1B, Z2B, Z1C, Z2C, Y)
    geedata <- geedata[order(geedata$ptid),]
    rm(ptid, Y, Z1A, Z2A, Z1B, Z2B, Z1C, Z2C, XA, XB, XC)
    try({
      mod1 <- geepack::geeglm(Y~XA+XB+XC+Z1A+Z2A+Z1B+Z2B+Z1C+Z2C-1, family=poisson(link="log"), data=geedata, id=ptid,corstr = "independence")
      beta_hat <- mod1$coefficients[1:3];
      sd_beta_hat <- summary(mod1)$coef[1:3,2];
      pi_hat <- rbind(pi_hat,exp(beta_hat));
      sd_pi_hat <- rbind(sd_pi_hat,exp(beta_hat)*sd_beta_hat);
      b_hat <- coef(mod1);
      grad <- exp(rbind(c(2*b_hat[1]+b_hat[4], b_hat[2]+b_hat[5], b_hat[1]+b_hat[2]+b_hat[5]),
                        c(2*b_hat[1]+b_hat[4], b_hat[3]+b_hat[5], b_hat[1]+b_hat[3]+b_hat[5]),
                        c(2*b_hat[2]+b_hat[6], b_hat[1]+b_hat[7], b_hat[2]+b_hat[1]+b_hat[7]),
                        c(2*b_hat[2]+b_hat[6], b_hat[3]+b_hat[7], b_hat[2]+b_hat[3]+b_hat[7]),
                        c(2*b_hat[3]+b_hat[8], b_hat[1]+b_hat[9], b_hat[3]+b_hat[1]+b_hat[9]),
                        c(2*b_hat[3]+b_hat[8], b_hat[2]+b_hat[9], b_hat[3]+b_hat[2]+b_hat[9])));
      pi_DTR_hat <- rbind(pi_DTR_hat, c(1, 1, -1) %*% t(grad))
      grad[,3] <- -grad[,3]
      sigma_b <- mod1$geese$vbeta
      L1 <- matrix(c(2, 0, 0, 1, 0, 0, 0, 0, 0,
                     0, 1, 0, 0, 1, 0, 0, 0, 0,
                     1, 1, 0, 0, 1, 0, 0, 0, 0), nrow=3, ncol=9, byrow=T)
      sigma_g <- L1 %*% sigma_b %*% t(L1)
      seAB <- sqrt(grad[1,] %*% sigma_g %*% grad[1,])

      L2 <- matrix(c(2, 0, 0, 1, 0, 0, 0, 0, 0,
                     0, 0, 1, 0, 1, 0, 0, 0, 0,
                     1, 0, 1, 0, 1, 0, 0, 0, 0), nrow=3, ncol=9, byrow=T)
      sigma_g <- L2 %*% sigma_b %*% t(L2)
      seAC <- sqrt(grad[2,] %*% sigma_g %*% grad[2,])

      L3 <- matrix(c(0, 2, 0, 0, 0, 1, 0, 0, 0,
                     1, 0, 0, 0, 0, 0, 1, 0, 0,
                     1, 1, 0, 0, 0, 0, 1, 0, 0), nrow=3, ncol=9, byrow=T)
      sigma_g <- L3 %*% sigma_b %*% t(L3)
      seBA <- sqrt(grad[3,] %*% sigma_g %*% grad[3,])

      L4 <- matrix(c(0, 2, 0, 0, 0, 1, 0, 0, 0,
                     0, 0, 1, 0, 0, 0, 1, 0, 0,
                     0, 1, 1, 0, 0, 0, 1, 0, 0), nrow=3, ncol=9, byrow=T)
      sigma_g <- L4 %*% sigma_b %*% t(L4)
      seBC <- sqrt(grad[4,] %*% sigma_g %*% grad[4,])

      L5 <- matrix(c(0, 0, 2, 0, 0, 0, 0, 1, 0,
                     1, 0, 0, 0, 0, 0, 0, 0, 1,
                     1, 0, 1, 0, 0, 0, 0, 0, 1), nrow=3, ncol=9, byrow=T)
      sigma_g <- L5 %*% sigma_b %*% t(L5)
      seCA <- sqrt(grad[5,] %*% sigma_g %*% grad[5,])

      L6 <- matrix(c(0, 0, 2, 0, 0, 0, 0, 1, 0,
                     0, 1, 0, 0, 0, 0, 0, 0, 1,
                     0, 1, 1, 0, 0, 0, 0, 0, 1), nrow=3, ncol=9, byrow=T)
      sigma_g <- L6 %*% sigma_b %*% t(L6)
      seCB <- sqrt(grad[6,] %*% sigma_g %*% grad[6,])
      pi_DTR_se <- rbind(pi_DTR_se, c(seAB, seAC, seBA, seBC, seCA, seCB))
    })

  } else {
    Z1 = Z1A + Z1B + Z1C
    Z2 = Z2A + Z2B + Z2C
    geedata <- data.frame(ptid, XA, XB, XC, Z1, Z2, Y)
    geedata <- geedata[order(geedata$ptid),]
    rm(ptid, Y, Z1A, Z2A, Z1B, Z2B, Z1C, Z2C, XA, XB, XC, Z1, Z2)
    try({
      mod1 <- geepack::geeglm(Y~XA+XB+XC+Z1+Z2-1, family=poisson(link="log"), data=geedata, id=ptid,corstr = "independence")
      beta_hat <- mod1$coefficients[1:3];
      sd_beta_hat <- summary(mod1)$coef[1:3,2];
      pi_hat <- rbind(pi_hat,exp(beta_hat));
      sd_pi_hat <- rbind(sd_pi_hat,exp(beta_hat)*sd_beta_hat);
      b_hat <- coef(mod1);
      grad <- exp(rbind(c(2*b_hat[1]+b_hat[4], b_hat[2]+b_hat[5], b_hat[1]+b_hat[2]+b_hat[5]),
                        c(2*b_hat[1]+b_hat[4], b_hat[3]+b_hat[5], b_hat[1]+b_hat[3]+b_hat[5]),
                        c(2*b_hat[2]+b_hat[4], b_hat[1]+b_hat[5], b_hat[2]+b_hat[1]+b_hat[5]),
                        c(2*b_hat[2]+b_hat[4], b_hat[3]+b_hat[5], b_hat[2]+b_hat[3]+b_hat[5]),
                        c(2*b_hat[3]+b_hat[4], b_hat[1]+b_hat[5], b_hat[3]+b_hat[1]+b_hat[5]),
                        c(2*b_hat[3]+b_hat[4], b_hat[2]+b_hat[5], b_hat[3]+b_hat[2]+b_hat[5])));
      pi_DTR_hat <- rbind(pi_DTR_hat, c(1, 1, -1) %*% t(grad))
      grad[,3] <- -grad[,3]
      sigma_b <- mod1$geese$vbeta
      L1 <- matrix(c(2, 0, 0, 1, 0,
                     0, 1, 0, 0, 1,
                     1, 1, 0, 0, 1), nrow=3, ncol=5, byrow=T)
      sigma_g <- L1 %*% sigma_b %*% t(L1)
      seAB <- sqrt(grad[1,] %*% sigma_g %*% grad[1,])

      L2 <- matrix(c(2, 0, 0, 1, 0,
                     0, 0, 1, 0, 1,
                     1, 0, 1, 0, 1), nrow=3, ncol=5, byrow=T)
      sigma_g <- L2 %*% sigma_b %*% t(L2)
      seAC <- sqrt(grad[2,] %*% sigma_g %*% grad[2,])

      L3 <- matrix(c(0, 2, 0, 1, 0,
                     1, 0, 0, 0, 1,
                     1, 1, 0, 0, 1), nrow=3, ncol=5, byrow=T)
      sigma_g <- L3 %*% sigma_b %*% t(L3)
      seBA <- sqrt(grad[3,] %*% sigma_g %*% grad[3,])

      L4 <- matrix(c(0, 2, 0, 1, 0,
                     0, 0, 1, 0, 1,
                     0, 1, 1, 0, 1), nrow=3, ncol=5, byrow=T)
      sigma_g <- L4 %*% sigma_b %*% t(L4)
      seBC <- sqrt(grad[4,] %*% sigma_g %*% grad[4,])

      L5 <- matrix(c(0, 0, 2, 1, 0,
                     1, 0, 0, 0, 1,
                     1, 0, 1, 0, 1), nrow=3, ncol=5, byrow=T)
      sigma_g <- L5 %*% sigma_b %*% t(L5)
      seCA <- sqrt(grad[5,] %*% sigma_g %*% grad[5,])

      L6 <- matrix(c(0, 0, 2, 1, 0,
                     0, 1, 0, 0, 1,
                     0, 1, 1, 0, 1), nrow=3, ncol=5, byrow=T)
      sigma_g <- L6 %*% sigma_b %*% t(L6)
      seCB <- sqrt(grad[6,] %*% sigma_g %*% grad[6,])
      pi_DTR_se <- rbind(pi_DTR_se, c(seAB, seAC, seBA, seBC, seCA, seCB))
    })
  }

  result = list("GEE_output" = mod1, "beta_hat" = beta_hat, "sd_beta_hat" = sd_beta_hat, "pi_hat" = pi_hat, "sd_pi_hat" = sd_pi_hat, "pi_DTR_hat" = pi_DTR_hat,
                "pi_DTR_se" = pi_DTR_se)


  return(result)
}
