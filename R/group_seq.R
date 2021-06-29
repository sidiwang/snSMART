#'
#'
#'
#'
#'
#'
#'
#'
#'
#' data generation function add final analysis (below functions should add to both data_gen_group_seq_1step.R and .....2step.R)
#' find a way to write functions on this page as a pure real trial analysis function
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
group_seq = function(data, rule.type = 1){

  mydata = data
  mydata$trt.2nd <- unlist(mydata$trt.2nd)
  stage1_count <- rbind(stage1_count,table(mydata$trt.1st))
  stage2_count <- rbind(stage2_count,table(mydata$trt.2nd))
  rand_prob_output <- rbind(rand_prob_output,outcome[[2]])
  dropped_arm <- c(dropped_arm,outcome[[3]])
  dropped_look <- c(dropped_look,outcome[[4]])
  mydata$disc <- 2 * mydata$trt.1st - (mydata$resp.1st == 0)
  error_ind <- 0
  tryCatch({
    jags <- jags.model(file.path(file_path,jags.model.name),
                       data=list(n = nrow(mydata),
                                 num_arms = NUM_ARMS,
                                 Y1 = mydata$resp.1st,
                                 Y2 = mydata$resp.2nd,
                                 treatment_stageI = mydata$trt.1st,
                                 treatment_stageII = mydata$trt.2nd,
                                 response_stageI_disc = mydata$disc,
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
    print(i)     # show the number of iterations run
    print(iter)
  }
  )
  out_post <- posterior_sample[[1]]
  response_rate_bayes <- rbind(response_rate_bayes,apply(out_post[,7:9],2,mean))
}

mean_rr <- apply(response_rate_bayes,2,mean)
se_rr <- apply(response_rate_bayes,2,sd)
expected <- c(pi_1A,pi_1B,pi_1C)
result_table <- data.frame(response_rate = mean_rr,
                           se_rr = se_rr,
                           expected = expected)
result_table$bias2 <- result_table$response_rate - result_table$expected
result_table$rMSE2 <- sqrt(result_table$bias2^2+result_table$se_rr^2)
result_table <- format(round(result_table,3),3)
assignment <- cbind(apply(stage1_count,2,mean),apply(stage2_count,2,mean))

CI_matrix <- matrix(c(paste0("(",paste(specify_decimal(quantile(stage1_count[,1],probs = c(0.025,0.975)),2),collapse=","),")"),
                      paste0("(",paste(specify_decimal(quantile(stage1_count[,2],probs = c(0.025,0.975)),2),collapse=","),")"),
                      paste0("(",paste(specify_decimal(quantile(stage1_count[,3],probs = c(0.025,0.975)),2),collapse=","),")"),
                      paste0("(",paste(specify_decimal(quantile(stage2_count[,1],probs = c(0.025,0.975)),2),collapse=","),")"),
                      paste0("(",paste(specify_decimal(quantile(stage2_count[,2],probs = c(0.025,0.975)),2),collapse=","),")"),
                      paste0("(",paste(specify_decimal(quantile(stage2_count[,3],probs = c(0.025,0.975)),2),collapse=","),")")),
                    nrow=3,ncol=2,byrow=F)
result_table$s1_count <- assignment[,1]
result_table$s1_count_CI <- CI_matrix[,1]
result_table$s2_count <- assignment[,2]
result_table$s2_count_CI <- CI_matrix[,2]
rownames(result_table) <- c("A","B","C")

result_table_tab <- rbind(result_table_tab,result_table)

not_largest_pi <- which(c(pi_1A,pi_1B,pi_1C)!=max(c(pi_1A,pi_1B,pi_1C)))
smallest_pi <- which(c(pi_1A,pi_1B,pi_1C)==min(c(pi_1A,pi_1B,pi_1C)))
if (length(smallest_pi)!=3){
  good_drop <- as.numeric(dropped_arm %in% not_largest_pi)
  best_drop <- as.numeric(dropped_arm %in% smallest_pi)
} else {
  good_drop <- NA
  best_drop <- NA
}
dropped_look_tab <- matrix(NA,nrow=n.sim,ncol=n.update)
for(j in 1:n.update){
  dropped_look_tab[,j] <- (dropped_look==j)
}
operating_char <- c(paste(drop_threshold.large.pool[[iter]],collapse=","),
                    paste(drop_threshold.small.pool[[iter]],collapse=","),
                    apply(dropped_look_tab,2,mean),mean(dropped_arm!=0),
                    specify_decimal(sum(good_drop)/sum(dropped_arm!=0),2),
                    specify_decimal(sum(best_drop)/sum(dropped_arm!=0),2))
operating_char_tab <- rbind(operating_char_tab,operating_char)
colnames(operating_char_tab) <- c("threshold_largest","threshold_smallest",sprintf("dropping_rate_look_%d",1:n.update),"total_dropping_rate","good_drop_rate","best_drop_rate")
}





