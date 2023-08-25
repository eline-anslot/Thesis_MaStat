
binom_gen <- function(x, n, py, px) {
  if (py <= px) {   
    y <- rbinom(n, 1, x*py/px)}
  else {
    y <- as.numeric( !rbinom(n, 1, (1-x)*(1-py)/(1-px)) )
  }  
  y
}

#############################################################
#########Stage 1########
#############################################################

#-----------------------------------------------------------#
#--------------------#Standard method#----------------------#
#-----------------------------------------------------------#
standard_method <- function(data, 
                            outcome,
                            recruitment_acc,
                            time_long,
                            n1, 
                            r1){
  
  #time at which participants have completed the long-term endpoint y
  time_ppt_long <- recruitment_acc + time_long
  
  #time at which the last n1 participant completed they long-term endpoint y
  time_n1_ppt <- time_ppt_long[n1]
  
  ###Cohort1####
  # Return true or false depending on whether participants have completed 
  # the long-term endpoint at stage 1 
  
  #indices_long <- time_ppt_long <= time_n1_ppt
  indices_long <- time_ppt_long <= time_n1_ppt & 1:length(recruitment_acc) <= n1
  
  #Subgroup of the sample at interim which have all information 
  #(baseline covariates, short_term and long-term endpoints)
  cohort1_l <- data[indices_long, outcome]
  
  data.frame(decision = ifelse(sum(cohort1_l) > r1, 
                               "continue", 
                               "stop"), 
             n_long = length(cohort1_l)
  )
}

#-----------------------------------------------------------#
#-----------------------#Proposed method#-------------------------#
#-----------------------------------------------------------#

to_predictors <- function(var){
  paste(var, collapse = "+")
}

to_formula <- function(outcome, predictor){
  as.formula(paste(outcome, to_predictors(predictor), sep = "~"))
}

VL_method <- function(data,
                      outcome,
                      n1,
                      recruitment_acc, 
                      time_long, 
                      covariate_vec = NULL,
                      short_vec = NULL,
                      time_short = NULL, 
                      adjustment = c("baseline&short",
                                     "short",
                                     "baseline",
                                     "nothing")) {
  
  #################################################
  #######Sample at interim#######
  #################################################
  #time at which participants have completed the long-term endpoint y
  time_ppt_long <- recruitment_acc + time_long
  
  #time at which the last n1 participant completed they long-term endpoint y
  time_n1_ppt <- time_ppt_long[n1]
  
  ###Cohort1####
  # Return true or false depending on whether participants have completed 
  # the long-term endpoint at stage 1 
  indices_long <- time_ppt_long <= time_n1_ppt & 1:length(recruitment_acc) <= n1
  #indices_long <- time_ppt_long <= time_n1_ppt
  
  #Subgroup of the sample at interim which have all information 
  #(baseline covariates, short_term and long-term endpoints)
  cohort1_l <- data[indices_long, ]
  
  if (adjustment == "baseline&short") {
    #The baseline covariates that we will take into account
    cohort_cov_var <- covariate_vec
    
    #time at which participants have completed the short-term endpoint x
    time_ppt_short <- recruitment_acc + time_short
    
    # Return true or false depending on whether we have
    # information on the short-term endpoint at stage 1 
    indices_short <- time_ppt_short <= time_n1_ppt
    
    # Return true or false depending on whether we have
    # information on the baseline covariates at stage 1 
    indices_bascov <- recruitment_acc <= time_n1_ppt
    
    # Return true or false depending on whether we have
    # information on the baseline covariates and short term but not
    # of the long term endpoint at stage 1 ==> Cohort 2
    indices_cohort2 <- indices_short &! indices_long
    
    # Return true or false depending on whether we have
    # information on the baseline covariates but not on the short-term and 
    # on the long term endpoint at stage 1 => Cohort 3
    indices_cohort3 <- indices_bascov &! indices_short
    
    #The baseline covariates and short term endpoint variables that we will take
    #into account
    cohort_short_var <- c(covariate_vec, short_vec)
    
    #####Subgroups######
    ###Cohort2####
    #Participants with baseline covariates and short-term endpoint but no 
    #long-term endpoint
    cohort2_s <- data[indices_cohort2, cohort_short_var, drop = FALSE]
    
    ###Cohort3###
    #Participants with only baseline covariates information
    cohort3_bc     <- data[indices_cohort3, cohort_cov_var, drop = FALSE]
    
    #################################################
    #######Estimate#######
    #################################################
    model   <- glm(to_formula(outcome, cohort_short_var),
                   family = binomial(link = "logit"),
                   data = cohort1_l)
    
    data_cohort1_2 <- rbind(cohort1_l[, cohort_short_var, drop = FALSE], 
                            cohort2_s)
    
    predict2 <- predict(model, 
                        type = "response",
                        newdata = data_cohort1_2)
    
    model2 <- glm(to_formula("predict2", cohort_cov_var),
                  family = binomial(link="logit"),
                  data = cbind(data_cohort1_2, predict2))
    
    predict <- predict(model2,
                       type = "response",
                       newdata = rbind(data_cohort1_2[, cohort_cov_var, drop = FALSE],
                                       cohort3_bc))
    
    estimate <- mean(predict)
    
    #################################################
    #######Variance#######
    #################################################
    new_indices_long <- indices_long[indices_bascov]
    new_indices_short <- indices_short[indices_bascov]
    cycx <- new_indices_long/(mean(new_indices_long))#*mean(new_indices_short))
    predict2_adj <- c(predict2, rep(0, times = sum(indices_bascov) - length(predict2)))
    yi_yihat <- data[indices_bascov, outcome] - predict2_adj
    predict_adj <- c(predict, rep(0, times = sum(indices_bascov) - length(predict)))
    cix_yihat_yihat2 <- (new_indices_short/mean(new_indices_short) ) *(predict2_adj - predict_adj)
    variance <- var( (cycx*yi_yihat + cix_yihat_yihat2 + predict_adj) )
    
    #################################################
    #######Results#######
    #################################################
    
    results <- list(estimate = estimate,
                    n_interim = sum(indices_cohort2) + 
                      sum(indices_cohort3) + 
                      sum(indices_long),
                    prediction = predict,
                    variance = variance)
    results
    
    ##CHECK##
    # predict2_check <- predict(model,
    #                           type = "response",
    #                           newdata=cohort2_s)
    # 
    # new_outcome <- do.call("c", list(cohort1_l[[outcome]], as.vector(predict2_check)))
    # 
    # model_2_check <- glm(to_formula("new_outcome", cohort_cov_var),
    #                      family = binomial(link="logit"),
    #                      data = cbind(data_cohort1_2, new_outcome))
    # 
    # predict3_check <- predict(model_2_check,
    #                           type = "response",
    #                           newdata=cohort3_bc)
    # 
    # estimate_check <- mean(c(cohort1_l[, outcome],
    #                          predict2_check,
    #                          predict3_check))
    #c(estimate, estimate_check)
    
  } else if (adjustment == "baseline") {
    #The baseline covariates that we will take into account
    cohort_cov_var <- covariate_vec
    
    #####Subgroups######
    # Return true or false depending on whether we have
    # information on the baseline covariates at stage 1
    indices_bascov <- recruitment_acc <= time_n1_ppt
    
    # Return true or false depending on whether we have
    # information on the baseline covariates but
    # not on the long term endpoint at stage 1 => Cohort 2
    indices_cohort2 <- indices_bascov &! indices_long
    
    ###Cohort2####
    #Subgroup of the sample at interim which only have baseline covariate information
    cohort2_bc     <- as.data.frame(data[indices_cohort2, cohort_cov_var])
    colnames(cohort2_bc) <- cohort_cov_var
    
    #################################################
    #######Estimate#######
    #################################################
    model   <- glm(to_formula(outcome, cohort_cov_var),
                   family = binomial(link="logit"),
                   data = cohort1_l)
    
    predict <- predict(model,
                       type = "response",
                       newdata = rbind(cohort1_l[, cohort_cov_var, drop = FALSE],
                                       cohort2_bc))
    
    estimate <- mean(predict)
    
    n <- dim(data)[1]
    predict <- c(predict, rep(0, times = sum(indices_bascov) - length(predict)))
    yi_yihat <- data[indices_bascov, outcome] - predict
    variance <- var( indices_long[indices_bascov]/mean(indices_long[indices_bascov])*yi_yihat + predict )
    
    results <- list(estimate = estimate,
                    n_interim = sum(indices_cohort2) + sum(indices_long),
                    prediction = predict,
                    variance = variance)
    results

    ##CHECK##
    # predict_cohort2 <- predict(model,
    #                            type = "response",
    #                            newdata=cohort2_bc)
    #
    # estimate_check <- mean(c(cohort1_l[, outcome],
    #                          predict_cohort2))
    #
    # c(estimate, estimate_check)
    } else if (adjustment == "short") {
    
    #time at which participants have completed the short-term endpoint x
    time_ppt_short <- recruitment_acc + time_short
    
    # Return true or false depending on whether we have
    # information on the short-term endpoint at stage 1 
    indices_short <- time_ppt_short <= time_n1_ppt
    
    # Return true or false depending on whether we have
    # information on the baseline covariates and short term but not
    # of the long term endpoint at stage 1 ==> Cohort 2
    indices_cohort2 <- indices_short &! indices_long
    
    #####Subgroups######
    ###Cohort2####
    #Participants with baseline covariates and short-term endpoint but no 
    #long-term endpoint
    cohort2_s <- data[indices_cohort2, short_vec, drop = FALSE]
    
    #################################################
    #######Estimate#######
    #################################################
    
    model <- glm(to_formula(outcome, short_vec), 
                 family = binomial(link = "logit"),
                 data=cohort1_l)
    
    predict <- predict(model, 
                       type = "response",
                       newdata=rbind(cohort1_l[, short_vec, drop=FALSE],
                                            cohort2_s))
    
    estimate <- mean(predict)   
    predict <- c(predict, rep(0, times = sum(indices_short) - length(predict)))
    yi_yihat <- data[indices_short, outcome] - predict
    variance <- var( indices_long[indices_short]/mean(indices_long[indices_short])*yi_yihat + predict )
    
    
    results <- list(estimate = estimate, 
                   variance = variance,
                   prediction = predict,
                   n_interim = nrow(cohort1_l) + nrow(cohort2_s))
    results
    
  } else {
    estimate <- mean(cohort1_l[, outcome])
    variance <- estimate*(1 - estimate)
    results <- list(estimate = estimate,
                    variance = variance,
                    n_interim = length(cohort1_l))
    results
  }
}

#### Adjustment ####
VL_deci <- function(variance,
                    estimate,
                    n_interim,
                    p0,
                    PET) {
  
  se <- sqrt(variance/n_interim)
  t1 <- (estimate - p0)/se
  ifelse(t1 <= qnorm(PET), "stop", "continue")
}

#-----------------------------------------------------------#
#-----------------------------------------------------------#

#############################################################
#########Stage 2########
#############################################################

decision_stage2 <- function(outcome, r2){
  ifelse(sum(outcome) > r2, 
         "reject H0", 
         "not reject H0")}
