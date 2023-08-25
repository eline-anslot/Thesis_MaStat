library(goldilocks)

source("Functions2.R")

sim <- function(
    c, b0_x, b0_y0, b1_y1_x, b1_y0, b2_y0, b2_y1,
    p0, p1, n, n1, r1, r2, PET, b0_y1, n_sim, directory_name
  ) {

  b1_y1 <- b1_y1_x
  b1_x <- b1_y1_x
  
  alpha <- .05
  power_start <- .80 

  lambda <- 2

  p_LTE_0 <- p0
  p_LTE_1 <- p1
  end_n149 <- data.frame(lambda = c(0.5, 1, 2, 3, 4),
                         mean = c(295.7055, 147.983, 73.88, 49.1437, 36.8025),
                         sd = c(24.28236, 12.05218, 6.049154, 4.071377, 3.098944))
  
  max_rec <- end_n149[end_n149$lambda == lambda, "mean"]
  sd_rec <- end_n149[end_n149$lambda == lambda, "sd"]

  #Time in month after entering the study for when the long-term endpoint is reached 
  months_long <- 12*4

  #Time in months after entering the study for when the short-term endpoint is reached
  months_short <- 6*4

  #Data from the simulation
  sim_data <- data.frame(
                        Interim_method = character(0),
                        n              = numeric(0),
                        n1             = numeric(0),
                        n_interim      = numeric(0),
                        Lambda         = numeric(0),
                        Decision_underH0_st1 = character(0),
                        Decision_underH0_st2 = character(0),
                        Decision_underH1_st1 = character(0),
                        Decision_underH1_st2 = character(0),
                        Alpha = numeric(0),
                        Power = numeric(0))

  #The different methods
  interim_methods <- c("Standard", "Proposed")
  n_methods <- length(interim_methods)

  index <- 1

  for (j in 1:n_sim) {
    #indexation for saving the data in sim_data 
    indexation <- index:(index + n_methods - 1)
    
    ####Baseline covariate#####
    # z <- sample(data_start$z_start, n, replace = TRUE) #baseline covariate: e.g., normalized age
    # 
    # ####Short-term endpoint####
    # x0 <- sample(data_start$x_0_start, n, replace = TRUE)
    # x1 <- sample(data_start$x_1_start, n, replace = TRUE)
    z <- rnorm(n) #baseline covariate: e.g., normalized age
    
    ####Short-term endpoint####
    # prog_score_x_0 <- -1 + 2*z
    # expit_x_0 <- 1/(1 + exp(-prog_score_x_0))
    # x0 <- rbinom(n, 1, expit_x_0)
    # 
    # prog_score_x_1 <- -1 + 2*z
    # expit_x_1 <- 1/(1 + exp(-prog_score_x_1))
    # x1 <- rbinom(n, 1, expit_x_1)
    prog_score_x <- b0_x + b1_x*z
    expit_x <- 1/(1 + exp(-prog_score_x))
    x <- rbinom(n, 1, expit_x)
    
    ####Long-term endpoint####
    prog_score_y0 <- b0_y0 + c*b1_y0*z + c*b2_y0*x
    expit_y0 <- 1/(1 + exp(-prog_score_y0))
    outcome_H0 <- rbinom(n, 1, expit_y0)
    
    prog_score_y1 <- b0_y1 + c*b1_y1*z + c*b2_y1*x
    expit_y1 <- 1/(1 + exp(-prog_score_y1))
    outcome_H1 <- rbinom(n, 1, expit_y1)


    data <- data.frame(outcome_H0, outcome_H1, x, z)
    
    #Patient accrual in months
    recruitment_acc <- enrollment(lambda = lambda, N_total = n, lambda_time = 0)
    
    generated_max_rec <- rnorm(1, mean=max_rec, sd=sd_rec)
    recruitment_acc <- generated_max_rec * recruitment_acc/max(recruitment_acc)

    #Saving the design
    sim_data[indexation, "Lambda"]         <- rep(lambda, times = n_methods)
    sim_data[indexation, "n"]              <- rep(n, times = n_methods)
    sim_data[indexation, "n1"]             <- rep(n1, times = n_methods)
    sim_data[indexation, "Interim_method"] <- interim_methods
    sim_data[indexation, "py_0"]           <- rep(mean(outcome_H0), 
                                                  time = n_methods)
    sim_data[indexation, "py_1"]           <- rep(mean(outcome_H1), 
                                                  time = n_methods)
    sim_data[indexation, "constant"]       <- rep(c,
                                                  time = n_methods)
    sim_data[indexation, "cor_H0_xy"]      <- rep(cor(outcome_H0, x),
                                                  time = n_methods)
    sim_data[indexation, "cor_H1_xy"]      <- rep(cor(outcome_H1, x),
                                                  time = n_methods)
    sim_data[indexation, "cor_H0_zy"]      <- rep(cor(outcome_H0, z),
                                                  time = n_methods)
    sim_data[indexation, "cor_H1_zy"]      <- rep(cor(outcome_H1, z),
                                                  time = n_methods)
    
    #######################################################################
    ###Decision stage1###
    #######################################################################
    
    ##Standard method##
    #Under H0#
    sm_s1 <- standard_method(data = data, 
                            outcome = "outcome_H0", 
                            recruitment_acc = recruitment_acc,
                            time_long = months_long, 
                            n1 = n1, 
                            r1 = r1)
    
    #Under H1#
    sm_s1_h1 <- standard_method(data = data, 
                                outcome = "outcome_H1",
                                recruitment_acc = recruitment_acc,
                                time_long = months_long,
                                n1 = n1,
                                r1 = r1)
    ##Adjusment baseline##
    vl_baseline_H0 <- VL_method(data = data,
                                outcome = "outcome_H0",
                                n1 = n1,
                                recruitment_acc = recruitment_acc, 
                                time_long = months_long, 
                                time_short = months_short,
                                covariate_vec = c("z"),
                                short_vec = c("x"),
                                adjustment = "baseline&short")
    
    VL_baseline_H0_deci <- VL_deci(variance = vl_baseline_H0$variance,
                                  estimate = vl_baseline_H0$estimate,
                                  n_interim = vl_baseline_H0$n_interim,
                                  p0 = p_LTE_0,
                                  PET = PET) 
    
    vl_baseline_H1 <- VL_method(data = data,
                                outcome = "outcome_H1",
                                n1 = n1,
                                recruitment_acc = recruitment_acc, 
                                time_long = months_long, 
                                time_short = months_short,
                                covariate_vec = c("z"),
                                short_vec = c("x"),
                                adjustment = "baseline&short")
    
    VL_baseline_H1_deci <- VL_deci(variance = vl_baseline_H1$variance,
                                  estimate = vl_baseline_H1$estimate,
                                  n_interim = vl_baseline_H1$n_interim,
                                  p0 = p_LTE_0,
                                  PET = PET) 
    
    #Saving the decisions into a data frame
    sim_data[indexation, "Decision_underH0_st1"] <- c(sm_s1$decision, 
                                                      VL_baseline_H0_deci)
    
    sim_data[indexation, "Decision_underH1_st1"] <- c(sm_s1_h1$decision,
                                                      VL_baseline_H1_deci)
    
    sim_data[indexation, "n_interim"] <- c(sm_s1$n_long,
                                          vl_baseline_H1$n_interim)
    
    sim_data[indexation, "n_extra"] <- c(0, 
                                        vl_baseline_H1$n_interim - sm_s1$n_long)
    
    sim_data[indexation, "n_cohort1"] <- c(0, vl_baseline_H1$n_cohort1)
    sim_data[indexation, "n_cohort2"] <- c(0, vl_baseline_H1$n_cohort2)
    sim_data[indexation, "n_cohort3"] <- c(0, vl_baseline_H1$n_cohort3)
    
    #######################################################################
    ###Decision stage2###
    #######################################################################
    
    dst2 <- decision_stage2(data[, "outcome_H0"], 
                            r2)
    
    dst2_h1 <- decision_stage2(data[, "outcome_H1"], 
                              r2)
    
    sim_data[indexation, "Decision_underH0_st2"] <- rep(dst2, 
                                                        times = n_methods)
    sim_data[indexation, "Decision_underH1_st2"] <- rep(dst2_h1,
                                                        times = n_methods)
    
    
    #######################################################################
    ###Checking continue in stage 1 = rejection  in  stage 2###
    #######################################################################
    coherence_reject <- data.frame(# type 1 error
      s_underH0 = 
        (sm_s1$decision == "continue") & 
        (dst2 == "reject H0"), 
      vl_baseline_underH0 = (VL_baseline_H0_deci == "continue") & 
        (dst2 == "reject H0"),
      
      # Power
      s_underH1 = 
        (sm_s1_h1$decision == "continue") & 
        (dst2_h1 == "reject H0"),
      vl_baseline_underH1 = (VL_baseline_H1_deci == "continue") & 
        (dst2_h1 == "reject H0")
    )
    
    sim_data[indexation, "Alpha"] <- c(coherence_reject$s_underH0, 
                                      coherence_reject$vl_baseline_underH0)
    
    sim_data[indexation, "Power"] <- c(coherence_reject$s_underH1, 
                                      coherence_reject$vl_baseline_underH1)
    #######################################################################
    #updating the indexation
    index <- index + n_methods
  }

  # if (b0 == -1) {
  #   directory_name <- "./data_expit_diff_p1"
  #   dir.create(directory_name, showWarnings = FALSE)
  #   write.csv(sim_data, paste(directory_name, "/sim_data_", c, ".csv", sep=""))
  # } else {
  #   directory_name <- "./data_expit_same_p1"
  #   dir.create(directory_name, showWarnings = FALSE)
  #   write.csv(sim_data, paste(directory_name, "/sim_data_", c, "_", b0, ".csv", sep=""))
  # }
  dir.create(directory_name, showWarnings = FALSE)
  write.csv(sim_data, paste(directory_name, "/sim_data_", c, ".csv", sep=""))
}
