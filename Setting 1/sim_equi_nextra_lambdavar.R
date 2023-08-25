library(goldilocks)

source("Functions2.R")


sim <- function(p_short, l, max_rec, sd_rec, n_sim, p_LTE_1, design_test) {
  p_LTE_0 <- .25
  alpha <- .05
  power_start <- .80 
  
  #Time in month after entering the study for when the long-term endpoint is reached 
  months_long <- 12*4
  
  #Time in months after entering the study for when the short-term endpoint is reached
  months_short <- 6*4
  
  #Data from the simulation
  sim_data <- data.frame(Design_method  = character(0),
                         Interim_method = character(0),
                         n              = numeric(0),
                         n1             = numeric(0),
                         n_interim      = numeric(0),
                         Lambda         = numeric(0),
                         p_short        = numeric(0),
                         Decision_underH0_st1 = character(0),
                         Decision_underH0_st2 = character(0),
                         Decision_underH1_st1 = character(0),
                         Decision_underH1_st2 = character(0),
                         Alpha = numeric(0),
                         Power = numeric(0))
  
  #The different methods
  interim_methods <- c("Standard", "VL_short")
  n_methods <- length(interim_methods)
  
  index <- 1
  
  for (i in seq_len(nrow(design_test))){
    
    design_method <- design_test[i, "method"]
    n <- design_test[i, "n"]
    n1 <- design_test[i, "n1"]
    r1 <- design_test[i, "r1"]
    r2 <- design_test[i, "r2"]
    PET <- design_test[i, "PET"]
    for (j in 1:n_sim) {
      
      #indexation for saving the data in sim_data 
      indexation <- index:(index + n_methods - 1)
      
      s <- rbinom(n, 1, p_short)
      
      short_vec <- c("s")
      
      outcome_H0 <- binom_gen(x = s, n = n, py = p_LTE_0, px = p_short)
      outcome_H1 <- binom_gen(x = s, n = n, py = p_LTE_1, px = p_short)
      
      data <- data.frame(outcome_H0, outcome_H1, s)
      
      #Patient accrual in months
      recruitment_acc <- enrollment(lambda = l, N_total = n, lambda_time = 0)
      
      generated_max_rec <- rnorm(1, mean=max_rec, sd=sd_rec)
      recruitment_acc <- generated_max_rec * recruitment_acc/max(recruitment_acc)
      #recruitment_acc <- 1:n
      
      #Saving the design
      sim_data[indexation, "Design_method"]  <- rep(design_method,
                                                    times = n_methods)
      sim_data[indexation, "Lambda"]         <- rep(l, times = n_methods)
      sim_data[indexation, "n"]              <- rep(n, times = n_methods)
      sim_data[indexation, "n1"]             <- rep(n1, times = n_methods)
      sim_data[indexation, "Interim_method"] <- interim_methods
      sim_data[indexation, "py_0"]           <- rep(mean(outcome_H0), 
                                                    time = n_methods)
      sim_data[indexation, "py_1"]           <- rep(mean(outcome_H1), 
                                                    time = n_methods)
      sim_data[indexation, "cor_H0"]         <- rep(cor(s, outcome_H0), 
                                                    time = n_methods)   
      sim_data[indexation, "cor_H1"]         <- rep(cor(s, outcome_H1), 
                                                    time = n_methods)                                     
      sim_data[indexation, "p_short"]        <- rep(p_short,
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
      ##Adjusment short##
      vl_short_H0 <- VL_method(data = data,
                               outcome = "outcome_H0",
                               n1 = n1,
                               recruitment_acc = recruitment_acc, 
                               time_long = months_long, 
                               short_vec = short_vec,
                               time_short = months_short, 
                               adjustment = "short")
      
      VL_short_H0_deci <- VL_deci(variance = vl_short_H0$variance,
                                  estimate = vl_short_H0$estimate,
                                  n_interim = vl_short_H0$n_interim,
                                  p0 = p_LTE_0,
                                  PET = PET) 
      
      vl_short_H1 <- VL_method(data = data,
                               outcome = "outcome_H1",
                               n1 = n1,
                               recruitment_acc = recruitment_acc, 
                               time_long = months_long, 
                               short_vec = short_vec,
                               time_short = months_short, 
                               adjustment = "short")
      
      VL_short_H1_deci <- VL_deci(variance = vl_short_H1$variance,
                                  estimate = vl_short_H1$estimate,
                                  n_interim = vl_short_H1$n_interim,
                                  p0 = p_LTE_0,
                                  PET = PET) 
      
      #Saving the decisions into a data frame
      sim_data[indexation, "Decision_underH0_st1"] <- c(sm_s1$decision, 
                                                        VL_short_H0_deci)
      
      sim_data[indexation, "Decision_underH1_st1"] <- c(sm_s1_h1$decision,
                                                        VL_short_H1_deci)
      
      sim_data[indexation, "n_interim"] <- c(sm_s1$n_long,
                                             vl_short_H1$n_interim)
      
      sim_data[indexation, "n_extra"] <- c(0, 
                                           vl_short_H1$n_interim - sm_s1$n_long)
      
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
        vl_short_underH0 = (VL_short_H0_deci == "continue") & 
          (dst2 == "reject H0"),
        
        # Power
        s_underH1 = 
          (sm_s1_h1$decision == "continue") & 
          (dst2_h1 == "reject H0"),
        vl_short_underH1 = (VL_short_H1_deci == "continue") & 
          (dst2_h1 == "reject H0")
      )
      
      sim_data[indexation, "Alpha"] <- c(coherence_reject$s_underH0, 
                                         coherence_reject$vl_short_underH0)
      
      sim_data[indexation, "Power"] <- c(coherence_reject$s_underH1, 
                                         coherence_reject$vl_short_underH1)
      #######################################################################
      #updating the indexation
      index <- index + n_methods
    }
  }
  directory_name <- paste("./data_nextra_lambda_n_", design_test$n, "_", n_sim, sep="")
  dir.create(directory_name, showWarnings = FALSE)
  write.csv(sim_data, paste(directory_name, "/sim_data", p_short, "_", l, ".csv", sep=""))
}