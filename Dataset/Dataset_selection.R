library(goldilocks)

set.seed(123)

data_url <-
  paste0("https://github.com/jbetz-jhu/CovariateAdjustmentTutorial",
         "/raw/main/Simulated_MISTIE_III_v1.2.csv")

sim_miii_full <- read.csv(file = url(data_url))

#based on https://cancer.unc.edu/biostatistics/program/ivanova/SimonsTwoStageDesign.aspx
#type I error = .025, power = .80, p0 = .25, p1 = .38
design <- data.frame(method = "Optimal",
                     p0 = .25,
                     p1 = .38,
                     n = 146,
                     n1 = 48,
                     r1 = 13,
                     r2 = 46,
                     PET = .6986)

#Selecting relevant variables
data <- sim_miii_full

#Collapsing the short and long-term endpoints into binary variablesS
data$mrs_30d_binary <- ifelse(data$mrs_30d_complete >= 4, 0, 1)
data$mrs_180d_binary <- ifelse(data$mrs_180d_complete >= 4, 0, 1)
data$mrs_365d_binary <- ifelse(data$mrs_365d_complete >= 4, 0, 1)

#Keeping only participants assigned to the treatment group
data <- data[data$arm == "surgical", ]

#Sampling n participant for the study
id_thesis <- sample(data$sim_participant_id, design$n)
data_thesis <- data[data$sim_participant_id %in% id_thesis, ]


recruitment_rate <- enrollment(lambda = .3, N_total = design$n, lambda_time = 0)
# recruitment_rate <- enrollment(lambda = .4, N_total = design$n, lambda_time = 0)
# recruitment_rate <- enrollment(lambda = .2, N_total = design$n, lambda_time = 0)
