library(goldilocks)

means <- c()
sds <- c()
lambda <- c(0.5, 1, 2, 3, 4)
n <- 149

for(l in lambda){
  rec_data <- c()
  for (i in 1:10000){
    recruitment_acc <- enrollment(lambda = l, N_total = n, lambda_time = 0)
    rec_data <- c(rec_data, max(recruitment_acc))
  }
  means <- c(means, mean(rec_data))
  sds <- c(sds, sd(rec_data))
  #hist(rec_data)
}

data.frame(lambda, means, sds)

# data.frame(lambda = c(0.5, 1, 2, 3, 4),
#            mean = c(295.7055, 147.983, 73.88, 49.1437, 36.8025),
#            sd = c(24.28236, 12.05218, 6.049154, 4.071377, 3.098944))
