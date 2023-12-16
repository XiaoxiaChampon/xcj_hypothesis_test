##################################
#
# Genetic Algorithm Based Optimization for finding mu1, mu2, and lp_intercept values.
#

library(foreach)
library(doParallel)
library(doRNG)

library(rmoo)

# library(cdata)
# library(reshape2)
# library(rgl)
# library(ecr)
# library(emoa)

begin_exp_time <- Sys.time()

source("source_code/R/data_generator.R")

timeseries_length = 180
timestamps01 <- seq(from = 0.01, to = 0.99, length=timeseries_length)

# mu1_coef=c(-6.67,-2.47,5.42)
# mu2_coef=c(-3.14,-0.99,3.91)

run_test <- function(x, flc){
  mu1_coef <- x[1:3]
  mu2_coef <- x[4:6]
  intercept <- x[7]
  results <- GenerateCategoricalFDTest(3, mu1_coef, mu2_coef, 500, timeseries_length, timestamps01, flc, intercept)

  tab_y <- table(results$true$yis)
  tab_y <- tab_y / sum(tab_y)

  tab_y_without <- table(results$true$yis_without)
  tab_y_without <- tab_y_without / sum(tab_y_without)

  # balance01 is a measure of balance between 0s and 1s in the Ys.
  # It is the absolute difference of the fraction of 0s and 1s.
  # 0 <= balance01 <= 1
  # We prefer when it is minimal. 0 is the best value. 1 is the worst.
  # We will be happy if balance01 <= 0.3
  balance01 <- 1.0
  if(length(tab_y) >= 2 && length(tab_y_without) >= 2){
    balance01 <- max(abs(tab_y[[1]] - tab_y[[2]]), abs(tab_y_without[[1]] - tab_y_without[[2]]))
  }

  # We give high penalty when balance01 > 0.3
  # 0 <= balance_penalty < 10
  balance_penalty <- balance01
  # if(balance01 > 0.3){
  #   balance_penalty <- balance01 * 10
  # }

  # lp_distance is the absolute distance between LP means.
  # Theoretically: 0 <= lp_distance < infinity
  # We prefer when lp_distance is maximum.
  # We will be happy if lp_distance >= 1
  lp_distance <- abs(mean(results$true$linear_predictor$linearw) - mean(results$true$linear_predictor$linearwo))

  # We take the negative of the lp_distance as the penalty since we are minimizing
  # We give high penalty when lp_distance < 1
  distance_penalty <- lp_distance * -1.0
  # if(lp_distance >= 1){
  #   distance_penalty <- lp_distance * -1.0 - 1.1
  # }

  minimizing_fitness <- cbind(balance01, balance_penalty, lp_distance, distance_penalty)
  return(list("yw"=tab_y, "ywo"=tab_y_without, "fit"=minimizing_fitness))
}

my_fitness <- function(mu1_coef, mu2_coef, intercept, flc){
  results <- GenerateCategoricalFDTest(3, mu1_coef, mu2_coef, 500, timeseries_length, timestamps01, flc, intercept)

  tab_y <- table(results$true$yis)
  tab_y <- tab_y / sum(tab_y)

  tab_y_without <- table(results$true$yis_without)
  tab_y_without <- tab_y_without / sum(tab_y_without)

  # balance01 is a measure of balance between 0s and 1s in the Ys.
  # It is the absolute difference of the fraction of 0s and 1s.
  # 0 <= balance01 <= 1
  # We prefer when it is minimal. 0 is the best value. 1 is the worst.
  # We will be happy if balance01 <= 0.3
  balance01 <- 1.0
  if(length(tab_y) >= 2 && length(tab_y_without) >= 2){
    balance01 <- max(abs(tab_y[[1]] - tab_y[[2]]), abs(tab_y_without[[1]] - tab_y_without[[2]]))
  }

  # We give high penalty when balance01 > 0.3
  # 0 <= balance_penalty < 10
  balance_penalty <- balance01
  # if(balance01 > 0.3){
  #   balance_penalty <- balance01 * 10
  # }

  # lp_distance is the absolute distance between LP means.
  # Theoretically: 0 <= lp_distance < infinity
  # We prefer when lp_distance is maximum.
  # We will be happy if lp_distance >= 1
  lp_distance <- abs(mean(results$true$linear_predictor$linearw) - mean(results$true$linear_predictor$linearwo))

  # We take the negative of the lp_distance as the penalty since we are minimizing
  # We give high penalty when lp_distance < 1
  distance_penalty <- lp_distance * -1.0
  # if(lp_distance >= 1){
  #   distance_penalty <- lp_distance * -1.0 - 1.1
  # }

  minimizing_fitness <- cbind(balance_penalty, distance_penalty)
  return(minimizing_fitness)
}

set.seed(123)

# my_fitness(mu1_coef, mu2_coef, 1, "6")

fitness_func <- function(x){
  if(is.null(dim(x))) {
    x <- matrix(x, ncol = 1)
  }
  flc_result <- list()
  for(flc in c("6", "8", "10")){
    res_list <- list()
    for (replica_num in c(1:3)) {
      res <- my_fitness(x[1:3], x[4:6], x[7], flc)
      res_list <- rbind(res_list, res)
    }
    median_results <- apply(matrix(unlist(res_list), ncol=2), 2, median)
    flc_result <- rbind(flc_result, median_results)
  }

  median_fitness <- apply(matrix(unlist(flc_result), ncol=2), 2, median)
  return(median_fitness)
}

ga <- nsga2(type = "real-valued",
            fitness = fitness_func,
            nObj = 2,
            lower = rep(-10.0,7),
            upper = rep(10.0,7),
            popSize = 10,
            summary = TRUE,
            monitor = TRUE,
            parallel = TRUE,
            maxiter = 10)

summary(ga)
plot(ga)

end_exp_time <- Sys.time()

cat("\n====================\n",
    "\tAll Experiemnts Took:", capture.output(end_exp_time - begin_exp_time),
    "\n====================\n")

