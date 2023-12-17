
# ---- For: parallelization ----
# For: foreach loop
library(foreach)

run_parallel <- FALSE
time_elapsed <- list()
if(run_parallel)
{
  print("RUNNING PARALLEL")
  
  # For: makeCluster
  library(doParallel)
  
  # For: %dorng% or registerDoRNG for reproducable parallel random number generation
  library(doRNG)
  
  if(exists("initialized_parallel") && initialized_parallel == TRUE)
  {
    parallel::stopCluster(cl = my.cluster)
  }
  n.cores <- parallel::detectCores() - 1
  my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
  doParallel::registerDoParallel(cl = my.cluster)
  cat("Parellel Registered: ", foreach::getDoParRegistered(), "\n")
  initialized_parallel <- TRUE
  
  # registerDoRNG(123) # ///<<<< THIS CREATES THE ERROR FOR FADPClust !!!
}


library(rmoo)

source("source_code/R/optimize_essentials.R")

timeseries_length = 180
timestamps01 <- seq(from = 0.01, to = 0.99, length=timeseries_length)

my_fitness <- function(mu1_coef, mu2_coef, intercept, flc){
  results <- GenerateCategoricalFDTest(3, mu1_coef, mu2_coef, 500, timeseries_length, timestamps01, flc, intercept)

  tab_y <- table(results$true$yis)
  tab_y <- tab_y / sum(tab_y)

  tab_y_without <- table(results$true$yis_without)
  tab_y_without <- tab_y_without / sum(tab_y_without)

  balance01 <- 1.0
  if(length(tab_y) >= 2 && length(tab_y_without) >= 2){
    balance01 <- max(abs(tab_y[[1]] - tab_y[[2]]), abs(tab_y_without[[1]] - tab_y_without[[2]]))
  }

  balance_penalty <- balance01
  
  lp_distance <- abs(mean(results$true$linear_predictor$linearw) - mean(results$true$linear_predictor$linearwo))

  distance_penalty <- -lp_distance

  minimizing_fitness <- cbind(balance_penalty, distance_penalty)
  return(minimizing_fitness)
}

set.seed(123)

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
    median_results <- apply(matrix(unlist(res_list), ncol=2), 2, max)
    flc_result <- rbind(flc_result, median_results)
  }

  median_fitness <- apply(matrix(unlist(flc_result), ncol=2), 2, max)
  return(median_fitness)
}

ga <- nsga2(type = "real-valued",
            fitness = fitness_func,
            nObj = 2,
            lower = rep(-10000.0,7),
            upper = rep(10000.0,7),
            popSize = 10,
            summary = FALSE,
            monitor = TRUE,
            parallel = FALSE,
            maxiter = 2)

summary(ga)
plot(ga)

if(run_parallel)
{
  parallel::stopCluster(cl = my.cluster)
  initialized_parallel <- FALSE
}

