library(devtools)
install_github("stchen3/glmmVCtest")
library("glmmVCtest")

devtools::install_github("Evolutionary-Optimization-Laboratory/rmoo")
library(ecr)
# ---- For: parallelization ----
# For: foreach loop
library(foreach)

run_parallel <- TRUE
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

test_multiple <- function(population, indexes){
  for(idx in indexes){
    cat("---", idx, "---\n")
    print(population[idx,])
    test_run(population[idx,])
  }
}

test_run <- function(x){
  mu1 <- c(16.114973, -85.344580, -87.923671)
  mu2 <- c(-97.980461,  50.390122,  58.193313)
  lpintercept <- 3.562851
  flparam <- x[1:3]
  timeseries_length = 90
  timestamps01 <- seq(from = 0.01, to = 0.99, length=timeseries_length)
  results <- NULL
  for (flc in c("30")) {
    cat("\tFL Choice:",flc,"\n")
    res <- GenerateCategoricalFDTest(3, mu1, mu2, 300, timeseries_length, timestamps01, flc, lpintercept, flparam)
    intpen <- calc_integral_penalty(timeseries_length, timestamps01, res, x[8:9])
    print(intpen)
    
    tab_y_raw <- table(res$yis)
    tab_y <- tab_y_raw / sum(tab_y_raw)
    cat("\t\tY:", tab_y_raw, "-->", tab_y,"\n")
    
    tab_y_without_raw <- table(res$yis_without)
    tab_y_without <- tab_y_without_raw / sum(tab_y_without_raw)
    cat("\t\tY_without:", tab_y_without_raw,"-->", tab_y_without, "\n")
    
    balance01 <- 1.0
    if(length(tab_y) >= 2 && length(tab_y_without) >= 2){
      balance01 <- max(abs(tab_y[[1]] - tab_y[[2]]), abs(tab_y_without[[1]] - tab_y_without[[2]]))
    }
    cat("\t\tMax 0-1 Balance in Y:", balance01, "\n")
    
    lp_distance <- abs(mean(res$linear_predictor$linearw) - mean(res$linear_predictor$linearwo))
    cat("\t\tLP mean Distance:", lp_distance, "\n")
    
    categories_in_w <- length(table(array(unlist(res$Truecatcurve))))
    cat("\t\tNum categories in W:", categories_in_w, "\n")
    
    results <- rbind(results, res)
  }
  return(results)
}


library(rmoo)

fitness_func <- function(x){
  if(is.null(dim(x))) {
    x <- matrix(x, ncol = 1)
  }
  print(t(x))
  
  flcs <- c("30")
  num_replications <- 5
  
  flc_result3 <- foreach(replica_num = 1:num_replications, .combine = cbind) %dorng% {
        source("source_code/R/optimize_essentials_03_fixed_mu.R")
        # my_fitness(x[1:3], x[4:6], x[7], flcs[1], x[8:9])
        my_fitness(c(16.114973, -85.344580, -87.923671), 
                   c(-97.980461,  50.390122,  58.193313), 
                   3.562851, 
                   flcs[1], 
                   x[1:3])
  }
  
  median_fitness <- apply(matrix(apply(flc_result3, 2, mean), ncol=num_replications), 1, max)
  
  return(median_fitness)
}

# known_candidates <- rbind(cbind(-40.2339963, -5.3899194, -3.1525938, 40.3651894, -38.1102743, 1.0417336, 0.2908304, 1.0, 1.0),
#                           cbind(-6.67, -2.47, 5.42, -3.14, -0.99, 3.91, 0.1, 1.0, 1.0),
#                           cbind(-6.67, -2.47, 5.42, -3.14, -0.99, 3.91, 1.0, 1.0, 1.0),
#                           cbind(1,2,3,1,2,3,1, 1.0, 1.0),
#                           cbind(1,2,3,1,2,3,0, 1.0, 1.0),
#                           cbind(-53.46193, 9.838961, -30.932283, 0.1438664, -60.11443, -8.940435, 0.9943165, 1.0, 1.0),
#                           cbind(-59.93567, 10.020723, 9.866337, 0.1438664, -69.19495, -14.787011, 1.3283204, 1.0, 1.0),
#                           cbind(-72.68827, 10.256293, 4.410807, 0.1438664, -63.97833, 2.120423, 1.3169213, 1.0, 1.0),
#                           cbind(-72.67520, -25.894426, 5.453164, 0.1438664, -63.97833, 2.120423, 1.3169213, 1.0, 1.0),
#                           cbind(-72.67520, -25.894426, 53.902214, 0.1438664, -63.97833, 2.120423, 1.3169213, 1.0, 1.0),
#                           cbind(-60.45411, 9.838961, 9.658552, 0.1438664, -69.19495, -14.787011, 1.3403095, 1.0, 1.0),
#                           cbind(-60.45411, 9.838961, 9.658552, 0.1438664, -69.19495, -14.787011, 1.3403095, 1.0, 1.0))

known_candidates <- NULL

begin_exp_time <- Sys.time()

set.seed(123)

ga <- nsga2(type = "real-valued",
             fitness = fitness_func,
             nObj = 2,
             lower = rep(-100.0, 3),
             upper = rep(100.0, 3),
             popSize = 200,
             summary = FALSE,
             parallel = FALSE,
             #monitor=FALSE,
             #suggestions = known_candidates,
             maxiter = 100)

summary(ga)
plot(ga)
ga_params = list("flcs"= c("30"), "popSize"=200, "num_indv"=300, "maxiter"=100, "method"="max of medians", "count_iter_indv"=100)
save(ga, ga_params, file="ga_run_02.RData")

end_exp_time <- Sys.time()

cat("\n====================\n",
    "\tAll Experiemnts Took:", capture.output(end_exp_time - begin_exp_time),
    "\n====================\n")


good_idxs <- intersect(which(ga@fitness[,1] < 0.4), which(ga@fitness[,2] < -0.8))
cat("\nGood ones at:", good_idxs, "\n")
print(data.frame(ga@population)[good_idxs,])

if(run_parallel)
{
  parallel::stopCluster(cl = my.cluster)
  initialized_parallel <- FALSE
}
