
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
  # library(doRNG)
  
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

fitness_func <- function(x){
  if(is.null(dim(x))) {
    x <- matrix(x, ncol = 1)
  }
  print(t(x))
  
  # ####################################################################
  # flc_result <- list()
  # for(flc in c("6", "8", "10")){
  #   res_list <- list()
  #   for (replica_num in c(1:3)) {
  #     # res <- my_fitness(x[1:3], x[4:6], x[7], flc)
  #     res <- testness(x[1:3], x[4:6], x[7], flc)
  #     res_list <- rbind(res_list, res)
  #   }
  #   median_results <- apply(matrix(unlist(res_list), ncol=2), 2, max)
  #   flc_result <- rbind(flc_result, median_results)
  # }
  # median_fitness <- apply(matrix(unlist(flc_result), ncol=2), 2, max)
  # ####################################################################
  # flc_result2 <- foreach(flc = c("6", "8", "10"), .combine = rbind, .init = NULL) %do% {
  #     res_list <- list()
  #     for (replica_num in c(1:3)) {
  #       # res <- my_fitness(x[1:3], x[4:6], x[7], flc)
  #       res <- testness(x[1:3], x[4:6], x[7], flc)
  #       res_list <- rbind(res_list, res)
  #     }
  #     median_results <- apply(matrix(unlist(res_list), ncol=2), 2, max)
  #   }
  # median_fitness2 <- apply(matrix(unlist(flc_result2), ncol=2), 2, max)
  # ####################################################################
  flcs <- c("6", "7", "8", "9", "10")
  rng <- rngtools::RNGseq( 5 * 3, 1234)
  
  flc_result3 <- foreach(flcidx = 1:5, .combine = rbind) %:%
    foreach(replica_num = 1:3, .combine = cbind, r=rng[(flcidx-1)*5 + 1:5]) %dopar% {
        rngtools::setRNG(r)
        source("source_code/R/optimize_essentials.R")
        my_fitness(x[1:3], x[4:6], x[7], flcs[flcidx])
      }
  
  median_fitness <- apply(matrix(apply(flc_result3, 2, mean), ncol=3), 1, mean)
  ####################################################################
  
  print(median_fitness)
  return(median_fitness)
}

known_candidates <- rbind(cbind(-40.2339963, -5.3899194, -3.1525938, 40.3651894, -38.1102743, 1.0417336, 0.2908304),
                          cbind(-6.67, -2.47, 5.42, -3.14, -0.99, 3.91, 0.1),
                          cbind(-6.67, -2.47, 5.42, -3.14, -0.99, 3.91, 1.0),
                          cbind(1,2,3,1,2,3,1),
                          cbind(1,2,3,1,2,3,0))

begin_exp_time <- Sys.time()

set.seed(123)

ga <- nsga2(type = "real-valued",
            fitness = fitness_func,
            nObj = 2,
            lower = rep(-100.0,7),
            upper = rep(100.0,7),
            popSize = 100,
            summary = FALSE,
            parallel = FALSE,
            suggestions = known_candidates,
            maxiter = 100)

summary(ga)
plot(ga)
save(ga, file = "ga_run_main.RData")

end_exp_time <- Sys.time()

cat("\n====================\n",
    "\tAll Experiemnts Took:", capture.output(end_exp_time - begin_exp_time),
    "\n====================\n")

if(run_parallel)
{
  parallel::stopCluster(cl = my.cluster)
  initialized_parallel <- FALSE
}
