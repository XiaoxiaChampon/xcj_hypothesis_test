

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
  flcs <- c("6", "7","8", "9","10")
  num_fls <- 5
  num_replications <- 3
  rng <- rngtools::RNGseq( num_fls * num_replications, 1234)
  
  flc_result3 <- foreach(flcidx = 1:num_fls, .combine = rbind) %:%
    foreach(replica_num = 1:num_replications, .combine = cbind, r=rng[(flcidx-1)*num_fls + 1:num_fls]) %dopar% {
        rngtools::setRNG(r)
        source("source_code/R/optimize_essentials.R")
        my_fitness(x[1:3], x[4:6], x[7], flcs[flcidx])
      }
  
  #median_fitness <- apply(matrix(apply(flc_result3, 2, mean), ncol=3), 1, mean)
  median_fitness <- apply(matrix(apply(flc_result3, 2, mean), ncol=3), 1, max)
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

##only 6, 7, 8
ga <- nsga2(type = "real-valued",
             fitness = fitness_func,
             nObj = 2,
             lower = rep(-100.0,7),
             upper = rep(100.0,7),
             popSize = 100,
             summary = FALSE,
             parallel = FALSE,
             #monitor=FALSE,
             suggestions = known_candidates,
             maxiter = 100)
save(ga,file="ga.RData")
summary(ga)
plot(ga)
#intersect(which(ga@fitness[,1]<0.3) ,which(ga@fitness[,2]< -1))
#[1]  6  9 90 91

#ga@population[c(6,9,90,91),]
#[,1]      [,2]     [,3]      [,4]      [,5]     [,6]     [,7]
#[1,] -22.07539 -2.881389 1.909785 -1.122699 -28.24956 3.959064 1.106851
#[2,] -15.44999 -2.881389 2.559723 -1.152536 -27.90200 3.903929 1.106851
#[3,] -21.57212 -2.881389 1.897356 -1.152536 -27.78895 3.952971 1.109069
#[4,] -22.07539 -2.881389 1.909785 -1.122699 -28.24956 3.959064 1.106851


#6, ,7, ,8, 9, 10
ga1 <- nsga2(type = "real-valued",
            fitness = fitness_func,
            nObj = 2,
            lower = rep(-100.0,7),
            upper = rep(100.0,7),
            popSize = 100,
            summary = FALSE,
            parallel = FALSE,
            #monitor=FALSE,
            suggestions = known_candidates,
            maxiter = 100)
save(ga1,file="ga1.RData")
summary(ga1)
plot(ga1)

load("ga1.RData")
intersect(which(ga1@fitness[,1]<0.3) ,which(ga1@fitness[,2]< -1))
# 7 22
#ga1@population[c(7,22),]

end_exp_time <- Sys.time()

cat("\n====================\n",
    "\tAll Experiemnts Took:", capture.output(end_exp_time - begin_exp_time),
    "\n====================\n")

if(run_parallel)
{
  parallel::stopCluster(cl = my.cluster)
  initialized_parallel <- FALSE
}
