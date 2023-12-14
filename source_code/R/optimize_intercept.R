
# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation
# files (the “Software”), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge,
# publish, distribute, sublicense, and/or sell copies of the Software,
# and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:

# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
# THE USE OR OTHER DEALINGS IN THE SOFTWARE.
######################################################################
#
# Purpose: Simulations for Categorical Functional Data Hypothesis Testing
#         
# Author:  Xiaoxia Champon
# Date: 10/26/2023
#
##############################################################




library(mgcv)
library(fda)
library(fda.usc)
library(devtools)
# install_github("stchen3/glmmVCtest")
library("glmmVCtest")
library(RLRsim)
library(MASS)
library(splines)

library(GA)


# ---- For: parallelization ----
# For: foreach loop
# library(foreach)
# 
# run_parallel <- TRUE
# time_elapsed <- list()
# if(run_parallel)
# {
#   print("RUNNING PARALLEL")
#   
#   # For: makeCluster
#   library(doParallel)
#   
#   # For: %dorng% or registerDoRNG for reproducable parallel random number generation
#   library(doRNG)
#   
#   if(exists("initialized_parallel") && initialized_parallel == TRUE)
#   {
#     parallel::stopCluster(cl = my.cluster)
#   }
#   n.cores <- parallel::detectCores() - 1
#   my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
#   doParallel::registerDoParallel(cl = my.cluster)
#   cat("Parellel Registered: ", foreach::getDoParRegistered(), "\n")
#   initialized_parallel <- TRUE
#   
#   # registerDoRNG(123) # ///<<<< THIS CREATES THE ERROR FOR FADPClust !!!
# }


source("source_code/R/data_generator.R")

begin_exp_time <- Sys.time()

timeseries_length = 180
timestamps01 <- seq(from = 0.01, to = 0.99, length=timeseries_length)

mu1_coef=c(-6.67,-2.47,5.42)
mu2_coef=c(-3.14,-0.99,3.91)

run_test <- function(mu1_coef, mu2_coef, intercept, flc){
  results <- GenerateCategoricalFDTest(3, mu1_coef, mu2_coef, 500, timeseries_length, timestamps01, flc, intercept)
  
  tab_y <- table(results$true$yis)
  tab_y <- tab_y / sum(tab_y)
  
  tab_y_without <- table(results$true$yis_without)
  tab_y_without <- tab_y_without / sum(tab_y_without)
  
  if(length(tab_y) < 2 || length(tab_y_without) < 2){
    # return(list("tab_y"=tab_y, "tab_y_without"=tab_y_without, "val"=1.0))
    return(1.0)
  }
  
  val <- max(abs(tab_y[[1]] - tab_y[[2]]), abs(tab_y_without[[1]] - tab_y_without[[2]]))
  
  # return(list("tab_y"=tab_y, "tab_y_without"=tab_y_without, "val"=val))
  return(val)
}

set.seed(123)

# run_test(mu1_coef, mu2_coef, 1, "6")

# ed_table <- expand.grid(seq(-5,5, length.out = 10), c("6", "7", "8", "9", "10", "200"))
# colnames(ed_table) <- c("intercept", "fl_choice")
# print(ed_table)
# 
# results <- foreach (row_index = 1:dim(ed_table)[1], .combine = cbind, .init = NULL) %dorng% {
#   res_list <- list()
#   for (replica_num in c(1:5)) {
#     res <- run_test(mu1_coef, mu2_coef, ed_table[row_index,]$intercept, as.character(ed_table[row_index,]$fl_choice))
#     res_list <- rbind(res_list, res)
#   }
#   return(median(unlist(res_list)))
# }
# 
# print(array(results[1,]))
# plot(array(results[1,]))

fitness_func <- function(x){
  flc_result <- list()
  for(flc in c("6", "8", "10")){
    res_list <- list()
    for (replica_num in c(1:3)) {
      res <- run_test(x[1:3], x[4:6], x[7], flc)
      res_list <- rbind(res_list, res)
    }
    meanval <- median(unlist(res_list))
    flc_result <- rbind(meanval, flc_result)
  }
  return(-median(unlist(flc_result)))
}


ga <- ga(type = "real-valued",
         fitness = fitness_func,
         lower = rep(-10, 7),
         upper = rep(10,7),
         popSize = 10,
         maxiter = 10,
         run = 10,
         parallel = TRUE,
         monitor = TRUE)

summary(ga)
plot(ga)

end_exp_time <- Sys.time()

cat("\n====================\n",
    "\tAll Experiemnts Took:", capture.output(end_exp_time - begin_exp_time), 
    "\n====================\n")

# if(run_parallel)
# {
#   parallel::stopCluster(cl = my.cluster)
#   initialized_parallel <- FALSE
# }