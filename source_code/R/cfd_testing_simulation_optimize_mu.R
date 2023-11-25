
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


# # ---- For: parallelization ----
# # For: foreach loop
# library(foreach)
# 
# run_parallel <- FALSE
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

library(GA)
  
fitness_mu <- function(x){
  source("source_code/R/data_generator.R")
  
  timeseries_length = 2000
  timestamps01 <- seq(from = 0.01, to = 0.99, length=timeseries_length)
  cfd_test_data <- GenerateCategoricalFDTest(klen=3,
                                             mu1_coef=c(x[1],x[2],x[3]),
                                             mu2_coef=c(x[4],x[5],x[6]),
                                             num_indvs=100,
                                             timeseries_length = timeseries_length,
                                             time_interval = timestamps01,
                                             fl_choice=3)
  
  linear_predictor=cfd_test_data$true$linear_predictor
  Y_indvs=cfd_test_data$true$yis
  prob_ind=cfd_test_data$true$prob_ind
  data_sample=data.frame(linear_predictor,as.factor(Y_indvs),prob_ind)
  colnames(data_sample)=c("linear_predictor","Group","Probability")
  mean_0 <- mean(data_sample[data_sample$Group=="0","linear_predictor"])
  mean_1=mean(data_sample[data_sample$Group=="1","linear_predictor"])
  
  abs_diff = abs(mean_0 - mean_1) # this should be high
  
  alloc_diff = abs(table(Y_indvs)[[1]] - table(Y_indvs)[[2]]) # this should be low
  
  penalty = abs_diff - alloc_diff
  
  return(penalty)
}

ga <- ga(type = "real-valued",
         fitness = fitness_mu,
         lower = rep(-100, 6),
         upper = rep(100,6),
         popSize = 5,
         maxiter = 10,
         run = 10,
         parallel = TRUE,
         monitor = TRUE)

summary(ga)
plot(ga)

# if(run_parallel)
# {
#   parallel::stopCluster(cl = my.cluster)
#   initialized_parallel <- FALSE
# }

