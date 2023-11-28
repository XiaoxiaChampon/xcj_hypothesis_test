
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

source("source_code/R/data_generator.R")


#record the time for one 
time_elapsed <<- list()
last_time <- 0
row_name <- NULL
timeKeeperStart <- function(rn)
{
  row_name <<- rn
  if(FALSE == row_name %in% names(time_elapsed))
  {
    time_elapsed[[row_name]] <<- NULL
  }
  last_time <<- Sys.time()
}
timeKeeperNext <- function()
{
  this_time <- Sys.time()
  this_section_time <- this_time - last_time
  cat(row_name, "calc time taken:", capture.output(this_section_time), "\n")
  time_elapsed[[row_name]] <<- append(time_elapsed[[row_name]], this_section_time)
  last_time <<- this_time
}
  
timeKeeperStart("GA Run")

timestamps01 <- seq(from = 0.01, to = 0.99, length=2000)

cfd_testing_part1 <- function(start_time, end_time, timeseries_length,
                        num_indvs, mu1_coef, mu2_coef, fl_choice,response_family,test_type,
                        klen=3){
  
  
  cfd_test_data <- GenerateCategoricalFDTest(klen=3,mu1_coef,mu2_coef,
                                             num_indvs=num_indvs,
                                             timeseries_length = timeseries_length,
                                             time_interval = timestamps01,
                                             fl_choice=fl_choice)
  return(cfd_test_data)
}
  
cfd_testing_part2 <- function(cfd_test_data){
  
  result <- cfd_hypothesis_test(cfd_test_data$true$yis,
                                cfd_test_data$true$Truecatcurve,
                                time_interval = timestamps01,
                                response_family='bernoulli', 
                                test_type='Functional')
  
  return(list("pvalue"=result$pvalue,"test_statistics"=result$statistics,
              "yis"=cfd_test_data$true$yis,"flt"=cfd_test_data$true$fl,
              "W"=cfd_test_data$true$Truecatcurve,
              "linear_predictor"=cfd_test_data$true$linear_predictor,
              "prob_ind"=cfd_test_data$true$prob_ind))
}

fitness_mufl <- function(x){
  
  num_replicas = 10
  
  p_value=NULL
  for (i in 1:num_replicas){
    print( i)
    cfd_test_data <- cfd_testing_part1(start_time=0.01, end_time=0.99, timeseries_length=2000,
                          num_indvs=100,
                          mu1_coef=c(x[1],x[2],x[3]),
                          mu2_coef=c(x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12]),
                          fl_choice=3,
                          response_family='bernoulli', test_type='Functional', klen=3)
    
    if(length(table(cfd_test_data$true$Truecatcurve)) >= 3){
      if(length(table(cfd_test_data$true$yis)) >= 2){
        result <- cfd_testing_part2(cfd_test_data)
        
        if(!is.null(result)){
          p_value[i]=result$pvalue
        }
      }
    }
  }
  
  if(is.null(p_value)){
    return(-1000000000)
  }
  
  power <- mean(p_value < 0.05)
  
  objective_value <- power * num_replicas
  
  return(objective_value)
}

ga <- ga(type = "real-valued",
         fitness = fitness_mufl,
         lower = rep(-10, 12),
         upper = rep(10,12),
         popSize = 10,
         maxiter = 20,
         run = 20,
         parallel = FALSE,
         monitor = TRUE)

summary(ga)
plot(ga)

timeKeeperNext()

# fitness_mu <- function(x){
#   source("source_code/R/data_generator.R")
#   
#   timeseries_length = 2000
#   timestamps01 <- seq(from = 0.01, to = 0.99, length=timeseries_length)
#   cfd_test_data <- GenerateCategoricalFDTest(klen=3,
#                                              mu1_coef=c(x[1],x[2],x[3]),
#                                              mu2_coef=c(x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12]),
#                                              num_indvs=100,
#                                              timeseries_length = timeseries_length,
#                                              time_interval = timestamps01,
#                                              fl_choice=3)
#   
#   if(length(table(cfd_test_data$true$Truecatcurve)) < 3){
#     return(-10000000000)
#   }
#   
#   linear_predictor=cfd_test_data$true$linear_predictor
#   Y_indvs=cfd_test_data$true$yis
#   if(length(table(Y_indvs)) < 2){
#     return(-10000000000)
#   }
#   prob_ind=cfd_test_data$true$prob_ind
#   data_sample=data.frame(linear_predictor,as.factor(Y_indvs),prob_ind)
#   colnames(data_sample)=c("linear_predictor","Group","Probability")
#   mean_0 <- mean(data_sample[data_sample$Group=="0","linear_predictor"])
#   mean_1 <- mean(data_sample[data_sample$Group=="1","linear_predictor"])
#   
#   # bonus1 = 0
#   abs_diff = abs(mean_0 - mean_1) # this should be high ( better to be more than 3)
#   # if(abs_diff > 1){
#   #   bonus1 = 10
#   # }
#   # if(abs_diff > 2){
#   #   bonus1 = 20
#   # }
#   # if(abs_diff > 3){
#   #   bonus1 = 30
#   # }
#   
#   # bonus2 = 0
#   # alloc_diff = abs(table(Y_indvs)[[1]] - table(Y_indvs)[[2]]) # this should be low
#   # if(alloc_diff < 10){
#   #   bonus2 = 1
#   # }
#   
#   # objective_value = bonus1 + abs_diff * 10 - alloc_diff + bonus2
#   objective_value = abs_diff
#   
#   return(objective_value)
# }
