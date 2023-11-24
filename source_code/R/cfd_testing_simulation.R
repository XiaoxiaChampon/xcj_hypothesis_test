
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
install_github("stchen3/glmmVCtest")
library("glmmVCtest")
library(RLRsim)
library(MASS)
library(splines)


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


cfd_testing_simulation <- function (num_replicas, start_time, end_time, timeseries_length,
                                    num_indvs,fl_choice,response_family,test_type,
                                    klen=3){
  cat("CFD Testing Simulation \nNum Replicas:\t", num_replicas)
  
  result_all <- foreach (number_simulation = 1:num_replicas, .combine = cbind, .init = NULL,
                         packages=c("splines","mgcv","fda","fda.usc","devtools","glmmVCtest","RLRsim","MASS")) %dorng% {
    #set.seed(123 + 8 * number_simulation)
    print(number_simulation)
    source("source_code/R/data_generator.R")
    result <- cfd_testing(start_time, end_time, timeseries_length,
                          num_indvs,fl_choice,response_family,test_type, klen=3)
    return(list("pvalue"=result$pvalue,"teststat"=result$test_statistics))
  } 
  return(result_all)
  
}

cfd_testing_simulation_no_paralel <- function (num_replicas, start_time, end_time, timeseries_length,
                                    num_indvs,fl_choice,response_family,test_type,
                                    klen=3){
  cat("CFD Testing Simulation \nNum Replicas:\t", num_replicas)
  source("source_code/R/data_generator.R")
  p_value=c(0)
  test_stats=c(0)
  #browser("mystop")
  for (i in 1:num_replicas){
    print( i)
    result <- cfd_testing(start_time, end_time, timeseries_length,
                          num_indvs,fl_choice,response_family,test_type, klen=3)
    p_value[i]=result$pvalue
    test_stats[i]=result$test_statistics
  }

  return(list("pvalue"=p_value,"teststat"=test_stats))

}

######test
timeKeeperStart("n100t2000")
set.seed(123456)
n100t2000_nopara <- cfd_testing_simulation_no_paralel(num_replicas=50, start_time=0.01, end_time=0.99,
                                    timeseries_length=2000,
                                    num_indvs=100,fl_choice=3,
                                    response_family='bernoulli',test_type='Functional',
                                    klen=3)
timeKeeperNext()
powern100t2000=mean(n100t2000_nopara$pvalue<0.05)


timeKeeperStart("n100t2000")
set.seed(123456)
n100t2000_nopara <- cfd_testing_simulation_no_paralel(num_replicas=50, start_time=0.01, end_time=0.99,
                                                      timeseries_length=2000,
                                                      num_indvs=100,fl_choice=4,
                                                      response_family='bernoulli',test_type='Functional',
                                                      klen=3)
timeKeeperNext()
powern100t2000=mean(n100t2000_nopara$pvalue<0.05)
powern100t2000
############


######test
timeKeeperStart("n100t2000")
set.seed(123456)
n100t2000 <- cfd_testing_simulation(num_replicas=50, start_time=0.01, end_time=0.99,
                                    timeseries_length=2000,
                                    num_indvs=100,fl_choice=3,
                                    response_family='bernoulli',test_type='Functional',
                                    klen=3)
timeKeeperNext()
n100t2000_data=matrix(n100t2000,nrow=2,ncol=50)
powern100t2000=mean(n100t2000_data[1,] < 0.05)
############




#record the time for one 
time_elapsed <<- list()
# "Xiaoxia"=NULL, "univfpca"=NULL, "kmeans"=NULL, "fadp"=NULL, "dbscan"=NULL, "cfd"=NULL)
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


#hyp3 expect to reject, it is rejecting , needs to see power 
# set.seed(1234)
# #power
timeKeeperStart("n100t300")
set.seed(1234)
n100t300 <- cfd_testing_simulation(num_replicas=5000, start_time=0.01, end_time=0.99,
                                                  timeseries_length=300,
                                                  num_indvs=100,fl_choice=3,
                                                  response_family='bernoulli',test_type='Functional',
                                                  klen=3)
timeKeeperNext()
n100t300_data=matrix(n100t300,nrow=2,ncol=5000)
powern100t300=mean(n100t300_data[1,] < 0.05)
test_statisticsn100t300=mean(unlist(n100t300_data[2,]))
test_statisticsn100t300sd=sd(unlist(n100t300_data[2,]))/sqrt(5000)


####
timeKeeperStart("n100t750")
set.seed(1234)
n100t750 <- cfd_testing_simulation(num_replicas=5000, start_time=0.01, end_time=0.99,
                                   timeseries_length=750,
                                   num_indvs=100,fl_choice=3,
                                   response_family='bernoulli',test_type='Functional',
                                   klen=3)
timeKeeperNext()
n100t750_data=matrix(n100t750,nrow=2,ncol=5000)
powern100t750=mean(n100t750_data[1,] < 0.05)
test_statisticsn100t750=mean(unlist(n100t750_data[2,]))
test_statisticsn100t750sd=sd(unlist(n100t750_data[2,]))/sqrt(5000)
####
timeKeeperStart("n100t1050")
set.seed(1234)
n100t1050 <- cfd_testing_simulation(num_replicas=5000, start_time=0.01, end_time=0.99,
                                   timeseries_length=1050,
                                   num_indvs=100,fl_choice=3,
                                   response_family='bernoulli',test_type='Functional',
                                   klen=3)
timeKeeperNext()
n100t1050_data=matrix(n100t1050,nrow=2,ncol=5000)
powern100t1050=mean(n100t1050_data[1,] < 0.05)
test_statisticsn100t1050=mean(unlist(n100t1050_data[2,]))
test_statisticsn100t1050sd=sd(unlist(n100t1050_data[2,]))/sqrt(5000)
####
timeKeeperStart("n100t1350")
set.seed(1234)
n100t1350 <- cfd_testing_simulation(num_replicas=5000, start_time=0.01, end_time=0.99,
                                    timeseries_length=1350,
                                    num_indvs=100,fl_choice=3,
                                    response_family='bernoulli',test_type='Functional',
                                    klen=3)
timeKeeperNext()
n100t1350_data=matrix(n100t1350,nrow=2,ncol=5000)
powern100t1350=mean(n100t1350_data[1,] < 0.05)
test_statisticsn100t1350=mean(unlist(n100t1350_data[2,]))
test_statisticsn100t1350sd=sd(unlist(n100t1350_data[2,]))/sqrt(5000)

####
timeKeeperStart("n100t2000")
set.seed(123)
n100t2000 <- cfd_testing_simulation(num_replicas=50, start_time=0.01, end_time=0.99,
                                    timeseries_length=2000,
                                    num_indvs=100,fl_choice=3,
                                    response_family='bernoulli',test_type='Functional',
                                    klen=3)
timeKeeperNext()
n100t2000_data=matrix(n100t2000,nrow=2,ncol=4)
powern100t2000=mean(n100t2000_data[1,] < 0.05)
test_statisticsn100t2000=mean(unlist(n100t2000_data[2,]))
test_statisticsn100t2000sd=sd(unlist(n100t2000_data[2,]))/sqrt(5000)

power_table=c(powern100t300,powern100t750,powern100t1050,
              powern100t1350,powern100t2000)

test_stat=c(test_statisticsn100t300,test_statisticsn100t750,
            test_statisticsn100t1050,test_statisticsn100t1350,
            test_statisticsn100t2000)
test_statsd=c(test_statisticsn100t300sd,test_statisticsn100t750sd,
            test_statisticsn100t1050sd,test_statisticsn100t1350sd,
            test_statisticsn100t2000sd)
save(power_table,test_stat,test_statsd,file="n100power.RData")
######
#hy4 is constant, test constant, expect to fail to reject
#type I error rate
# timeKeeperStart("n100t300")
# set.seed(1234)
# #check which is less 0.05 , type I error rate
# n100t300p <- cfd_testing_simulation(num_replicas=5000, start_time=0.01, end_time=0.99,
#                                    timeseries_length=300,
#                                    num_indvs=100,fl_choice=4,
#                                    response_family='bernoulli',test_type='Functional',
#                                    klen=3)
# timeKeeperNext()
# n100t300_datap=matrix(n100t300p,nrow=2,ncol=5000)
# powern100t300p=mean(n100t300_datap[1,] < 0.05)
# test_statisticsn100t300p=mean(unlist(n100t300_datap[2,]))
# test_statisticsn100t300sdp=sd(unlist(n100t300_datap[2,]))/sqrt(5000)
# 
# 
# ####
# timeKeeperStart("n100t750")
# set.seed(1234)
# n100t750p <- cfd_testing_simulation(num_replicas=5000, start_time=0.01, end_time=0.99,
#                                    timeseries_length=750,
#                                    num_indvs=100,fl_choice=4,
#                                    response_family='bernoulli',test_type='Functional',
#                                    klen=3)
# timeKeeperNext()
# n100t750_datap=matrix(n100t750p,nrow=2,ncol=5000)
# powern100t750p=mean(n100t750_datap[1,] < 0.05)
# test_statisticsn100t750p=mean(unlist(n100t750_datap[2,]))
# test_statisticsn100t750sdp=sd(unlist(n100t750_datap[2,]))/sqrt(5000)
# ####
# timeKeeperStart("n100t1050")
# set.seed(1234)
# n100t1050p <- cfd_testing_simulation(num_replicas=5000, start_time=0.01, end_time=0.99,
#                                     timeseries_length=1050,
#                                     num_indvs=100,fl_choice=4,
#                                     response_family='bernoulli',test_type='Functional',
#                                     klen=3)
# timeKeeperNext()
# n100t1050_datap=matrix(n100t1050p,nrow=2,ncol=5000)
# powern100t1050p=mean(n100t1050_datap[1,] < 0.05)
# test_statisticsn100t1050p=mean(unlist(n100t1050_datap[2,]))
# test_statisticsn100t1050sdp=sd(unlist(n100t1050_datap[2,]))/sqrt(5000)
# ####
# timeKeeperStart("n100t1350")
# set.seed(1234)
# n100t1350p <- cfd_testing_simulation(num_replicas=5000, start_time=0.01, end_time=0.99,
#                                     timeseries_length=1350,
#                                     num_indvs=100,fl_choice=4,
#                                     response_family='bernoulli',test_type='Functional',
#                                     klen=3)
# timeKeeperNext()
# n100t1350_datap=matrix(n100t1350p,nrow=2,ncol=5000)
# powern100t1350p=mean(n100t1350_datap[1,] < 0.05)
# test_statisticsn100t1350p=mean(unlist(n100t1350_datap[2,]))
# test_statisticsn100t1350sdp=sd(unlist(n100t1350_datap[2,]))/sqrt(5000)
# 
# ####
# timeKeeperStart("n100t2000")
# set.seed(1234)
# n100t2000p <- cfd_testing_simulation(num_replicas=5000, start_time=0.01, end_time=0.99,
#                                     timeseries_length=2000,
#                                     num_indvs=100,fl_choice=4,
#                                     response_family='bernoulli',test_type='Functional',
#                                     klen=3)
# timeKeeperNext()
# n100t2000_datap=matrix(n100t2000p,nrow=2,ncol=5000)
# powern100t2000p=mean(n100t2000_datap[1,] < 0.05)
# test_statisticsn100t2000p=mean(unlist(n100t2000_datap[2,]))
# test_statisticsn100t2000sdp=sd(unlist(n100t2000_datap[2,]))/sqrt(5000)
# 
# power_tablep=c(powern100t300p,powern100t750p,powern100t1050p,
#               powern100t1350p,powern100t2000p)
# 
# test_statp=c(test_statisticsn100t300p,test_statisticsn100t750p,
#             test_statisticsn100t1050p,test_statisticsn100t1350p,
#             test_statisticsn100t2000p)
# test_statsdp=c(test_statisticsn100t300sdp,test_statisticsn100t750sdp,
#               test_statisticsn100t1050sdp,test_statisticsn100t1350sdp,
#               test_statisticsn100t2000sdp)
# save(power_tablep,test_statp,test_statsdp,file="n100typeIerror.RData")
###########

# set.seed(123)
# hy3test=cfd_testing (start_time=0.01, end_time=0.99, timeseries_length=2000,
#                         num_indvs=100,fl_choice=3,response_family="'bernoulli'",test_type='Functional',
#                         klen=3)

if(run_parallel)
{
  parallel::stopCluster(cl = my.cluster)
  initialized_parallel <- FALSE
}
