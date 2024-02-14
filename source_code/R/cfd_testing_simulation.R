
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
                                    mu1_coef, mu2_coef,
                                    num_indvs,fl_choice,response_family,test_type,gam_choice,
                                    klen=3){
  cat("CFD Testing Simulation \nNum Replicas:\t", num_replicas)
  #num_replicas=2
  result_all <- foreach (number_simulation = 1:num_replicas, .combine = cbind, .init = NULL,
                        .packages=c("glmmVCtest","splines","mgcv","fda","fda.usc","devtools","RLRsim","MASS")) %dorng% {
  # result_all <- foreach (number_simulation = 1:num_replicas, .combine = cbind, .init = NULL) %dorng% {
    source("./source_code/R/data_generator.R")
    result <- cfd_testing(start_time, end_time, timeseries_length,
                          num_indvs,mu1_coef,mu2_coef,fl_choice,response_family,
                          test_type, gam_choice,klen=3)
    
    
    #return(list("pvalue"=result$pvalue,"teststat"=result$test_statistics,"fl"=result$flt))
    #return(list("pvalue"=result$pvalue,"yip"=result$yip,"yip_wo"=result$yip_wo,"pvalue2"=result$pvalue2))
    return(list("pvalue"=result$pvalue,"yip"=result$yip,"yip_wo"=result$yip_wo))
  } 
  return(result_all)
  
}



# cfd_testing_simulation_no_paralel <- function (num_replicas, start_time, end_time, timeseries_length,
#                                                mu1_coef,mu2_coef,num_indvs,fl_choice,response_family,test_type,
#                                     klen=3){
#   # cat("CFD Testing Simulation \nNum Replicas:\t", num_replicas)
#   source("source_code/R/data_generator.R")
#   p_value=c(0)
#   test_stats=c(0)
#   # browser("mystop")
#   for (i in 1:num_replicas){
#     print( i)
#     result <- cfd_testing(start_time, end_time, timeseries_length,
#                           num_indvs,mu1_coef,mu2_coef,fl_choice,response_family,test_type, klen=3)
#     p_value[i]=result$pvalue
#     test_stats[i]=result$test_statistics
#   }
# 
#   return(list("pvalue"=p_value,"teststat"=test_stats))
# }


source("./source_code/R/time_track_function.R")
run_experiment_hypothesis <- function(exp_idx,
                                      num_indvs,
                                      timeseries_length,
                                      fl_choice,
                                      test_type,
                                      gam_choice,
                                      num_replicas = 5000,
                                      alpha = 0.05){
  
  # mu1_coef=c(-6.67,-2.47,5.42)
  # mu2_coef=c(-3.14,-0.99,3.91)
  
  mu1_coef=c(-1.8270644 ,-2.4700275,  5.4299181)
  mu2_coef=c(-2.9990822, -0.8243365,  3.9100000  )
  exp_str <- paste("Track time for \nNum Subjects:\t", num_indvs,
        "\n timeserires_length:\t",timeseries_length,
        "\n fl_choice:\t",fl_choice,
        "\n test_type:\t",test_type)
  writeLines(exp_str)
  timeKeeperStart(exp_str)
  #timeseries_length=90
  #num_indvs=100
  # fl_choice=6
  # test_type="Inclusion"
  # gam_choice=0
  #num_replicas=5
  simulation_scenarios <- cfd_testing_simulation(num_replicas=num_replicas, start_time=0.01, end_time=0.99,
                                           timeseries_length=timeseries_length,
                                           mu1_coef=mu1_coef,
                                           mu2_coef=mu2_coef,
                                           num_indvs=num_indvs,fl_choice=fl_choice,
                                           response_family='bernoulli',test_type=test_type,gam_choice=gam_choice,
                                           klen=3)

# n100t2000_mu1mu2 <- cfd_testing_simulation_no_paralel(num_replicas=4, start_time=0.01, end_time=0.99,
#                                            timeseries_length=2000,
#                                            mu1_coef=c(1,2,3),
#                                            mu2_coef=c(4,5,6,7,8,9,10,11,12),
#                                            num_indvs=100,fl_choice=3,
#                                            response_family='bernoulli',test_type='Functional',
#                                            klen=3)
  #simulation_pvalues <- matrix(unlist(simulation_scenarios), nrow=4) #when gam has z1, z2 pvalue, it's 4 rows
  #simulation_pvalues <- matrix(unlist(simulation_scenarios), nrow=3)
  #add x2fl2-not adding, since it's a vector itself
  simulation_pvalues <- matrix(unlist(simulation_scenarios), nrow=3)
  save(simulation_pvalues, file = paste0("./outputsnongampower/simpvals3",
                                                    "_i", exp_idx,
                                                    "_fl", fl_choice,
                                                    "_ttype", test_type,
                                                    "_n", num_indvs,
                                                    "_tlen", timeseries_length,
                                                    ".RData"))
  non_null_count <- dim(simulation_pvalues)[2]
  power <- mean(simulation_pvalues[1,] < alpha)
  power_se <- sqrt(power*(1-power)/non_null_count)
  ############
  power_01 <- mean(simulation_pvalues[1,] < 0.1)
  power_se01 <- sqrt(power_01*(1-power_01)/non_null_count)
  ##############
  yip_mean= mean(simulation_pvalues[2,])
  yip_sd= sd(simulation_pvalues[2,])/sqrt(non_null_count)
  yip_wo_mean= mean(simulation_pvalues[3,])
  yip_wo_sd= sd(simulation_pvalues[3,])/sqrt(non_null_count)
  ################
  # x2fl2 <- mean(simulation_pvalues[4,])
  # x2fl2_se <- sd(simulation_pvalues[4,])/non_null_count
  ##############
  # power2 <- mean(simulation_pvalues[4,] < alpha)
  # power2_se <- sqrt(power2*(1-power2)/non_null_count)
  # power_012 <- mean(simulation_pvalues[4,] < 0.1)
  # power_se012 <- sqrt(power_012*(1-power_012)/non_null_count)
  ############
  # cat("\npower:", power,"\n", "power_se:", power_se, "\n")
  timeKeeperNext()
  # return(list("power"=power,"se"=power_se,"power_01"=power_01 ,"se01"=power_se01,
  #             "yip_mean"=yip_mean,"yip_wo_mean"=yip_wo_mean,"yip_sd"=yip_sd,
  #             "yip_wo_sd"=yip_wo_sd,"power2"=power2,"se2"=power2_se,"power_012"=power_012 ,
  #             "se012"=power_se012,"NAs"=num_replicas - non_null_count))
  
  return(list("power"=power,"se"=power_se,"power_01"=power_01 ,"se01"=power_se01,
              "yip_mean"=yip_mean,"yip_wo_mean"=yip_wo_mean,"yip_sd"=yip_sd,
              "yip_wo_sd"=yip_wo_sd,"NAs"=num_replicas - non_null_count))
}
# 
 # run_experiment_hypothesis( 0,
 #                            100,
 #                            300,
 #                           "11",
 #                            "Functional",
 #                             gam_choice=0,
 #                           num_replicas = 5,
 #                            alpha = 0.05 )

begin_exp_time <- Sys.time()

set.seed(123456)


generate_ed_table <- function(subjects_vector = c(500,300,100),
                               time_length_vector = c(180,90),
                               fl_choice_vector = c("6"),
                               test_type_vector = c("Inclusion", "Functional")){
  ed_table_ret <- expand.grid(fl_choice_vector, test_type_vector, subjects_vector, time_length_vector)
  return(ed_table_ret)
}

########
#type I error rate
# ed_table1=generate_ed_table()
# ed_table2=generate_ed_table(fl_choice_vector = c("200","7","21"),
#                                                          test_type_vector = c("Functional"))
# 
# 
# ed_table <- rbind(ed_table1,ed_table2)

###################
#power
ed_table1 <- generate_ed_table(subjects_vector = c(1000,500,300,100),
                                fl_choice_vector = c("6","7", "8","9","10"),
                               test_type_vector = c("Inclusion"))
ed_table2 <- generate_ed_table(subjects_vector = c(1000,500,300,100),
                             fl_choice_vector = c("21","22","23","24","25"))
ed_table <- rbind(ed_table1,ed_table2)
###################

colnames(ed_table) <- c("fl_choice", "test_type", "num_subjects", "num_timepoints")
########################
gam_choice=0
##########################
all_experiment_outputs <- list()
for (row_index in 1:dim(ed_table)[1]){
  num_indvs <- ed_table[row_index,]$num_subjects
  timeseries_length <- ed_table[row_index,]$num_timepoints
  fl_choice <- as.character(ed_table[row_index,]$fl_choice)
  test_type <- as.character(ed_table[row_index,]$test_type)
  experiment_output <- run_experiment_hypothesis( row_index,
                                                  num_indvs , 
                                                  timeseries_length,
                                                  fl_choice,
                                                  test_type,
                                                  gam_choice)
  save(experiment_output, file = paste0("./outputsnongampower/exp3_", 
                                       "_i", row_index, 
                                       "_fl", fl_choice, 
                                       "_ttype", test_type, 
                                       "_n", num_indvs, 
                                       "_tlen", timeseries_length,
                                       ".RData"))
  all_experiment_outputs <- rbind(all_experiment_outputs, experiment_output)
}

final_table <- cbind(ed_table, all_experiment_outputs)

mu1_coef=c(-1.8270644 ,-2.4700275,  5.4299181)
mu2_coef=c(-2.9990822, -0.8243365,  3.9100000  )
save(final_table,mu1_coef,mu2_coef,gam_choice,file = "EXP3_r5000_cfdanongampower.RData")

end_exp_time <- Sys.time()

cat("\n====================\n",
    "\tAll Experiemnts Took:", capture.output(end_exp_time - begin_exp_time), 
    "\n====================\n")

if(run_parallel)
{
  parallel::stopCluster(cl = my.cluster)
  initialized_parallel <- FALSE
}
