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
# Date: 02/28/2023
#
##############################################################


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

###test to get T


# start_time=0.01
# end_time=0.99
# timeseries_length=90
# time_interval=seq(start_time,end_time,length=timeseries_length)
# num_indvs=500
# #fl_choice=3 #not constant, expect to reject
# #fl_choice=8 #not constant, expect to reject
# 
# fl_choice=6
# klen=3
# mu1_coef=c(-1.8270644 ,-2.4700275,  5.4299181)
# mu2_coef=c(-2.9990822, -0.8243365,  3.9100000  )
# 
# 
# WY_sample=GenerateCategoricalFDTest(klen, mu1_coef,mu2_coef,num_indvs, timeseries_length,
#                                       time_interval, fl_choice, lp_intercept=0.9998364)
# 
# num_indvs=1000
# T_example=get_T(WY_sample$true$Truecatcurve, WY_sample$true$yis,time_interval, 
#                 number_basis =30,est_choice="binomial" )
# 
# 
# library(MASS)
# fit_data=c(T_example$T_vector)
# fit_result <- fitdistr(fit_data, "log-normal")
# estimated_df <- fit_result$estimate
# print(estimated_df)
# 
# num_indvs=500
# 
# get_T_for_hist=function(klen, mu1_coef,mu2_coef,num_indvs, timeseries_length,
#                         time_interval, fl_choice, lp_intercept=0.9998364){
#     WY_sample=GenerateCategoricalFDTest(klen, mu1_coef,mu2_coef,num_indvs, timeseries_length,
#                                         time_interval, fl_choice, lp_intercept=0.9998364)
#     
#     T_example=get_T(WY_sample$true$Truecatcurve, WY_sample$true$yis,time_interval, 
#                     number_basis =30,est_choice="binomial" )
#     fit_data=c(T_example$T_vector)
#     #fit_result <- fitdistr(fit_data, "log-normal")
#     fit_result <- fitdistr(fit_data, "chi-squared", list(df=3), lower = 0.001)
#     estimated_df <- fit_result$estimate
#     print(estimated_df)
#     
#     #####
#     fit_datap=c(T_example$T_vectorp)
#     #fit_result <- fitdistr(fit_data, "log-normal")
#     fit_resultp <- fitdistr(fit_datap, "chi-squared", list(df=3), lower = 0.001)
#     estimated_dfp <- fit_resultp$estimate
#     print(estimated_dfp)
#     ##
#     two_sample_result <- t.test(T_example$rv_XF,T_example$rv_E_PF, var.equal = FALSE)
#     print(two_sample_result)
#     
#     par(mfrow=c(1,2))
#     hist(T_example$T_vector,xlab="T",main=paste0("Historgram of T for", num_indvs," Subjects"))
#     hist(T_example$T_vectorp,xlab="T",main=paste0("Historgram of Tp for", num_indvs," Subjects"))
#     
# }
# num_indvs=100
# get_T_for_hist(klen, mu1_coef,mu2_coef,num_indvs, timeseries_length,
#                         time_interval, fl_choice, lp_intercept=0.9998364)
# num_indvs=500
# get_T_for_hist(klen, mu1_coef,mu2_coef,num_indvs, timeseries_length,
#                time_interval, fl_choice, lp_intercept=0.9998364)
# num_indvs=1000
# get_T_for_hist(klen, mu1_coef,mu2_coef,num_indvs, timeseries_length,
#                time_interval, fl_choice, lp_intercept=0.9998364)




cfd_T_testing_simulation <- function (num_replicas, start_time, end_time, timeseries_length,
                                    mu1_coef, mu2_coef,
                                    num_indvs,fl_choice,
                                    klen=3){
    cat("CFD Testing Simulation \nNum Replicas:\t", num_replicas, "\nNum indvs:\t", num_indvs)
    #num_replicas=2
    time_interval=seq(start_time,end_time,length.out=timeseries_length)
    
    result_all <- foreach (number_simulation = 1:num_replicas, .combine = cbind, .init = NULL,
                           .packages=c("splines","mgcv","fda","fda.usc","MASS","stats")) %dorng% {
                               # result_all <- foreach (number_simulation = 1:num_replicas, .combine = cbind, .init = NULL) %dorng% {
                               source("./source_code/R/T_testing_functions.R")
                               result <- cfd_T_testing(klen, mu1_coef,mu2_coef,num_indvs, timeseries_length,
                                                       time_interval, fl_choice, lp_intercept=0.9998364)
                               
                               #return(list("pvalue"=result$pvalue,"teststat"=result$test_statistics,"fl"=result$flt))
                               #return(list("pvalue"=result$pvalue,"yip"=result$yip,"yip_wo"=result$yip_wo,"pvalue2"=result$pvalue2))
                               return(list("pvalue"=result$pvalue,"rvmean"=result$rvmean,"rvemean"=result$rvemean))
                           } 
    return(result_all)
    
}


source("./source_code/R/time_track_function.R")
run_experiment_hypothesis <- function(exp_idx,
                                      num_indvs,
                                      timeseries_length,
                                      fl_choice,
                                      num_replicas = 5000,
                                      alpha = 0.05, 
                                      start_time=0.01,
                                      end_time=0.99){
    
    mu1_coef=c(-1.8270644 ,-2.4700275,  5.4299181)
    mu2_coef=c(-2.9990822, -0.8243365,  3.9100000  )
    exp_str <- paste("Track time for \nNum Subjects:\t", num_indvs,
                     "\n timeserires_length:\t",timeseries_length
                    
                     )
    writeLines(exp_str)
    timeKeeperStart(exp_str)

    simulation_scenarios <- cfd_T_testing_simulation (num_replicas, start_time, end_time, timeseries_length,
                                                      mu1_coef, mu2_coef,
                                                      num_indvs,fl_choice,
                                                      klen=3)
    simulation_pvalues <- matrix(unlist(simulation_scenarios), nrow=3)
    save(simulation_pvalues, file = paste0("./outputsT/simpvals3",
                                           "_i", exp_idx,
                                           "_fl", fl_choice,
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
    rv_mean= mean(simulation_pvalues[2,])
    rv_sd= sd(simulation_pvalues[2,])/sqrt(non_null_count)
    rve_mean= mean(simulation_pvalues[3,])
    rve_sd= sd(simulation_pvalues[3,])/sqrt(non_null_count)
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
                "rv_mean"=rv_mean,"rv_sd"=rv_sd,"rve_mean"=rve_mean,
                "rve_sd"=rve_sd,"NAs"=num_replicas - non_null_count))
}
# 
# run_experiment_hypothesis (0,
#                                      100,
#                                       90,
#                                       6,
#                                       num_replicas = 5,
#                                       alpha = 0.05)

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
ed_table1 <- generate_ed_table(
                               fl_choice_vector = c("6"),
                               time_length_vector = c(90),
                               test_type_vector = c("Inclusion", "Functional"))
ed_table2=generate_ed_table(fl_choice_vector = c("200","7","21"),time_length_vector = c(90),
                                                         test_type_vector = c("Functional"))


ed_table <- rbind(ed_table1,ed_table2)

###################
#power
# ed_table1 <- generate_ed_table(subjects_vector = c(500,300,100),
#                                fl_choice_vector = c("6","7", "8","9","10"),
#                                time_length_vector = c(90),
#                                test_type_vector = c("Inclusion"))
# ed_table2 <- generate_ed_table(subjects_vector = c(1000,500,300,100),
#                                fl_choice_vector = c("21","22","23","24","25"))
# ed_table <- rbind(ed_table1,ed_table2)
# ###################

colnames(ed_table) <- c("fl_choice", "test_type", "num_subjects", "num_timepoints")
########################
#gam_choice=0
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
                                                    test_type)
    save(experiment_output, file = paste0("./outputsT/exp3_", 
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
save(final_table,mu1_coef,mu2_coef,file = "EXP3_outputsT.RData")

end_exp_time <- Sys.time()

cat("\n====================\n",
    "\tAll Experiemnts Took:", capture.output(end_exp_time - begin_exp_time), 
    "\n====================\n")


# @param W 2D array, t*n: t is the timestamp and n is the number of the observation
# @return X 3D array, n*t*Q, Q: the total number of the category
# GetXFromW <- function(W)




if(run_parallel)
{
    parallel::stopCluster(cl = my.cluster)
    initialized_parallel <- FALSE
}
