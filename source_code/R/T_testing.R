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
library(mgcv)
library(fda)
library(fda.usc)
#library(devtools)
# install_github("stchen3/glmmVCtest")
#library("glmmVCtest")
#library(RLRsim)
library(MASS)
library(splines)
library(parallel)
library(stats)

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





# cfd_T_testing_simulation <- function (num_replicas, start_time, end_time, timeseries_length,
#                                     mu1_coef, mu2_coef,
#                                     num_indvs,fl_choice,
#                                     klen=3){
#     cat("CFD Testing Simulation \nNum Replicas:\t", num_replicas, "\nNum indvs:\t", num_indvs)
#     #num_replicas=2
#     time_interval=seq(start_time,end_time,length.out=timeseries_length)
#     
#     result_all <- foreach (number_simulation = 1:num_replicas, .combine = cbind, .init = NULL,
#                            .packages=c("splines","mgcv","fda","fda.usc","MASS","stats")) %dorng% {
#                                # result_all <- foreach (number_simulation = 1:num_replicas, .combine = cbind, .init = NULL) %dorng% {
#                                source("./source_code/R/T_testing_functions.R")
#                                result <- cfd_T_testing(klen, mu1_coef,mu2_coef,num_indvs, timeseries_length,
#                                                        time_interval, fl_choice, lp_intercept=0.9998364)
#                                
#                                #return(list("pvalue"=result$pvalue,"teststat"=result$test_statistics,"fl"=result$flt))
#                                #return(list("pvalue"=result$pvalue,"yip"=result$yip,"yip_wo"=result$yip_wo,"pvalue2"=result$pvalue2))
#                                return(list("pvalue"=result$pvalue,"rvmean"=result$rvmean,"rvemean"=result$rvemean))
#                            } 
#     return(result_all)
#     
# }

cfd_T_testing_simulation=function(klen, mu1_coef,mu2_coef,num_indvs, timeseries_length,
                           time_interval, fl_choice,num_replicas, 
                           lp_intercept=0.9998364,boot_number=99){
    T_rep <- foreach(this_row = 1:num_replicas ) %dorng%
        { source("./source_code/R/data_generator.R")
            source("./source_code/R/integral_penalty_function.R")
            source("./source_code/R/T_testing_functions.R")
            
            number_basis =30
            
            
            #T_rv_erv <- list()
            WY_sample=GenerateCategoricalFDTest(klen, mu1_coef,mu2_coef,num_indvs, timeseries_length,
                                                time_interval, fl_choice, lp_intercept=0.9998364)
            
           
            #W is t*n
            #categFD_est <- EstimateCategFuncDataX(est_choice, time_interval, WY_sample$true$Truecatcurve)
        
            temp=get_T(WY_sample$true$TrueX1,WY_sample$true$TrueX2,WY_sample$true$TrueX3, 
                       WY_sample$true$yis,time_interval,
                       number_basis =number_basis,est_choice="binomial" )
            T_stat=array(0,3)
            T_stat[1]=temp$T_statistics #scalar
            betals=temp$betals
            ####################
            #T_star_series=c(0)
            ####################
            #########
            #get Y from X, and betals, betals 1: intercept, 2:31, 32:62
            get_Y_star=function( X_2t,X_3t,betals,time_interval,num_indvs,number_basis){
                
                vec <- matrix(1:num_indvs, nrow=num_indvs, ncol=1)
                beta0=betals[1]
                betal=betals[2:(number_basis+1)]*0
                betal3=betals[(number_basis+2):(2*number_basis+1)]
                
                knots <- construct.knots(time_interval,knots=(number_basis-3),knots.option='equally-spaced')
                bspline <- splineDesign(knots=knots,x=time_interval,ord=4)
                
                x1fl1 <- rep(beta0,num_indvs)
                x2fl2 <- apply(vec, 1, function(x) {fda.usc::int.simpson2(time_interval, X_2t[x,]*(bspline%*%betal), equi = TRUE, method = "TRAPZ")})
                x3fl3 <- apply(vec, 1, function(x) {fda.usc::int.simpson2(time_interval, X_3t[x,]*(bspline%*%betal3), equi = TRUE, method = "TRAPZ")})
                
                linear_predictor <- matrix(x1fl1 + x2fl2+ x3fl3 )
                y_star <- apply(linear_predictor, 1, function(x){ rbinom(1, 1, 1/(1+exp(-x))) })
                return(y_star)
            }
            # y_star=get_Y_star(WY_sample$true$TrueX2,WY_sample$true$TrueX3,
            #                   betals,time_interval,num_indvs,number_basis)
            
            
            
            #bootstrap
            ################
            #####################################################################
            #temp_series <- foreach(this_col = 1:num_replicas ) %do%
            #     temp_series <- foreach(this_col = 1:10 ) %do%
            #     {
            #         source("./source_code/R/integral_penalty_function.R")
            #         source("./source_code/R/T_testing_functions.R")
            #         boot_index=sample(1:num_indvs, num_indvs,replace=T)
            # 
            #         temp[this_col ]=get_T(WY_sample$true$Truecatcurve[,boot_index],
            #                               WY_sample$true$yis[boot_index],time_interval,
            #                               number_basis =30,est_choice="binomial")$T_statistics
            # 
            #         return(temp[this_col ])
            #     }
            # temp_series  <- do.call(rbind, temp_series)
            # T_stat[2]=(T_stat<=quantile(unlist(temp_series), .05))[[1]]
            # T_stat[3]=(T_stat<=quantile(unlist(temp_series), .10))[[1]]
            # start_time_boot=Sys.time()
             temp_series=c(0)
            for (this_col in 1:boot_number){
                
                boot_index=sample(1:num_indvs, num_indvs,replace=T)
                y_star=get_Y_star(WY_sample$true$TrueX2[boot_index,],
                                  WY_sample$true$TrueX3[boot_index,],
                                  betals,time_interval,num_indvs,number_basis)
                temp_series[this_col ]=get_T(WY_sample$true$TrueX1[boot_index,],
                                             WY_sample$true$TrueX2[boot_index,],
                                             WY_sample$true$TrueX3[boot_index,],
                                             y_star,time_interval,
                                      number_basis =number_basis,est_choice="binomial")$T_statistics
            }
            # end_time_boot=Sys.time()
            # cat("boot 1000 for 500 useres takes", end_time_boot-start_time_boot)
            # boot 1000 for 500 useres takes 50.5146
            ###############
            # T_stat[2]=(T_stat<=quantile(temp_series, .05))[[1]]
            # T_stat[3]=(T_stat<=quantile(temp_series, .10))[[1]]
            
            T_stat[2]=(T_stat>=quantile(temp_series, .95))[[1]]
            T_stat[3]=(T_stat>=quantile(temp_series, .90))[[1]]
            ################
            #T_star_series=temp_series
            ################
            # T_rv_erv[2]=temp$rv_XF #1D vector
            # T_rv_erv[3]=temp$rv_E_PF #scalar
            ######save T star series as well
            return(T_stat)
            ##############################
        }
    T_rep <- do.call(rbind, T_rep)
    #three columns, T, and T_binary, T_binary0.1
    return(T_rep)
}

source("./source_code/R/time_track_function.R")
run_experiment_hypothesis <- function(exp_idx,
                                      num_indvs,
                                      timeseries_length,
                                      fl_choice,
                                      num_replicas = 100,
                                      alpha = 0.05, 
                                      start_time=0.01,
                                      end_time=0.99,
                                      klen=3){
    
    mu1_coef=c(-1.8270644 ,-2.4700275,  5.4299181)
    mu2_coef=c(-2.9990822, -0.8243365,  3.9100000  )
    exp_str <- paste("Track time for \nNum Subjects:\t", num_indvs,
                     "\n timeserires_length:\t",timeseries_length,
                     "\n fl_choice:\t",fl_choice
                    
                     )
    writeLines(exp_str)
    timeKeeperStart(exp_str)
    time_interval=seq(start_time,end_time,length.out=timeseries_length)
    simulation_scenarios <- cfd_T_testing_simulation (klen, mu1_coef,mu2_coef,num_indvs, timeseries_length,
                                                      time_interval, fl_choice,num_replicas, lp_intercept=0.9998364)
    #simulation_pvalues <- matrix(unlist(simulation_scenarios), nrow=3)
    save(simulation_scenarios, file = paste0("./outputsTbootstrap/simpvals3",
                                           "_i", exp_idx,
                                           "_fl", fl_choice,
                                           "_n", num_indvs,
                                           "_tlen", timeseries_length,
                                           ".RData"))
    #non_null_count <- dim(simulation_pvalues)[2]
    power <- mean(simulation_scenarios[,2] )
    power_se <- sqrt(power*(1-power)/num_replicas)
    ############
    power_01 <- mean(simulation_scenarios[,3] )
    power_se01 <- sqrt(power_01*(1-power_01)/num_replicas)
    ##############
    T_rv= simulation_scenarios[,1]
    #rv_sd= sd(simulation_pvalues[2,])/sqrt(non_null_count)
    # rve_mean= mean(simulation_pvalues[3,])
    # rve_sd= sd(simulation_pvalues[3,])/sqrt(non_null_count)
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
    
    # return(list("power"=power,"se"=power_se,"power_01"=power_01 ,"se01"=power_se01,
    #             "rv_mean"=rv_mean,"rv_sd"=rv_sd,"rve_mean"=rve_mean,
    #             "rve_sd"=rve_sd,"NAs"=num_replicas - non_null_count))
    return(list("power"=power,"se"=power_se,"power_01"=power_01 ,"se01"=power_se01,
                "T_rv"=T_rv))
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
ed_table1 <- generate_ed_table(subjects_vector = c(500),
                               fl_choice_vector = c("6"),
                               time_length_vector = c(90),
                               test_type_vector = c("Inclusion"))
# ed_table2=generate_ed_table(subjects_vector = c(500),fl_choice_vector = c("200","7","21"),time_length_vector = c(90),
#                                                          test_type_vector = c("Functional"))


#ed_table <- rbind(ed_table1,ed_table2)
ed_table <- ed_table1
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
                                                    fl_choice
                                                    )
    save(experiment_output, file = paste0("./outputsTbootstrap/exp3_", 
                                          "_i", row_index, 
                                          "_fl", fl_choice, 
                            
                                          "_n", num_indvs, 
                                          "_tlen", timeseries_length,
                                          ".RData"))
    all_experiment_outputs <- rbind(all_experiment_outputs, experiment_output)
}

final_table <- cbind(ed_table, all_experiment_outputs)
final_table_pvalue=final_table[1:8]
final_table_rv=final_table[9]
#hist(final_table_rv$T_rv[[1]])
mu1_coef=c(-1.8270644 ,-2.4700275,  5.4299181)
mu2_coef=c(-2.9990822, -0.8243365,  3.9100000  )
save(final_table_pvalue,final_table_rv,mu1_coef,mu2_coef,file = "EXP3_outputsTbootstrap.RData")

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
