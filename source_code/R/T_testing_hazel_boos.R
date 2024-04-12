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
#library(pracma)


###########
library(optparse)

# Define options
option_list <- list(
    make_option(c("-j", "--jobid"), type="integer", default=123,
                help="Job Index", metavar="JOBID"),
    make_option(c("-n", "--numcpus"), type="integer", default=32,
                help="Num CPUs", metavar="NUMCPUS"),
    make_option(c("-r", "--replicas"), type="integer", default=100,
                help="Num Replicas", metavar="NUMREPLICAS"),
    make_option(c("-b", "--bootscov"), type="integer", default=100,
                help="Num Bootstraps", metavar="NUMBOOTS"),
    make_option(c("-m", "--bootststar"), type="integer", default=99,
                help="Num Bootstraps", metavar="NUMBOOTS"),
    make_option(c("-s", "--subjects"), type="integer", default=100,
                help="Num Subjects/Individuals", metavar="NUMSUBJECTS")
)

# Create parser and parse options
parser <- OptionParser(option_list=option_list)
options <- parse_args(parser)

options_jobid <- options$jobid
options_numcpus <- options$numcpus
options_replicas <- options$replicas
options_bootscov <- options$bootscov
options_bootststar <- options$bootststar
options_subjects <- options$subjects

# options_jobid <- 1
# options_numcpus <- 3
# options_replicas <- 2
# options_bootscov <- 5
# options_bootststar <- 4
# options_subjects <- 100

# Use the options
cat("Job Idx:", options_jobid, "\n")
cat("Num CPUs:", options_numcpus, "\n")
cat("Num Replicas:", options_replicas, "\n")
cat("Num Bootstraps for Covariance:", options_bootscov, "\n")
cat("Num Bootstraps for T star:", options_bootststar, "\n")
cat("Num Subjects:", options_subjects, "\n")

###########
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
    # n.cores <- parallel::detectCores()
    n.cores <- options_numcpus
    my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
    doParallel::registerDoParallel(cl = my.cluster)
    cat("Parellel Registered: ", foreach::getDoParRegistered(), " (num cores=", n.cores, ")\n")
    initialized_parallel <- TRUE
    
    # registerDoRNG(123) # ///<<<< THIS CREATES THE ERROR FOR FADPClust !!!
}

ensure_dir_exist <- function(directory_path){
    # Check if the directory exists
    if(!dir.exists(directory_path)) {
        # Directory doesn't exist, so create it
        dir.create(directory_path, recursive = TRUE)
        cat("Directory created:", directory_path, "\n")
    } else {
        cat("Directory already exists:", directory_path, "\n")
    }
}

scenario_folder = "outputsTbootstrap_p_boos"
ensure_dir_exist(scenario_folder)

final_table_folder = paste0("final_table_output_p_boos_n", options_subjects)
ensure_dir_exist(final_table_folder)


cfd_T_testing_simulation=function(klen, mu1_coef,mu2_coef,num_indvs, timeseries_length,
                                  time_interval, fl_choice,num_replicas, 
                                  lp_intercept=0.9998364,boot_1=options_bootscov,boot_2=options_bootststar,
                                  shuffle_option=TRUE){
    #num_replicas=6
    p_values <- foreach(pval_idx = 1:num_replicas, .combine = 'c') %dorng%
        
        #T_rep <- foreach(this_row = 1:5) %dorng%
        { source("./source_code/R/data_generator.R")
            source("./source_code/R/integral_penalty_function.R")
            source("./source_code/R/T_testing_functions.R")
            
            number_basis =30

            WY_sample=GenerateCategoricalFDTest(klen, mu1_coef,mu2_coef,num_indvs, timeseries_length,
                                                time_interval, fl_choice, lp_intercept=0.9998364)
            X_1t=WY_sample$true$TrueX1
            X_2t=WY_sample$true$TrueX2
            X_3t=WY_sample$true$TrueX3
            Y=WY_sample$true$yis #time_interval
            pval <- calculate_double_boot_pvalue(X_2t, X_3t, Y, 
                                                 time_interval, 
                                                 boot_1, boot_2, 
                                                 number_basis =30, 
                                                 category_count=3)
            return(pval)
            
        }
    return(p_values)
}


# calculate_new_T <- function(X_1t, X_2t, X_3t, Y, time_interval, number_basis=30, 
#                             est_choice, category_count=3, 
#                             replicas=1000, boot_1=100,  boot_2=99){
#     p_values <- foreach(pval_idx = 1:replicas, .combine = 'c') %do% 
#         {
#             pval <- calculate_double_boot_pvalue(X_2t, X_3t, Y, 
#                                                  time_interval, 
#                                                  boot_1, boot_2, 
#                                                  number_basis =30, 
#                                                  category_count=3)
#             return(pval)
#         }
#     
#     return(p_values)
#     
#     # power_005 <- mean(p_values < 0.05)
#     # stderr_005 <- sqrt(power * (1-power) / replicas)
#     # power_01 <- mean(p_values < 0.1)
#     # stderr_01 <- sqrt(power_01 * (1-power_01) / replicas)
#     # 
#     # return(c(power_005, stderr_005, power_01, stderr_01))
# }


source("./source_code/R/time_track_function.R")
run_experiment_hypothesis <- function(exp_idx,
                                      num_indvs,
                                      timeseries_length,
                                      fl_choice,
                                      num_replicas = options_replicas,
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
    save(simulation_scenarios, file = paste0("./", scenario_folder, "/simpvals3",
                                             "_i", exp_idx,
                                             "_fl", fl_choice,
                                             "_n", num_indvs,
                                             "_tlen", timeseries_length,
                                             "_",options_numcpus,
                                             "_",options_jobid,
                                             ".RData"))
    
    
    power_005 <- mean(simulation_scenarios < 0.05)
    power_01 <- mean(simulation_scenarios < 0.1) 
    
    # power2 <- simulation_scenarios[,5] 
    # power_012 <- simulation_scenarios[,6] 
    stderr_005 <- sqrt(power_005 * (1-power_005) / num_replicas)
     
   stderr_01 <- sqrt(power_01 * (1-power_01) / num_replicas)
    
    ##############
    #T_rv= simulation_scenarios[,1]
    #T_rv2= simulation_scenarios[,4]
    
    ############
    # cat("\npower:", power,"\n", "power_se:", power_se, "\n")
    timeKeeperNext()
    
    # return(list("power"=power,"power_01"=power_01,
    #             "power2"=power2,"power_012"=power_012,
    #             "T_rv"=T_rv,"T_rv2"=T_rv2))
    
    # 
    # return(list("power"=power,"power_01"=power_01,
    #             
    #             "T_rv"=T_rv))
    
    return(list("power_005"=power_005 ,"stderr_005"= stderr_005,
                "power_01"= power_01,"stderr_01"=stderr_01
                
                ))
}
# 
# run_experiment_hypothesis (0,
#                                      100,
#                                       90,
#                                       6,
#                                       num_replicas = 5,
#                                       alpha = 0.05)

begin_exp_time <- Sys.time()

set.seed(123456 + 10 * options_jobid)


generate_ed_table <- function(subjects_vector = c(500,300,100),
                              time_length_vector = c(180,90),
                              fl_choice_vector = c("6"),
                              test_type_vector = c("Inclusion", "Functional")){
    ed_table_ret <- expand.grid(fl_choice_vector, test_type_vector, subjects_vector, time_length_vector)
    return(ed_table_ret)
}

########
#type I error rate
ed_table1 <- generate_ed_table(subjects_vector = c(options_subjects),
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
    save(experiment_output, file = paste0("./", scenario_folder, "/exp3_", 
                                          "_i", row_index, 
                                          "_fl", fl_choice, 
                                          
                                          "_n", num_indvs, 
                                          "_tlen", timeseries_length,
                                          "_", options_numcpus,
                                          "_", options_jobid,
                                          ".RData"))
    all_experiment_outputs <- rbind(all_experiment_outputs, experiment_output)
}

final_table <- cbind(ed_table, all_experiment_outputs)
# power_vector=final_table$power
# final_rv=final_table$T_rv
#hist(final_table_rv$T_rv[[1]])
mu1_coef=c(-1.8270644 ,-2.4700275,  5.4299181)
mu2_coef=c(-2.9990822, -0.8243365,  3.9100000  )
save(final_table,file =paste0("./", final_table_folder, "/Hazel_outputsTbootstrap_",
                              "_", options_subjects,
                              "_", options_replicas,
                              "_", options_bootscov,
                              "_", options_bootststar,
                              "_", options_numcpus,
                              "_", options_jobid, ".RData"))

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
