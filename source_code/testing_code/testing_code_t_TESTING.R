#testing code for two sample T testing

###test to get T
# quantile(temp_series,0.95)
# 95% 
# 0.0002349054 
####example code
#source("./source_code/R/T_testing_functions.R")
start_time=0.01
end_time=0.99
timeseries_length=90
time_interval=seq(start_time,end_time,length.out=timeseries_length)
#num_indvs=500
#fl_choice=3 #not constant, expect to reject
#fl_choice=8 #not constant, expect to reject

#fl_choice=6
#fl_choice="25"
klen=3
mu1_coef=c(-1.8270644 ,-2.4700275,  5.4299181)
mu2_coef=c(-2.9990822, -0.8243365,  3.9100000  )

#est_choice="binomial"
get_T_simulations=function(klen, mu1_coef,mu2_coef,num_indvs, timeseries_length,
                           time_interval, fl_choice,num_replications,boot_number, lp_intercept=0.9998364){
    #T_rep <- foreach(this_row = 1:num_replications ) %dorng%
        T_rep <- foreach(this_row = 1:num_replications ) %do%
        { source("./source_code/R/data_generator.R")
            source("./source_code/R/integral_penalty_function.R")
            source("./source_code/R/T_testing_functions.R")
            #T_rv_erv <- list()
            number_basis=30
            WY_sample=GenerateCategoricalFDTest(klen, mu1_coef,mu2_coef,num_indvs, timeseries_length,
                                                time_interval, fl_choice, lp_intercept=0.9998364)
            
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
               betal=betals[2:(number_basis+1)]
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
           #bootstrap
           ################
           # temp_series=c(0)
           # start_time_boot=Sys.time()
           # for (this_col in 1:boot_number){
           #     
           #     boot_index=sample(1:num_indvs, num_indvs,replace=T)
           #     y_star=get_Y_star(WY_sample$true$TrueX2[boot_index,],
           #                       WY_sample$true$TrueX3[boot_index,],
           #                       betals,time_interval,num_indvs,number_basis)
           #     temp_series[this_col ]=get_T(WY_sample$true$TrueX1[boot_index,],
           #                                  WY_sample$true$TrueX2[boot_index,],
           #                                  WY_sample$true$TrueX3[boot_index,],
           #                                  y_star,time_interval,
           #                                  number_basis =number_basis,est_choice="binomial")$T_statistics
           # }
           # end_time_boot=Sys.time()
           # cat("boot 1000 for 500 useres takes", end_time_boot-start_time_boot)
           #boot 1000 for 500 useres takes 1.598265
           
           
           start_time_boot=Sys.time()
           temp_series=foreach(this_col = 1:boot_number) %dorng%{
               
               source("./source_code/R/data_generator.R")
               source("./source_code/R/integral_penalty_function.R")
               source("./source_code/R/T_testing_functions.R")
               
               get_Y_star=function( X_2t,X_3t,betals,time_interval,num_indvs,number_basis){
                   
                   vec <- matrix(1:num_indvs, nrow=num_indvs, ncol=1)
                   beta0=betals[1]
                   betal=betals[2:(number_basis+1)]
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
               
               boot_index=sample(1:num_indvs, num_indvs,replace=T)
               y_star=get_Y_star(WY_sample$true$TrueX2[boot_index,],
                                 WY_sample$true$TrueX3[boot_index,],
                                 betals,time_interval,num_indvs,number_basis)
               temp_series_this_col=get_T(WY_sample$true$TrueX1[boot_index,],
                                            WY_sample$true$TrueX2[boot_index,],
                                            WY_sample$true$TrueX3[boot_index,],
                                            y_star,time_interval,
                                            number_basis =number_basis,est_choice="binomial")$T_statistics
               return(temp_series_this_col)
           }
           
           end_time_boot=Sys.time()
           cat("boot 1000 for 500 useres takes", end_time_boot-start_time_boot)
           #boot 1000 for 500 useres takes 8.015576 
           #hist(unlist(temp_series))
           
           # quantile(temp_series,0.95)
           # 95% 
           # 0.0002743217 
           
           # quantile(unlist(temp_series),0.95)
           # 95% 
           # 0.0002993409 
           
           T_stat[2]=(T_stat>=quantile(temp_series, .95))[[1]]
           T_stat[3]=(T_stat>=quantile(temp_series, .90))[[1]]
           # T_rv_erv[2]=temp$rv_XF #1D vector
           # T_rv_erv[3]=temp$rv_E_PF #scalar
           return(T_stat)
        }
    T_rep <- do.call(rbind, T_rep)
    #three columns, T, and T_binary, T_binary0.1
    return(T_rep)
}
    
fl_choice="6"
num_indvs=500
number_basis=30
boot_number=1000
num_replications=1
source("./source_code/R/time_track_function.R")
exp_str <- paste("Track time for \nNum Subjects:\t", num_indvs,
                 "\n timeserires_length:\t",timeseries_length,
                 "\n fl_choice:\t",fl_choice
)
writeLines(exp_str)
timeKeeperStart(exp_str)
time_interval=seq(start_time,end_time,length.out=timeseries_length)
n100_rep5=get_T_simulations(klen, mu1_coef,mu2_coef,num_indvs, timeseries_length,
                                     time_interval=time_interval, fl_choice,num_replications, boot_number,lp_intercept=0.9998364)
# mean(n100_rep5[,2])
# hist(n100_rep5[,1])
# quantile(n100_rep5[,1],0.95)
timeKeeperNext()
#fl_choice="25"


get_T_distribution=function(klen, mu1_coef,mu2_coef,num_indvs, timeseries_length,
                        time_interval, fl_choice,num_replications,lp_intercept=0.9998364){
    T_rep <- foreach(this_row = 1:num_replications ) %dorng%
        { source("./source_code/R/data_generator.R")
            source("./source_code/R/integral_penalty_function.R")
            source("./source_code/R/T_testing_functions.R")
            #T_rv_erv <- list()
            number_basis=30
            WY_sample=GenerateCategoricalFDTest(klen, mu1_coef,mu2_coef,num_indvs, timeseries_length,
                                                time_interval, fl_choice, lp_intercept=0.9998364)
            
            temp=get_T(WY_sample$true$TrueX1,WY_sample$true$TrueX2,WY_sample$true$TrueX3, 
                       WY_sample$true$yis,time_interval,
                       number_basis =number_basis,est_choice="binomial" )
            
            T_stat=temp$T_statistics #scalar
            return(T_stat)
        }
    T_rep <- do.call(rbind, T_rep)
    #three columns, T, and T_binary, T_binary0.1
    return(T_rep)
    }

source("./source_code/R/time_track_function.R")
exp_str <- paste("Track time for \nNum Subjects:\t", num_indvs,
                 "\n timeserires_length:\t",timeseries_length,
                 "\n fl_choice:\t",fl_choice
)
writeLines(exp_str)
timeKeeperStart(exp_str)
num_replications=1000
time_interval=seq(start_time,end_time,length.out=timeseries_length)
n500_rep_justT=get_T_distribution(klen, mu1_coef,mu2_coef,num_indvs, timeseries_length,
                            time_interval=time_interval, fl_choice,num_replications,
                            lp_intercept=0.9998364)

timeKeeperNext()

# --------------------
#     Track time for 
# Num Subjects:	 500 
# timeserires_length:	 90 
# fl_choice:	 6 
# took: Time difference of 7.039114 mins 
# ====================

# hist_T_simulation=function(fl_choice,num_replications=1000){
#     par(mfrow=c(1,3))
#     num_indvs=100
#     source("./source_code/R/time_track_function.R")
#     exp_str <- paste("Track time for \nNum Subjects:\t", num_indvs,
#                      "\n timeserires_length:\t",timeseries_length,
#                      "\n fl_choice:\t",fl_choice
#                     )
#     writeLines(exp_str)
#     timeKeeperStart(exp_str)
#     null_100=get_T_simulations(klen, mu1_coef,mu2_coef,num_indvs, timeseries_length,
#                             time_interval, fl_choice, num_replications,lp_intercept=0.9998364)
#     timeKeeperNext()
#     num_indvs=500
#     exp_str <- paste("Track time for \nNum Subjects:\t", num_indvs,
#                      "\n timeserires_length:\t",timeseries_length,
#                      "\n fl_choice:\t",fl_choice)
#     writeLines(exp_str)
#     timeKeeperStart(exp_str)
#     null_500=get_T_simulations(klen, mu1_coef,mu2_coef,num_indvs, timeseries_length,
#                             time_interval, fl_choice,num_replications, lp_intercept=0.9998364)
#     timeKeeperNext()
#     num_indvs=1000
#     exp_str <- paste("Track time for \nNum Subjects:\t", num_indvs,
#                      "\n timeserires_length:\t",timeseries_length,
#                      "\n fl_choice:\t",fl_choice)
#     null_1000=get_T_simulations(klen, mu1_coef,mu2_coef,num_indvs, timeseries_length,
#                              time_interval, fl_choice,num_replications, lp_intercept=0.9998364)
#     timeKeeperNext()
#     par(mfrow=c(1,3))
#     num_indvs=100
#     hist(null_100,xlab="T",main=paste0("T for", num_indvs," Subjects", " with", num_replications," simulations"))
#     num_indvs=500
#     hist(null_500,xlab="T",main=paste0("T for", num_indvs," Subjects", " with", num_replications," simulations"))
#     num_indvs=1000
#     hist(null_1000,xlab="T",main=paste0("T for", num_indvs," Subjects", " with", num_replications," simulations"))
# }
# 
# 
# par(mfrow=c(1,2))
# num_indvs=100
# hist(null_100,xlab="T",main=paste0("T for", num_indvs," Subjects", " with", num_replications," simulations"))
# num_indvs=500
# hist(null_500,xlab="T",main=paste0("T for", num_indvs," Subjects", " with", num_replications," simulations"))
# 
# 
# hist_T_simulation("6",num_replications=1000)
# hist_T_simulation("25",num_replications=1000)

#
#
# WY_sample=GenerateCategoricalFDTest(klen, mu1_coef,mu2_coef,num_indvs, timeseries_length,
#                                       time_interval, fl_choice, lp_intercept=0.9998364)
# #
# # num_indvs=1000
# T_example=get_T(WY_sample$true$Truecatcurve, WY_sample$true$yis,time_interval,
#                 number_basis =30,est_choice="binomial" )
# 
# par(mfrow=c(1,2))
#      hist(T_example$rv_XF,xlab="RV",main=paste0("Historgram of rv for", num_indvs," Subjects"))
#      hist(T_example$rv_E_PF,xlab="E of RV",main=paste0("Historgram of E rv for", num_indvs," Subjects"))
# two_sample_result <- t.test(T_example$rv_XF,T_example$rv_E_PF, var.equal = FALSE)
#     print(two_sample_result)
# 
# 
# library(MASS)
# fit_data=c(T_example$T_vector)
# fit_result <- fitdistr(fit_data, "log-normal")
# estimated_df <- fit_result$estimate
# print(estimated_df)

#num_indvs=500

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
#     # fit_datap=c(T_example$T_vectorp)
#     # #fit_result <- fitdistr(fit_data, "log-normal")
#     # fit_resultp <- fitdistr(fit_datap, "chi-squared", list(df=3), lower = 0.001)
#     # estimated_dfp <- fit_resultp$estimate
#     # print(estimated_dfp)
#     ##
#     two_sample_result <- t.test(T_example$rv_XF,T_example$rv_E_PF, var.equal = FALSE)
#     print(two_sample_result)
# 
#     #par(mfrow=c(1,2))
#     hist(T_example$T_vector,xlab="T",main=paste0("Historgram of T for", num_indvs," Subjects"))
#     #hist(T_example$T_vectorp,xlab="T",main=paste0("Historgram of Tp for", num_indvs," Subjects"))
#     
#     return(list("T_example$rv_XF"=T_example$rv_XF,"T_example$rv_E_PF"=T_example$rv_E_PF))
# }
# 
# fl_choice="6"
# #fl_choice="25"
# par(mfrow=c(1,3))
# num_indvs=100
# null_100=get_T_for_hist(klen, mu1_coef,mu2_coef,num_indvs, timeseries_length,
#                         time_interval, fl_choice, lp_intercept=0.9998364)
# num_indvs=500
# null_500=get_T_for_hist(klen, mu1_coef,mu2_coef,num_indvs, timeseries_length,
#                time_interval, fl_choice, lp_intercept=0.9998364)
# num_indvs=1000
# null_1000=get_T_for_hist(klen, mu1_coef,mu2_coef,num_indvs, timeseries_length,
#                time_interval, fl_choice, lp_intercept=0.9998364)
# 
# 
# 
# fl_choice="25"
# par(mfrow=c(1,3))
# num_indvs=100
# alt_100=get_T_for_hist(klen, mu1_coef,mu2_coef,num_indvs, timeseries_length,
#                         time_interval, fl_choice, lp_intercept=0.9998364)
# num_indvs=500
# alt_500=get_T_for_hist(klen, mu1_coef,mu2_coef,num_indvs, timeseries_length,
#                         time_interval, fl_choice, lp_intercept=0.9998364)
# num_indvs=1000
# alt_1000=get_T_for_hist(klen, mu1_coef,mu2_coef,num_indvs, timeseries_length,
#                          time_interval, fl_choice, lp_intercept=0.9998364)
# 
# two_sample_test_function=function(null_data,alt_data){
#     two_sample_result <- t.test(null_data,alt_data, var.equal = FALSE)
#     print(two_sample_result)
# }
# 
# two_sample_test_function(null_100$`T_example$rv_XF`,alt_100$`T_example$rv_XF`)
# two_sample_test_function(null_500$`T_example$rv_XF`,alt_500$`T_example$rv_XF`)
# two_sample_test_function(null_1000$`T_example$rv_XF`,alt_1000$`T_example$rv_XF`)
