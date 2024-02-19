#sol = ? #

source("source_code/R/data_generator.R")

timeseries_length = 2000
# timestamps01 <- seq(from = 0.01, to = 0.99, length=timeseries_length)

test_results <- NULL
for (test_replica in 1:100) {
  test_results[test_replica] <- cfd_testing(start_time=0.01, end_time=0.99, timeseries_length=timeseries_length,
              num_indvs=100, 
              mu1_coef=c(sol_a[1],sol_a[2],sol_a[3]),
              mu2_coef=c(sol_a[4],sol_a[5],sol_a[6],sol_a[7],sol_a[8],sol_a[9],sol_a[10],sol_a[11],sol_a[12]), 
              fl_choice=3,
              response_family='bernoulli', test_type='Functional',
              klen=3)
}

# sol <- ga@solution
# 
# 
# cfd_test_data <- GenerateCategoricalFDTest(klen=3,
#                                            mu1_coef=c(sol[1],sol[2],sol[3]),
#                                            mu2_coef=c(sol[4],sol[5],sol[6],sol[7],sol[8],sol[9],sol[10],sol[11],sol[12]),
#                                            num_indvs=100,
#                                            timeseries_length = timeseries_length,
#                                            time_interval = timestamps01,
#                                            fl_choice=3)
# 
# if(length(table(cfd_test_data$true$Truecatcurve)) < 3){
#   print("MYERROR: Less Than 3 categories!!!")
#   
# }else{
# 
#   linear_predictor=cfd_test_data$true$linear_predictor
#   Y_indvs=cfd_test_data$true$yis
#   if(length(table(Y_indvs)) < 2){
#     print("MYERROR: Only contain eitehr 0s or 1s!!!")
#     
#   } else {
#     
#     prob_ind=cfd_test_data$true$prob_ind
#     data_sample=data.frame(linear_predictor,as.factor(Y_indvs),prob_ind)
#     colnames(data_sample)=c("linear_predictor","Group","Probability")
#     mean_0 <- mean(data_sample[data_sample$Group=="0","linear_predictor"])
#     mean_1 <- mean(data_sample[data_sample$Group=="1","linear_predictor"])
#     
#     abs_diff = abs(mean_0 - mean_1) # this should be high
#     print(abs_diff)
#     
#     alloc_diff = abs(table(Y_indvs)[[1]] - table(Y_indvs)[[2]]) # this should be low
#     print(table(Y_indvs))
#     print(alloc_diff)
#     
#     penalty = abs_diff - alloc_diff
#     print(penalty)
#   }
# }
