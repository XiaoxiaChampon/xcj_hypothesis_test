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
source("./source_code/R/data_generator.R")

start_time=0.01
end_time=0.99
timeseries_length=90
time_interval=seq(start_time,end_time,length=timeseries_length)
num_indvs=500
#fl_choice=3 #not constant, expect to reject
#fl_choice=8 #not constant, expect to reject

fl_choice=6
klen=3
mu1_coef=c(-1.8270644 ,-2.4700275,  5.4299181)
mu2_coef=c(-2.9990822, -0.8243365,  3.9100000  )


WY_sample=GenerateCategoricalFDTest(klen, mu1_coef,mu2_coef,num_indvs, timeseries_length,
                                      time_interval, fl_choice, lp_intercept=0.9998364)

T_example=get_T(WY_sample$true$Truecatcurve, WY_sample$true$yis,time_interval, 
                number_basis =30,est_choice="binomial" )


library(MASS)
fit_data=c(T_example$T_vector)
fit_result <- fitdistr(fit_data, "log-normal")
estimated_df <- fit_result$estimate
print(estimated_df)

num_indvs=500

get_T_for_hist=function(klen, mu1_coef,mu2_coef,num_indvs, timeseries_length,
                        time_interval, fl_choice, lp_intercept=0.9998364){
    WY_sample=GenerateCategoricalFDTest(klen, mu1_coef,mu2_coef,num_indvs, timeseries_length,
                                        time_interval, fl_choice, lp_intercept=0.9998364)
    
    T_example=get_T(WY_sample$true$Truecatcurve, WY_sample$true$yis,time_interval, 
                    number_basis =30,est_choice="binomial" )
    fit_data=c(T_example$T_vector)
    fit_result <- fitdistr(fit_data, "log-normal")
    estimated_df <- fit_result$estimate
    print(estimated_df)
    hist(T_example$T_vector,xlab="T",main=paste0("Historgram of T for", num_indvs," Subjects"))
}
num_indvs=100
get_T_for_hist(klen, mu1_coef,mu2_coef,num_indvs, timeseries_length,
                        time_interval, fl_choice, lp_intercept=0.9998364)
num_indvs=500
get_T_for_hist(klen, mu1_coef,mu2_coef,num_indvs, timeseries_length,
               time_interval, fl_choice, lp_intercept=0.9998364)
num_indvs=1000
get_T_for_hist(klen, mu1_coef,mu2_coef,num_indvs, timeseries_length,
               time_interval, fl_choice, lp_intercept=0.9998364)
# @param W 2D array, t*n: t is the timestamp and n is the number of the observation
# @return X 3D array, n*t*Q, Q: the total number of the category
# GetXFromW <- function(W)

#for getxfromW function
source("./source_code/R/cfd_hypothesis_test.R")

get_T <- function(W, Y,time_interval, number_basis =30,est_choice ){
    # W=WY_sample$true$Truecatcurve
    # Y=WY_sample$true$yis
    #est_choice="binomial"
    #Estimation
    
    categFD_est <- EstimateCategFuncDataX(est_choice, time_interval, W)
    X_array=categFD_est$X_array
    num_indv <- nrow(X_array[,,1])
    category_count <- dim(X_array)[3]
    timeseries_length<- length(time_interval)
    # @return list of 2D array: True Z curves , Est Z curves, True p curves, Est p curves
    #                          all have dimension t*n
    
    #p<-array(0, c( timeseries_length,num_indv, category_count))
    
    pl_matrix=array(c(categFD_est$pl$p1_est,categFD_est$pl$p2_est,categFD_est$pl$p3_est),dim=c(timeseries_length,num_indvs,category_count))
    number_col <- number_basis
    knots <- construct.knots(time_interval,knots=(number_basis-3),knots.option='equally-spaced')
    bspline <- splineDesign(knots=knots,x=time_interval,ord=4)
    mub_matrix <- foreach(this_row = 1:num_indv ) %do%
        {
            source("./source_code/R/integral_penalty_function.R")
            
            temp <- array(-123, number_col)
            for(this_col in 1:number_col){
                temp[this_col] <- integral_penalty(time_interval,integral_function(time_interval,pl_matrix[,this_row,category_count-1]*bspline[,this_col]))$value
            }
            return(temp)
        }
    mub_matrix <- do.call(rbind, mub_matrix)
    
    
    time_interval_matrix=do.call("rbind", replicate(length(Y), time_interval, simplify = FALSE)) 
    
    ##parallel use each row of X to multiple each column of bspline
    
    # J_matrix <- matrix(0,nrow=nrow(X_matrix),ncol=number_basis) #empty
    # for(row in 1:number_row){
    #   for(col in 1:number_col){
    #     J_matrix[row,col] <- integral_penalty(time_interval,X_matrix[row,]*bspline[,col])$value
    #   }
    # }
    
    
    # J_matrix <- foreach(this_row = 1:number_row) %do%
    #     {
    #         source("./source_code/R/integral_penalty_function.R")
    #         
    #         temp <- array(-123, number_col)
    #         for(this_col in 1:number_col){
    #             temp[this_col] <- integral_penalty(time_interval,X_matrix[this_row,]*bspline[,this_col])$value
    #         }
    #         return(temp)
    #     }
    # J_matrix <- do.call(rbind, J_matrix)
    
    
    DBB_matrix <- foreach(this_row = 1:number_col) %do%
        {
            source("./source_code/R/integral_penalty_function.R")
            
            temp <- array(-123, number_col)
            for(this_col in 1:number_col){
                temp[this_col] <- (integral_penalty(time_interval,bspline[,this_row]*bspline[,this_col])$value)*(cov(X_array[,,2])[this_row,this_col])
            }
            return(temp)
        }
    DBB_matrix <- do.call(rbind, DBB_matrix)
    
    
    logit_model=gam(Y~s(time_interval_matrix,by=X_array[,,2],k = number_basis)+s(time_interval_matrix,by=X_array[,,3],k = number_basis),family = 'binomial')
    betal=logit_model$coefficients[2:(number_basis+1)]
    
    
    
    ## 
    T_vector <- foreach(this_row = 1:num_indv ) %do%
    {
        source("./source_code/R/integral_penalty_function.R")

        
            temp <-t(betal)%*% (mub_matrix[this_row,]%*%t(mub_matrix[this_row,])+DBB_matrix)%*%(betal)
        
        return(temp)
    }
    T_vector <- do.call(cbind, T_vector)
    
    
    #R*R, R*1, n*R
    return(list("DBB_matrix"=DBB_matrix, "betal"=betal,"mub_matrix"=mub_matrix,"T_vector"=T_vector))}





#' Function to select 
#' @param choice "probit", "binomial",  or "multinormial"
#' @param timestamps01, 1D array, time interval that cfd is observed
#' @param  W: 2D array, t*n, t: the number of time points, n: the number of individuals
#' @param  basis_size=25, the number of basis function used 
#' @param  method="ML"
#' @return list of 2D array: True Z curves , Est Z curves, True p curves, Est p curves
#'                           all have dimension t*n
EstimateCategFuncDataX <- function(choice, timestamps01, W, basis_size=25, method="ML")
{
    if(choice == "probit"){
        X <- GetXFromW(W)
        pl=EstimateCategFuncData_probit(timestamps01, X, basis_size, method, 1/150)
        return(list("pl"=pl,"X_array"=X))
    }else if(choice == "binomial"){
        X <- GetXFromW(W)
        #timestamps01=time_interval
        # basis_size=25
        # method="ML"
        pl=EstimateCategFuncData_binorm(timestamps01, X, basis_size, method)
        return(list("pl"=pl,"X_array"=X))
    }
}




#'Function to estimate z and p using probit
#' @param timestamps01, 1D array, time interval that cfd is observed
#' @param  X: 3D array, t*n*Q, t: the number of time points, n: the number of individuals, Q: the number of categories
#' @param  basis_size=25, the number of basis function used 
#' @param  method="ML"
#' @return list of 2D array: True Z curves , Est Z curves, True p curves, Est p curves
#'                           all have dimension t*n
EstimateCategFuncData_probit <- function(timestamps01, X, basis_size=25, method="ML", threshold_probability=0.004)
{
    num_indv<- dim(X)[1]
    timeseries_length<- dim(X)[2]
    category_count <- dim(X)[3]
    
    Z<-NULL
    p<-array(0, c(num_indv, timeseries_length, category_count))
    # for (indv in 1:num_indv){
    #   #i=93
    #   #basis_size=25
    #   x1<- X[indv,,1]
    #   x2<- X[indv,,2]
    #   x3<- X[indv,,3]
    # 
    #   if (timeseries_length<=301 && sum(x1)/timeseries_length<threshold_probability){
    # 
    #     gam_result_1 <- RunGam(timestamps01, x1, "probit", basis_size, method)
    #     p1 <- gam_result_1$prob
    #     p1_linpred <- gam_result_1$linpred
    # 
    #     gam_result_2 <- RunGam(timestamps01, x2, "probit", basis_size, method)
    #     p2 <- gam_result_2$prob
    #     p2_linpred <- gam_result_2$linpred
    # 
    #     gam_result_3 <- RunGam(timestamps01, x3, "probit", basis_size, method)
    #     p3 <- gam_result_3$prob
    #     p3_linpred <- gam_result_3$linpred
    #     denominator_p <- 1+exp(p3_linpred)
    #     z1<- (p1_linpred-p3_linpred)-log( (1+exp(p1_linpred))/(denominator_p))
    #     z2<- (p2_linpred-p3_linpred)-log( (1+exp(p2_linpred))/(denominator_p))
    #   }else{
    #     gam_result_1 <- RunGam(timestamps01, x1, "binomial", basis_size, method)
    #     p1 <- gam_result_1$prob
    #     p1_linpred <- gam_result_1$linpred
    # 
    #     gam_result_2 <- RunGam(timestamps01, x2, "binomial", basis_size, method)
    #     p2 <- gam_result_2$prob
    #     p2_linpred <- gam_result_2$linpred
    # 
    #     gam_result_3 <- RunGam(timestamps01, x3, "binomial", basis_size, method)
    #     p3 <- gam_result_3$prob
    #     p3_linpred <- gam_result_3$linpred
    # 
    #     # estimate the latent curves Z
    #     exp_p3_linepred <- exp(p3_linpred)
    #     z1 <- (p1_linpred-p3_linpred)-log( (1+exp(p1_linpred))/(1+exp_p3_linepred))
    #     z2 <- (p2_linpred-p3_linpred)-log( (1+exp(p2_linpred))/(1+exp_p3_linepred))
    # 
    # 
    #   } # end if special case for probit
    # 
    #   Z <- cbind(Z, c(z1,z2))
    #   psum <- (p1+p2+p3)
    #   p[indv,,] <- cbind(p1/psum, p2/psum, p3/psum)
    # }
    # return(list(Z1_est=Z[1:timeseries_length,],
    #             Z2_est=Z[1:timeseries_length+timeseries_length,],
    #             p1_est=t(p[,,1]),
    #             p2_est=t(p[,,2]),
    #             p3_est=t(p[,,3]) ))
    
    zp <- foreach (indv = 1:num_indv, .combine = cbind, .init = NULL, .packages = c("mgcv")) %dorng%
        {
            source("./source_code/R/run_gam_function.R")
            
            x1<- X[indv,,1]
            x2<- X[indv,,2]
            x3<- X[indv,,3]
            
            probit_binom <- function(x_binary){
                if (sum(x_binary)/timeseries_length < threshold_probability){
                    gam_result_binary <- RunGam(timestamps01, x_binary, "probit", basis_size, method)
                    p_binary <- gam_result_binary$prob
                    p_binary_linpred <- gam_result_binary$linpred
                }else{
                    gam_result_binary <- RunGam(timestamps01, x_binary, "binomial", basis_size, method)
                    p_binary <- gam_result_binary$prob
                    p_binary_linpred <- gam_result_binary$linpred
                }
                return(list("p_binary"=p_binary,"p_binary_linpred"=p_binary_linpred))
            }
            
            r_1 <- probit_binom(x1)
            p1 <- r_1$p_binary
            p1_linpred <- r_1$p_binary_linpred
            
            r_2 <- probit_binom(x2)
            p2 <- r_2$p_binary
            p2_linpred <- r_2$p_binary_linpred
            
            r_3 <- probit_binom(x3)
            p3 <- r_3$p_binary
            p3_linpred <- r_3$p_binary_linpred
            
            # estimate the latent curves Z
            denominator_p <- 1 + exp(p3_linpred)
            z1 <- (p1_linpred-p3_linpred)-log( (1+exp(p1_linpred))/(denominator_p))
            z2 <- (p2_linpred-p3_linpred)-log( (1+exp(p2_linpred))/(denominator_p))
            
            psum <- p1 + p2 + p3
            return(c(c(z1,z2), cbind(p1/psum, p2/psum, p3/psum)))
        }
    # Unravel the two variables from zp
    z_rows_count <- timeseries_length * 2
    Z <- array(zp[1:z_rows_count, ], c(z_rows_count, num_indv))
    p <- array(t(matrix(zp[(z_rows_count + 1):dim(zp)[1], ], ncol=num_indv)), c(num_indv, timeseries_length, category_count))
    
    
    return(list(Z1_est=Z[1:timeseries_length,],
                Z2_est=Z[1:timeseries_length+timeseries_length,],
                p1_est=t(p[,,1]),
                p2_est=t(p[,,2]),
                p3_est=t(p[,,3]) ))
}

#'Function to estimate z and p using binom
#' @param timestamps01, 1D array, time interval that cfd is observed
#' @param  X: 3D array, t*n*Q, t: the number of time points, n: the number of individuals, Q: the number of categories
#' @param  basis_size=25, the number of basis function used 
#' @param  method="ML"
#' @return list of 2D array: True Z curves , Est Z curves, True p curves, Est p curves
#'                           all have dimension t*n

EstimateCategFuncData_binorm <- function(timestamps01, X, basis_size=25, method="ML")
{
    num_indv<- dim(X)[1]
    timeseries_length<- dim(X)[2]
    category_count <- dim(X)[3]
    
    Z<-NULL
    p<-array(0, c(num_indv, timeseries_length, category_count))
    ##########################
    ###########################
    # num_indv is the subject and this step is done by subject level
    # can parallel
    #############################
    #############################
    #   for (indv in 1:num_indv)
    #   {
    #     x1<- X[indv,,1]
    #     x2<- X[indv,,2]
    #     x3<- X[indv,,3]
    #
    #     # fit the Binom model
    #     ###################################################
    #     ###################################################
    #     ##updated estimation
    #
    #     gam_result_1 <- RunGam(timestamps01, x1, "binomial", basis_size, method)
    #     p1 <- gam_result_1$prob
    #     p1_linpred <- gam_result_1$linpred
    #
    #     gam_result_2 <- RunGam(timestamps01, x2, "binomial", basis_size, method)
    #     p2 <- gam_result_2$prob
    #     p2_linpred <- gam_result_2$linpred
    #
    #     gam_result_3 <- RunGam(timestamps01, x3, "binomial", basis_size, method)
    #     p3 <- gam_result_3$prob
    #     p3_linpred <- gam_result_3$linpred
    #
    #     # estimate the latent tranjecotries Z
    #     exp_p3_linepred <- exp(p3_linpred)
    #     z1 <- (p1_linpred-p3_linpred)-log( (1+exp(p1_linpred))/(1+exp_p3_linepred))
    #     z2 <- (p2_linpred-p3_linpred)-log( (1+exp(p2_linpred))/(1+exp_p3_linepred))
    #
    #     Z <- cbind(Z, c(z1,z2))
    #     psum <- (p1+p2+p3)
    #     p[indv,,] <- cbind(p1/psum, p2/psum, p3/psum)
    #   }
    # return(p)
    # return(list(Z=Z,p=p))
    
    
    zp <- foreach (indv = 1:num_indv, .combine = cbind, .init = NULL, .packages = c("mgcv")) %dorng%
        {
            source("./source_code/R/run_gam_function.R")
            
            x1<- X[indv,,1]
            x2<- X[indv,,2]
            x3<- X[indv,,3]
            
            gam_result_1 <- RunGam(timestamps01, x1, "binomial", basis_size, method)
            p1 <- gam_result_1$prob
            p1_linpred <- gam_result_1$linpred
            
            gam_result_2 <- RunGam(timestamps01, x2, "binomial", basis_size, method)
            p2 <- gam_result_2$prob
            p2_linpred <- gam_result_2$linpred
            
            gam_result_3 <- RunGam(timestamps01, x3, "binomial", basis_size, method)
            p3 <- gam_result_3$prob
            p3_linpred <- gam_result_3$linpred
            
            # estimate the latent tranjecotries Z
            denominator_p <- 1 + exp(p3_linpred)
            z1 <- (p1_linpred-p3_linpred)-log( (1+exp(p1_linpred))/(denominator_p))
            z2 <- (p2_linpred-p3_linpred)-log( (1+exp(p2_linpred))/(denominator_p))
            
            psum <- p1 + p2 + p3
            return(c(c(z1,z2), cbind(p1/psum, p2/psum, p3/psum)))
        }
    # Unravel the two variables from zp
    z_rows_count <- timeseries_length * 2
    Z <- array(zp[1:z_rows_count, ], c(z_rows_count, num_indv))
    p <- array(t(matrix(zp[(z_rows_count + 1):dim(zp)[1], ], ncol=num_indv)), c(num_indv, timeseries_length, category_count))
    
    return(list(Z1_est=Z[1:timeseries_length,],
                Z2_est=Z[1:timeseries_length+timeseries_length,],
                p1_est=t(p[,,1]),
                p2_est=t(p[,,2]),
                p3_est=t(p[,,3]) ))
}




if(run_parallel)
{
    parallel::stopCluster(cl = my.cluster)
    initialized_parallel <- FALSE
}
