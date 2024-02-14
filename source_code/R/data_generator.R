
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
# Purpose: Functions to generate Categorical Functional Data Hypothesis Testing 
#         
# Author:  Xiaoxia Champon
# Date: 10/26/2023
#
##############################################################
library(mgcv)
library(fda)
library(fda.usc)
library(devtools)
#install_github("stchen3/glmmVCtest")
library("glmmVCtest")
library(RLRsim)
library(MASS)
#for splineDesign
library(splines)

source("./source_code/R/cfd_hypothesis_test.R")


#' Get mu_1, mu_2 functions, and score_vals objects for a given context.
#' @param klen number of points along the score decay axis
#' @param mu1_coef 1D array, length 3
#' @return A list that contains mu_1, mu_2, score_vals
#'
GetMuAndScore_2 <- function(klen,mu1_coef,mu2_coef)
{
    # mu_1 <- function(t){ -0.64+4*t }
    # 
    # mu_2 <- function(t){ 0.97+6*t^2 }
  
    mu_1 <- function(t){ mu1_coef[1] + mu1_coef[2] * t + mu1_coef[3] * t^2 }

    mu_2 <- function(t){ mu2_coef[1] + mu2_coef[2] * t + mu2_coef[3] * t^2 }

    all_score_values <- rep(0, klen)

    score_vector <- cbind(1,sqrt(1/2),1/2)

    for(idx in 1:length(score_vector))
    {
        all_score_values[idx] <- score_vector[idx]
    }

    return(list("mu_1" = mu_1, "mu_2" = mu_2, "score_vals" = all_score_values))
}

#' Generate cluster data for a given scenario
#' @param num_indvs number of individuals
#' @param timeseries_length length of time-series as an integer
#' @param k description?? number of psi functions
#' @param mu_1 description??
#' @param mu_2 description??
#' @param score_vals description??
#'


GenerateDataTest <- function(num_indvs,
                            timeseries_length,
                            mu_1,
                            mu_2,
                            score_vals,
                            start_time = 0.01,
                            end_time = 0.99,
                            k = 3)
{
    timestamps01 <- seq(from = start_time, to = end_time, length=timeseries_length)

    # noise octaves
    # cat("octave", num_indvs, k, num_indvs * k, "\n")
    scores_standard <- matrix(rnorm(num_indvs * k), ncol = k)
    scores <- scores_standard %*% diag(sqrt(score_vals))

    #
    BIG_mu <- c(mu_1(timestamps01), mu_2(timestamps01))
    BIG_phi <- PsiFunc(k, timestamps01)

    Z <- BIG_phi %*% t(scores) + BIG_mu
    expZ <- exp(Z)
    Z1 <- Z[1:timeseries_length, ]
    Z2 <- Z[1:timeseries_length + timeseries_length, ]
    expZ1 <- expZ[1:timeseries_length, ]
    expZ2 <- expZ[1:timeseries_length + timeseries_length, ]
    denom <- 1 + expZ1 + expZ2
    p1 <- expZ1 / denom
    p2 <- expZ2 / denom
    p3 <- 1 / denom

    # vectorize for future work!!!
    return(list(Z1 = Z1, Z2 = Z2,
                p1 = p1, p2 = p2, p3 = p3,
                MEAN = BIG_mu, PHI = BIG_phi, MFPC = scores))
}

#' Psi function
#'
PsiFunc <- function(klen, timestamps01)
{
    psi_k1 <- sapply(c(1:klen), function(i) sin((2 * i + 1) * pi * timestamps01))
    psi_k2 <- sapply(c(1:klen), function(i) cos(2 * i * pi * timestamps01))
    return(rbind(psi_k1, psi_k2))
}

GenerateCategFuncData <- function(prob_curves)
{
  curve_count <- length(prob_curves);
  
  num_indvs <- ncol(prob_curves$p1)
  timeseries_length <- nrow(prob_curves$p1)
  # cat("n:", num_indvs, "\tt:", timeseries_length, "\n")
  
  
  W <- matrix(0, ncol=num_indvs, nrow=timeseries_length)
  X_array <- array(0, c(num_indvs, timeseries_length, curve_count))
  
  for(indv in c(1:num_indvs))
  {
    X <- sapply(c(1:timeseries_length),
                function(this_time) rmultinom(n=1,
                                              size=1,
                                              prob = c(prob_curves$p1[this_time,indv],
                                                       prob_curves$p2[this_time,indv],
                                                       prob_curves$p3[this_time,indv]) ))
    W[,indv] <- apply(X, 2, which.max)
    X_array[indv,,] <- t(X)
  }
  
  return(list(X=X_array, W=W)) # X_binary W_catfd
}


GenerateCategFuncDataUpdate <- function(prob_curves,mu1_coef,mu2_coef)
{
    categ_func_data_list <- GenerateCategFuncData(prob_curves)
    W=categ_func_data_list$W
     num_indvs <- ncol(prob_curves$p1)
     timeseries_length <- nrow(prob_curves$p1)
     # cat("n:", num_indvs, "\tt:", timeseries_length, "\n")
    ##########add re generate W if one of the category is missing
    Q_vals <- unique(c(categ_func_data_list$W))
    if(is.numeric(Q_vals))
    {
      Q_vals <- sort(Q_vals)
    }
    #####################################
    
    
    for(indv in c(1:num_indvs))
    {##########add re generate W if one of the category is missing
        tolcat <- table(categ_func_data_list$W[,indv])
        catorder <- order(tolcat, decreasing = TRUE)
        numcat <- length(catorder)
        refcat <- catorder[numcat]
        count_iter <- 0
        while (count_iter < 100 && 
               ( (numcat < length(Q_vals))
                 ||(timeseries_length==300  && min(as.numeric(tolcat)) < 4)
                 ||(timeseries_length==750  && min(as.numeric(tolcat)) < 10)
               )
        )
        {
          count_iter <- count_iter + 1
          
          mns <- GetMuAndScore_2(klen=3,mu1_coef,mu2_coef)
          klen=3
          time_interval=seq(0.01,0.99,length=timeseries_length)
          generated_data <- GenerateDataTest(num_indvs = 5,
                                             timeseries_length = timeseries_length,
                                             mu_1 = mns$mu_1,
                                             mu_2 = mns$mu_2,
                                             score_vals = mns$score_vals,
                                             start_time = time_interval[1],
                                             end_time = tail(time_interval,1),
                                             k = klen)
          
          new_prob_curves <-  list(p1 = generated_data$p1, p2 = generated_data$p2, p3 = generated_data$p3)
          new_categ_func_data_list <- GenerateCategFuncData( new_prob_curves )
          
          categ_func_data_list$W[, indv] <- new_categ_func_data_list$W[, 3]
          #Z1[, indv] <- generated_data$Z1[, 3] # latent curves Z1 and Z2
          #create empty series
          categ_func_data_list$X[indv, , ] <- 0
          #Z2[, indv] <- generated_data$Z2[, 3]
          
          for (this_time in 1:timeseries_length)
          {
            categ_func_data_list$X[indv, this_time, which(Q_vals == categ_func_data_list$W[, indv][this_time])] <- 1
          }
          
          tolcat <- table(categ_func_data_list$W[, indv])
          catorder <- order(tolcat, decreasing = TRUE)
          numcat <- length(catorder)
          refcat <- catorder[numcat]
        } # end while
        #############################
        
    }
    #W: t*n, X: n*t*Q
    return(list(X=categ_func_data_list$X, W=categ_func_data_list$W)) # X_binary W_catfd
}



expit <- function(x){1/(1+exp(-x))}

fl3_tilda = function(t){
  return(-20*sin((2*pi/25)*(t-1))-6)
}

fl2_tilda = function(t){
  return(-10*sin((2*pi/25)*(t-1)))
}

fl2f = function(t) {
    return(t - 8/9)
}

fl3f = function(t) {
    return(0.5*t)
}

fl3fn = function(t) {
    return(-3*t^2 + 2*t - 0.9)
}

flx789 <- function(t, x){
  return(x[7] + x[8]*t + x[9]*t^2)
}

flx456 <- function(t, x){
  return(x[4] + x[5]*t + x[6]*t^2)
}

GenerateCategoricalFDTest <- function(klen, mu1_coef,mu2_coef,num_indvs, timeseries_length,
                                      time_interval, fl_choice, lp_intercept=0.9998364){
  
    mns <- GetMuAndScore_2(klen,mu1_coef,mu2_coef)

    generated_data <- GenerateDataTest(num_indvs = num_indvs,
                                      timeseries_length = timeseries_length,
                                      mu_1 = mns$mu_1,
                                      mu_2 = mns$mu_2,
                                      score_vals = mns$score_vals,
                                      start_time = time_interval[1],
                                      end_time = tail(time_interval,1),
                                      k = klen)
    prob_curves <- list(p1 = generated_data$p1, p2 = generated_data$p2, p3 = generated_data$p3)
    cat_data <- GenerateCategFuncDataUpdate(prob_curves,mu1_coef,mu2_coef)

    flfn <- switch(fl_choice,
                   "1"=list("fl1"=rep(0.45,timeseries_length),
                            "fl2"=rep(0.5,timeseries_length),
                            "fl3"=rep(-0.51,timeseries_length)),
                   
                   "2"=list("fl1"=rep(-0.1,timeseries_length),
                            "fl2"=matrix(-0.1*fl2f(time_interval),nrow=timeseries_length,ncol=1),
                            "fl3"=matrix(fl3f(time_interval),nrow=timeseries_length,ncol=1)),
                   #not constant
                   # fll2=-10*sin((2*pi/25)*(time_interval-1))
                   # fll3=-20*sin((2*pi/25)*(time_interval-1))-6
                   "3"=list("fl1"=rep(-0.2,timeseries_length),
                            # "fl2"=matrix(-0.15*fl2f(time_interval),nrow=timeseries_length,ncol=1),
                            # "fl3"=matrix(fl3f(time_interval),nrow=timeseries_length,ncol=1)),
                            "fl2"=matrix(-10*sin((2*pi/25)*(time_interval-1)),nrow=timeseries_length,ncol=1),
                            "fl3"=matrix(-20*sin((2*pi/25)*(time_interval-1))-6,nrow=timeseries_length,ncol=1)),

                   "100"=list("fl1"=rep(-0.2,timeseries_length),
                            # "fl3"=matrix(fl3f(time_interval,),nrow=timeseries_length,ncol=1)),
                            "fl2"=matrix(flx456(time_interval,mu2_coef),nrow=timeseries_length,ncol=1),
                            "fl3"=matrix(flx789(time_interval,mu2_coef),nrow=timeseries_length,ncol=1)),                  
                   "200"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
                            "fl2"=matrix(rep(2.5,timeseries_length),nrow=timeseries_length,ncol=1),
                            "fl3"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)),
                   
                   "4"=list("fl1"=rep(-0.2,timeseries_length),
                            "fl2"=matrix(1.23+1.56*time_interval+0.58*time_interval^2,nrow=timeseries_length,ncol=1),
                            "fl3"=matrix(-1.86-5.03*time_interval+3.68*time_interval^2,nrow=timeseries_length,ncol=1)),
                   "200"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
                            "fl2"=matrix(rep(2.5,timeseries_length),nrow=timeseries_length,ncol=1),
                            "fl3"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)),
                   
                   
                   "6"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
                             "fl2"=matrix(rep(0,timeseries_length),nrow=timeseries_length,ncol=1),
                             "fl3"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)),
                   "7"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
                             "fl2"=matrix(rep(5,timeseries_length),nrow=timeseries_length,ncol=1),
                             "fl3"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)),
                   "8"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
                             "fl2"=matrix(rep(10,timeseries_length),nrow=timeseries_length,ncol=1),
                             "fl3"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)),
                   "9"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
                             "fl2"=matrix(rep(15,timeseries_length),nrow=timeseries_length,ncol=1),
                             "fl3"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)),
                   "10"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
                             "fl2"=matrix(rep(20,timeseries_length),nrow=timeseries_length,ncol=1),
                             "fl3"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)),
                   
                   
                   "11"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
                            "fl2"=matrix(1+0*time_interval,nrow=timeseries_length,ncol=1),
                            "fl3"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)),
                   "12"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
                            "fl2"=matrix(1+20*time_interval,nrow=timeseries_length,ncol=1),
                            "fl3"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)),
                   "13"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
                            "fl2"=matrix(1+40*time_interval,nrow=timeseries_length,ncol=1),
                            "fl3"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)),
                   "14"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
                            "fl2"=matrix(1+60*time_interval,nrow=timeseries_length,ncol=1),
                            "fl3"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)),
                  "15"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
                            "fl2"=matrix(1+80*time_interval,nrow=timeseries_length,ncol=1),
                            "fl3"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)),
                  
                  "21"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
                            "fl2"=matrix(1+0*time_interval,nrow=timeseries_length,ncol=1),
                            "fl3"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)),
                  "22"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
                            "fl2"=matrix(1+10*time_interval,nrow=timeseries_length,ncol=1),
                            "fl3"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)),
                  "23"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
                            "fl2"=matrix(1+20*time_interval,nrow=timeseries_length,ncol=1),
                            "fl3"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)),
                  "24"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
                            "fl2"=matrix(1+30*time_interval,nrow=timeseries_length,ncol=1),
                            "fl3"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)),
                  "25"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
                            "fl2"=matrix(1+40*time_interval,nrow=timeseries_length,ncol=1),
                            "fl3"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)),
                  "26"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
                            "fl2"=matrix(1+5*time_interval,nrow=timeseries_length,ncol=1),
                            "fl3"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1))
                   )

    vec <- matrix(1:num_indvs, nrow=num_indvs, ncol=1)
   
    x1fl1 <- apply(vec, 1, function(x) {fda.usc::int.simpson2(time_interval, flfn$fl1, equi = TRUE, method = "TRAPZ")})
    x2fl2 <- apply(vec, 1, function(x) {fda.usc::int.simpson2(time_interval, cat_data$X[x,,2]*(flfn$fl2), equi = TRUE, method = "TRAPZ")})
    x3fl3 <- apply(vec, 1, function(x) {fda.usc::int.simpson2(time_interval, cat_data$X[x,,3]*(flfn$fl3), equi = TRUE, method = "TRAPZ")})
    #lp_intercept= 0.6206897
    linear_predictor <- matrix(x1fl1 + x2fl2+ x3fl3 + lp_intercept )
    linear_predictor_without <- matrix(x1fl1 + x3fl3+ lp_intercept )
   
    generate_y_indvs <- function(lp){
      ys <- apply(lp, 1, function(x){ rbinom(1, 1, 1/(1+exp(-x))) })
      leastOccr <- min(sum(ys), length(ys) - sum(ys))
      count_iter <- 1
      min_occurrence <- round(num_indvs * 0.2)
      max_iterations <- 500
      while (count_iter < max_iterations && (length(ys) - sum(ys) < min_occurrence || sum(ys) < min_occurrence)){
        count_iter <- count_iter + 1
        candidate <- apply(lp, 1, function(x){ rbinom(1, 1, 1/(1+exp(-x))) })
        num1s <- sum(candidate)
        num0s <- length(ys) - num1s
        if(leastOccr < min(num0s, num1s)){
          ys <- candidate
          leastOccr <- min(num0s, num1s)
        }
      }
      # cat("Y generation count:", count_iter, "\n")
      return(ys)
    }
    Y_indvs <- generate_y_indvs(linear_predictor)
    Y_indvs_without <- generate_y_indvs(linear_predictor_without)
    
    prob_ind=1/(1+exp(-linear_predictor))
    
    truelist=list("TrueX1"=cat_data$X[,,1],
                  "TrueX2"=cat_data$X[,,2],
                  "TrueX3"=cat_data$X[,,3],
                  "Truecatcurve"=cat_data$W,
                  "fl"=flfn,
                  "yis"=Y_indvs,
                  "x2fl2"=x2fl2,
                  "yis_without" = Y_indvs_without,
                  "linear_predictor"=list("linearw"=linear_predictor,"linearwo"=linear_predictor_without),
                  "prob_ind"=prob_ind)
   
    return(list("true"=truelist))
}

GenerateCategoricalFDTestIntercept <- function(klen, mu1_coef,mu2_coef,num_indvs, timeseries_length,
                                      time_interval, fl_choice, intercept_opt){
  
  mns <- GetMuAndScore_2(klen,mu1_coef,mu2_coef)
  
  generated_data <- GenerateDataTest(num_indvs = num_indvs,
                                     timeseries_length = timeseries_length,
                                     mu_1 = mns$mu_1,
                                     mu_2 = mns$mu_2,
                                     score_vals = mns$score_vals,
                                     start_time = time_interval[1],
                                     end_time = tail(time_interval,1),
                                     k = klen)
  prob_curves <- list(p1 = generated_data$p1, p2 = generated_data$p2, p3 = generated_data$p3)
  cat_data <- GenerateCategFuncDataUpdate(prob_curves,mu1_coef,mu2_coef)
  
  flfn <- switch(fl_choice,
                 "1"=list("fl1"=rep(0.45,timeseries_length),
                          "fl2"=rep(0.5,timeseries_length),
                          "fl3"=rep(-0.51,timeseries_length)),
                 
                 "2"=list("fl1"=rep(-0.1,timeseries_length),
                          "fl2"=matrix(-0.1*fl2f(time_interval),nrow=timeseries_length,ncol=1),
                          "fl3"=matrix(fl3f(time_interval),nrow=timeseries_length,ncol=1)),
                 #not constant
                 # fll2=-10*sin((2*pi/25)*(time_interval-1))
                 # fll3=-20*sin((2*pi/25)*(time_interval-1))-6
                 "3"=list("fl1"=rep(-0.2,timeseries_length),
                          # "fl2"=matrix(-0.15*fl2f(time_interval),nrow=timeseries_length,ncol=1),
                          # "fl3"=matrix(fl3f(time_interval),nrow=timeseries_length,ncol=1)),
                          "fl2"=matrix(-10*sin((2*pi/25)*(time_interval-1)),nrow=timeseries_length,ncol=1),
                          "fl3"=matrix(-20*sin((2*pi/25)*(time_interval-1))-6,nrow=timeseries_length,ncol=1)),
                 
                 "100"=list("fl1"=rep(-0.2,timeseries_length),
                            # "fl3"=matrix(fl3f(time_interval,),nrow=timeseries_length,ncol=1)),
                            "fl2"=matrix(flx456(time_interval,mu2_coef),nrow=timeseries_length,ncol=1),
                            "fl3"=matrix(flx789(time_interval,mu2_coef),nrow=timeseries_length,ncol=1)),                  
                 "200"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
                            "fl2"=matrix(rep(2.5,timeseries_length),nrow=timeseries_length,ncol=1),
                            "fl3"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)),
                 
                 "4"=list("fl1"=rep(-0.2,timeseries_length),
                          "fl2"=matrix(1.23+1.56*time_interval+0.58*time_interval^2,nrow=timeseries_length,ncol=1),
                          "fl3"=matrix(-1.86-5.03*time_interval+3.68*time_interval^2,nrow=timeseries_length,ncol=1)),
                 "200"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
                            "fl2"=matrix(rep(2.5,timeseries_length),nrow=timeseries_length,ncol=1),
                            "fl3"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)),
                 
                 
                 "6"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
                          "fl2"=matrix(rep(0,timeseries_length),nrow=timeseries_length,ncol=1),
                          "fl3"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)),
                 "7"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
                          "fl2"=matrix(rep(5,timeseries_length),nrow=timeseries_length,ncol=1),
                          "fl3"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)),
                 "8"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
                          "fl2"=matrix(rep(10,timeseries_length),nrow=timeseries_length,ncol=1),
                          "fl3"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)),
                 "9"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
                          "fl2"=matrix(rep(15,timeseries_length),nrow=timeseries_length,ncol=1),
                          "fl3"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)),
                 "10"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
                           "fl2"=matrix(rep(20,timeseries_length),nrow=timeseries_length,ncol=1),
                           "fl3"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)),
                 
                 
                 "11"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
                           "fl2"=matrix(1+0*time_interval,nrow=timeseries_length,ncol=1),
                           "fl3"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)),
                 "12"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
                           "fl2"=matrix(1+20*time_interval,nrow=timeseries_length,ncol=1),
                           "fl3"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)),
                 "13"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
                           "fl2"=matrix(1+40*time_interval,nrow=timeseries_length,ncol=1),
                           "fl3"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)),
                 "14"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
                           "fl2"=matrix(1+60*time_interval,nrow=timeseries_length,ncol=1),
                           "fl3"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)),
                 "15"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
                           "fl2"=matrix(1+80*time_interval,nrow=timeseries_length,ncol=1),
                           "fl3"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)),
                 
                 "21"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
                           "fl2"=matrix(1+0*time_interval,nrow=timeseries_length,ncol=1),
                           "fl3"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)),
                 "22"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
                           "fl2"=matrix(1+10*time_interval,nrow=timeseries_length,ncol=1),
                           "fl3"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)),
                 "23"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
                           "fl2"=matrix(1+20*time_interval,nrow=timeseries_length,ncol=1),
                           "fl3"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)),
                 "24"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
                           "fl2"=matrix(1+30*time_interval,nrow=timeseries_length,ncol=1),
                           "fl3"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)),
                 "25"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
                           "fl2"=matrix(1+40*time_interval,nrow=timeseries_length,ncol=1),
                           "fl3"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)),
                 "26"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
                           "fl2"=matrix(1+5*time_interval,nrow=timeseries_length,ncol=1),
                           "fl3"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1))
  )
  
  vec <- matrix(1:num_indvs, nrow=num_indvs, ncol=1)
  
  x1fl1 <- apply(vec, 1, function(x) {fda.usc::int.simpson2(time_interval, flfn$fl1, equi = TRUE, method = "TRAPZ")})
  x2fl2 <- apply(vec, 1, function(x) {fda.usc::int.simpson2(time_interval, cat_data$X[x,,2]*(flfn$fl2), equi = TRUE, method = "TRAPZ")})
  x3fl3 <- apply(vec, 1, function(x) {fda.usc::int.simpson2(time_interval, cat_data$X[x,,3]*(flfn$fl3), equi = TRUE, method = "TRAPZ")})
  
  linear_predictor <- matrix(x1fl1 + x2fl2+ x3fl3 + intercept_opt )
  linear_predictor_without <- matrix(x1fl1 + x3fl3+ intercept_opt )
 
  
  generate_y_indvs <- function(lp){
    ys <- apply(lp, 1, function(x){ rbinom(1, 1, 1/(1+exp(-x))) })
    count_iter <- 1
    min_occurrence <- round(num_indvs * 0.2)
    mc_step_size <- min_occurrence * 0.2
    max_iterations <- 500
    increment_every <- round(max_iterations * 0.2)
    while (count_iter < max_iterations && (length(ys) - min_occurrence < sum(ys) || sum(ys) < min_occurrence)){
      count_iter <- count_iter + 1
      ys <- apply(lp, 1, function(x){ rbinom(1, 1, 1/(1+exp(-x))) })
      if(count_iter %% increment_every){
        min_occurrence <- min_occurrence - mc_step_size
      }
    }
    # cat("Y generation count:", count_iter, "\n")
    return(ys)
  }
  Y_indvs <- generate_y_indvs(linear_predictor)
  Y_indvs_without <- generate_y_indvs(linear_predictor_without)
  
  prob_ind=1/(1+exp(-linear_predictor))
  
  truelist=list("TrueX1"=cat_data$X[,,1],
                "TrueX2"=cat_data$X[,,2],
                "TrueX3"=cat_data$X[,,3],
                "Truecatcurve"=cat_data$W,
                "fl"=flfn,
                "yis"=Y_indvs,
                "yis_without" = Y_indvs_without,
                "linear_predictor"=list("linearw"=linear_predictor,"linearwo"=linear_predictor_without),
                "prob_ind"=prob_ind)
  
  return(list("true"=truelist))
}

cfd_testing <- function(start_time, end_time, timeseries_length,
                        num_indvs, mu1_coef, mu2_coef, fl_choice,
                        response_family,test_type,gam_choice,
                        klen=3){
  # cat("CFD Testing \nNum Indvs:\t", num_indvs,
  #     "\nTimeseries Len:\t", timeseries_length,
  #     "\nfl_choice:\t", fl_choice,
  #     "\ntest_type:\t", test_type, "\n")
  
  timestamps01 <- seq(from = start_time, to = end_time, length=timeseries_length)
  cfd_test_data <- GenerateCategoricalFDTest(klen=3,mu1_coef,mu2_coef,
                                             num_indvs=num_indvs,
                                             timeseries_length = timeseries_length,
                                             time_interval = timestamps01,
                                             fl_choice=fl_choice)
  yip=sum(cfd_test_data$true$yis)/num_indvs
  yip_wo=sum(cfd_test_data$true$yis_without)/num_indvs
  result <- cfd_hypothesis_test(cfd_test_data$true$yis,
                                cfd_test_data$true$Truecatcurve,
                                time_interval = timestamps01,
                                response_family=response_family,
                                test_type=test_type,gam_choice)
  #save(result, file = "result2.RData")
  return(list("pvalue"=result$pvalue,"yip"=yip,"yip_wo"=yip_wo,"pvalue2"=result$pvalue2))

  # return(list("pvalue"=result$pvalue,"test_statistics"=result$statistics,
  #             "yis"=cfd_test_data$true$yis,"flt"=cfd_test_data$true$fl,
  #             "W"=cfd_test_data$true$Truecatcurve,
  #             "linear_predictor"=cfd_test_data$true$linear_predictor,
  #             "prob_ind"=cfd_test_data$true$prob_ind))
}



