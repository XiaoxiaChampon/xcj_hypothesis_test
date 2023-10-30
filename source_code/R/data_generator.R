
library(mgcv)
library(fda)
library(fda.usc)
library(devtools)
install_github("stchen3/glmmVCtest")
library("glmmVCtest")
library(RLRsim)
library(MASS)
library(splines)

source("./R/cfd_hypothesis_test.R")


#' Get mu_1, mu_2 functions, and score_vals objects for a given context.
#' @param klen number of points along the score decay axis
#' @return A list that contains mu_1, mu_2, score_vals
#'
GetMuAndScore_2 <- function(klen)
{
    mu_1 <- function(t){ -0.64+4*t }

    mu_2 <- function(t){ 0.97+6*t^2 }

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
    cat("octave", num_indvs, k, num_indvs * k, "\n")
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

    # we could have just passed these arguments ???
    num_indvs <- ncol(prob_curves$p1)
    timeseries_length <- nrow(prob_curves$p1)
    cat("n:", num_indvs, "\tt:", timeseries_length, "\n")

    # better names for W and X ???
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



expit <- function(x){1/(1+exp(-x))}

fl2f = function(t) {
    return(t - 8/9)
}

fl3f = function(t) {
    return(0.5*t)
}

fl3fn = function(t) {
    return(-3*t^2 + 2*t - 0.9)
}

GenerateCategoricalFDTest <- function(klen, num_indvs, timeseries_length,
                                      time_interval, fl_choice){

    mns <- GetMuAndScore_2(klen)

    generated_data <- GenerateDataTest(num_indvs = num_indvs,
                                      timeseries_length = timeseries_length,
                                      mu_1 = mns$mu_1,
                                      mu_2 = mns$mu_2,
                                      score_vals = mns$score_vals,
                                      start_time = time_interval[1],
                                      end_time = tail(time_interval,1),
                                      k = klen)

    prob_curves <- list(p1 = generated_data$p1, p2 = generated_data$p2, p3 = generated_data$p3)
    cat_data <- GenerateCategFuncData(prob_curves)

    flfn <- switch(fl_choice,

                   "1"=list("fl1"=rep(0.45,timeseries_length),
                            "fl2"=rep(0.5,timeseries_length),
                            "fl3"=rep(-0.51,timeseries_length)),

                   "2"=list("fl1"=rep(-0.1,timeseries_length),
                            "fl2"=matrix(-0.1*fl2f(time_interval),nrow=timeseries_length,ncol=1),
                            "fl3"=matrix(fl3f(time_interval),nrow=timeseries_length,ncol=1)),

                   "3"=list("fl1"=rep(-0.2,timeseries_length),
                            "fl2"=matrix(-0.15*fl2f(time_interval),nrow=timeseries_length,ncol=1),
                            "fl3"=matrix(fl3f(time_interval),nrow=timeseries_length,ncol=1)),

                   "4"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
                            "fl2"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)+1.3145,
                            "fl3"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)),
                   )

    vec <- matrix(1:num_indvs, nrow=num_indvs, ncol=1)



    # integral a function on a interval, returns a scalar
    # x1fl1 <- parApply(my.cluster, vec, 1, function(x) {fda.usc::int.simpson2(time_interval, cat_data$X[x,,1]*flfn$fl1, equi = TRUE, method = "TRAPZ")})
    # x2fl2 <- parApply(my.cluster, vec, 1, function(x) {fda.usc::int.simpson2(time_interval, cat_data$X[x,,2]*flfn$fl2, equi = TRUE, method = "TRAPZ")})
    # x3fl3 <- parApply(my.cluster, vec, 1, function(x) {fda.usc::int.simpson2(time_interval, cat_data$X[x,,3]*flfn$fl3, equi = TRUE, method = "TRAPZ")})
    # 
    x1fl1 <- apply(vec, 1, function(x) {fda.usc::int.simpson2(time_interval, cat_data$X[x,,1]*flfn$fl1, equi = TRUE, method = "TRAPZ")})
    x2fl2 <- apply( vec, 1, function(x) {fda.usc::int.simpson2(time_interval, cat_data$X[x,,2]*flfn$fl2, equi = TRUE, method = "TRAPZ")})
    x3fl3 <- apply( vec, 1, function(x) {fda.usc::int.simpson2(time_interval, cat_data$X[x,,3]*flfn$fl3, equi = TRUE, method = "TRAPZ")})
    
    sum_int_xtft <- matrix(x1fl1 + x2fl2+ x3fl3 + 0.02)

    #Y_indvs <- parApply(my.cluster, sum_int_xtft, 1, function(x){ rbinom(1,1, 1/(1+exp(-x))) })
    
    Y_indvs <- apply(sum_int_xtft, 1, function(x){ rbinom(1,1, 1/(1+exp(-x))) })


    truelist=list("TrueX1"=cat_data$X[,,1],
                  "TrueX2"=cat_data$X[,,2],
                  "TrueX3"=cat_data$X[,,3],
                  "Truecatcurve"=cat_data$W,
                  "fl"=flfn,
                  "yis"=Y_indvs)

    ########get zistart
    #recover Z_i1 hat using X_i[1,all j, all num_indvs] only related to p1
    #Z_i1hat=Z_ihat(X_i1,t)
    #recover Z_i2 hat using X_i[2,all j, all num_indvs] only related to p2
    ##Z_i2hat=Z_ihat(X_i2,t)
    #Z_i3hat=Z_ihat(X_i3,t)

    #Z_i1hatstar=Z_i1hat+log(1+exp(Z_i3hat))-log(1+exp(Z_i1hat))-Z_i3hat
    #Z_i2hatstar=Z_i2hat+log(1+exp(Z_i3hat))-log(1+exp(Z_i2hat))-Z_i3hat


    #truel=list("TrueZ_i1"=Z_i1,"TrueZ_i2"=Z_i2)
    #est=list("EstimateZ_i1"=Z_i1hatstar,"EstimateZ_i2"=Z_i2hatstar)
    #return(list("true"=truelist,"est"=est))
    return(list("true"=truelist))
}





cfd_testing <- function(start_time, end_time, timeseries_length,
                        num_indvs,fl_choice,response_family,test_type,
                        klen=3){
  cat("CFD Testing Simulation\nNum Indvs:\t", num_indvs,
      "\nTimeseries Len:\t", timeseries_length,
      "\nfl_choice:\t", fl_choice,
      "\nNtest_type:\t", test_type)
  
  timestamps01 <- seq(from = start_time, to = end_time, length=timeseries_length)
  cfd_test_data <- GenerateCategoricalFDTest(klen=3,
                                    num_indvs=num_indvs,
                                    timeseries_length = timeseries_length,
                                    time_interval = timestamps01,
                                    fl_choice=fl_choice)
  
  result <- cfd_hypothesis_test(cfd_test_data$true$yis,
                                cfd_test_data$true$Truecatcurve,
                                time_interval = timestamps01,
                                response_family=response_family,
                                test_type=test_type)
  return(list("pvalue"=result$pvalue,"test_statistics"=result$statistics,
              "yis"=cfd_test_data$true$yis,"flt"=cfd_test_data$true$fl,
              "W"=cfd_test_data$true$Truecatcurve))
}



