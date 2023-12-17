
library(foreach)
library(doRNG)
  
# library(profvis)
# 
# profvis({
  
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
  
  ############################################################
  # set.seed(123)
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
  ############################################################
  # set.seed(123)
  # retval <- foreach(indv=1:num_indvs, .combine = cbind, .init = NULL, .packages = c("base","stats")) %dorng% {
  #   X_indv <- sapply(c(1:timeseries_length),
  #                    function(this_time) rmultinom(n=1,
  #                                                  size=1,
  #                                                  prob = c(prob_curves$p1[this_time,indv],
  #                                                           prob_curves$p2[this_time,indv],
  #                                                           prob_curves$p3[this_time,indv]) ))
  #   W_indv <- apply(X_indv, 2, which.max)
  #   list("x"=X_indv, "w"=W_indv)
  # }
  # X_array <- aperm( array(as.double(unlist(retval[1,])), c(curve_count,timeseries_length,num_indvs)), c(3,2,1))
  # W <- matrix(as.double(unlist(retval[2,])), ncol=num_indvs)
  ############################################################
  
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
                                      time_interval, fl_choice, lp_intercept=0.6206897){
  
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
  
  linear_predictor <- matrix(x1fl1 + x2fl2+ x3fl3 + lp_intercept )
  linear_predictor_without <- matrix(x1fl1 + x3fl3+ lp_intercept )
  
  
  generate_y_indvs <- function(lp){
    ys <- apply(lp, 1, function(x){ rbinom(1, 1, 1/(1+exp(-x))) })
    leastOccr <- min(sum(ys), length(ys) - sum(ys))
    count_iter <- 1
    min_occurrence <- round(num_indvs * 0.2)
    max_iterations <- 100
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
                "yis_without" = Y_indvs_without,
                "linear_predictor"=list("linearw"=linear_predictor,"linearwo"=linear_predictor_without),
                "prob_ind"=prob_ind)
  
  return(list("true"=truelist))
}

# begin_exp_time <- Sys.time()
# 
# timeseries_length = 180
# timestamps01 <- seq(from = 0.01, to = 0.99, length=timeseries_length)
# mu1_coef=c(-6.67,-2.47,5.42)
# mu2_coef=c(-3.14,-0.99,3.91)
# for (idx in 1:20) {
#   results <- GenerateCategoricalFDTest(3, mu1_coef, mu2_coef, 500, timeseries_length, timestamps01, "8", 1.5)
# }
# 
# end_exp_time <- Sys.time()
# 
# cat("\n====================\n",
#     "\tAll Experiemnts Took:", capture.output(end_exp_time - begin_exp_time),
#     "\n====================\n")


# })