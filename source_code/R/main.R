############################################################
# Purpose: Functions Needed for Categorical Functional Data Hypothesis Testing
#           File that contains all the functions necessary to generate data 
# Author:  Xiaoxia Champon
# Date: 10/26/2023
##############################################################

#load thel library
# For: profiling and visualization of profiling
library(profvis)

# For: gam 
library(mgcv)

# For: cubicspline
library(pracma)

library(fda)
library(fda.usc)
#library needed to perform the hypothesis testing for continuous functional data
library(RLRsim)

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

#' Create directories 
if (!dir.exists("outputs")){
  dir.create("outputs")
}
if (!dir.exists("outputs/clustersims")){
  dir.create("outputs/clustersims")
}

#' function of category 2 effect
#' @param t: 1D vector, time interval where the categorical functional data observed
#' @return fl2: 1D vector, effect of category 2 over time
fl2f <- function(t) {
  fl2 <- t-8/9
  return(fl2)
}


#' function of category 3 effect
#' @param t: 1D vector, time interval where the categorical functional data observed
#' @return fl3: 1D vector, effect of category 3 over time
fl3f <- function(t) {
  fl3=0.5*t
  return(fl3)
}


#' function of category 3 effect
#' @param t: 1D vector, time interval where the categorical functional data observed
#' @return fl3: 1D vector, effect of category 3 over time
fl3fn = function(t) {
  fl3=3*t^2+2*t-0.9
  return(fl3)
}


p_ihat=function(Z_i1app,Z_i2app){
  denom=(1+exp(Z_i1app)+exp(Z_i2app))
  p_i1h=exp(Z_i1app)/denom
  p_i2h=exp(Z_i2app)/denom
  p_i3h=1/denom
  return(list("p_1hatmatrix"=p_i1h,"p_2hatmatrix"=p_i2h,"p_3hatmatrix"=p_i3h))
}

#' function to generate categorical funcitonal data and response
#' @param num_indvs: scalar- number of subjects
#' @param timeseries_length: scalar -number of time points
#' @param sparse=1 yes   sparse=0 no
#' @param scorevar=2 bigger var , scorevar=1 smaller var
#' @param k, scalar: number of eigen functions
#' @param q, scalar: level of the categorical level

generate_cfd_test=function(k,num_indvs,timeseries_length,sparse,scorevar,seed=123,st,et,fl,q=3){
  cat("Cluster Simulation\nNum Indvs:\t", num_indvs,
      "\nTimeseries Len:\t", timeseries_length,
      "\nScenario:\t", scenario,
      "\nNum Replicas:\t", num_replicas)
  
  
  time_elapsed <<- list()
  # "Xiaoxia"=NULL, "univfpca"=NULL, "kmeans"=NULL, "fadp"=NULL, "dbscan"=NULL, "cfd"=NULL)
  last_time <- 0
  row_name <- NULL
  timeKeeperStart <- function(rn)
  {
    row_name <<- rn
    if(FALSE == row_name %in% names(time_elapsed))
    {
      time_elapsed[[row_name]] <<- NULL
    }
    last_time <<- Sys.time()
  }
  timeKeeperNext <- function()
  {
    this_time <- Sys.time()
    this_section_time <- this_time - last_time
    cat(row_name, "calc time taken:", capture.output(this_section_time), "\n")
    time_elapsed[[row_name]] <<- append(time_elapsed[[row_name]], this_section_time)
    last_time <<- this_time
  }
  
  
  
  if(sparse==1){
    mu_1=function(t){
      3.8+4*t  #-0.64+4*t    #4.5+4*t   #sparse 3.8+4*t
      
    }
    mu_2=function(t){
      1.5+4*t^2    #0.97+6*t^2
      
    }
  
  }
  
  
  if (sparse==0){
    mu_1=function(t){
      #3.8+4*t  
      -0.64+4*t    #4.5+4*t   #sparse 3.8+4*t
      
    }
    mu_2=function(t){
      #1.5+4*t^2    
      0.97+6*t^2
    }
    
  }
  
  #####10/25/2022
  if (sparse==3){
    mu_1=function(t){
      #3.8+4*t  
      #-0.64+4*t    #4.5+4*t   #sparse 3.8+4*t
      t+1
      
    }
    mu_2=function(t){
      #1.5+4*t^2    
      #0.97+6*t^2
      0.3*t
      
    }
    
    p_ihat=function(Z_i1app,Z_i2app){
      denom=(1+exp(Z_i1app)+exp(Z_i2app))
      p_i1h=exp(Z_i1app)/denom
      p_i2h=exp(Z_i2app)/denom
      p_i3h=1- p_i1h- p_i2h
      return(list("p_1hatmatrix"=p_i1h,"p_2hatmatrix"=p_i2h,"p_3hatmatrix"=p_i3h))
    }
    
  }
  
  
  
  if (sparse==4){
    mu_1=function(t){
      #3.8+4*t  
      #-0.64+4*t    #4.5+4*t   #sparse 3.8+4*t
      #10*t+1
      t+1
      
    }
    mu_2=function(t){
      #1.5+4*t^2    
      #0.97+6*t^2
      #10*t+3
      #t+3
      t-1
      
    }
    
    p_ihat=function(Z_i1app,Z_i2app){
      denom=(1+exp(Z_i1app)+exp(Z_i2app))
      p_i1h=exp(Z_i1app)/denom
      p_i2h=exp(Z_i2app)/denom
      p_i3h=1- p_i1h- p_i2h
      #p_i3h=1/denom
      return(list("p_1hatmatrix"=p_i1h,"p_2hatmatrix"=p_i2h,"p_3hatmatrix"=p_i3h))
    }
    
  }
  
  mu_vec=rep(0,k)
  
  
  psi_fn=function(k){
    
    psi_k1=matrix(rep(1,length(t)*k),ncol=k)  
    psi_k2=matrix(rep(1,length(t)*k),ncol=k) 
    for (i in 1:k) {
      psi_k1[,i]=sin(2*i*pi*t )
      psi_k2[,i]=cos(2*i*pi*t )
    }
    list("psi_k1"=psi_k1,"psi_k2"=psi_k2)
  }
  
  
  t=seq(from = st,to = et, length=datapoints)
  
  X_i=array(0,dim=c(q,datapoints,n))  #multinormial results: row is level q, column is time points, n is the number of subjects, each column only has one row of 1 and every other rows are 0
  X_nt=matrix(rep(1,n*length(t)),nrow=n,ncol=length(t))  #true observations of categorical-valued outcome, each row represent one subject, columns represent time points
  score_matrix=matrix(rep(1,n*k),nrow=n,ncol=k)  #row is number of subjects, column is the number of eigen functions
  psi_score_matrix_1=matrix(rep(1,n*length(t)),ncol=n)  #dim: length(t)*nsubjects
  psi_score_matrix_2=matrix(rep(1,n*length(t)),ncol=n)
  Z_i1=matrix(rep(1,n*length(t)),nrow=n)  #True latent curves1:row is n subjects, col is t time points
  Z_i2=matrix(rep(1,n*length(t)),nrow=n) #True latent curve 2
  p_i1=matrix(rep(0,n*length(t)),nrow=n)  #True p_i1
  p_i2=matrix(rep(0,n*length(t)),nrow=n)  #True p_i2
  p_i3=matrix(rep(0,n*length(t)),nrow=n)  #True p_i3
  for (i in 1:n){
    set.seed(seed+i)
    
    if (k==3){
      if (scorevar==1){
        #score varies based on i
        score_1=score(0,1)
        score_2=score(0,sqrt(1/2))
        score_3=score(0,1/2)
        
      }
      
      if (scorevar==2){
        #score varies based on i
        score_1=score(0,1)
        score_2=score(0,sqrt(1/3))
        score_3=score(0,1/3)
      }
      
      if (scorevar==3){
        #score varies based on i
        score_1=score(0,1)
        score_2=score(0,sqrt(1/4))
        score_3=score(0,1/4)
      }
      
      
      if (scorevar==4){
        #score varies based on i
        score_1=score(-0.5,1)
        score_2=score(1,sqrt(1/2))
        score_3=score(0.25,1/2)
      }
      
      score_vector=cbind(score_1,score_2,score_3)
    }
    
    
    
    if (k==4){
      
      if (scorevar==1){
        score_1=score(0,1)
        score_2=score(0,sqrt(1/2))
        score_3=score(0,1/2)
        score_4=score(0,sqrt(1/8))
        # cpve=cumsum(c(1,sqrt(1/2),1/2,sqrt(1/8)))/sum(c(1,sqrt(1/2),1/2,sqrt(1/8)))
        # cvar=c(1,sqrt(1/2),1/2,sqrt(1/8))
      }
      if (scorevar==2){
        score_1=score(0,1)
        score_2=score(0,sqrt(1/3))
        score_3=score(0,1/3)
        score_4=score(0,sqrt(1/27))
        # cpve=cumsum(c(1,sqrt(1/3),1/3,sqrt(1/27)))/sum(c(1,sqrt(1/3),1/3,sqrt(1/27)))
        # cvar=c(1,sqrt(1/3),1/3,sqrt(1/27))
      }
      
      
      
      if (scorevar==3){
        score_1=score(0,1)
        score_2=score(0,sqrt(1/4))
        score_3=score(0,1/4)
        score_4=score(0,sqrt(1/64))
        # cpve=cumsum(c(1,sqrt(1/4),1/4,sqrt(1/64)))/sum(c(1,sqrt(1/4),1/4,sqrt(1/64)))
        # cvar=c(1,sqrt(1/4),1/4,sqrt(1/64))
      }
      score_vector=cbind(score_1,score_2,score_3,score_4)
      
    }
    
    
    
    
    
    
    
    psi_k1=psi_fn(k)$psi_k1
    psi_k2=psi_fn(k)$psi_k2
    
    #Z varies based on i
    #psi t*k, score: t*k,  psi%*%t(score)
    psi_score_matrix_1[,i]=psi_k1%*%t(score_vector)
    Z_i1[i,]=mu_1(t)+psi_score_matrix_1[,i]
    
    psi_score_matrix_2[,i]=psi_k2%*%t(score_vector)
    Z_i2[i,]=mu_2(t)+psi_score_matrix_2[,i]
    
    
    #p varies based on i
    denominator=(1+exp(as.vector(Z_i1[i,]))+exp(as.vector(Z_i2[i,])))
    p_i1[i,]=(exp(as.vector(Z_i1[i,])))/denominator
    p_i2[i,]=(exp(as.vector(Z_i2[i,])))/denominator
    p_i3[i,]=1-p_i1[i,]-p_i2[i,]
    
    
    #X_i varies based on i
    #X_i=matrix(rep(1,k*length(t)),nrow=k,ncol=length(t))
    
    for (j in 1:length(t)){
      X_i[,j,i]=rmultinom(n=1, size=1, prob=c(p_i1[i,j],p_i2[i,j],p_i3[i,j]))
    }
    
    #X_it varies based on i
    X_it=c(1)
    for (j in 1:length(t)){
      X_it[j]=as.vector(which(X_i[,j,i] == 1))
    }
    X_nt[i,]=X_it
    
    #collect score matrix
    score_matrix[i,]=score_vector
  }
  
  #collect value and graph
  #collect first two rows of observed binary curves
  X_i1=t(X_i[1,,])  #all n row subjects , t columns values related to p1
  X_i2=t(X_i[2,,]) #all n row subjects , t columns values related to p2
  X_i3=t(X_i[3,,]) #all n row subjects , t columns values related to p3
  
  #generate Fl functions
  #Formof delta0(t ): (a) Scalar: , (b) Linear: 1+delta1t , (c) Trigonometric: 1+ t +delta2 cos(2pit )
  ##################################################  
  if(fl==1){
    fl1=rep(0.45,datapoints)
    fl2=rep(0.5,datapoints)
    #fl2=matrix(fl2f(t),nrow=datapoints,ncol=1)
    fl3=rep(-0.51,datapoints)
  }
  
  if (fl==2){
    # fl1=rep(-0.5,datapoints)
    # fl2=matrix(fl2f(t)+0.1,nrow=datapoints,ncol=1)
    # #fl2=rep(-0.3,datapoints)
    # fl3=matrix(fl3f(t)-0.02,nrow=datapoints,ncol=1)
    
    
    fl1=rep(-0.1,datapoints)
    fl2=matrix(-0.1*fl2f(t),nrow=datapoints,ncol=1)
    #fl2=rep(-0.3,datapoints)
    fl3=matrix(fl3f(t),nrow=datapoints,ncol=1)
    
  }
  #############################
  #73, 27
  # if (fl==3){
  #  fl1=rep(-0.2,datapoints)
  #  fl2=matrix(-0.1*fl2f(t),nrow=datapoints,ncol=1)
  #  #fl2=rep(-0.3,datapoints)
  #  fl3=matrix(fl3f(t),nrow=datapoints,ncol=1)
  # }
  
  ##############################
  #
  if (fl==3){
    fl1=rep(-0.2,datapoints)
    fl2=matrix(-0.15*fl2f(t),nrow=datapoints,ncol=1)
    #fl2=rep(-0.3,datapoints)
    fl3=matrix(fl3f(t),nrow=datapoints,ncol=1)
  }
  
  ####################
  #
  if (fl==4){
    fl3=matrix(fl3fn(t),nrow=datapoints,ncol=1)
    fl1=fl3-0.09
    fl2=fl3+1.3145
  }
  
  ####################
  
  flfn=list("fl1"=fl1,"fl2"=fl2,"fl3"=fl3)
  
  ###############################
  #generate pi
  vec=matrix(1:n,nrow=n,ncol=1)
  #integral a function on a interval, returns a scalar
  x1fl1=apply(vec,1, function (x) {int.simpson2(t, X_i1[x,]*fl1, equi = TRUE, method = "TRAPZ")})
  x2fl2=apply(vec,1, function(x) {int.simpson2(t, X_i2[x,]*fl2, equi = TRUE, method = "TRAPZ")})
  
  x3fl3=apply(vec,1, function(x) {int.simpson2(t, X_i3[x,]*fl3, equi = TRUE, method = "TRAPZ")})
  
  pis=c(0)
  z_lin=c(0)
  for (i in 1:n){
    #pi=beta0+xqfl1+x2fl2
    z_lin[i]=0.02+sum(x1fl1,x2fl2,x3fl3)
    pis[i]=expit(z_lin[i])
    
  }
  
  #############
  #generate Yi
  yis=c(0)
  for (i in 1:n){
    yis[i]=rbinom(1,1,pis[i])
  }
  
  
  truelist=list("TrueX1"=X_i1,"TrueX2"=X_i2,"TrueX3"=X_i3,"Truecatcurve"=X_nt,"fl"=flfn,"yis"=yis,"z_lin"=z_lin)
  
  ########get zistart
  #recover Z_i1 hat using X_i[1,all j, all n] only related to p1
  #Z_i1hat=Z_ihat(X_i1,t)
  #recover Z_i2 hat using X_i[2,all j, all n] only related to p2 
  ##Z_i2hat=Z_ihat(X_i2,t)
  #Z_i3hat=Z_ihat(X_i3,t)
  
  #Z_i1hatstar=Z_i1hat+log(1+exp(Z_i3hat))-log(1+exp(Z_i1hat))-Z_i3hat
  #Z_i2hatstar=Z_i2hat+log(1+exp(Z_i3hat))-log(1+exp(Z_i2hat))-Z_i3hat
  
  
  #truel=list("TrueZ_i1"=Z_i1,"TrueZ_i2"=Z_i2)
  #est=list("EstimateZ_i1"=Z_i1hatstar,"EstimateZ_i2"=Z_i2hatstar)
  #return(list("true"=truelist,"est"=est))
  return(list("true"=truelist))
}




#' Function to produce functional dummy variables X from categorical functional data W
#' @param W 2D array, t*n: t is the timestamp and n is the number of the observation
#' @return X 3D array, n*t*Q, Q: the total number of the category
GetXFromW <- function(W)
{
  num_indv <- ncol(W)
  timeseries_length <-nrow(W)
  category_count<- length(unique(c(W)))
  Q_vals <- unique(c(W))
  if(is.numeric(Q_vals)) Q_vals<- sort(Q_vals)
  
  X<- array(0, c(num_indv,timeseries_length,category_count))
  for(indv in 1:num_indv)
  {
    for(timestamps01 in 1:timeseries_length)
    {
      X[indv, timestamps01, which(Q_vals==W[, indv][timestamps01])] <- 1
    }
  }
  return(X)
}

#' Function to select 
#' @param choice "probit", "binomial",  or "multinormial"
#' @param timestamps01, 1D array, time interval that cfd is observed
#' @param  W: 2D array, t*n, t: the number of time points, n: the number of individuals
#' @param  basis_size=25, the number of basis function used 
#' @param  method="ML"
#' @return list of 2D array: True Z curves , Est Z curves, True p curves, Est p curves
#'                           all have dimension t*n
EstimateCategFuncData <- function(choice, timestamps01, W, basis_size=25, method="ML")
{
  if(choice == "probit"){
    X <- GetXFromW(W)
    return(EstimateCategFuncData_probit(timestamps01, X, basis_size, method, 1/150))
  }else if(choice == "binomial"){
    X <- GetXFromW(W)
    return(EstimateCategFuncData_binorm(timestamps01, X, basis_size, method))
  }else if(choice == "multinomial"){
    return(EstimateCategFuncData_multinormial(timestamps01, W, basis_size, method))
  }
}

#'Function to estimate z and p using wood_multinormial
#' @param timestamps01, 1D array, time interval that cfd is observed
#' @param  W: 2D array, t*n, t: the number of time points, n: the number of individuals
#' @param  basis_size=25, the number of basis function used 
#' @param  method="ML"
#' @return list of 2D array: True Z curves , Est Z curves, True p curves, Est p curves
#'                           all have dimension t*n

EstimateCategFuncData_multinormial <- function(timestamps01, W, basis_size=25, method="ML")
{
  
  num_indv<- ncol(W)
  timeseries_length <-nrow(W)
  
  Z<-NULL
  prob<-array(0, c(num_indv, timeseries_length , 3))
  for (i in 1:num_indv){
    fit_binom<-gam(list(W[,i]-1~s(timestamps01,bs = "cr", m=2, k = basis_size),
                        ~s(timestamps01,bs = "cr", m=2, k = basis_size)),
                   family=multinom(K=2), method = method,
                   control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
                   optimizer=c("outer","bfgs")) 
    
    z1<- fit_binom$linear.predictors[,1]
    z2<- fit_binom$linear.predictors[,2]
    Z<- cbind(Z, c(z1,z2))
    ##find probability
    Z_cbind=cbind(z1,z2)
    exp_z=exp(Z_cbind)
    denominator_p=1+exp_z[,1]+exp_z[,2]
    p1 <- exp_z[,1]/denominator_p
    p2 <- exp_z[,2]/denominator_p
    p3=1/denominator_p
    prob[i,,] <- cbind(p1, p2, p3)
    
  }
  
  return(list(Z1_est=Z[1:timeseries_length ,], Z2_est=Z[1:timeseries_length +timeseries_length ,], 
              p1_est=t(prob[,,1]), p2_est=t(prob[,,2]), p3_est=t(prob[,,3]) ))
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
      source("R/acj/run_gam_function.R")
      
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
      source("R/acj/run_gam_function.R")
      
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



GenerateCategFuncData <- function(prob_curves)
{
  curve_count <- length(prob_curves);
  
  # we could have just passed these arguments ???
  num_indvs <- ncol(prob_curves$p1)
  timeseries_length <- nrow(prob_curves$p1)
  
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

#' Get clustered data
#'
#'
GenerateClusterData <- function(setting, scenario, k, num_indvs, timeseries_length)
{
  setting_object <- GetMuAndScore(setting, scenario, k)
  cluster_f <- GenerateClusterDataScenario(num_indvs,
                                           timeseries_length,
                                           k,
                                           mu_1 = setting_object$mu_1,
                                           mu_2 = setting_object$mu_2,
                                           score_vals = setting_object$score_vals)
  return (cluster_f)
}

#' Get fraction of occurrence of each class for a given scenario
#' @param scenario scenario name as a string "A", "B", "C"
#' @return a vector containing the fractions
#'
GetOccurrenceFractions <- function(scenario)
{
  occur_fraction <- switch (scenario,
                            "A" = c(0.75, 0.22, 0.03),
                            "B" = c(0.5, 0.3, 0.2),
                            "C" = c(0.1, 0.6, 0.3)
  )
  
  return (occur_fraction)
}

#' Get mu_1, mu_2 functions, and score_vals objects for a given context.
#' @param setting setting identified as an integer 1,2,3
#' @param scenario scenario name as a string "A", "B", "C"
#' @param k number of points along the score decay axis
#' @return A list that contains mu_1, mu_2, score_vals
#'
GetMuAndScore <- function(setting, scenario, k)
{
  all_score_values = rep(0, k)
  
  if(1 == setting)
  {
    mu_1 <- function(t) -1 + 2 * t + 2 * t^2
    
    mu_2 <- switch(scenario,
                   "A" = function(t) -2.5 + exp(t * 2),
                   "B" = function(t) -0.5 + exp(t * 2),
                   "C" = function(t) -2.5 + exp(t * 2)
    )
    score_front <- switch(scenario,
                          "A" = c(1, 1/2, 1/4),
                          "B" = c(1, 1/2, 1/4),
                          "C" = c(50, 25, 5)
    )
  } else if(2 == setting)
  {
    mu_1 <- function(t) 4 * t^2 - 1.2
    
    mu_2 <- function(t) 4 * t^2 - 3.5
    
    score_front <- c(1, 1/2, 1/4)
  } else if(3 == setting)
  {
    mu_1 <- function(t) -2.2 + 4 * t^2
    
    mu_2 <- function(t) -7 + 6 * t^2
    
    score_front <- c(1, 1/4, 1/16)
  }
  
  for(idx in 1:length(score_front))
  {
    all_score_values[idx] <- score_front[idx]
  }
  
  return(list("mu_1" = mu_1, "mu_2" = mu_2, "score_vals" = all_score_values))
}


#' Generate cluster data for a given scenario
#' @param num_indvs number of individuals
#' @param timeseries_length length of time-series as an integer
#' @param k  number of eigen(psi) functions
#' @param mu_1 mean function for the first latent curve
#' @param mu_2 mean function for the second latent curve
#' @param score_vals the variance of the principal component scores
#'
GenerateClusterDataScenario <- function(num_indvs,
                                        timeseries_length,
                                        k = 3,
                                        mu_1,
                                        mu_2,
                                        score_vals)
{
  timestamps01 <- seq(from = 0.0001, to = 1, length=timeseries_length)
  
  # noise octaves
  # cat("octave", num_indvs, k, num_indvs * k, "\n")
  scores_standard <- matrix(rnorm(num_indvs * k), ncol = k)
  scores <- scores_standard %*% diag(sqrt(score_vals))
  
  #
  BIG_mu <- c(mu_1(timestamps01), mu_2(timestamps01))
  BIG_phi <- PsiFunc(k, timestamps01)
  
  Z <- BIG_phi %*% t(scores) + BIG_mu
  Z1 <- Z[1:timeseries_length, ]
  Z2 <- Z[1:timeseries_length + timeseries_length, ]
  expZ1 <- exp(Z1)
  expZ2 <- exp(Z2)
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

  