############################################################
# Copyright 2023 Xiaoxia Champon
  
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
# Purpose: Functions to perform Categorical Functional Data Hypothesis Testing
#           File that contains all the functions necessary to generate data 
# Author:  Xiaoxia Champon
# Date: 10/26/2023
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

#' Create directories 
if (!dir.exists("outputs")){
  dir.create("outputs")
}

#' Funciton to do the hypothesis testing on categorical functional data
#' @param Y: 1D array, length num_indvs, class membership for num_indvs individuals
#' @param cfd: 2d array, num_indvs* timeseries_length, categorical functional data 
#'              observed for num_indvs individuals and timeseries_length time points
#' @param time_interval, 1D array, observational time
#' @param response_family, currently only support 'bernoulli'           
#' @param test_type, "Function" or "constant"
#' @return list: statistics is test statistics, pvalue          
cfd_hypothesis_test <- function (Y, cfd, time_interval, response_family, test_type){
  
  # test.mat2<-data.frame(Y=Y,X=Xmat_Func2,Z.test=Zmat_Func.mat2,Z.test3=Zmat_Func.mat32,ones=rep(1,nsub))
  # names(test.mat2)<-c('Y','X1','X2',"X3",paste0('Z.test',1:ncol(Zmat_Func.mat2)),paste0('Z.test3',1:ncol(Zmat_Func.mat32)),"ones")
  # 
  
  alternative_fit <- fit.glmmPQL(test_matrix, response_family, num_indivs, test_type)
  
  result <- try(test.aRLRT(alternative_fit), silent=T)$aRLRT # Functional only
  
  return(list(statistics=result$statistic, pvalue=result$p.value))
}


#' Funciton to get the Zmatrix for hyppothesis testing
#' @param X_matrix: 2d array, num_indvs* timeseries_length, bernoulli functional data 
#'              observed for num_indvs individuals and timeseries_length time points
#' @param time_interval, 1D array, observational time
#' @param test_type, "Function" or "constant"
#' @param number_basis, scalar, default value is 30
#' @return list: Zmatrix, X_matrix_part2,
#'               J_matrix,D_difference_matrix,
#'               phi=bspline,Q=Q,Q2=Q2,Lambda1.inv.half=Lambda1.inv.half
get_Zmatrix <- function(X_matrix, time_interval, test_type, number_basis =30){
  
  knots<-construct.knots(time_interval,knots=(number_basis-3),knots.option='equally-spaced')
  bspline<-splineDesign(knots=knots,x=time_interval,ord=4)
  
  number_row <- nrow(X_matrix)
  number_col <- number_basis
  
  ##parallel use each row of X to multiple each column of bspline
  
  # J_matrix<-matrix(0,nrow=nrow(X_matrix),ncol=number_basis) #empty
  # for(row in 1:number_row){
  #   for(col in 1:number_col){
  #     J_matrix[row,col]<-integral_penalty(time_interval,X_matrix[row,]*bspline[,col])$value
  #   }
  # }
  
  J_list <- foreach(row = 1:number_row) %doRNG%
    {
      source("./R/integral_penalty_function.R")
      
      temp <- array(-123,col)
      for(col in 1:number_col){
        temp[col] <- integral_penalty(time_interval,X_matrix[row,]*bspline[,col])$value
      }
      return(temp)
    }
  J_matrix <- do.call(rbind, J_list)
  
  constant<-sqrt(1/number_basis)*rep(1,number_basis) #Q2_1
  range_time_interval<-max(time_interval)-min(time_interval) #what is full scale
  constant<-rep(1,number_basis)/range_time_interval # unscaled
  lin<-seq(min(time_interval),max(time_interval),length.out=number_basis)/range_time_interval #unscaled
  line<-sqrt(c(1/t(lin)%*%lin))*lin #Q2_2
  
  D<-diag(ncol(J_matrix))
  if(test=='Inclusion'){
    difference_penalty=0
    #return(list(Zmat=(ximat %*% J_matrix), X.g2=NULL,J_matrix=J_matrix,D=D))
    return(list(Zmat=( J_matrix), X.g2=NULL,J=J_matrix,D=D))
  }
  # if(test=='Linearity'){ #Test for linearity
  #   d=2
  #   Q2=as.matrix(cbind(constant,lin))
  # }
  if(test=='Functional'){ #Test for functional form
    difference_penalty=1
    Q2=as.matrix(constant)
  }
  D<-diff(D,differences=difference_penalty)
  P<- t(D)%*%D #penalty matrix
  P2<-1/2*(P+t(P))
  P.eigen<-eigen(P2)
  evalues<-P.eigen$values[1:nrow(D)]
  Q<-P.eigen$vectors
  Lambda1.inv.half<-diag(sqrt(1/evalues))
  #Q2<-Q[,(number_basis-d+1):number_basis]
  # 
  # Ztilde<- ximat %*% J_matrix %*% Q[,1:(number_basis-d)]
  # X.g2<- ximat %*% J_matrix %*% Q2
  Ztilde<- J_matrix %*% Q[,1:(number_basis-d)]
  X.g2<- J_matrix %*% Q2
  list(Zmat=Ztilde%*%Lambda1.inv.half, X.g2=X.g2,
       J=J_matrix,D=D,phi=bspline,Q=Q,Q2=Q2,Lambda1.inv.half=Lambda1.inv.half)
  
}


#use fpca.face to generate eigen 5functions
construct.knots <- function(argvals,knots,knots.option,p){
  if(length(knots)==1){
    allknots <- select.knots(argvals,knots,option=knots.option)
  }
  if(length(knots)>1){
    K = length(knots)-1 
    knots_left <- 2*knots[1]-knots[p:1+1]
    knots_right <- 2*knots[K] - knots[K-(1:p)]
    allknots <- c(knots_left,knots,knots_right)
  }
  return(allknots)
}
select.knots <- function(t,knots=27,p=3,option="equally-spaced"){
  qs <- seq(0,1,length=knots+1)
  if(option=="equally-spaced"){
    knots <- (max(t)-min(t))*qs + min(t)
  }
  if(option=="quantile"){
    knots <- as.vector(quantile(t,qs))
  }
  K <- length(knots)
  knots_left <- 2*knots[1]-knots[p:1+1]
  knots_right <- 2*knots[K] - knots[K-(1:p)]
  return(c(knots_left,knots,knots_right))
}

if(run_parallel)
{
  parallel::stopCluster(cl = my.cluster)
  initialized_parallel <- FALSE
}