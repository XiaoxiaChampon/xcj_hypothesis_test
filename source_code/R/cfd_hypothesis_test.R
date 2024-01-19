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

#' Create directories
if (!dir.exists("outputs")){
  dir.create("outputs")
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

  X <- array(0, c(num_indv,timeseries_length,category_count))
  for(indv in 1:num_indv)
  {
    for(timestamps01 in 1:timeseries_length)
    {
      X[indv, timestamps01, which(Q_vals==W[, indv][timestamps01])] <- 1
    }
  }
  return(X)
}

#' Funciton to do the hypothesis testing on categorical functional data
#' @param Y: 1D array, length num_indvs, class membership for num_indvs individuals
#' @param cfd: 2d array,timeseries_length* num_indvs, categorical functional data
#'              observed for num_indvs individuals and timeseries_length time points
#' @param time_interval, 1D array, observational time
#' @param response_family, currently only support 'bernoulli'
#' @param test_type, "Function" or "constant"
#' @return list: statistics is test statistics, pvalue
cfd_hypothesis_test <- function(Y, cfd, time_interval, response_family, test_type){
  #X will be n*t*Q
  X_cfd <- GetXFromW(cfd)
  
  num_indvs <- length(Y)

  Zmat_test_type_2 <- get_Zmatrix(X_cfd[,,2], time_interval, test_type)
  Zmat_test_type_3 <- get_Zmatrix(X_cfd[,,3], time_interval, test_type)

  Xmat_test_type <- matrix(rep(1, num_indvs), ncol=1)
  
  if (test_type=="Functional"){
    Xmat_test_type <- cbind(Xmat_test_type, Zmat_test_type_2$X.g2, Zmat_test_type_3$X.g2)
  }
  
  test_matrix <- data.frame(Y=Y,
                            X=Xmat_test_type,
                            Z.test=Zmat_test_type_2$Zmat,
                            Z.test3=Zmat_test_type_3$Zmat,
                            ones=rep(1,num_indvs))
 
  if(test_type=="Inclusion"){
    names(test_matrix) <- c('Y','X1',
                            paste0('Z.test',1:ncol(Zmat_test_type_2$Zmat)),
                            paste0('Z.test3',1:ncol(Zmat_test_type_3$Zmat)),
                            "ones")
  } else if (test_type=="Functional"){
    names(test_matrix) <- c('Y','X1','X2',"X3",
                            paste0('Z.test',1:ncol(Zmat_test_type_2$Zmat)),
                            paste0('Z.test3',1:ncol(Zmat_test_type_3$Zmat)),
                            "ones")
  }
  
  
  #For testing in models with multiple variance
  #' components, the fitted model \code{m} must contain \bold{only} the random
  #' effect set to zero under the null hypothesis, while \code{mA} and \code{m0}
  #' are the models under the alternative and the null, respectively. 
  
  # alternative_fit <- fit.glmmPQL(test_matrix, response_family, num_indvs, test_type)
  # 
  # result_try <- try(test.aRLRT(alternative_fit), silent=T)
  # 
  # if(is.atomic(result_try)){
  #   return(list(statistics=NULL, pvalue=NULL))
  # }
  # 
  # result <- result_try$aRLRT 
  
  gam_test <- gam(cbind(Y, num_indvs - Y) ~ 0 + Xmat_test_type + 
                     s(Zmat_test_type_2$Zmat, bs = 're')+ s(Zmat_test_type_3$Zmat, bs = 're'), family = 'binomial')
  
  
 
  #return(list(statistics=result$statistic, pvalue=result$p.value))
  
  return(list(statistics=gam_test[1,3], pvalue=gam_test[1,4]))
}


#' Function to get the Zmatrix for hyppothesis testing
#' @param X_matrix: 2d array, num_indvs* timeseries_length, bernoulli functional data
#'              observed for num_indvs individuals and timeseries_length time points
#' @param time_interval, 1D array, observational time
#' @param test_type, "Function" or "constant"
#' @param number_basis, scalar, default value is 30
#' @return list: Zmatrix, X_matrix_part2,
#'               J_matrix,D_difference_matrix,
#'               phi=bspline,Q=Q,Q2=Q2,Lambda1.inv.half=Lambda1.inv.half
get_Zmatrix <- function(X_matrix, time_interval, test_type, number_basis =30){
  #test
  ######
  #test_type="Functional"
  #X_matrix=X_cfd[,,2]
  #setequal(X_matrix,X_matrix_not)
  ######
  knots <- construct.knots(time_interval,knots=(number_basis-3),knots.option='equally-spaced')
  bspline <- splineDesign(knots=knots,x=time_interval,ord=4)
  
  number_row <- nrow(X_matrix)
  number_col <- number_basis
  
  ##parallel use each row of X to multiple each column of bspline

  # J_matrix <- matrix(0,nrow=nrow(X_matrix),ncol=number_basis) #empty
  # for(row in 1:number_row){
  #   for(col in 1:number_col){
  #     J_matrix[row,col] <- integral_penalty(time_interval,X_matrix[row,]*bspline[,col])$value
  #   }
  # }

  J_matrix <- foreach(this_row = 1:number_row) %do%
    {
      source("source_code/R/integral_penalty_function.R")

      temp <- array(-123, number_col)
      for(this_col in 1:number_col){
        temp[this_col] <- integral_penalty(time_interval,X_matrix[this_row,]*bspline[,this_col])$value
      }
      return(temp)
    }
  J_matrix <- do.call(rbind, J_matrix)

  #constant <- sqrt(1/number_basis)*rep(1,number_basis) #Q2_1
  range_time_interval <- max(time_interval)-min(time_interval) #what is full scale
  constant <- rep(1,number_basis)/range_time_interval # unscaled
  lin <- seq(min(time_interval),max(time_interval),length.out=number_basis)/range_time_interval #unscaled
  line <- sqrt(c(1/t(lin)%*%lin))*lin #Q2_2

  D <- diag(ncol(J_matrix))
  if(test_type=='Inclusion'){
    difference_penalty=0
    #return(list(Zmat=(ximat %*% J_matrix), X.g2=NULL,J_matrix=J_matrix,D=D))
    return(list(Zmat=( J_matrix), X.g2=NULL,J=J_matrix,D=D))
  }
  # if(test_type=='Linearity'){ #Test for linearity
  #   d=2
  #   Q2=as.matrix(cbind(constant,lin))
  # }
  if(test_type=='Functional'){ #Test for functional form
    difference_penalty=1
    Q2=as.matrix(constant)
  }
  D <- diff(D,differences=difference_penalty)
  P <-  t(D)%*%D #penalty matrix
  P2 <- 1/2*(P+t(P))
  P.eigen <- eigen(P2)
  evalues <- P.eigen$values[1:nrow(D)]
  Q <- P.eigen$vectors
  Lambda1.inv.half <- diag(sqrt(1/evalues))
  #Q2 <- Q[,(number_basis-difference_penalty+1):number_basis]
  #
  # Ztilde <-  ximat %*% J_matrix %*% Q[,1:(number_basis-d)]
  # X.g2 <-  ximat %*% J_matrix %*% Q2
  Ztilde <-  J_matrix %*% Q[,1:(number_basis-difference_penalty)]
  X.g2 <-  J_matrix %*% Q2
  
  

  return(list(Zmat=Ztilde%*%Lambda1.inv.half, X.g2=X.g2,
       J=J_matrix,D=D,phi=bspline,Q=Q,Q2=Q2,Lambda1.inv.half=Lambda1.inv.half))
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

fit.glmmPQL<-function(test.mat,family,n,test.type,ku=30){
  # can't automate the fixed effects :(
  family.glmm=family
  if(family=='bernoulli'){
    family.glmm='binomial'
  }
  if(family.glmm=='binomial'){ # binomial, needs success and failures
    test.mat$prop<-test.mat$Y/n
    test.mat$n<-n
    if(test.type=='Inclusion'){
      Z.test.names<-c("0",paste0('Z.test',1:ku))
      Z.test.formula<-as.formula(paste("~",paste(Z.test.names,collapse="+")))
      glmm.fit<-try(glmmPQL.mod(prop~0+X1,
                                random=list(ones=pdIdent(Z.test.formula)),family=family.glmm,
                                data=test.mat,weights=n),silent=T)
    }
    if(test.type=='Functional'){
      Z.test.names<-c("0",paste0('Z.test',1:(ku-1)))
      Z.test.formula<-as.formula(paste("~",paste(Z.test.names,collapse="+")))
      glmm.fit<-try(glmmPQL.mod(prop~0+X1+X2,
                                random=list(ones=pdIdent(Z.test.formula)),family=family.glmm,
                                data=test.mat,weights=n),silent=T)
    }
    # if(test.type=='Linearity'){
    #   Z.test.names<-c("0",paste0('Z.test',1:(ku-2)))
    #   Z.test.formula<-as.formula(paste("~",paste(Z.test.names,collapse="+")))
    #   glmm.fit<-try(glmmPQL.mod(prop~0+X1+X2+X3,
    #                             random=list(ones=pdIdent(Z.test.formula)),family=family.glmm,
    #                             data=test.mat,weights=n),silent=T)
    # }
    return(glmm.fit)
  } else { # non-binomial (poisson or bernoulli)
    if(test.type=='Inclusion'){
      Z.test.names<-c("0",paste0('Z.test',1:ku))
      Z.test.formula<-as.formula(paste("~",paste(Z.test.names,collapse="+")))
      glmm.fit<-try(glmmPQL.mod(Y~0+X1,
                                random=list(ones=pdIdent(Z.test.formula)),family=family.glmm,
                                data=test.mat),silent=T)
    }
    if(test.type=='Functional'){
      Z.test.names<-c("0",paste0('Z.test',1:(ku-1)))
      Z.test.formula<-as.formula(paste("~",paste(Z.test.names,collapse="+")))
      glmm.fit<-try(glmmPQL.mod(Y~0+X1+X2,
                                random=list(ones=pdIdent(Z.test.formula)),family=family.glmm,
                                data=test.mat),silent=T)
    }
    if(test.type=='Linearity'){
      Z.test.names<-c("0",paste0('Z.test',1:(ku-2)))
      Z.test.formula<-as.formula(paste("~",paste(Z.test.names,collapse="+")))
      glmm.fit<-try(glmmPQL.mod(Y~0+X1+X2+X3,
                                random=list(ones=pdIdent(Z.test.formula)),family=family.glmm,
                                data=test.mat),silent=T)
    }
    return(glmm.fit)
  }
}

#
# if(run_parallel)
# {
#   parallel::stopCluster(cl = my.cluster)
#   initialized_parallel <- FALSE
# }
