#cfda n100 t500
#rejection rate
library(mgcv)
library(fda)
library(fda.usc)
library(devtools)
install_github("stchen3/glmmVCtest")
library("glmmVCtest")
library(RLRsim)
#######
options(warn = -1) 
#######
#Function to generate data
logit <- function(x){
  return(log(x/(1-x)))
}

expit <- function(x){1/(1+exp(-x))}

#####
#Function: wrapper function for gam() which outputs the fitted values
#
#Inputs: 
# z : index z = 1,...,N 
# Curves : N x D matrix of observed binary series
# tt : grid of timepoints going from 0 to 1 with D observations
# k : number of basis functions
# method: method used to evaluate the gam
#
#Output: 
# Fitted values from the game function for subject z 
#
#####
#get smoothed curves
regression_g = function(z, Curves, tt, k=10, method="ML"){   #changed from 10 to 25
  z1 = Curves[z,]
  gam1 <- gam(z1~s(tt, bs = "cr", m=2, k = k),
              family="binomial", method = method,control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),optimizer=c("outer","bfgs"))
  return(gam1$fitted.values)
}


Z_ihat=function(Curves_train,tt){
  N_train=dim(Curves_train)[1]
  vec = matrix(1:(N_train), ncol = 1)
  smoothed_x = logit(t(apply(vec, 1, function(x) regression_g(x, Curves_train, tt))))
  smoothed_x
}



score=function(mu,sd){
  rnorm(1,mean=mu,sd=sd)
}

#> x1 = rnorm(1000)           # some continuous variables 
# x2 = rnorm(1000)
#z = 1 + 2*x1 + 3*x2

fl2f = function(t) {
  fl2=c(0)
  for (i in 1:length(t)){
    #x1 = rnorm(1)
    #x2 = rnorm(1)
    # if (t[i] < 0.5) {
    #fl2[i]=0.5-8/9
    #}
    #else {
    #fl2[i]= t[i]-8/9
    #}
    #fl2[i]= 1+5*t[i]
    #fl2[i]= 1 + 2*x1 + 3*x2
    fl2[i]= t[i]-8/9
  }
  
  fl2
}

fl3f = function(t) {
  fl3=c(0)
  for (i in 1:length(t)){
    # if (t[i] < 0.5) {
    #fl2[i]=0.5-8/9
    #}
    #else {
    #fl2[i]= t[i]-8/9
    #}
    #x1 = rnorm(1)
    #x2 = rnorm(1)
    #fl3[i]= 1+t[i]-1.6**cos(2*pi*t[i])
    #fl3[i]= 1 - 3*x1 + 0.5*x2
    fl3[i]=0.5*t[i]
  }
  
  fl3
}

fl3fn = function(t) {
  fl3=c(0)
  for (i in 1:length(t)){
    # if (t[i] < 0.5) {
    #fl2[i]=0.5-8/9
    #}
    #else {
    #fl2[i]= t[i]-8/9
    #}
    #x1 = rnorm(1)
    #x2 = rnorm(1)
    #fl3[i]= 1+t[i]-1.6**cos(2*pi*t[i])
    #fl3[i]= 1 - 3*x1 + 0.5*x2
    fl3[i]=-3*t[i]^2+2*t[i]-0.9
  }
  
  fl3
}
#n number of subjects
#datapoints
#sparse=1 yes   sparse=0 no
#scorevar=2 bigger var , scorevar=1 smaller var
#ps=1 find z1,z2,z3, find p1,p2,p3, logp1-logp3
#ps=2 find z1,z2,z3, find p1,p2,p3=1-p1-p2 logp1-logp3
#ps=3 find z1,z2 staicu find z1hat z2hat
#ps=4 find z1,z2 staicu find z1hat z2hat but only use numerator
#k  #number of eigen functions
#q  #level of the categorical level

generate_cfd_test=function(k,n,datapoints,sparse,scorevar,ps,seed=123,st,et,fl){
  k=k
  seed=seed
  st=st
  et=et
  scorevar=scorevar
  #k=3  #number of eigen functions
  q=3  #level of the categorical level
  
  if(sparse==1){
    mu_1=function(t){
      3.8+4*t  #-0.64+4*t    #4.5+4*t   #sparse 3.8+4*t
      
    }
    mu_2=function(t){
      1.5+4*t^2    #0.97+6*t^2
      
    }
    
    
    p_ihat=function(Z_i1app,Z_i2app){
      denom=(1+exp(Z_i1app)+exp(Z_i2app))
      p_i1h=exp(Z_i1app)/denom
      p_i2h=exp(Z_i2app)/denom
      p_i3h=1/denom
      # p_i3h=1-p_i1h- p_i2h
      return(list("p_1hatmatrix"=p_i1h,"p_2hatmatrix"=p_i2h,"p_3hatmatrix"=p_i3h))
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
    
    p_ihat=function(Z_i1app,Z_i2app){
      denom=(1+exp(Z_i1app)+exp(Z_i2app))
      p_i1h=exp(Z_i1app)/denom
      p_i2h=exp(Z_i2app)/denom
      p_i3h=1- p_i1h- p_i2h
      return(list("p_1hatmatrix"=p_i1h,"p_2hatmatrix"=p_i2h,"p_3hatmatrix"=p_i3h))
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
  for (i in 1:n){
    #pi=beta0+xqfl1+x2fl2
    pis[i]=expit(0.02+sum(x1fl1,x2fl2,x3fl3))
    
  }
  
  #############
  #generate Yi
  yis=c(0)
  for (i in 1:n){
    yis[i]=rbinom(1,1,pis[i])
  }
  
  
  truelist=list("TrueX1"=X_i1,"TrueX2"=X_i2,"TrueX3"=X_i3,"Truecatcurve"=X_nt,"fl"=flfn,"yis"=yis)
  
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
##############
#https://rdrr.io/cran/fdasrvf/src/R/utils.R
#derivative of a function, returns to a curve
gradientfn=function (f, binsize) 
{
  n = nrow(f)
  if (is.null(n)) {
    f = as.vector(f)
    n = length(f)
    g = rep(0, n)
    h = binsize * (1:n)
    g[1] = (f[2] - f[1])/(h[2] - h[1])
    g[n] = (f[n] - f[(n - 1)])/(h[length(h)] - h[(length(h) - 
                                                    1)])
    h = h[3:n] - h[1:(n - 2)]
    g[2:(n - 1)] = (f[3:n] - f[1:(n - 2)])/h[1]
  }
  else {
    f = as.matrix(f)
    p = ncol(f)
    g = matrix(0, n, p)
    h = binsize * (1:n)
    g[1, ] = (f[2, ] - f[1, ])/(h[2] - h[1])
    g[n, ] = (f[n, ] - f[(n - 1), ])/(h[length(h)] - h[(length(h) - 
                                                          1)])
    h = h[3:n] - h[1:(n - 2)]
    g[2:(n - 1), ] = (f[3:n, ] - f[1:(n - 2), ])/h[rep(1, 
                                                       p)]
  }
  return(g)
}

#returns to a curve
ndims <- function(x){
  return(length(dim(x)))
}

cumtrapzfn <- function(x,y,dims=1){
  if ((dims-1)>0){
    perm = c(dims:max(ndims(y),dims), 1:(dims-1))
  } else {
    perm = c(dims:max(ndims(y),dims))
  }
  
  if (ndims(y) == 0){
    n = 1
    m = length(y)
  } else {
    if (length(x) != dim(y)[dims])
      stop('Dimension Mismatch')
    y = aperm(y, perm)
    m = nrow(y)
    n = ncol(y)
  }
  
  if (n==1){
    dt = diff(x)/2.0
    z = c(0, cumsum(dt*(y[1:(m-1)] + y[2:m])))
    dim(z) = c(m,1)
  } else {
    tmp = diff(x)
    dim(tmp) = c(m-1,1)
    dt = repmat(tmp/2.0,1,n)
    z = rbind(rep(0,n), apply(dt*(y[1:(m-1),] + y[2:m,]),2,cumsum))
    perm2 = rep(0, length(perm))
    perm2[perm] = 1:length(perm)
    z = aperm(z, perm2)
  }
  
  return(z)
}
#gamhat=cumtrapz(seq(0, 1, length.out = M), psihat * psihat)
#returns to a number
trapzfnum <- function(x,y,dims=1){
  if ((dims-1)>0){
    perm = c(dims:max(ndims(y),dims), 1:(dims-1))
  } else {
    perm = c(dims:max(ndims(y),dims))
  }
  
  if (ndims(y) == 0){
    m = 1
  } else {
    if (length(x) != dim(y)[dims])
      stop('Dimension Mismatch')
    y = aperm(y, perm)
    m = nrow(y)
  }
  
  if (m==1){
    M = length(y)
    out = sum(diff(x)*(y[-M]+y[-1])/2)
  } else {
    slice1 = y[as.vector(outer(1:(m-1), dim(y)[1]*( 1:prod(dim(y)[-1])-1 ), '+')) ]
    dim(slice1) = c(m-1, length(slice1)/(m-1))
    slice2 = y[as.vector(outer(2:m, dim(y)[1]*( 1:prod(dim(y)[-1])-1 ), '+'))]
    dim(slice2) = c(m-1, length(slice2)/(m-1))
    out = t(diff(x)) %*% (slice1+slice2)/2.
    siz = dim(y)
    siz[1] = 1
    out = array(out, siz)
    perm2 = rep(0, length(perm))
    perm2[perm] = 1:length(perm)
    out = aperm(out, perm2)
    ind = which(dim(out) != 1)
    out = array(out, dim(out)[ind])
  }
  
  return(out)
}
#trapz(time, (tmp - q0[, i])^2)
##sintegral integral of penalty
sintegral=function (x, fx, n.pts = max(256, length(x))) 
{
  if (class(fx) == "function") 
    fx = fx(x)
  n.x = length(x)
  if (n.x != length(fx)) 
    stop("Unequal input vector lengths")
  ap = approx(x, fx, n = 2 * n.pts + 1)
  h = diff(ap$x)[1]
  integral = h * (ap$y[2 * (1:n.pts) - 1] + 4 * ap$y[2 * (1:n.pts)] + 
                    ap$y[2 * (1:n.pts) + 1])/3
  results = list(value = sum(integral), cdf = list(x = ap$x[2 * 
                                                              (1:n.pts)], y = cumsum(integral)))
  class(results) = "sintegral"
  return(results)
}
#############
#Function to test
#use construct bassi functions with equal knots
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


#gen.Zs<-function(efuncs,times,ximat,ku=30,test){
gen.Zscfd<-function(Xil,times,ku=30,test){
  knots<-construct.knots(times,knots=(ku-3),knots.option='equally-spaced')
  bspline<-splineDesign(knots=knots,x=times,ord=4)
  J<-matrix(0,nrow=nrow(Xil),ncol=ku) #empty
  for(row in 1:nrow(J)){
    for(col in 1:ncol(J)){
      J[row,col]<-sintegral(times,Xil[row,]*bspline[,col])$value
    }
  }
  constant<-sqrt(1/ku)*rep(1,ku) #Q2_1
  range<-max(times)-min(times) #what is full scale
  constant<-rep(1,ku)/range # unscaled
  lin<-seq(min(times),max(times),length.out=ku)/range #unscaled
  line<-sqrt(c(1/t(lin)%*%lin))*lin #Q2_2
  
  D<-diag(ncol(J))
  if(test=='Inclusion'){
    d=0
    #return(list(Zmat=(ximat %*% J), X.g2=NULL,J=J,D=D))
    return(list(Zmat=( J), X.g2=NULL,J=J,D=D))
  }
  # if(test=='Linearity'){ #Test for linearity
  #   d=2
  #   Q2=as.matrix(cbind(constant,lin))
  # }
  if(test=='Functional'){ #Test for functional form
    d=1
    Q2=as.matrix(constant)
  }
  D<-diff(D,differences=d)
  P<- t(D)%*%D #penalty matrix
  P2<-1/2*(P+t(P))
  P.eigen<-eigen(P2)
  evalues<-P.eigen$values[1:nrow(D)]
  Q<-P.eigen$vectors
  Lambda1.inv.half<-diag(sqrt(1/evalues))
  #Q2<-Q[,(ku-d+1):ku]
  # 
  # Ztilde<- ximat %*% J %*% Q[,1:(ku-d)]
  # X.g2<- ximat %*% J %*% Q2
  Ztilde<- J %*% Q[,1:(ku-d)]
  X.g2<- J %*% Q2
  list(Zmat=Ztilde%*%Lambda1.inv.half, X.g2=X.g2,
       J=J,D=D,phi=bspline,Q=Q,Q2=Q2,Lambda1.inv.half=Lambda1.inv.half)
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
    if(test.type=='Linearity'){
      Z.test.names<-c("0",paste0('Z.test',1:(ku-2)))
      Z.test.formula<-as.formula(paste("~",paste(Z.test.names,collapse="+")))
      glmm.fit<-try(glmmPQL.mod(prop~0+X1+X2+X3,
                                random=list(ones=pdIdent(Z.test.formula)),family=family.glmm,
                                data=test.mat,weights=n),silent=T)
    }
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
#####

#################
#generate data
# hyp4=generate_cfd_test(3,nsub,datapoints,0,1,1,seed=123,0.01,0.99,4) 
# cfd2=hyp4$true$TrueX2
# cfd3=hyp4$true$TrueX3
# 
# Xmat_Inc<-matrix(rep(1,nsub),ncol=1)
# Zmat_Func<-gen.Zscfd(cfd2,seq(0.01,0.99,length=numt),ku=30,test='Functional')
# Zmat_Func3<-gen.Zscfd(cfd3,seq(0.01,0.99,length=numt),ku=30,test='Functional')
# Zmat_Func.mat <- Zmat_Func$Zmat
# Zmat_Func.mat3 <- Zmat_Func3$Zmat
# Xmat_Func<-cbind(Xmat_Inc,Zmat_Func$X.g2,Zmat_Func3$X.g2)
# 
# test.mat<-data.frame(Y=hyp4$true$yis,X=Xmat_Func,Z.test=Zmat_Func.mat,Z.test3=Zmat_Func.mat3,ones=rep(1,n4))
# names(test.mat)<-c('Y','X1','X2',"X3",paste0('Z.test',1:ncol(Zmat_Func.mat)),paste0('Z.test3',1:ncol(Zmat_Func.mat3)),"ones")
# 
# 
# #Proposed
# pro.alt.fit <-fit.glmmPQL(test.mat,'bernoulli',1,'Functional') 
# pro.test.Func <- try(test.aRLRT(pro.alt.fit),silent=T) # Functional only
# datatest=c(pro.test.Func$aRLRT$statistic,pro.test.Func$aRLRT$p.value)
# names(datatest)=c("RLRT","p-value")
# datatest


#(k,n,datapoints,sparse,scorevar,ps,seed=123,st,et,fl){

sim_cfdatest=function(k,nsub,datapoints,sparse,scorevar,ps,seedlength,st,et,fl){
  datatest=matrix(0,ncol=2,nrow=seedlength)
  for (i in 1:seedlength){
    seed=123+i
    hyp4=generate_cfd_test(k,nsub,datapoints,sparse,scorevar,ps,seed=seed,st,et,fl) 
    cfd2=hyp4$true$TrueX2
    cfd3=hyp4$true$TrueX3
    numt=datapoints
    Xmat_Inc<-matrix(rep(1,nsub),ncol=1)
    Zmat_Func<-gen.Zscfd(cfd2,seq(st,et,length=numt),ku=30,test='Functional')
    Zmat_Func3<-gen.Zscfd(cfd3,seq(st,et,length=numt),ku=30,test='Functional')
    Zmat_Func.mat <- Zmat_Func$Zmat
    Zmat_Func.mat3 <- Zmat_Func3$Zmat
    Xmat_Func<-cbind(Xmat_Inc,Zmat_Func$X.g2,Zmat_Func3$X.g2)
    
    test.mat<-data.frame(Y=hyp4$true$yis,X=Xmat_Func,Z.test=Zmat_Func.mat,Z.test3=Zmat_Func.mat3,ones=rep(1,nsub))
    names(test.mat)<-c('Y','X1','X2',"X3",paste0('Z.test',1:ncol(Zmat_Func.mat)),paste0('Z.test3',1:ncol(Zmat_Func.mat3)),"ones")
    
    
    #Proposed
    pro.alt.fit <-fit.glmmPQL(test.mat,'bernoulli',1,'Functional') 
    pro.test.Func <- try(test.aRLRT(pro.alt.fit),silent=T) # Functional only
    pvalue=pro.test.Func$aRLRT$p.value<0.05
    datatest[i,]=c(pro.test.Func$aRLRT$statistic,pvalue)
    names(datatest)=c("RLRT","p-value")
  }
  
  
  testst=round(mean(datatest[,1]),3)
  pvalue= round(mean(datatest[,2]),3)
  
  final=c(testst,pvalue,round(nsub),round(datapoints))
  names(final)=c("RLRT","p-value","n","datapoints")
  final
}

#n100t250=sim_cfdatest(3,100,250,0,1,1,5000,0.01,0.99,4)
n100t500=sim_cfdatest(3,100,500,0,1,1,5000,0.01,0.99,4)
#n100t2000=sim_cfdatest(3,100,2000,0,1,1,5000,0.01,0.99,4)
#save(n100t2000,file="typeIerror1002000.RData")
save(n100t500,file="typeIerror100500.RData")
#load("typeIerror1002000.RData")
