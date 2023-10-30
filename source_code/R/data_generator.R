
library(fda.usc)

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

#n number of subjects
#timeseries_length
#sparse=1 yes   sparse=0 no
#scorevar=2 bigger var , scorevar=1 smaller var
#ps=1 find z1,z2,z3, find p1,p2,p3, logp1-logp3
#ps=2 find z1,z2,z3, find p1,p2,p3=1-p1-p2 logp1-logp3
#ps=3 find z1,z2 staicu find z1hat z2hat
#ps=4 find z1,z2 staicu find z1hat z2hat but only use numerator
#k  #number of eigen functions
#q  #level of the categorical level

generate_cfd_test <- function(klen, num_indvs, timeseries_length, time_interval, fl){
    seed=123
    start_time=time_interval[1]
    end_time=tail(time_interval,1)
    #k=3  #number of eigen functions
    q_val=3  #level of the categorical level

    mu_1=function(t){ -0.64+4*t }
    mu_2=function(t){ 0.97+6*t^2 }

    p_ihat=function(Z_i1app,Z_i2app){
        denom=(1+exp(Z_i1app)+exp(Z_i2app))
        p_i1h=exp(Z_i1app)/denom
        p_i2h=exp(Z_i2app)/denom
        p_i3h=1- p_i1h- p_i2h
        return(list("p_1hatmatrix"=p_i1h,"p_2hatmatrix"=p_i2h,"p_3hatmatrix"=p_i3h))
    }

    mu_vec=rep(0,klen)


    psi_fn=function(klen){

        psi_k1=matrix(rep(1,length(t)*klen),ncol=klen)
        psi_k2=matrix(rep(1,length(t)*klen),ncol=klen)
        for (i in 1:klen) {
            psi_k1[,i]=sin(2*i*pi*t )
            psi_k2[,i]=cos(2*i*pi*t )
        }
        list("psi_k1"=psi_k1,"psi_k2"=psi_k2)
    }


    t=seq(from = start_time,to = end_time, length=timeseries_length)

    X_i=array(0,dim=c(q_val,timeseries_length,num_indvs))  #multinormial results: row is level q_val, column is time points, n is the number of subjects, each column only has one row of 1 and every other rows are 0
    X_nt=matrix(rep(1,num_indvs*length(t)),nrow=num_indvs,ncol=length(t))  #true observations of categorical-valued outcome, each row represent one subject, columns represent time points
    score_matrix=matrix(rep(1,num_indvs*klen),nrow=num_indvs,ncol=klen)  #row is number of subjects, column is the number of eigen functions
    psi_score_matrix_1=matrix(rep(1,num_indvs*length(t)),ncol=num_indvs)  #dim: length(t)*nsubjects
    psi_score_matrix_2=matrix(rep(1,num_indvs*length(t)),ncol=num_indvs)
    Z_i1=matrix(rep(1,num_indvs*length(t)),nrow=num_indvs)  #True latent curves1:row is n subjects, col is t time points
    Z_i2=matrix(rep(1,num_indvs*length(t)),nrow=num_indvs) #True latent curve 2
    p_i1=matrix(rep(0,num_indvs*length(t)),nrow=num_indvs)  #True p_i1
    p_i2=matrix(rep(0,num_indvs*length(t)),nrow=num_indvs)  #True p_i2
    p_i3=matrix(rep(0,num_indvs*length(t)),nrow=num_indvs)  #True p_i3
    for (i in 1:num_indvs){
        set.seed(seed+i)

        if (klen==3){

            #score varies based on i
            score_1=score(0,1)
            score_2=score(0,sqrt(1/2))
            score_3=score(0,1/2)

            score_vector=cbind(score_1,score_2,score_3)
        }

        if (klen==4){

            score_1=score(0,1)
            score_2=score(0,sqrt(1/2))
            score_3=score(0,1/2)
            score_4=score(0,sqrt(1/8))
            # cpve=cumsum(c(1,sqrt(1/2),1/2,sqrt(1/8)))/sum(c(1,sqrt(1/2),1/2,sqrt(1/8)))
            # cvar=c(1,sqrt(1/2),1/2,sqrt(1/8))

            score_vector=cbind(score_1,score_2,score_3,score_4)

        }

        psi_k = psi_fn(klen)

        #Z varies based on i
        #psi t*klen, score: t*klen,  psi%*%t(score)
        psi_score_matrix_1[,i]=psi_k$psi_k1%*%t(score_vector)
        Z_i1[i,]=mu_1(t)+psi_score_matrix_1[,i]

        psi_score_matrix_2[,i]=psi_k$psi_k2%*%t(score_vector)
        Z_i2[i,]=mu_2(t)+psi_score_matrix_2[,i]


        #p varies based on i
        denominator=(1+exp(as.vector(Z_i1[i,]))+exp(as.vector(Z_i2[i,])))
        p_i1[i,]=(exp(as.vector(Z_i1[i,])))/denominator
        p_i2[i,]=(exp(as.vector(Z_i2[i,])))/denominator
        p_i3[i,]=1-p_i1[i,]-p_i2[i,]


        #X_i varies based on i
        #X_i=matrix(rep(1,klen*length(t)),nrow=klen,ncol=length(t))

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
    }# end for i in 1:num indvs

    #collect value and graph
    #collect first two rows of observed binary curves
    X_i1=t(X_i[1,,])  #all num_indvs row subjects , t columns values related to p1
    X_i2=t(X_i[2,,]) #all num_indvs row subjects , t columns values related to p2
    X_i3=t(X_i[3,,]) #all num_indvs row subjects , t columns values related to p3

    #generate Fl functions
    #Formof delta0(t ): (a) Scalar: , (b) Linear: 1+delta1t , (c) Trigonometric: 1+ t +delta2 cos(2pit )
    ##################################################
    if(fl==1){
        fl1=rep(0.45,timeseries_length)
        fl2=rep(0.5,timeseries_length)
        #fl2=matrix(fl2f(t),nrow=timeseries_length,ncol=1)
        fl3=rep(-0.51,timeseries_length)
    }

    if (fl==2){
        # fl1=rep(-0.5,timeseries_length)
        # fl2=matrix(fl2f(t)+0.1,nrow=timeseries_length,ncol=1)
        # #fl2=rep(-0.3,timeseries_length)
        # fl3=matrix(fl3f(t)-0.02,nrow=timeseries_length,ncol=1)


        fl1=rep(-0.1,timeseries_length)
        fl2=matrix(-0.1*fl2f(t),nrow=timeseries_length,ncol=1)
        #fl2=rep(-0.3,timeseries_length)
        fl3=matrix(fl3f(t),nrow=timeseries_length,ncol=1)

    }

    ##############################
    #
    if (fl==3){
        fl1=rep(-0.2,timeseries_length)
        fl2=matrix(-0.15*fl2f(t),nrow=timeseries_length,ncol=1)
        #fl2=rep(-0.3,timeseries_length)
        fl3=matrix(fl3f(t),nrow=timeseries_length,ncol=1)
    }

    ####################
    #
    if (fl==4){
        fl3=matrix(fl3fn(t),nrow=timeseries_length,ncol=1)
        fl1=fl3-0.09
        fl2=fl3+1.3145
    }

    ####################

    flfn=list("fl1"=fl1,"fl2"=fl2,"fl3"=fl3)

    ###############################
    #generate pi
    vec=matrix(1:num_indvs,nrow=num_indvs,ncol=1)
    #integral a function on a interval, returns a scalar
    x1fl1=apply(vec,1, function(x) {int.simpson2(t, X_i1[x,]*fl1, equi = TRUE, method = "TRAPZ")})
    x2fl2=apply(vec,1, function(x) {int.simpson2(t, X_i2[x,]*fl2, equi = TRUE, method = "TRAPZ")})
    x3fl3=apply(vec,1, function(x) {int.simpson2(t, X_i3[x,]*fl3, equi = TRUE, method = "TRAPZ")})

    prob_for_indv=c(0)
    for (i in 1:num_indvs){
        #pi=beta0+xqfl1+x2fl2
        prob_for_indv[i]=expit(0.02+sum(x1fl1,x2fl2,x3fl3))

    }

    #############
    #generate Yi
    yis=c(0)
    for (i in 1:num_indvs){
        yis[i]=rbinom(1,1,prob_for_indv[i])
    }


    truelist=list("TrueX1"=X_i1,"TrueX2"=X_i2,"TrueX3"=X_i3,"Truecatcurve"=X_nt,"fl"=flfn,"yis"=yis)

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

expit <- function(x){1/(1+exp(-x))}

fl2f = function(t) {
    return(t - 8/9)
}

fl3f = function(t) {
    return(0.5*t)
}

fl3fn = function(t) {
    return(-3*t[i]^2 + 2*t[i] - 0.9)
}

GenerateCategoricalFDTest <- function(klen, num_indvs, timeseries_length,
                                      time_interval, fl_choice){

    mns <- GetMuAndScore_2(3)

    generated_data <- GenerateDataTest(num_indvs = 100,
                                      timeseries_length = 100,
                                      mu_1 = mns$mu_1,
                                      mu_2 = mns$mu_2,
                                      score_vals = mns$score_vals,
                                      start_time = time_interval[1],
                                      end_time = tail(time_interval,1),
                                      k = 3)

    prob_curves <- list(p1 = generated_data$p1, p2 = generated_data$p2, p3 = generated_data$p3)
    cat_data <- GenerateCategFuncData(prob_curves)

    flfn <- switch(fl_choice,

                   "1"=list("fl1"=rep(0.45,timeseries_length),
                            "fl2"=rep(0.5,timeseries_length),
                            "fl3"=rep(-0.51,timeseries_length)),

                   "2"=list("fl1"=rep(-0.1,timeseries_length),
                            "fl2"=matrix(-0.1*fl2f(t),nrow=timeseries_length,ncol=1),
                            "fl3"=matrix(fl3f(t),nrow=timeseries_length,ncol=1)),

                   "3"=list("fl1"=rep(-0.2,timeseries_length),
                            "fl2"=matrix(-0.15*fl2f(t),nrow=timeseries_length,ncol=1),
                            "fl3"=matrix(fl3f(t),nrow=timeseries_length,ncol=1)),

                   "4"=list("fl1"=matrix(fl3fn(t),nrow=timeseries_length,ncol=1)-0.09,
                            "fl2"=matrix(fl3fn(t),nrow=timeseries_length,ncol=1)+1.3145,
                            "fl3"=matrix(fl3fn(t),nrow=timeseries_length,ncol=1)),
                   )

    vec <- matrix(1:num_indvs, nrow=num_indvs, ncol=1)

    # integral a function on a interval, returns a scalar
    x1fl1 <- parApply(my.cluster, vec,1, function(x) {int.simpson2(t, X_i1[x,]*flfn$fl1, equi = TRUE, method = "TRAPZ")})
    x2fl2 <- parApply(my.cluster, vec,1, function(x) {int.simpson2(t, X_i2[x,]*flfn$fl2, equi = TRUE, method = "TRAPZ")})
    x3fl3 <- parApply(my.cluster, vec,1, function(x) {int.simpson2(t, X_i3[x,]*flfn$fl3, equi = TRUE, method = "TRAPZ")})

    sum_int_xtft <- matrix(x1fl1 + x2fl2+ x3fl3 + 0.02)
    Y_indvs <- parApply(my.cluster, sum_int_xtft, 1, function(x){ rbinom(1,1, expit(x)) })


    truelist=list("TrueX1"=cat_data$X_array[,,1],
                  "TrueX2"=cat_data$X_array[,,2],
                  "TrueX3"=cat_data$X_array[,,3],
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

start_time <- 0.01
end_time <- 0.99
timeseries_length <- 250
timestamps01 <- seq(from = start_time, to = end_time, length=timeseries_length)

cfdt <- GenerateCategoricalFDTest(klen=3,
                                  num_indvs=100,
                                  timeseries_length = timeseries_length,
                                  time_interval = timestamps01,
                                  fl_choice=2)

result <- cfd_hypothesis_test(cfdt$yis,
                            cfdt$Truecatcurve,
                            time_interval = timestamps01,
                            response_family='bernoulli',
                            test_type='Functional')

if(run_parallel)
{
  parallel::stopCluster(cl = my.cluster)
  initialized_parallel <- FALSE
}
