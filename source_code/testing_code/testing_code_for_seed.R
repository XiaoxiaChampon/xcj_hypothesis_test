#extra testing code
#check
start_time=0.01
end_time=0.99
timeseries_length=180
time_interval=seq(start_time,end_time,length=timeseries_length)
num_indvs=500
#fl_choice=3 #not constant, expect to reject
fl_choice=6 #not constant, expect to reject
response_family='bernoulli'
test_type='Functional'
klen=3
mu1_coef=c(-6.67,-2.47,5.42)
mu2_coef=c(-3.14,-0.99,3.91)

###############################################
MAGIC_NUM_DUM_DUM <<- 0.1 #0.6206897
source("source_code/R/data_generator.R")

set.seed(123456) #working

#GenerateCategoricalFDTest <- function(klen, mu1_coef,mu2_coef,num_indvs, timeseries_length,
#                                      lp_intercept=0.6206897)

test_noconstant=GenerateCategoricalFDTest(klen=3,mu1_coef, mu2_coef,timeseries_length,
                                          time_interval, fl_choice,MAGIC_NUM_DUM_DUM)
test_noconstant$pvalue
linear_predictor_w=test_noconstant$linear_predictor$linearw
linear_predictor_wo=test_noconstant$linear_predictor$linearwo
linear_predictor=c(linear_predictor_w,linear_predictor_wo)
indicator <- c(rep("With",length(linear_predictor_w)), rep("Without", length(linear_predictor_wo)))
Y_indvs=test_noconstant$yis

Y_without = apply(linear_predictor_wo, 1, function(x){ rbinom(1, 1, 1/(1+exp(-x))) })
#table(Y_without)
Y_perc_wo = table(Y_without)/sum(table(Y_without))
Y_perc_wo

#table(Y_indvs)
Y_perc_w = table(Y_indvs)/sum(table(Y_indvs))
Y_perc_w


#####################
set.seed(123456) #working
intercept_values <- seq(-2,2,length.out=100)
result_y <- foreach (intercept_index = 1:length(intercept_values), .combine = cbind, .init = NULL) %dorng% {
  source("source_code/R/data_generator.R")
  test_noconstant=cfd_testing(start_time, end_time, timeseries_length,
                              num_indvs,mu1_coef, mu2_coef,fl_choice,response_family,test_type,
                              klen=3, intercept_values[intercept_index])
  test_noconstant$pvalue
  linear_predictor_w=test_noconstant$linear_predictor$linearw
  linear_predictor_wo=test_noconstant$linear_predictor$linearwo
  linear_predictor=c(linear_predictor_w,linear_predictor_wo)
  indicator <- c(rep("With",length(linear_predictor_w)), rep("Without", length(linear_predictor_wo)))
  Y_indvs=test_noconstant$yis
  
  Y_without = test_noconstant$yis_without
  #table(Y_without)
  Y_perc_wo = table(Y_without)/sum(table(Y_without))
  
  #table(Y_indvs)
  Y_perc_w = table(Y_indvs)/sum(table(Y_indvs))
  return(list("with" = Y_perc_w, "without" = Y_perc_wo ))
}
result_compare = cbind(unlist(result_y[1,]),
                       unlist(result_y[2,]))
result_compare 

result_compare[,1]-0.5<0.2

intercept_values[20]



#old use intercept -1
#Y_indvs
#0   1 
#383 117 
#Y_indvs
#0     1 
#0.766 0.234 
prob_ind=test_noconstant$prob_ind
data_sample=data.frame(linear_predictor,as.factor(indicator),prob_ind)
colnames(data_sample)=c("linear_predictor","Group","Probability")

save(data_sample,file="data_sample_lineargroup_pie_beta0is1.RData")

save(data_sample,file="data_sample_lineargroup_pie.RData")
load("data_sample_lineargroup_pie.RData")
#odd beta0=-1
data_sample_non_linear =matrix(data_sample[data_sample$Group=="Without",]$linear_predictor,ncol=1)
Y_without_old = apply(data_sample_non_linear, 1, function(x){ rbinom(1, 1, 1/(1+exp(-x))) })
table(Y_without_old)
table(Y_without_old)/sum(table(Y_without_old))


library(ggplot2)
linear_group=ggplot( data_sample, aes(x = linear_predictor, fill=Group, colour = Group)) +
    geom_histogram(position = "dodge")+
    #theme(axis.title.x=element_blank())+
    theme(text = element_text(size = 20))+
    labs(x = "Linear Predictor")+
    theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
    )
linear_group
ggsave("linear_group.png")

linear_group_square
ggsave("linear_group_square.png")

linear_group_square_1=ggplot( data_sample, aes(x = linear_predictor, fill=Group, colour = Group)) +
  geom_histogram(position = "dodge")+
  #theme(axis.title.x=element_blank())+
  theme(text = element_text(size = 20))+
  labs(x = "Linear Predictor")+
  theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  )
linear_group_square_1
ggsave("linear_group_square1.png")
cfd=test_noconstant$W

##pie chart for   W
pie_data_cc=table(cfd)
cc_pie_zero=as.data.frame(pie_data_cc)
cc_pie_zero$pct=round(cc_pie_zero$Freq/(sum(cc_pie_zero$Freq)),2)
cc_pie_zero$labels= scales::percent(cc_pie_zero$pct)
colnames(cc_pie_zero)[colnames(cc_pie_zero) == "cfd"] <- "Category"

pie_w=ggplot(cc_pie_zero, aes(x = "", y = pct, fill = Category)) +
  geom_col() +
  geom_text(aes(label = labels),size=8,
            position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y")+ggtitle("")+theme(plot.title = element_text(hjust = 0.5))+
  theme(text=element_text(size = 20))
pie_w
ggsave("pie_w.png")


pie_w1=ggplot(cc_pie_zero, aes(x = "", y = pct, fill = Category)) +
  geom_col() +
  geom_text(aes(label = labels),size=8,
            position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y")+ggtitle("")+theme(plot.title = element_text(hjust = 0.5))+
  theme(text=element_text(size = 20))
pie_w1
ggsave("pie_w1.png")
library(gridExtra)
#grid.arrange(linear_group, pie_w, ncol = 2)
grid.arrange(linear_group_square_1, pie_w1, ncol = 2)
################################
set.seed(123) #notworking
test_noconstant_not=cfd_testing(start_time, end_time, timeseries_length,
                                num_indvs,fl_choice,response_family,test_type,
                                klen=3)
cfd_not=test_noconstant_not$W
test_noconstant_not$pvalue

linear_predictor_not=test_noconstant_not$linear_predictor
Y_indvs_not=test_noconstant_not$yis
table(Y_indvs_not)
prob_ind_not=test_noconstant_not$prob_ind
data_sample_not=data.frame(linear_predictor_not,as.factor(Y_indvs_not),prob_ind_not)
colnames(data_sample_not)=c("linear_predictor","Group","Probability")
linear_group2 = ggplot( data_sample_not, aes(x = linear_predictor, fill=Group, colour = Group)) +
    geom_histogram(position = "dodge")+
    #theme(axis.title.x=element_blank())+
    theme(text = element_text(size = 20))+
    labs(x = "Linear Predictor 2")

library(gridExtra)
grid.arrange(linear_group,linear_group2,ncol=2)


###############
test.aRLRT<-function(fit){
    if(class(fit)=='glmmPQL'){
        stop('Please use glmmPQL.mod to estimate the alternative hypothesis')
    }
    if(class(fit)=='lme'){ # normal responses
        if(length(fit$call$random)==2){ # only 1 random effect
            RLRT <- exactRLRT(fit)
            return(list(aRLRT=RLRT,std.data=fit$data,fit.alt=fit,fit.null=NULL,fit.test=NULL))
        } else { # multiple randome effects
            null.call <- test.call <- fit$call
            null.call$random<- fit$call$random[-2] # remove effect being tested
            fit.null<-eval(null.call)
            test.call$random <- fit$call$random[1:2] # effect being tested only
            fit.test <- eval(test.call)
            RLRT <- exactRLRT(m=fit.test,mA=fit,m0=fit.null)
            return(list(aRLRT=RLRT,std.data=fit$data,fit.alt=fit,fit.null=fit.null,fit.test=fit.test))
        }
    } else if(class(fit)=='list'){ # glmmPQL.mod output
        if(class(fit$fit)[1]=='glmmPQL'){ # generalized responses
            # Extract and standardize to Ytilde
            mcall.orig<-fit$mcall # lme call (X and Z haven't been adjusted to iid)
            mcall.std<-std.glmmPQL(fit) # updated mcall with standardized Ytilde, Xtilde, Ztilde
            std.data<-mcall.std$data
            
            # refit under null and alt hypothesis
            fit.alt.std<-try(eval(mcall.std),silent=T)
            if('try-error' %in% class(fit.alt.std)){
                stop('Error in lme model estimation under alternative. Consider simplifying
             or rescaling variables.')
            }
            
            if(length(mcall.std$random)==1){ # no nuisance random effects
                fit.null.std<-try(lm(mcall.std$fixed,data=std.data)) # fixed effects only
                if('try-error' %in% class(fit.null.std)){
                    stop('Error in lm model estimation under null. Consider simplifying
               or rescaling variables.')
                }
                fit.test.std<-fit.alt.std # same as alternative model
            } else { # nuisance r.effect
                mcall.null<-mcall.test<-mcall.std # update from alt model fit
                mcall.null$random<-mcall.std$random[2:length(mcall.std$random)] #only nuis
                mcall.test$random<-mcall.std$random[1] # only test
                fit.null.std<-try(eval(mcall.null),silent=T)
                if('try-error' %in% class(fit.null.std)){
                    stop('Error in lme model estimation under null. Consider simplifying
                 or rescaling variables.')
                }
                fit.test.std<-try(eval(mcall.test),silent=T)
                if('try-error' %in% class(fit.test.std)){
                    stop('Error in lme model estimation for testing variable. Consider rescaling.')
                }
            }
            # testing
            #For testing in models with multiple variance
            #' components, the fitted model \code{m} must contain \bold{only} the random
            #' effect set to zero under the null hypothesis, while \code{mA} and \code{m0}
            #' are the models under the alternative and the null, respectively. 
            aRLRT<-exactRLRT(m=fit.test.std,mA=fit.alt.std,m0=fit.null.std)
            return(list(aRLRT=aRLRT,std.data=std.data,fit.alt=fit.alt.std,fit.null=fit.null.std,fit.test=fit.test.std))
        } else {
            stop('Only lme and glmmPQL.mod class models are supported!')
        }
    } else {
        stop('Only lme and glmmPQL.mod class models are supported!')
    }
}




#########fl=4
fl_choice=4 # constant, expect to not reject
set.seed(1234) #working
test_constant=cfd_testing(start_time, end_time, timeseries_length,
                            num_indvs,fl_choice,response_family,test_type,
                            klen=3)
test_constant$pvalue
linear_predictor=test_constant$linear_predictor
Y_indvs=test_constant$yis
table(Y_indvs)
prob_ind=test_noconstant$prob_ind
data_sample=data.frame(linear_predictor,as.factor(Y_indvs),prob_ind)
colnames(data_sample)=c("linear_predictor","Group","Probability")
linear_group=ggplot( data_sample, aes(x = linear_predictor, fill=Group, colour = Group)) +
    geom_histogram(position = "dodge")+
    #theme(axis.title.x=element_blank())+
    theme(text = element_text(size = 20))+
    labs(x = "Linear Predictor")
linear_group
cfd=test_noconstant$W

#########fl=3 inclusion
fl_choice=3 # not constant, expect to reject, and it rejects after 5 iterations
test_type='Inclusion'
set.seed(1234) #working
test_constant=cfd_testing(start_time, end_time, timeseries_length,
                          num_indvs,fl_choice,response_family,test_type,
                          klen=3)
test_constant$pvalue
library(ggplot2)
linear_predictor=test_constant$linear_predictor
Y_indvs=test_constant$yis
table(Y_indvs)
prob_ind=test_constant$prob_ind
data_sample=data.frame(linear_predictor,as.factor(Y_indvs),prob_ind)
colnames(data_sample)=c("linear_predictor","Group","Probability")
linear_group=ggplot( data_sample, aes(x = linear_predictor, fill=Group, colour = Group)) +
    geom_histogram(position = "dodge")+
    #theme(axis.title.x=element_blank())+
    theme(text = element_text(size = 20))+
    labs(x = "Linear Predictor")
linear_group
cfd=test_constant$W

#########fl=1 inclusion
fl_choice=1 # 0, expect to not reject
test_type='Inclusion'
set.seed(1234) #working
test_constant=cfd_testing(start_time, end_time, timeseries_length,
                          num_indvs,fl_choice,response_family,test_type,
                          klen=3)
test_constant$pvalue
library(ggplot2)
linear_predictor=test_constant$linear_predictor
Y_indvs=test_constant$yis
table(Y_indvs)
prob_ind=test_constant$prob_ind
data_sample=data.frame(linear_predictor,as.factor(Y_indvs),prob_ind)
colnames(data_sample)=c("linear_predictor","Group","Probability")
linear_group=ggplot( data_sample, aes(x = linear_predictor, fill=Group, colour = Group)) +
    geom_histogram(position = "dodge")+
    #theme(axis.title.x=element_blank())+
    theme(text = element_text(size = 20))+
    labs(x = "Linear Predictor")
linear_group
cfd=test_constant$W

#############after add mu1_coef, and mu2_coef
fl_choice=1 # 0, expect to not reject
test_type='Inclusion'
set.seed(1234) #working
test_constant=cfd_testing(start_time, end_time, timeseries_length,
                          num_indvs,mu1_coef,mu2_coef,fl_choice,response_family,test_type,
                          klen=3)
test_constant$pvalue
library(ggplot2)
linear_predictor=test_constant$linear_predictor
Y_indvs=test_constant$yis
table(Y_indvs)
prob_ind=test_constant$prob_ind
data_sample=data.frame(linear_predictor,as.factor(Y_indvs),prob_ind)
colnames(data_sample)=c("linear_predictor","Group","Probability")
linear_group=ggplot( data_sample, aes(x = linear_predictor, fill=Group, colour = Group)) +
    geom_histogram(position = "dodge")+
    #theme(axis.title.x=element_blank())+
    theme(text = element_text(size = 20))+
    labs(x = "Linear Predictor")
linear_group
cfd=test_constant$W



#####################################
library(foreach)
library(doRNG)
library(doParallel)
#RNGkind("L'Ecuyer-CMRG")
RNGkind("Mersenne-Twister")


set.seed(123)
rn3 <- foreach(i=1:10, .combine = 'c') %dopar%{ 
    rs <- .Random.seed
    return(list("a"=rnorm(3,0,1), "s"=rs))
}

rn1 <- foreach(i=1:10, .combine = 'c', .options.RNG=123) %dorng%{ 
    rs <- .Random.seed
    return(list("a"=rnorm(3,0,1), "s"=rs))
}

set.seed(123)
rn2 <- foreach(i=1:10, .combine = 'c') %dorng%{ 
    rs <- .Random.seed
    return(list("a"=rnorm(3,0,1), "s"=rs))
}


rn4 <- rep(0,10)
set.seed(123)
for(i in 1:20){
    rn4[i] <- (rnorm(1,0,1))
}

identical(rn1, rn2) 
identical(rn1, rn3)
identical(rn1, rn4)
identical(rn3, rn4)

set.seed(123)
rn6 <- rep(0,20)
for(i in 1:20){
    .Random.seed <- attr(rn1,"rng")[[i]] #using seeds from rn1 from question
    rn6[i] <- (rnorm(1,0,1))
}

#######
# linear_predictor=test_noconstant$linear_predictor
# Y_indvs=test_noconstant$yis
#table(Y_indvs)
# prob_ind=test_noconstant$prob_ind
# data_sample=data.frame(linear_predictor,as.factor(Y_indvs),prob_ind)
# colnames(data_sample)=c("linear_predictor","Group","Probability")
# ggplot( data_sample, aes(x = linear_predictor, fill=Group, colour = Group)) +
#   geom_histogram(position = "dodge")+
#   #theme(axis.title.x=element_blank())+
#   theme(text = element_text(size = 20))+
#   labs(x = "Linear Predictor")
# 
# 
# #fl=4, constant, expect to not reject
# fl_choice=4
# cfd_testing(start_time, end_time, timeseries_length,
#             num_indvs,fl_choice,response_family,test_type,
#             klen=3)$pvalue


# linear_predictor=cfd_test_data$true$linear_predictor
# Y_indvs=cfd_test_data$true$yis
# 
# table(Y_indvs)
# prob_ind=cfd_test_data$true$prob_ind
# data_sample=data.frame(linear_predictor,as.factor(Y_indvs),prob_ind)
# colnames(data_sample)=c("linear_predictor","Group","Probability")
# ggplot( data_sample, aes(x = linear_predictor, fill=Group, colour = Group)) +
#   geom_histogram(position = "dodge")+
#   #theme(axis.title.x=element_blank())+
#   theme(text = element_text(size = 20))+
#   labs(x = "Linear Predictor")



# start_time=0.01
# end_time=0.99
# timeseries_length=2000
#time_interval=seq(start_time,end_time,length=timeseries_length)
# num_indvs=100
# fl_choice=3
# response_family='bernoulli'
# test_type='Functional'
# klen=3
#table(cat_data$W)/sum(cat_data$W)
#table(Y_indvs)
#plot( flfn$fl2)
#plot( flfn$fl3)
#####################
# linear_predictor=sum_int_xtft
# prob_ind=1/(1+exp(-sum_int_xtft))
# data_sample=data.frame(linear_predictor,as.factor(Y_indvs),prob_ind)
# library(ggplot2)
# ggplot( data_sample, aes(x = linear_predictor, fill=as.factor(Y_indvs), colour = as.factor(Y_indvs))) +
#   geom_histogram(position = "dodge")
# table(Y_indvs)
##################
# a=50
# b=50
# 
# set.seed(666)
#  x1 = rnorm(a)           # some continuous variables 
#  x2 = rnorm(b)
# linear_predictor=1 + 2*x1 + 3*x2
# prob_ind=1/(1+exp(-linear_predictor))
# Y_indvs=rbinom(a+b,1,prob_ind)
# table(Y_indvs)
# #data_sample=data.frame(linear_predictor,as.factor(Y_indvs),prob_ind)
# data_sample=data.frame(linear_predictor,as.factor(Y_indvs))
# library(ggplot2)
# ggplot( data_sample, aes(x = linear_predictor, fill=as.factor(Y_indvs), colour = as.factor(Y_indvs))) +
#   geom_histogram(position = "dodge")




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

# data_sample=data.frame(linear_predictor,as.factor(Y_indvs),prob_ind)
# colnames(data_sample)=c("linear_predictor","Group","Probability")
# library(ggplot2)
# ggplot( data_sample, aes(x = linear_predictor, fill=Group, colour = Group)) +
#   geom_histogram(position = "dodge")+
#   theme(axis.title.x=element_blank())+
#   theme(text = element_text(size = 20))
# table(Y_indvs)
######

#range(sum_int_xtft)
###

#Y_indvs <- parApply(my.cluster, sum_int_xtft, 1, function(x){ rbinom(1,1, 1/(1+exp(-x))) })




# integral a function on a interval, returns a scalar
# x1fl1 <- parApply(my.cluster, vec, 1, function(x) {fda.usc::int.simpson2(time_interval, cat_data$X[x,,1]*flfn$fl1, equi = TRUE, method = "TRAPZ")})
# x2fl2 <- parApply(my.cluster, vec, 1, function(x) {fda.usc::int.simpson2(time_interval, cat_data$X[x,,2]*flfn$fl2, equi = TRUE, method = "TRAPZ")})
# x3fl3 <- parApply(my.cluster, vec, 1, function(x) {fda.usc::int.simpson2(time_interval, cat_data$X[x,,3]*flfn$fl3, equi = TRUE, method = "TRAPZ")})
# 
# x1fl1 <- apply(vec, 1, function(x) {fda.usc::int.simpson2(time_interval, cat_data$X[x,,1]*flfn$fl1, equi = TRUE, method = "TRAPZ")})
# x2fl2 <- apply( vec, 1, function(x) {fda.usc::int.simpson2(time_interval, cat_data$X[x,,2]*flfn$fl2, equi = TRUE, method = "TRAPZ")})
# x3fl3 <- apply( vec, 1, function(x) {fda.usc::int.simpson2(time_interval, cat_data$X[x,,3]*flfn$fl3, equi = TRUE, method = "TRAPZ")})
# 
#############
# x1fl1 <- apply(vec, 1, function(x) {fda.usc::int.simpson2(time_interval, flfn$fl1, equi = TRUE, method = "TRAPZ")})
# x2fl2 <- apply( vec, 1, function(x) {fda.usc::int.simpson2(time_interval, cat_data$X[x,,2]*10*(flfn$fl2-flfn$fl1), equi = TRUE, method = "TRAPZ")})
# x3fl3 <- apply( vec, 1, function(x) {fda.usc::int.simpson2(time_interval, cat_data$X[x,,3]*10*(flfn$fl3-flfn$fl1), equi = TRUE, method = "TRAPZ")})

# fll2=-10*sin((2*pi/25)*(time_interval-1))
# fll3=-20*sin((2*pi/25)*(time_interval-1))-6


# num_indvs=100
# timeseries_length=2000
# start_time=0.01
# end_time=0.99
# fl_choice=3
# response_family='bernoulli'
# test_type='Functional'
# klen=3
#time_interval=seq(start_time,end_time,length=timeseries_length)


#########################
#test
#######
# X_matrix_not=X_cfd_not[,,2]
# ######
# knots_not <- construct.knots(time_interval,knots=(number_basis-3),knots.option='equally-spaced')
# bspline_not <- splineDesign(knots=knots,x=time_interval,ord=4)
# 
# number_row_not <- nrow(X_matrix_not)
# number_col_not <- number_basis

# J_matrix_not<- matrix(0,nrow=nrow(X_matrix_not),ncol=number_basis) #empty
# for(row in 1:number_row_not){
#   for(col in 1:number_col_not){
#     J_matrix_not[row,col] <- integral_penalty(time_interval,X_matrix_not[row,]*bspline_not[,col])$value
#   }
# }
# J_matrix_not <- foreach(this_row = 1:number_row_not) %do%
#   {
#     source("source_code/R/integral_penalty_function.R")
# 
#     temp <- array(-123, number_col_not)
#     for(this_col in 1:number_col_not){
#       temp[this_col] <- integral_penalty(time_interval,X_matrix_not[this_row,]*bspline_not[,this_col])$value
#     }
#     return(temp)
#   }
# J_matrix_not <- do.call(rbind, J_matrix_not)
# 
# constant_not <- sqrt(1/number_basis)*rep(1,number_basis) #Q2_1
# range_time_interval <- max(time_interval)-min(time_interval) #what is full scale
# constant_not <- rep(1,number_basis)/range_time_interval # unscaled
# lin_not <- seq(min(time_interval),max(time_interval),length.out=number_basis)/range_time_interval #unscaled
# line_not <- sqrt(c(1/t(lin_not)%*%lin_not))*lin_not #Q2_2
# 
# D_not <- diag(ncol(J_matrix_not))
# if(test_type=='Functional'){ #Test for functional form
#   difference_penalty=1
#   Q2_not=as.matrix(constant_not)
# }
# D_not <- diff(D_not,differences=difference_penalty)
# P_not <-  t(D_not)%*%D_not #penalty matrix
# P2_not <- 1/2*(P_not+t(P_not))
# P.eigen_not <- eigen(P2_not)
# evalues_not <- P.eigen_not$values[1:nrow(D_not)]
# Q_not <- P.eigen_not$vectors
# Lambda1.inv.half_not <- diag(sqrt(1/evalues_not))
# Q2_not <- Q_not[,(number_basis-difference_penalty+1):number_basis]
# #
# # Ztilde <-  ximat %*% J_matrix %*% Q[,1:(number_basis-d)]
# # X.g2 <-  ximat %*% J_matrix %*% Q2
# Ztilde_not <-  J_matrix_not %*% Q_not[,1:(number_basis-difference_penalty)]
# X.g2_not <-  J_matrix_not %*% Q2_not
# Zmat_not=Ztilde_not%*%Lambda1.inv.half_not


# alternative_fit_not <- fit.glmmPQL(test_matrix_not, response_family, num_indvs, test_type)
# 
# result_not <- try(test.aRLRT(alternative_fit_not), silent=T)$aRLRT # Functional only
# result_not$p.value


###################
#test both if equal
# X_cfd_not <- GetXFromW(cfd_not)
# Y_not=Y_indvs_not

# Zmat_Func2_not <- get_Zmatrix(X_cfd_not[,,2], time_interval, test_type)
# Zmat_Func3_not <- get_Zmatrix(X_cfd_not[,,3], time_interval, test_type)
# 
# num_indvs_not <- length(Y_not)
# Xmat_Inc_not<-matrix(rep(1, num_indvs_not),ncol=1)
# Xmat_Func_not <- cbind(Xmat_Inc_not, Zmat_Func2_not$X.g2, Zmat_Func3_not$X.g2)
# 
# test_matrix_not <- data.frame(Y=Y_not,
#                           X=Xmat_Func_not,
#                           Z.test=Zmat_Func2_not$Zmat,
#                           Z.test3=Zmat_Func3_not$Zmat,
#                           ones=rep(1,num_indvs_not))
# names(test_matrix_not) <- c('Y','X1','X2',"X3",
#                         paste0('Z.test',1:ncol(Zmat_Func2$Zmat)),
#                         paste0('Z.test3',1:ncol(Zmat_Func3$Zmat)),
#                         "ones")



###################
#old fls that works for rejecting, fl choice=3
# flfn <- switch(fl_choice,
#                
#                "1"=list("fl1"=rep(0.45,timeseries_length),
#                         "fl2"=rep(0.5,timeseries_length),
#                         "fl3"=rep(-0.51,timeseries_length)),
#                
#                "2"=list("fl1"=rep(-0.1,timeseries_length),
#                         "fl2"=matrix(-0.1*fl2f(time_interval),nrow=timeseries_length,ncol=1),
#                         "fl3"=matrix(fl3f(time_interval),nrow=timeseries_length,ncol=1)),
#                #not constant
#                # fll2=-10*sin((2*pi/25)*(time_interval-1))
#                # fll3=-20*sin((2*pi/25)*(time_interval-1))-6
#                "3"=list("fl1"=rep(-0.2,timeseries_length),
#                         # "fl2"=matrix(-0.15*fl2f(time_interval),nrow=timeseries_length,ncol=1),
#                         # "fl3"=matrix(fl3f(time_interval),nrow=timeseries_length,ncol=1)),
#                         "fl2"=matrix(-10*sin((2*pi/25)*(time_interval-1))+10,nrow=timeseries_length,ncol=1),
#                         "fl3"=matrix(-20*sin((2*pi/25)*(time_interval-1))-6,nrow=timeseries_length,ncol=1)),
#                
#                "4"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
#                         "fl2"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)+1.3145,
#                         "fl3"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)),
# )
# 
# vec <- matrix(1:num_indvs, nrow=num_indvs, ncol=1)
# 
# 
# x1fl1 <- apply(vec, 1, function(x) {fda.usc::int.simpson2(time_interval, flfn$fl1-40, equi = TRUE, method = "TRAPZ")})
# x2fl2 <- apply( vec, 1, function(x) {fda.usc::int.simpson2(time_interval, cat_data$X[x,,2]*(5)*(flfn$fl2), equi = TRUE, method = "TRAPZ")})
# x3fl3 <- apply( vec, 1, function(x) {fda.usc::int.simpson2(time_interval, cat_data$X[x,,3]*(-4)*(flfn$fl3), equi = TRUE, method = "TRAPZ")})






GenerateCategoricalFDTest <- function(klen, mu1_coef,mu2_coef,num_indvs, timeseries_length,
                                      time_interval, fl_choice){
    
    #mns <- GetMuAndScore_2(klen)
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
                   # 
                   #                    "1"=list("fl1"=rep(0.45,timeseries_length),
                   #                             "fl2"=rep(0.5,timeseries_length),
                   #                             "fl3"=rep(-0.51,timeseries_length)),
                   #                      
                   "1"=list("fl1"=rep(0.45,timeseries_length),
                            "fl2"=rep(0,timeseries_length),
                            "fl3"=rep(-0.51,timeseries_length)),
                   
                   
                   
                   "2"=list("fl1"=rep(-0.1,timeseries_length),
                            "fl2"=matrix(-0.1*fl2f(time_interval),nrow=timeseries_length,ncol=1),
                            "fl3"=matrix(fl3f(time_interval),nrow=timeseries_length,ncol=1)),
                   #not constant
                   # fll2=-10*sin((2*pi/25)*(time_interval-1))
                   # fll3=-20*sin((2*pi/25)*(time_interval-1))-6
                   
                   "3"=list("fl1"=rep(sample(c(mu1_coef,mu2_coef),1),timeseries_length),
                            #"fl1"=rep(-0.2,timeseries_length),
                            # "fl2"=matrix(-0.15*fl2f(time_interval),nrow=timeseries_length,ncol=1),
                            # "fl3"=matrix(fl3f(time_interval),nrow=timeseries_length,ncol=1)),
                            ###same as 2.495-2.49*time_interval
                            #works when needs to reject
                            # "fl2"=matrix(-10*sin((2*pi/25)*(time_interval-1))+10,nrow=timeseries_length,ncol=1),
                            # "fl3"=matrix(-20*sin((2*pi/25)*(time_interval-1))-6,nrow=timeseries_length,ncol=1)),
                            "fl2"=matrix(mu1_coef[1] + mu1_coef[2] * time_interval + mu1_coef[3] * time_interval^2,nrow=timeseries_length,ncol=1),
                            "fl3"=matrix(mu2_coef[1] + mu2_coef[2] * time_interval + mu2_coef[3] * time_interval^2,nrow=timeseries_length,ncol=1)),
                   #constant        
                   # "4"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
                   #          "fl2"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)+1.3145,
                   #          "fl3"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)),
                   
                   "4"=list("fl1"=matrix(fl3fn(time_interval),nrow=timeseries_length,ncol=1)-0.09,
                            "fl2"=matrix(sample(c(mu1_coef,mu2_coef),1),nrow=timeseries_length,ncol=1),
                            "fl3"=matrix(sample(c(mu1_coef,mu2_coef),1),nrow=timeseries_length,ncol=1)),
    )
    
    vec <- matrix(1:num_indvs, nrow=num_indvs, ncol=1)
    
    
    x1fl1 <- apply(vec, 1, function(x) {fda.usc::int.simpson2(time_interval, flfn$fl1, equi = TRUE, method = "TRAPZ")})
    x2fl2 <- apply( vec, 1, function(x) {fda.usc::int.simpson2(time_interval, cat_data$X[x,,2]*(flfn$fl2), equi = TRUE, method = "TRAPZ")})
    x3fl3 <- apply( vec, 1, function(x) {fda.usc::int.simpson2(time_interval, cat_data$X[x,,3]*(flfn$fl3), equi = TRUE, method = "TRAPZ")})
    
    #########
    linear_predictor <- matrix(x1fl1 + x2fl2+ x3fl3 -1)
    
    ######
    Y_indvs <- apply(linear_predictor, 1, function(x){ rbinom(1,1, 1/(1+exp(-x))) })
    #table(Y_indvs)
    
    prob_ind=1/(1+exp(-linear_predictor))
    
    
    truelist=list("TrueX1"=cat_data$X[,,1],
                  "TrueX2"=cat_data$X[,,2],
                  "TrueX3"=cat_data$X[,,3],
                  "Truecatcurve"=cat_data$W,
                  "fl"=flfn,
                  "yis"=Y_indvs,
                  "linear_predictor"=linear_predictor,
                  "prob_ind"=prob_ind)
    
    return(list("true"=truelist))
}



cfd_testing_simulation_no_paralel <- function (num_replicas, start_time, end_time, timeseries_length,
                                               mu1_coef,mu2_coef,num_indvs,fl_choice,response_family,test_type,
                                    klen=3){
  cat("CFD Testing Simulation \nNum Replicas:\t", num_replicas)
  source("source_code/R/data_generator.R")
  p_value=c(0)
  test_stats=c(0)
  #browser("mystop")
  for (i in 1:num_replicas){
    print( i)
    result <- cfd_testing(start_time, end_time, timeseries_length,
                          num_indvs,mu1_coef,mu2_coef,fl_choice,response_family,test_type, klen=3)
    p_value[i]=result$pvalue
    test_stats[i]=result$test_statistics
  }

  return(list("pvalue"=p_value,"teststat"=test_stats))

}



# source("source_code/R/time_track_function.R")
# ################
# mu1_coef=c(-6.67,-2.47,5.42)
# mu2_coef=c(-3.14,-0.99,3.91)
# 


# mu_1 <- function(t){ -0.64+4*t }
# 
# mu_2 <- function(t){ 0.97+6*t^2 }

# mu1_coef=c(-0.64,4,0)
# mu2_coef=c(0.97,0,6)
# n100t2000_mu1mu2 <- cfd_testing_simulation(num_replicas=4, start_time=0.01, end_time=0.99,
#                                                       timeseries_length=2000,
#                                                       mu1_coef=mu1_coef,
#                                                       mu2_coef=mu2_coef,
#                                                       num_indvs=100,fl_choice=3,
#                                                       response_family='bernoulli',test_type='Functional',
#                                                       klen=3)
# mean(n100t2000_mu1mu2[1,]<0.05)
# 
# n100t2000_mu1mu2 <- cfd_testing_simulation(num_replicas=100, start_time=0.01, end_time=0.99,
#                                            timeseries_length=2000,
#                                            mu1_coef=mu1_coef,
#                                            mu2_coef=mu2_coef,
#                                            num_indvs=100,fl_choice=3,
#                                            response_family='bernoulli',test_type='Inclusion',
#                                            klen=3)
# mean(n100t2000_mu1mu2[1,]<0.05)
# 
# 
# n100t2000_mu1mu2 <- cfd_testing_simulation(num_replicas=4, start_time=0.01, end_time=0.99,
#                                            timeseries_length=2000,
#                                            mu1_coef=mu1_coef,
#                                            mu2_coef=mu2_coef,
#                                            num_indvs=100,fl_choice=4,
#                                            response_family='bernoulli',test_type='Inclusion',
#                                            klen=3)
# mean(n100t2000_mu1mu2[1,]<0.05)

# n100t2000_mu1mu2 <- cfd_testing_simulation_no_paralel(num_replicas=4, start_time=0.01, end_time=0.99,
#                                            timeseries_length=2000,
#                                            mu1_coef=mu1_coef,
#                                            mu2_coef=mu2_coef,
#                                            num_indvs=100,fl_choice=3,
#                                            response_family='bernoulli',test_type='Functional',
#                                            klen=3)

# n100t2000_mu1mu2 <- cfd_testing_simulation_no_paralel(num_replicas=4, start_time=0.01, end_time=0.99,
#                                                       timeseries_length=2000,
#                                                       mu1_coef=mu1_coef,
#                                                       mu2_coef=mu2_coef,
#                                                       num_indvs=100,fl_choice=3,
#                                                       response_family='bernoulli',test_type='Inclusion',
#                                                       klen=3)
# ######test
# timeKeeperStart("n100t2000")
# set.seed(123456)
# n100t2000_nopara <- cfd_testing_simulation_no_paralel(num_replicas=5, start_time=0.01, end_time=0.99,
#                                     timeseries_length=2000,
#                                     num_indvs=100,fl_choice=3,
#                                     response_family='bernoulli',test_type='Functional',
#                                     klen=3)
# timeKeeperNext()
# powern100t2000=mean(n100t2000_nopara$pvalue<0.05)
# ################inclusion hypothesis test fl(t)=0, use fl3, which is not 0
# timeKeeperStart("n100t2000")
# set.seed(123456)
# n100t2000_nopara <- cfd_testing_simulation_no_paralel(num_replicas=5, start_time=0.01, end_time=0.99,
#                                                       timeseries_length=2000,
#                                                       num_indvs=100,fl_choice=3,
#                                                       response_family='bernoulli',test_type='Inclusion',
#                                                       klen=3)
# timeKeeperNext()
# powern100t2000=mean(n100t2000_nopara$pvalue<0.05)

################inclusion hypothesis test fl(t)=0, use fl1, which is 0, expect to fail to reject
# timeKeeperStart("n100t2000")
# set.seed(123456)
# n100t2000_nopara <- cfd_testing_simulation_no_paralel(num_replicas=5, start_time=0.01, end_time=0.99,
#                                                       timeseries_length=2000,
#                                                       mu1_coef=mu1_coef,
#                                                       mu2_coef=mu2_coef,
#                                                       num_indvs=100,fl_choice=1,
#                                                       response_family='bernoulli',test_type='Inclusion',
#                                                       klen=3)
# timeKeeperNext()
# powern100t2000=mean(n100t2000_nopara$pvalue<0.05)

##############################################
# timeKeeperStart("n100t2000")
# set.seed(123456)
# n100t2000_nopara <- cfd_testing_simulation_no_paralel(num_replicas=50, start_time=0.01, end_time=0.99,
#                                                       timeseries_length=2000,
#                                                       num_indvs=100,fl_choice=4,
#                                                       response_family='bernoulli',test_type='Functional',
#                                                       klen=3)
# timeKeeperNext()
# powern100t2000=mean(n100t2000_nopara[1,]<0.05)
# powern100t2000
############


######test
# timeKeeperStart("n100t2000")
# set.seed(123456)
# n100t2000 <- cfd_testing_simulation(num_replicas=10, start_time=0.01, end_time=0.99,
#                                     timeseries_length=2000,
#                                     num_indvs=100,fl_choice=3,
#                                     response_family='bernoulli',test_type='Functional',
#                                     klen=3)
# timeKeeperNext()
# n100t2000_data=matrix(n100t2000,nrow=2,ncol=50)
#powern100t2000=mean(n100t2000_data[1,] < 0.05)
############





# Chathura :
########fl4. =, fl6 5  fl7 10
timeKeeperStart("ch_n10000t2000")
set.seed(123456)
#need rejection rate
ch_n10000t2000 <- cfd_testing_simulation(num_replicas=20, start_time=0.01, end_time=0.99,
                                         timeseries_length=300,
                                         mu1_coef=mu1_coef,
                                         mu2_coef=mu2_coef,
                                         num_indvs=500,fl_choice=6,
                                         response_family='bernoulli',test_type='Inclusion',
                                         klen=3)
timeKeeperNext()
ch_n10000t2000_data=matrix(ch_n10000t2000,nrow=2,ncol=20)
ch_power_n10000t2000=mean(ch_n10000t2000_data[1,] < 0.05)
ch_power_n10000t2000_se=sd(unlist(ch_n10000t2000_data[1,]))/sqrt(num_replicas)

# ch_testsstat_n10000t2000=mean(unlist(ch_n10000t2000_data[2,]))
# ch_testsstat_n10000t2000sd=sd(unlist(ch_n10000t2000_data[2,]))/sqrt(50)

#######needs power 0.3
timeKeeperStart("ch_n10000t2000")
set.seed(123456)
ch_n10000t2000 <- cfd_testing_simulation(num_replicas=10, start_time=0.01, end_time=0.99,
                                         timeseries_length=300,
                                         mu1_coef=mu1_coef,
                                         mu2_coef=mu2_coef,
                                         num_indvs=1000,fl_choice=6,
                                         response_family='bernoulli',test_type='Inclusion',
                                         klen=3)
timeKeeperNext()
ch_n10000t2000_data=matrix(ch_n10000t2000,nrow=2,ncol=10)
ch_power_n10000t2000=mean(ch_n10000t2000_data[1,] < 0.05)
ch_testsstat_n10000t2000=mean(unlist(ch_n10000t2000_data[2,]))
ch_testsstat_n10000t2000sd=sd(unlist(ch_n10000t2000_data[2,]))/sqrt(50)


###needs power 1000,300, 1.  , 500, 100, 0.8
timeKeeperStart("ch_n10000t2000")
set.seed(123456)
ch_n10000t2000 <- cfd_testing_simulation(num_replicas=10, start_time=0.01, end_time=0.99,
                                         timeseries_length=100,
                                         mu1_coef=mu1_coef,
                                         mu2_coef=mu2_coef,
                                         num_indvs=500,fl_choice=7,
                                         response_family='bernoulli',test_type='Inclusion',
                                         klen=3)
timeKeeperNext()
ch_n10000t2000_data=matrix(ch_n10000t2000,nrow=2,ncol=10)
ch_power_n10000t2000=mean(ch_n10000t2000_data[1,] < 0.05)
ch_testsstat_n10000t2000=mean(unlist(ch_n10000t2000_data[2,]))
ch_testsstat_n10000t2000sd=sd(unlist(ch_n10000t2000_data[2,]))/sqrt(50)



timeKeeperStart("ch_n10000t2000")
set.seed(123456) #500, 100, 1
ch_n10000t2000 <- cfd_testing_simulation(num_replicas=10, start_time=0.01, end_time=0.99,
                                         timeseries_length=300,
                                         mu1_coef=mu1_coef,
                                         mu2_coef=mu2_coef,
                                         num_indvs=1000,fl_choice=8,
                                         response_family='bernoulli',test_type='Inclusion',
                                         klen=3)
timeKeeperNext()
ch_n10000t2000_data=matrix(ch_n10000t2000,nrow=2,ncol=10)
ch_power_n10000t2000=mean(ch_n10000t2000_data[1,] < 0.05)
ch_testsstat_n10000t2000=mean(unlist(ch_n10000t2000_data[2,]))
ch_testsstat_n10000t2000sd=sd(unlist(ch_n10000t2000_data[2,]))/sqrt(50)




timeKeeperStart("ch_n10000t2000")
set.seed(123456) #500
ch_n10000t2000 <- cfd_testing_simulation(num_replicas=10, start_time=0.01, end_time=0.99,
                                         timeseries_length=100,
                                         mu1_coef=mu1_coef,
                                         mu2_coef=mu2_coef,
                                         num_indvs=500,fl_choice=9,
                                         response_family='bernoulli',test_type='Functional',
                                         klen=3)
timeKeeperNext()
ch_n10000t2000_data=matrix(ch_n10000t2000,nrow=2,ncol=10)
ch_power_n10000t2000=mean(ch_n10000t2000_data[1,] < 0.05)
ch_testsstat_n10000t2000=mean(unlist(ch_n10000t2000_data[2,]))
ch_testsstat_n10000t2000sd=sd(unlist(ch_n10000t2000_data[2,]))/sqrt(50)
# ------------------------------
# Chathura :
# timeKeeperStart("ch_n10000t2000")
# set.seed(1234)
# ch_n10000t2000 <- cfd_testing_simulation(num_replicas=10, start_time=0.01, end_time=0.99,
#                                          timeseries_length=300,
#                                          mu1_coef=mu1_coef,
#                                          mu2_coef=mu2_coef,
#                                          num_indvs=1000,fl_choice=3,
#                                          response_family='bernoulli',test_type='Inclusion',
#                                          klen=3)
# timeKeeperNext()
# ch_n10000t2000_data=matrix(ch_n10000t2000,nrow=2,ncol=10)
# ch_power_n10000t2000=mean(ch_n10000t2000_data[1,] < 0.05)
# ch_testsstat_n10000t2000=mean(unlist(ch_n10000t2000_data[2,]))
# ch_testsstat_n10000t2000sd=sd(unlist(ch_n10000t2000_data[2,]))/sqrt(50)

#hyp3 expect to reject, it is rejecting , needs to see power 
# set.seed(1234)
# #power
timeKeeperStart("n100t300")
set.seed(1234)
n100t300 <- cfd_testing_simulation(num_replicas=10, start_time=0.01, end_time=0.99,
                                   timeseries_length=300,
                                   mu1_coef=mu1_coef,
                                   mu2_coef=mu2_coef,
                                   num_indvs=1000,fl_choice=3,
                                   response_family='bernoulli',test_type='Functional',
                                   klen=3)
timeKeeperNext()
n100t300_data=matrix(n100t300,nrow=2,ncol=10)
powern100t300=mean(n100t300_data[1,] < 0.1)
test_statisticsn100t300=mean(unlist(n100t300_data[2,]))
test_statisticsn100t300sd=sd(unlist(n100t300_data[2,]))/sqrt(5000)


####
timeKeeperStart("n100t750")
set.seed(1234)
n100t750 <- cfd_testing_simulation(num_replicas=5000, start_time=0.01, end_time=0.99,
                                   timeseries_length=750,
                                   num_indvs=100,fl_choice=3,
                                   response_family='bernoulli',test_type='Functional',
                                   klen=3)
timeKeeperNext()
n100t750_data=matrix(n100t750,nrow=2,ncol=5000)
powern100t750=mean(n100t750_data[1,] < 0.05)
test_statisticsn100t750=mean(unlist(n100t750_data[2,]))
test_statisticsn100t750sd=sd(unlist(n100t750_data[2,]))/sqrt(5000)
####
timeKeeperStart("n100t1050")
set.seed(1234)
n100t1050 <- cfd_testing_simulation(num_replicas=5000, start_time=0.01, end_time=0.99,
                                    timeseries_length=1050,
                                    num_indvs=100,fl_choice=3,
                                    response_family='bernoulli',test_type='Functional',
                                    klen=3)
timeKeeperNext()
n100t1050_data=matrix(n100t1050,nrow=2,ncol=5000)
powern100t1050=mean(n100t1050_data[1,] < 0.05)
test_statisticsn100t1050=mean(unlist(n100t1050_data[2,]))
test_statisticsn100t1050sd=sd(unlist(n100t1050_data[2,]))/sqrt(5000)
####
timeKeeperStart("n100t1350")
set.seed(1234)
n100t1350 <- cfd_testing_simulation(num_replicas=5000, start_time=0.01, end_time=0.99,
                                    timeseries_length=1350,
                                    num_indvs=100,fl_choice=3,
                                    response_family='bernoulli',test_type='Functional',
                                    klen=3)
timeKeeperNext()
n100t1350_data=matrix(n100t1350,nrow=2,ncol=5000)
powern100t1350=mean(n100t1350_data[1,] < 0.05)
test_statisticsn100t1350=mean(unlist(n100t1350_data[2,]))
test_statisticsn100t1350sd=sd(unlist(n100t1350_data[2,]))/sqrt(5000)

####
timeKeeperStart("n100t2000")
set.seed(123)
n100t2000 <- cfd_testing_simulation(num_replicas=50, start_time=0.01, end_time=0.99,
                                    timeseries_length=2000,
                                    num_indvs=100,fl_choice=3,
                                    response_family='bernoulli',test_type='Functional',
                                    klen=3)
timeKeeperNext()
n100t2000_data=matrix(n100t2000,nrow=2,ncol=4)
powern100t2000=mean(n100t2000_data[1,] < 0.05)
test_statisticsn100t2000=mean(unlist(n100t2000_data[2,]))
test_statisticsn100t2000sd=sd(unlist(n100t2000_data[2,]))/sqrt(5000)

power_table=c(powern100t300,powern100t750,powern100t1050,
              powern100t1350,powern100t2000)

test_stat=c(test_statisticsn100t300,test_statisticsn100t750,
            test_statisticsn100t1050,test_statisticsn100t1350,
            test_statisticsn100t2000)
test_statsd=c(test_statisticsn100t300sd,test_statisticsn100t750sd,
              test_statisticsn100t1050sd,test_statisticsn100t1350sd,
              test_statisticsn100t2000sd)
save(power_table,test_stat,test_statsd,file="n100power.RData")
######
#hy4 is constant, test constant, expect to fail to reject
#type I error rate
# timeKeeperStart("n100t300")
# set.seed(1234)
# #check which is less 0.05 , type I error rate
# n100t300p <- cfd_testing_simulation(num_replicas=5000, start_time=0.01, end_time=0.99,
#                                    timeseries_length=300,
#                                    num_indvs=100,fl_choice=4,
#                                    response_family='bernoulli',test_type='Functional',
#                                    klen=3)
# timeKeeperNext()
# n100t300_datap=matrix(n100t300p,nrow=2,ncol=5000)
# powern100t300p=mean(n100t300_datap[1,] < 0.05)
# test_statisticsn100t300p=mean(unlist(n100t300_datap[2,]))
# test_statisticsn100t300sdp=sd(unlist(n100t300_datap[2,]))/sqrt(5000)
# 
# 
# ####
# timeKeeperStart("n100t750")
# set.seed(1234)
# n100t750p <- cfd_testing_simulation(num_replicas=5000, start_time=0.01, end_time=0.99,
#                                    timeseries_length=750,
#                                    num_indvs=100,fl_choice=4,
#                                    response_family='bernoulli',test_type='Functional',
#                                    klen=3)
# timeKeeperNext()
# n100t750_datap=matrix(n100t750p,nrow=2,ncol=5000)
# powern100t750p=mean(n100t750_datap[1,] < 0.05)
# test_statisticsn100t750p=mean(unlist(n100t750_datap[2,]))
# test_statisticsn100t750sdp=sd(unlist(n100t750_datap[2,]))/sqrt(5000)
# ####
# timeKeeperStart("n100t1050")
# set.seed(1234)
# n100t1050p <- cfd_testing_simulation(num_replicas=5000, start_time=0.01, end_time=0.99,
#                                     timeseries_length=1050,
#                                     num_indvs=100,fl_choice=4,
#                                     response_family='bernoulli',test_type='Functional',
#                                     klen=3)
# timeKeeperNext()
# n100t1050_datap=matrix(n100t1050p,nrow=2,ncol=5000)
# powern100t1050p=mean(n100t1050_datap[1,] < 0.05)
# test_statisticsn100t1050p=mean(unlist(n100t1050_datap[2,]))
# test_statisticsn100t1050sdp=sd(unlist(n100t1050_datap[2,]))/sqrt(5000)
# ####
# timeKeeperStart("n100t1350")
# set.seed(1234)
# n100t1350p <- cfd_testing_simulation(num_replicas=5000, start_time=0.01, end_time=0.99,
#                                     timeseries_length=1350,
#                                     num_indvs=100,fl_choice=4,
#                                     response_family='bernoulli',test_type='Functional',
#                                     klen=3)
# timeKeeperNext()
# n100t1350_datap=matrix(n100t1350p,nrow=2,ncol=5000)
# powern100t1350p=mean(n100t1350_datap[1,] < 0.05)
# test_statisticsn100t1350p=mean(unlist(n100t1350_datap[2,]))
# test_statisticsn100t1350sdp=sd(unlist(n100t1350_datap[2,]))/sqrt(5000)
# 
# ####
# timeKeeperStart("n100t2000")
# set.seed(1234)
# n100t2000p <- cfd_testing_simulation(num_replicas=5000, start_time=0.01, end_time=0.99,
#                                     timeseries_length=2000,
#                                     num_indvs=100,fl_choice=4,
#                                     response_family='bernoulli',test_type='Functional',
#                                     klen=3)
# timeKeeperNext()
# n100t2000_datap=matrix(n100t2000p,nrow=2,ncol=5000)
# powern100t2000p=mean(n100t2000_datap[1,] < 0.05)
# test_statisticsn100t2000p=mean(unlist(n100t2000_datap[2,]))
# test_statisticsn100t2000sdp=sd(unlist(n100t2000_datap[2,]))/sqrt(5000)
# 
# power_tablep=c(powern100t300p,powern100t750p,powern100t1050p,
#               powern100t1350p,powern100t2000p)
# 
# test_statp=c(test_statisticsn100t300p,test_statisticsn100t750p,
#             test_statisticsn100t1050p,test_statisticsn100t1350p,
#             test_statisticsn100t2000p)
# test_statsdp=c(test_statisticsn100t300sdp,test_statisticsn100t750sdp,
#               test_statisticsn100t1050sdp,test_statisticsn100t1350sdp,
#               test_statisticsn100t2000sdp)
# save(power_tablep,test_statp,test_statsdp,file="n100typeIerror.RData")
###########

# set.seed(123)
# hy3test=cfd_testing (start_time=0.01, end_time=0.99, timeseries_length=2000,
#                         num_indvs=100,fl_choice=3,response_family="'bernoulli'",test_type='Functional',
#                         klen=3)
