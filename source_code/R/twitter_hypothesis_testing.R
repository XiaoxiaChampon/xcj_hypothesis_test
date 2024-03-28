
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
# Date: 12/08/2023
#
##############################################################

library(mgcv)
library(fda)
library(fda.usc)
library(devtools)
# install_github("stchen3/glmmVCtest")
library("glmmVCtest")
library(RLRsim)
library(MASS)
library(splines)
library(parallel)

##W needs t*n
##X needs n * t *3

source("/Users/xzhao17/Documents/GitHub/xcj_hypothesis_test_cfd/source_code/R/data_generator.R")
source("/Users/xzhao17/Documents/GitHub/xcj_hypothesis_test_cfd/source_code/R/cfd_hypothesis_test.R")

#' Funciton to do the hypothesis testing on categorical functional data
#' @param Y: 1D array, length num_indvs, class membership for num_indvs individuals
#' @param cfd: 2d array, num_indvs* timeseries_length, categorical functional data
#'              observed for num_indvs individuals and timeseries_length time points
#' @param time_interval, 1D array, observational time
#' @param response_family, currently only support 'bernoulli'
#' @param test_type, "Function" or "constant"
#' @return list: statistics is test statistics, pvalue
#' @param  test_option: 1-test.aRLRT 2-gam
cfd_hypothesis_test_twitter <- function(Y, cfd, time_interval, response_family, test_type,test_option){
    X_cfd <- GetXFromW(cfd)
    # 
     num_indvs <- length(Y)
    
    if (test_option==1){
        
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
        
        # 
        #For testing in models with multiple variance
        #' components, the fitted model \code{m} must contain \bold{only} the random
        #' effect set to zero under the null hypothesis, while \code{mA} and \code{m0}
        #' are the models under the alternative and the null, respectively. 
        # 
        alternative_fit <- fit.glmmPQL(test_matrix, response_family, num_indvs, test_type)

        result_try <- try(test.aRLRT(alternative_fit), silent=T)

        if(is.atomic(result_try)){
            return(list(statistics=NULL, pvalue=NULL))
        }
        # 
         result <- result_try$aRLRT 
        # 
         return(list(statistics=result$statistic, pvalue=result$p.value,fit=alternative_fit))
    }else if (test_option==2){
        
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
        gam_test <- try(gam(cbind(Y, num_indvs - Y) ~ 0 + Xmat_test_type + 
                                s(Zmat_test_type_2$Zmat, bs = 're')+ 
                                s(Zmat_test_type_3$Zmat, bs = 're'), family = 'binomial'))
        
        if('try-error' %in% class(gam_test)){return(list(statistics=NULL, pvalue=NULL))}
        return(list(statistics=gam_test[1,3], pvalue=gam_test[1,4]))
    }

}
#W label no tweets: 1, positive 2, neutral 3, negative 4.
##load the data
#Y_label
#0         1 
#0.7843137 0.2156863  
#positive negative neutral
# W_real= read.csv("/Users/xzhao17/Documents/GitHub/xcj_hypothesis_test_cfd/Real_data/output_wdf.csv")
# Y_real= read.csv("/Users/xzhao17/Documents/GitHub/xcj_hypothesis_test_cfd/Real_data/Y_data114days.csv")

#mixed, consistent, 114 days, June 1, Sep 23 (collection 9/12-10/27)
W_real= read.csv("/Users/xzhao17/Documents/GitHub/xcj_hypothesis_test_cfd/Real_data/W_sentiment_mix_type_matrix_f1D_mt10.csv")
Y_real= read.csv("/Users/xzhao17/Documents/GitHub/xcj_hypothesis_test_cfd/Real_data/Y_sentiment_mix_type_diff_f1D_mt10.csv")



#####try 102 days, June 1, to Sep 11


# > dim(Y_real)
# [1] 529  11
# > dim(W_real)
# [1] 546 117

user_index_intersect = which(W_real$Author %in% Y_real$Author)
Y_label=Y_real$Y[user_index_intersect]
# Y_label
# 0        1 
# 0.502924 0.497076
W_final=t(W_real[user_index_intersect,-1])
#######
dim(W_final)
table(W_final)/sum(table(W_final))
#W_final
#0          1          2 
#0.48155270 0.43999739 0.07844991 
#######
W_merge=W_final
#W_merge[W_merge==4]=3
#table(W_merge)/sum(table(W_merge))

table(W_merge)
#W_merge
#0     1     2 
#29550 27000  4814
####

W_merge=W_merge[1:104,]
table(W_merge)
# W_merge
# 0     1     2 
# 27041 23758  4217 

# 0          1          2 
# 0.49151156 0.43183801 0.07665043 





X_cfd_twitter <- GetXFromW(W_merge)
#sum(X_cfd_twitter[,,1])
#[1] 29550
#> sum(X_cfd_twitter[,,2])
#[1] 27000
#> sum(X_cfd_twitter[,,3])
#[1] 4814

#time_interval=seq(0,1,length.out=dim(W_final)[1])
time_interval=seq(0,1,length=dim(W_merge)[1])
num_indvs <- length(Y_label)
# num_indvs
# [1] 529
#Y_label
#0         1 
#0.8050682 0.1949318 

#table(Y_label)/sum(table(Y_label))
#Y_label
#0        1 
#0.502924 0.497076 



#######################Dr. Xiao suggested Feb 22
#estimate Fl
time_interval_matrix=do.call("rbind", replicate( 
    length(Y_label), time_interval, simplify = FALSE)) 

fl_gam30=gam(Y_label~s(time_interval_matrix,by=X_cfd_twitter[,,2],k = 30)+s(time_interval_matrix,by=X_cfd_twitter[,,3],k = 30),family = 'binomial')

#summary(fl_gam)
# Family: binomial 
# Link function: logit 
# 
# Formula:
#     Y_label ~ s(time_interval_matrix, by = X_cfd_twitter[, , 2]) + 
#     s(time_interval_matrix, by = X_cfd_twitter[, , 3])
# 
# Parametric coefficients:
#     Estimate Std. Error z value Pr(>|z|)
# (Intercept)  0.07601    0.19266   0.395    0.693
# 
# Approximate significance of smooth terms:
#     edf Ref.df Chi.sq
# s(time_interval_matrix):X_cfd_twitter[, , 2] 2.272  2.489  3.437
# s(time_interval_matrix):X_cfd_twitter[, , 3] 2.000  2.000  1.632
# p-value
# s(time_interval_matrix):X_cfd_twitter[, , 2]   0.186
# s(time_interval_matrix):X_cfd_twitter[, , 3]   0.442
# 
# R-sq.(adj) =  0.00495   Deviance explained = 0.977%
# UBRE = 0.39327  Scale est. = 1         n = 513

##############
library(mgcViz)
fl_gam_model=getViz(fl_gam)
print(plot(fl_gam_model, allTerms = T,xlab="Time",ylab="Value"), pages = 1)
##########
check(fl_gam_model,
      a.qq = list(method = "tnorm", 
                  a.cipoly = list(fill = "light blue")), 
      a.respoi = list(size = 0.5), 
      a.hist = list(bins = 10))

# Method: UBRE   Optimizer: outer newton
# full convergence after 7 iterations.
# Gradient range [-3.441343e-07,3.181516e-09]
# (score 0.3932711 & scale 1).
# Hessian positive definite, eigenvalue range [3.440704e-07,0.0002321639].
# Model rank =  21 / 21 
# 
# Basis dimension (k) checking results. Low p-value (k-index<1) may
# indicate that k is too low, especially if edf is close to k'.
# 
#                                                 k'   edf k-index p-value
# s(time_interval_matrix):X_cfd_twitter[, , 2] 10.00  2.27      NA      NA
# s(time_interval_matrix):X_cfd_twitter[, , 3] 10.00  2.00      NA      NA
#########
xlfl2plot <- plot( sm(fl_gam_model, 1) )
xlfl2plot + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
    l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
    xlab("Time")+
    ylab("Value")+
    theme(text = element_text(size = 20)) +
    l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()
########################################




Zmat_Inc2<-get_Zmatrix(X_cfd_twitter[,1:104,2],time_interval,test_type='Inclusion')
Zmat_Inc.mat2 <- Zmat_Inc2$Zmat

Zmat_Inc3<-get_Zmatrix(X_cfd_twitter[,1:104,3],time_interval,test_type='Inclusion')
Zmat_Inc.mat3 <- Zmat_Inc3$Zmat


Xmat_Inc<-matrix(rep(1, num_indvs),ncol=1)

test_matrix <- data.frame(Y=Y_label,
                          X=Xmat_Inc,
                          Z.test=Zmat_Inc2$Zmat,
                          Z.test3=Zmat_Inc3$Zmat,
                          #Z.test4=Zmat_Func4$Zmat,
                          ones=rep(1,num_indvs))
names(test_matrix) <- c('Y','X1',
                        paste0('Z.test',1:ncol(Zmat_Inc2$Zmat)),
                        paste0('Z.test3',1:ncol(Zmat_Inc3$Zmat)),
                        #paste0('Z.test4',1:ncol(Zmat_Func4$Zmat)),
                        "ones")

Zmat_Func2 <- get_Zmatrix(X_cfd_twitter[,1:104,2], time_interval, test_type="Functional")
#Zmat_Func3 <- get_Zmatrix(X_cfd_twitter[,,3], time_interval, test_type="Functional")
Zmat_Func3 <- get_Zmatrix(X_cfd_twitter[,1:104,3], time_interval, test_type="Linearity")

Zmat_Func.mat2=Zmat_Func2$Zmat
Zmat_Func.mat3=Zmat_Func3$Zmat

Xmat_Func <- cbind(Xmat_Inc, Zmat_Func2$X.g2, Zmat_Func3$X.g2)

num_indvs <- length(Y_label)
gam_Inc <- gam(cbind(Y_label, 1 - Y_label) ~ 0 + Xmat_Inc + 
                   s(Zmat_Inc.mat2, bs = 're')+ s(Zmat_Inc.mat3, bs = 're'), family = 'binomial')

# gam_m=gam(cbind(Y_label, num_indvs - Y_label) ~ 0 +s(Zmat_Inc.mat2, bs = 're'), family = 'binomial')
# gam_m0 <- gam(cbind(Y_label, num_indvs - Y_label) ~ 0 + Xmat_Inc + s(Zmat_Inc.mat3, bs = 're'), family = 'binomial')
# 
# exactRLRT(gam_m, gam_Inc, gam_m0)
gam.vcomp(gam_Inc )
# gam.vcomp(gam_Inc )
# s(Zmat_Inc.mat2) s(Zmat_Inc.mat3) 
# 0.003202755      0.394893930

summary(gam_Inc)

# Family: binomial 
# Link function: logit 
# 
# Formula:
#     cbind(Y_label, 1 - Y_label) ~ 0 + Xmat_Inc + s(Zmat_Inc.mat2, 
#                                                    bs = "re") + s(Zmat_Inc.mat3, bs = "re")
# 
# Parametric coefficients:
#     Estimate Std. Error z value Pr(>|z|)
# Xmat_Inc -0.02822    0.09192  -0.307    0.759
# 
# Approximate significance of smooth terms:
#     edf Ref.df Chi.sq p-value
# s(Zmat_Inc.mat2) 9.445e-05      1  0.000   0.946
# s(Zmat_Inc.mat3) 2.972e-01      1  0.422   0.234
# 
# R-sq.(adj) =  0.000819   Deviance explained = 0.101%
#UBRE = 0.38992  Scale est. = 1         n = 513

# gam_Inc$coefficients
# Xmat_Inc s(Zmat_Inc.mat2).1 s(Zmat_Inc.mat3).1 
# -2.821905e-02      -2.106267e-06       2.149345e-01


# 
# gam_Func <- try(gam(cbind(Y_label, 1- Y_label) ~ 0 + Xmat_Func + s(Zmat_Func.mat2, bs = 're')+
#                         s(Zmat_Func.mat3, bs = 're'), family = 'binomial'), silent = F)
gam_Func <- gam(cbind(Y_label, 1- Y_label) ~ 0 + Xmat_Func + s(Zmat_Func.mat2, bs = 're')+
                        s(Zmat_Func.mat3, bs = 're'), family = 'binomial')

library(mgcv)
gam_Func_check <- gam(cbind(Y_label, 1- Y_label) ~ 0 + Xmat_Func + s(Zmat_Func.mat2, bs = 'cr')+
                    s(Zmat_Func.mat3, bs = 'cr'), family = 'binomial')
# gam.vcomp(gam_Func )
# s(Zmat_Func.mat2) s(Zmat_Func.mat3) 
# 1.2942934         0.0152765 

summary(gam_Func)

# Family: binomial 
# Link function: logit 
# 
# Formula:
#     cbind(Y_label, 1 - Y_label) ~ 0 + Xmat_Func + s(Zmat_Func.mat2, 
#                                                     bs = "re") + s(Zmat_Func.mat3, bs = "re")
# 
# Parametric coefficients:
#     Estimate Std. Error z value Pr(>|z|)
# Xmat_Func           0.1534     0.2153   0.713    0.476
# Xmat_Func          -0.3357     0.3732  -0.900    0.368
# Xmat_Funcconstant   2.5890     3.5346   0.732    0.464
# Xmat_Funclin       -4.2712     6.9470  -0.615    0.539
# 
# Approximate significance of smooth terms:
#     edf Ref.df Chi.sq p-value  
# s(Zmat_Func.mat2) 0.8247220      1  4.648  0.0176 *
#     s(Zmat_Func.mat3) 0.0009601      1  0.001  0.3562  
# ---
#     Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.00637   Deviance explained = 1.02%
# UBRE = 0.3909  Scale est. = 1         n=513

gam_Func$coefficients
# Xmat_Func           Xmat_Func   Xmat_Funcconstant 
# 0.1534393106       -0.3357014655        2.5889744841 
# Xmat_Funclin s(Zmat_Func.mat2).1 s(Zmat_Func.mat3).1 
# -4.2711668762       -1.1682006636        0.0004364785

summary(gam_Inc)$s.table

summary(gam_Func)$s.table
# c(summary(gam_Inc)$s.table[1,4])
# 
# c(summary(gam_Func)$s.table[1,4])
# summary(gam_Inc)$s.table
# edf Ref.df       Chi.sq   p-value
# s(Zmat_Inc.mat2) 2.912717e-05      1 2.607989e-07 0.9246134
# s(Zmat_Inc.mat3) 9.819180e-05      1 7.001697e-05 0.3984303
# > 
#     > summary(gam_Func)$s.table
# edf Ref.df       Chi.sq    p-value
# s(Zmat_Func.mat2) 0.6509104295      1 1.860148e+00 0.09093586
# s(Zmat_Func.mat3) 0.0000578697      1 3.869628e-09 0.99347580
###L2 penalty for non test item
# summary(gam_Inc)$s.table
# edf Ref.df       Chi.sq   p-value
# s(Zmat_Inc.mat2) 9.444705e-05      1 4.325350e-07 0.9460461
# s(Zmat_Inc.mat3) 2.972242e-01      1 4.215366e-01 0.2336949
# > 
#     > summary(gam_Func)$s.table
# edf Ref.df       Chi.sq    p-value
# s(Zmat_Func.mat2) 0.8247220002      1 4.6477717842 0.01760119
# s(Zmat_Func.mat3) 0.0009601127      1 0.0008171459 0.35624727

###########################new days
# summary(gam_Inc)
# 
# Family: binomial 
# Link function: logit 
# 
# Formula:
#     cbind(Y_label, 1 - Y_label) ~ 0 + Xmat_Inc + s(Zmat_Inc.mat2, 
#                                                    bs = "re") + s(Zmat_Inc.mat3, bs = "re")
# 
# Parametric coefficients:
#     Estimate Std. Error z value Pr(>|z|)
# Xmat_Inc -0.02190    0.09062  -0.242    0.809
# 
# Approximate significance of smooth terms:
#     edf Ref.df Chi.sq p-value
# s(Zmat_Inc.mat2) 8.769e-05      1  0.000   0.938
# s(Zmat_Inc.mat3) 2.016e-01      1  0.252   0.263
# 
# R-sq.(adj) =  0.00049   Deviance explained = 0.0638%
# UBRE = 0.39006  Scale est. = 1         n = 513

# gam.vcomp(gam_Inc )
# s(Zmat_Inc.mat2) s(Zmat_Inc.mat3) 
# 0.002969493      0.304379611 

###############
#summary(gam_Func)
# Family: binomial 
# Link function: logit 
# 
# Formula:
#     cbind(Y_label, 1 - Y_label) ~ 0 + Xmat_Func + s(Zmat_Func.mat2, 
#                                                     bs = "re") + s(Zmat_Func.mat3, bs = "re")
# 
# Parametric coefficients:
#     Estimate Std. Error z value Pr(>|z|)
# Xmat_Func          0.05449    0.18937   0.288    0.774
# Xmat_Func         -0.18318    0.34083  -0.537    0.591
# Xmat_Funcconstant  4.43569    3.68187   1.205    0.228
# Xmat_Funclin      -7.73363    7.24927  -1.067    0.286
# 
# Approximate significance of smooth terms:
#     edf Ref.df Chi.sq p-value  
# s(Zmat_Func.mat2) 7.951e-01      1  3.836  0.0281 *
#     s(Zmat_Func.mat3) 8.834e-05      1  0.000  0.9939  
# ---
#     Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.00601   Deviance explained = 0.994%
# UBRE = 0.39117  Scale est. = 1         n = 513


# gam.vcomp(gam_Func )
# s(Zmat_Func.mat2) s(Zmat_Func.mat3) 
# 1.218281383       0.004983316






###############
library(gratia)
comp <- compare_smooths(gam_Inc, gam_Func)
draw(comp)
################
###
plot(gam_Func,pages=1,scheme=2)
#plot(gam_Func)
#####


library(mgcViz)
gam_Inc_model=getViz(gam_Inc)
print(plot(gam_Inc_model, allTerms = T), pages = 1)

gam_Func_model=getViz(gam_Func)
print(plot(gam_Func_model, allTerms = T), pages = 1)


############
wftest_inclusion =cfd_hypothesis_test_twitter(Y_label, W_merge, seq(0,1,length.out=dim(W_final)[1]),
                                     response_family='bernoulli',test_type='Inclusion',1)
wftest_inclusion$pvalue
#wftest_inclusion$fit$fit$coefficients

#[1] 1
#consistent
#wftest_inclusion$pvalue
#[1] 0.2723
##########
wftest_functional =cfd_hypothesis_test_twitter(Y_label, W_merge, seq(0,1,length.out=dim(W_final)[1]),
                                               response_family='bernoulli',test_type='Functional',1)
wftest_functional$pvalue
wftest_functional$fit$fit$coefficients
#wftest_functional$pvalue
#[1] 0.0205
# $fixed
# X1         X2 
# -6.8800119 -0.1581793 
# 
# $random
# $random$ones
# Z.test1     Z.test2       Z.test3     Z.test4     Z.test5     Z.test6
# 1 -0.0009789216 0.007053474 -0.0004366179 0.001589926 0.002239443 0.006296457
# Z.test7    Z.test8     Z.test9    Z.test10    Z.test11    Z.test12
# 1 -0.005454433 0.01641336 -0.01770479 0.002909879 -0.00135741 -0.02320621
# Z.test13   Z.test14   Z.test15   Z.test16  Z.test17    Z.test18    Z.test19
# 1 -0.01793816 0.01870753 0.01544055 0.03212505 0.0275855 -0.03265636 -0.02706912
# Z.test20   Z.test21   Z.test22  Z.test23   Z.test24   Z.test25   Z.test26
# 1 -0.05824782 0.01003825 0.08738908 0.0790089 0.01063501 0.02378343 -0.1391511
# Z.test27   Z.test28   Z.test29
# 1 0.04930252 -0.3321262 -0.5838096
##########
#march 8 , 2024, create a figfure for tweets
brand_controvercial_users=read.csv("/Users/xzhao17/Documents/GitHub/xcj_hypothesis_test_cfd/data/2041500801/catFDA_paper2_finaltesting.csv",
                                   skip=5)
save(brand_controvercial_users,file="brand_controvercial_users.RData")
#########
library(stringr)
library(dplyr)

find_unique_authors_brand_cond=function(key_brand){
  subset_brand_controvercial_users=brand_controvercial_users %>%
    filter(str_detect(Title, key_brand))
  unique(subset_brand_controvercial_users$Author)
}

#########
find_unique_authors_brand_cond("MacDonald")#1
#########
starbucks_user=find_unique_authors_brand_cond("Starbucks") #129
########
starbucks_positive=brand_controvercial_users[(brand_controvercial_users$Author%in% starbucks_user)&(brand_controvercial_users$Sentiment) %in% c("positive","negative"),]
unique(starbucks_positive$Author) #positive 77, neg, posi 98
#######nypost
find_sentiment_table=function(author_name_vector,index_start,index_end){
    
       for (i in index_start:index_end){
           cat("----------------","\n","User : ", author_name_vector[i], "\n",
                         table(brand_controvercial_users[brand_controvercial_users$Author==author_name_vector[i],c("Sentiment")]),
                        "----------------","\n")
           
          }
    
     }
 ####
  find_sentiment_table(starbucks_positive$Author,21,40)
##
  # User :  SBWorkersUnited 
  # 12 32 8 ---------------- 9/25-10/27
brand_controvercial_users[brand_controvercial_users$Author=="SBWorkersUnited",c("Title","Sentiment","Date")]
