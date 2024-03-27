
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
# Purpose: Power Curves for Categorical Functional Data Hypothesis Testing
#         
# Author:  Xiaoxia Champon
# Date: 12/04/2023
#
##############################################################
###
load("./catfda2hypothesis/EXP1_r5000_cfda2n100t90.RData")
n100t90deltapower=final_table[1:5,1:5]
n100t90deltapower
###
load("./catfda2hypothesis/EXP1_n500_tlen90_180_runs.RData")
n500t90deltapower=final_table[1:5,1:5]
n500t90deltapower
###
load("./catfda2hypothesis/EXP1_r5000_cfda2fulln500t300.RData")
n100n500t300deltapower=final_table[c(1:5,11:15),1:5]
n100n500t300deltapower
#####
n100n500t90t300deltapower=rbind(n100t90deltapower,n500t90deltapower,n100n500t300deltapower)
n100n500t90t300deltapower$fl_choice=as.numeric(n100n500t90t300deltapower$fl_choice)
str(n100n500t90t300deltapower)
# Define the conditions and replacement values
conditions <- c(1, 2,3, 4,5)
replacement_values <- c(0, 5,10,15,20)

replace_all <- function(old_vector, old_values, new_values){
    c(new_values, old_vector)[match(old_vector, c(old_values, old_vector))]
}

# Use replace() to replace the names in the 'Names' column
n100n500t90t300deltapower$fl_choice <- replace_all(n100n500t90t300deltapower$fl_choice, 
                                                   conditions, 
                                                   replacement_values)

save(n100n500t90t300deltapower,file="n100n500t90t300deltapower.RData")
####################################################
library(ggplot2)

n100n500t90t300deltapower_plot=ggplot(n100n500t90t300deltapower,
    aes(x=fl_choice,y=unlist(power),color=as.factor(num_subjects),shape=as.factor(num_subjects)))+
    geom_line() +
    facet_grid(. ~ paste0("Time Points: ",as.factor(num_timepoints)))+
    ylab("Power")+
    xlab(expression(~delta~": Deivation from 0"))+
    guides(color = guide_legend(title = "Subjects")) +
    theme(text = element_text(size = 20))  +
    theme(
        legend.position = c(.02, .95),
        legend.justification = c("left", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)
    )
n100n500t90t300deltapower_plot
ggsave("n100n500t90t300deltapower_plot.png")

##########
load("./EXP2_r5000_cfda2.RData")
n100t90t180final_table=final_table[,1:5]
n100t90t180final_table[,4:5]


load("./EXP3_r5000_cfda2n500t90180.RData")
n500t90t180final_table=final_table[,1:5]
# library(xtable)
# xtable(final_table,digits=4)

load("EXP3_r5000_cfda2n200300400t90fl11fl15.RData")
# library(xtable)
# xtable(final_table,digits=4)
n200ton400t90_table=final_table[,1:5]
n200ton400t90_table
#######
#get just 90 points, n=100 to 500 data
n100ton500_t90 = rbind(n100t90t180final_table[n100t90t180final_table$num_timepoints!=180,],
                 n500t90t180final_table[n500t90t180final_table$num_timepoints!=180,],
                 n200ton400t90_table)
not_list = c("14", "15")
n100ton500_t90_final= n100ton500_t90[!n100ton500_t90$fl_choice %in% not_list,]
n100ton500_t90_final

n100ton500_t90_final$fl_choice=as.numeric(n100ton500_t90_final$fl_choice)
# Define the conditions and replacement values
conditions_fl11fl13 <- c(1, 2,3)
replacement_values_fl11fl13 <- c(0, 20,40)

# Use replace() to replace the names in the 'Names' column
n100ton500_t90_final$fl_choice <- replace(n100ton500_t90_final$fl_choice, 
                                          n100ton500_t90_final$fl_choice %in% conditions_fl11fl13, 
                                               replacement_values_fl11fl13)
save(n100ton500_t90_final,file="n100ton500_t90_finalfl11fl13.RData")


n100ton500_t90_final
####graph the power

n100ton500t90deltapower_plot=ggplot(n100ton500_t90_final,
                                      aes(x=fl_choice,y=unlist(power),
                                          color=as.factor(num_subjects),
                                          shape=as.factor(num_subjects)))+
    geom_line() +
    facet_grid(. ~ paste0("Test Type: ",test_type))+
    ylab("Power")+
    xlab(expression(~1~"+"~delta~": Slope deviation from 0"))+
    guides(color = guide_legend(title = "Subjects")) +
    theme(text = element_text(size = 20))  +
    theme(
        legend.position = c(.02, .95),
        legend.justification = c("left", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)
    )
n100ton500t90deltapower_plot
ggsave("n100ton500t90deltapower_plot.png")

##########
# load("EXP2_r5000_cfda2n100n500t300fl11fl15.RData")
# library(xtable)
# xtable(final_table,digits=4)
load("EXP5_r5000_cfda2n100n500t90fl21fl25.RData")
xtable(final_table,digits=4)
n100n500t90fl21tofl25=final_table[,1:5]
n100n500t90fl21tofl25
save(n100n500t90fl21tofl25,file="n100n500t90fl21tofl25.RData")
#####power curve for fl21-fl25 just n=100 500 t90
n100n500t90fl21tofl25$fl_choice=as.numeric(n100n500t90fl21tofl25$fl_choice)
# Define the conditions and replacement values
conditions_fl21fl25 <- c(1, 2,3,4,5)
replacement_values_fl21fl25 <- c(0, 10,20,30,40)

# Use replace() to replace the names in the 'Names' column
n100n500t90fl21tofl25$fl_choice <- replace(n100n500t90fl21tofl25$fl_choice, 
                                           n100n500t90fl21tofl25$fl_choice %in% conditions_fl21fl25, 
                                           replacement_values_fl21fl25)
save(n100n500t90fl21tofl25,file="n100n500t90fl21tofl25.RData")

load("n100n500t90fl21tofl25.RData")
n100n500t90fl21fl25_plot=ggplot(n100n500t90fl21tofl25,
                                    aes(x=fl_choice,y=unlist(power),
                                        color=as.factor(num_subjects),
                                        shape=as.factor(num_subjects)))+
    geom_line() +
    facet_grid(. ~ paste0("Test Type: ",test_type))+
    ylab("Power")+
    xlab(expression(~delta~": Slope"))+
    guides(color = guide_legend(title = "Subjects")) +
    theme(text = element_text(size = 20))  +
    theme(
        legend.position = c(.02, .95),
        legend.justification = c("left", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)
    )
n100n500t90fl21fl25_plot
ggsave("n100n500t90fl21fl25_plot.png")

########100to500 t90 fl21tofl25
load("EXP4_r5000_cfda2t90.RData")
library(xtable)
n100ton500t90fl21tofl25_deltaneg1 = xtable(final_table,digits=4)
save(n100ton500t90fl21tofl25_deltaneg1,file="n100ton500t90fl21tofl25_deltaneg1.RData")

##############new final table, new slope with every combination
load("EXP2NewIntercept_r5000_cfda2.RData")
xtable(final_table,digits=4)

###
my_files <- list.files(path = "EXP2NewIntercept_r5000_cfda2_outputs", pattern = "*.RData", full.names = T)

# all_data <- lapply(my_files, load, .GlobalEnv)
# library(dplyr)
# data_newslope = bind_rows(mget(unlist(all_data)))
load(my_files[1])
mean(simulation_pvalues[1,] < 0.05)
mean(simulation_pvalues[1,] < 0.1)

load(my_files[2])
mean(simulation_pvalues[1,] < 0.05)
mean(simulation_pvalues[1,] < 0.1)

load(my_files[3])
mean(simulation_pvalues[1,] < 0.05)
mean(simulation_pvalues[1,] < 0.1)

###########not working
list_number =  length(my_files )
power_for_01=c(0)
power_for_005=c(0)
for (i in 1:list_number){
    file_sample = load(my_files[i])
    power_for_005[i] = mean(simulation_pvalues[1,] < 0.05)
    power_for_01[i] = mean(simulation_pvalues[1,] < 0.1)
    return (list(power_for_01,power_for_005))
}

#################fl5tofl10 new slope last trial including funciontal for fl=5
load("EXP2_r5000_cfda2outputstypeI.RData")
library(xtable)
xtable(final_table,digits=4)
n100n500fl678212226 =final_table[,1:6]
n100n500fl678212226
save(n100n500fl678212226,file="n100n500fl678212226.RData")
################ 180 points
load("EXP2_r5000_cfda2outputstypeI180.RData")
xtable(final_table,digits=4)
n100n500t180fl678212226 =final_table[,1:6]
n100n500t180fl678212226
save(n100n500t180fl678212226,file="n100n500t180fl678212226.RData")
###########################fl 6, 200, 7
load("EXP2_r5000_cfda2outputstypeI180fl200.RData")
xtable(final_table,digits=4)
typeI180fl200=final_table[,1:7]
typeI180fl200
xtable(typeI180fl200,digits=4)
save(typeI180fl200,file="typeI180fl200.RData")
###########################
load("EXP3_r5000_cfda2n300n1000t180typeI.RData")
xtable(final_table,digits=4)
n300n1000t180fl7fl200=final_table[,1:7]
save(n300n1000t180fl7fl200,file="n300n1000t180fl7fl200.RData")
#######################
load("EXP6_newintercept_t90fl21fl25power.RData")
xtable(final_table,digits=4)
#########################
load("EXP2_r5000_cfda2outputstypeI180fl200n300.RData")
xtable(final_table,digits=4)
n300t180_newintercept_fl6200721 =final_table[,1:7]
save(n300t180_newintercept_fl6200721,file="n300t180_newintercept_fl6200721.RData")
#############
load("EXP4_r5000_cfda2n1000t180typeI.RData")
xtable(final_table[,1:7],digits=4)

##############
load("EXP6_newintercept_fl6fl10t90power.RData")
xtable(final_table[,1:7],digits=4)

################
#final mu1, mu2
load("EXP7_outputfl6fl10.RData")
library(xtable)
xtable(final_table[,1:7],digits=4)

n100n500t180fl6fl10finalmu = final_table[,1:8]
save(n100n500t180fl6fl10finalmu,file="n100n500t180fl6fl10finalmu.RData")


#################
load("EXP8_outputfl21fl25.RData")
xtable(final_table[,1:8],digits=4)
n100n500t180fl21fl25finalmu = final_table[,1:8]
save(n100n500t180fl21fl25finalmu,file="n100n500t180fl21fl25finalmu.RData")
#########
load("EXP9_outputfl6200fl7.RData")
xtable(final_table[,1:8],digits=4)
n100n500t180fl6fl200fl7finalmu = final_table[,1:8]
save(n100n500t180fl6fl200fl7finalmu,file="n100n500t180fl6fl200fl7finalmu.RData")
#############################
load("EXP9_outputfl11fl15.RData")
xtable(final_table[,1:8],digits=4)
n100n500t180fl11fl15finalmu = final_table[,1:8]
save(n100n500t180fl11fl15finalmu,file="n100n500t180fl11fl15finalmu.RData")

###
# load("typeI_censor_balance.RData")
# xtable(typeI_censor_balance)
# ###3
# load("typeI_censor_unbalance.RData")
# xtable(typeI_censor_unbalance)
# 
# ###
# load("power_censor_balance.RData")
# xtable(power_censor_balance)
# 
# #
# load("power_censor_unbalance.RData")
# xtable(power_censor_unbalance)

##################################################
####################################################
####################################################
##manuscript table
load("EXP3_typeIJAN.RData")
typeIjan9=final_table
tpyeIjan9sub=typeIjan9[(typeIjan9$num_subjects!=1000 &typeIjan9$fl_choice %in% c(6,200,7,21)),]
tpyeIjan9sub

load("EXP3_typeIJAN3002.RData")
EXP3_typeIJAN300=final_table
EXP3_typeIJAN300_sub=EXP3_typeIJAN300[EXP3_typeIJAN300$fl_choice!=8,]
table1_manuscrip=rbind(tpyeIjan9sub,EXP3_typeIJAN300_sub)
table1_manuscrip
#save(table1_manuscrip,file="./simulation_data/table1_manuscrip.RData")


table_list_column=table1_manuscrip[,1:8]
table_part1_unlist=data.frame(matrix(unlist(table_list_column[,5:8]),nrow=15))
colnames(table_part1_unlist)=c("power","se","power_01","se_01")
table_part2=table1_manuscrip[,1:4]
table_1_final=cbind(table_part1_unlist,table_part2)
table_1_final
save(table_1_final,file="./simulation_data/table_1_final.RData")

########supplement for time points 180
#from simulation results before
load("EXP9_outputfl6200fl7.RData")
typeIjan9=final_table[c(1,5:9,13:16),1:8]
levels(typeIjan9$fl_choice) <- c("6", "200","7","21")

load("EXP3_aRLRTFeb5.RData")
EXP3_typeIJAN300=final_table[(final_table$num_subjects==300 & final_table$num_timepoints==180),1:8]
EXP3_typeIJAN300
table1_manuscrip=rbind(typeIjan9,EXP3_typeIJAN300)
table1_manuscrip
#save(table1_manuscrip,file="./simulation_data/table1_manuscrip.RData")


table_list_column=table1_manuscrip[,1:8]
table_part1_unlist=data.frame(matrix(unlist(table_list_column[,5:8]),nrow=15))
colnames(table_part1_unlist)=c("power","se","power_01","se_01")
table_part2=table1_manuscrip[,1:4]
table_1_final_180=cbind(table_part1_unlist,table_part2)
table_1_final_180
save(table_1_final_180,file="./simulation_data/table_1_final_180.RData")

#################################
####################################################
####################################################
######Jan 9, 2024
load("EXP3_r5000_cfda2jan.RData")
library(xtable)
xtable(final_table,digits = 4)
power_curve_final=final_table
save(power_curve_final,file="power_curve.RData")
#########
load("power_curve.RData")
# Define the conditions and replacement values
conditions <- c(6, 7,8, 9,10)

power_by_time=power_curve_final[power_curve_final$fl_choice %in% conditions,]

replacement_values <- c(0, 5,10,15,20)

# Use replace() to replace the names in the 'Names' column
power_by_time$fl_choice <- replace_all(power_by_time$fl_choice, 
                                       conditions, 
                                       replacement_values)
#########
save(power_by_time,file="power_by_time.RData")
####################################################
library(ggplot2)

power_by_time_plot=ggplot(power_by_time,
                                      aes(x=fl_choice,y=unlist(power),color=as.factor(num_subjects),shape=as.factor(num_subjects)))+
    geom_line() +
    facet_grid(. ~ paste0("Time Points: ",as.factor(num_timepoints)))+
    ylab("Power")+
    xlab(expression(~delta~": Deivation from 0"))+
    guides(color = guide_legend(title = "Subjects")) +
    theme(text = element_text(size = 20))  +
    theme(
        legend.position = c(.02, .95),
        legend.justification = c("left", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)
    )+
    theme(strip.text.x = element_blank())+
    annotate(geom="text", x=9, y=1, label=expression(~alpha~" =0.05"),
             color="black",size=8)
    #geom_text(x=3, y=0.6, label=expression(~alpha~" =0.05"))
power_by_time_plot
ggsave("power_by_time_plot.png")
power_by_time_new=power_by_time[,c(1,3:5,7)]

load("power_by_time.RData")
library(data.table)
power_by_time_new_long <- melt(setDT(power_by_time_new), id.vars = c("fl_choice", "num_subjects", "num_timepoints"),
                               variable.name = "power_both")
power_by_time_new_long

power_by_time_new_long_subset=power_by_time_new_long[power_by_time_new_long$num_timepoints==90,]
power_by_time_new_long_subset

power_by_time_new_long_subset$power_both_new <- c('alpha*" = 0.05"', 'alpha*" = 0.10"')[power_by_time_new_long_subset$power_both]
power_by_time_plot_new_sub=ggplot(power_by_time_new_long_subset,
                              aes(x=fl_choice,y=unlist(value),color=as.factor(num_subjects),shape=as.factor(num_subjects)))+
    geom_line() +
    facet_grid(.~power_both_new ,
               labeller = label_parsed)+
    ylab("Power")+
    #xlab(expression(~delta~"("~H[0]*~": "~tilde(F[l])*(t)~"=0 )"))+
    xlab(expression(~delta))+
    guides(color = guide_legend(title = "Subjects")) +
    theme(text = element_text(size = 20))  +
    theme(
        legend.position = c(.02, .95),
        legend.justification = c("left", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)
    )#+
#theme(strip.text.x = element_blank())
power_by_time_plot_new_sub
ggsave("power_by_time_plot_new_sub.png")

# New facet label names for dose variable
# power.labs <- c("alpha =0.05","alpha =0.10")
# names(power.labs) <- c("power", "power_01")
# time.labs=c("Time Points: 90", "Time Points: 180")
# names(time.labs) <- c("90", "180")
########################################################
########################################################
power_by_time_new_long$num_timepoints_new <- c('"Timepoints = 90"', '"Timepoints = 180"')[as.factor(power_by_time_new_long$num_timepoints)]
power_by_time_new_long$power_both_new <- c('alpha*" = 0.05"', 'alpha*" = 0.10"')[power_by_time_new_long$power_both]
##Feb7 power correct time points
save(power_by_time_new_long,file="./simulation_data/power_by_time_new_long.RData")
power_by_time_plot_new=ggplot(power_by_time_new_long,
                          aes(x=fl_choice,y=unlist(value),color=as.factor(num_subjects),shape=as.factor(num_subjects)))+
    geom_line() +
    facet_grid(as.factor(num_timepoints_new)~power_both_new,
             labeller = label_parsed)+
    ylab("Power")+
    xlab(expression(~delta))+
    guides(color = guide_legend(title = "Subjects")) +
    theme(text = element_text(size = 20))  +
    theme(
        legend.position = c(.02, .95),
        legend.justification = c("left", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)
    )#+
    #theme(strip.text.x = element_blank())
power_by_time_plot_new

########
power_by_time_plot_01=ggplot(power_by_time,
                          aes(x=fl_choice,y=unlist(power_01),color=as.factor(num_subjects),shape=as.factor(num_subjects)))+
    geom_line() +
    facet_grid(. ~ paste0("Time Points: ",as.factor(num_timepoints)))+
    ylab("Power")+
    xlab(expression(~delta~": Deivation from 0"))+
    guides(color = guide_legend(title = "Subjects")) +
    theme(text = element_text(size = 20))  +
    theme(
        legend.position = c(.02, .95),
        legend.justification = c("left", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)
    )
power_by_time_plot_01
ggsave("power_by_time_plot_01.png")
###############################################
conditions_more <- c(6, 7,8, 9,10,14,15)
power_by_test_type= power_curve_final[!power_curve_final$fl_choice %in% conditions_more,]
power_by_test_type

# Define the conditions and replacement values
conditions_rep=c(21,22,23,24,25)
replacement_values_fl <- c(0,10, 20,30,40)

# Use replace() to replace the names in the 'Names' column
power_by_test_type$fl_choice <- replace_all(power_by_test_type$fl_choice, 
                                          conditions_rep,
                                          replacement_values_fl)
####graph the power
#power_by_test_type$test_type <- c('"Test Type"*beta*" (t)= c"', '"Test Type"*beta*" (t)= 0"')[power_by_test_type$test_type]

save(power_by_test_type,file="power_by_test_type.RData")

load("power_by_test_type.RData")
power_by_test_type_plot=ggplot(power_by_test_type,
                                    aes(x=fl_choice,y=unlist(power),
                                        color=as.factor(num_subjects),
                                        shape=as.factor(num_subjects)))+
    geom_line() +
    facet_grid(. ~ paste0("Test Type: ",test_type))+
    #facet_grid(. ~ test_type,labeller = label_parsed)+
    ylab("Power")+
    xlab(expression(~1~"+"~delta~"t"))+
    guides(color = guide_legend(title = "Subjects")) +
    theme(text = element_text(size = 20))  +
    theme(
        legend.position = c(.02, .95),
        legend.justification = c("left", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)
    )
power_by_test_type_plot
ggsave("power_by_test_type_plot.png")
#######################################################
power_by_test_type_subset=power_by_test_type[,c(1:3,5,7)]
power_by_test_type_subset
power_by_test_type_subset_long <- melt(setDT(power_by_test_type_subset), 
                            id.vars = c("fl_choice", "num_subjects", "test_type"),
                               variable.name = "power_both")
power_by_test_type_subset_long

power_by_test_type_subset_long$power_both_new <- c('alpha*" = 0.05"', 'alpha*" = 0.10"')[power_by_test_type_subset_long$power_both]
power_by_test_type_subset_long$test_type_new <- c('tilde(F[l])(t)*" = 0"', 'tilde(F[l])(t)*" = "*c')[power_by_test_type_subset_long$test_type]


power_by_test_type_subset_plot_new=ggplot(power_by_test_type_subset_long,
                              aes(x=fl_choice,y=unlist(value),color=as.factor(num_subjects),shape=as.factor(num_subjects)))+
    geom_line() +
    facet_grid(  test_type_new~power_both_new,
               labeller = label_parsed)+
    ylab("Power")+
    #xlab(expression(~delta~": Deivation from 0"))+
    xlab(expression(~1~"+"~delta~"t"))+
    guides(color = guide_legend(title = "Subjects")) +
    theme(text = element_text(size = 20))  +
    theme(
        legend.position = c(.02, .95),
        legend.justification = c("left", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)
    )#+
#theme(strip.text.x = element_blank())
power_by_test_type_subset_plot_new
ggsave("power_by_test_type_subset_plot_new.png")

######################################################
#power 01
power_by_test_type_plot_01=ggplot(power_by_test_type,
                               aes(x=fl_choice,y=unlist(power_01),
                                   color=as.factor(num_subjects),
                                   shape=as.factor(num_subjects)))+
    geom_line() +
    facet_grid(. ~ paste0("Test Type: ",test_type))+
    #facet_grid(. ~ test_type,labeller = label_parsed)+
    ylab("Power")+
    xlab(expression(~1~"+"~delta~"t"))+
    guides(color = guide_legend(title = "Subjects")) +
    theme(text = element_text(size = 20))  +
    theme(
        legend.position = c(.02, .95),
        legend.justification = c("left", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)
    )
power_by_test_type_plot_01
ggsave("power_by_test_type_plot_01.png")
#############################
#typeI
load("EXP3_typeIJAN.RData")
xtable(final_table,digits = 4)
#######################
#type I 300
load("EXP3_r5000_cfda2jan2.RData")
xtable(final_table,digits = 4)
######################
load("EXP3_typeIJAN300.RData")
xtable(final_table,digits = 4)
##############
load("EXP3_typeIJAN3002.RData")
xtable(final_table,digits = 4)
#############
load("EXP3_typeIJAN1000.RData")
xtable(final_table,digits = 4)
###########################
load("EXP3_typeIJANFULL.RData")
library(xtable)
xtable(final_table,digits = 4)
check_100=final_table[final_table$num_subjects==100,]
check_300=final_table[final_table$num_subjects==300,]
check_500=final_table[final_table$num_subjects==500,]
check_700=final_table[final_table$num_subjects==700,]


load("EXP3_typeIJANFULLpower.RData")
#library(xtable)
xtable(final_table,digits = 4)
check_100p=final_table[final_table$num_subjects==100,]
check_300p=final_table[final_table$num_subjects==300,]
check_500p=final_table[final_table$num_subjects==500,]
check_700=final_table[final_table$num_subjects==700,]
check_1000=final_table[final_table$num_subjects==1000,]

load("EXP3_typeIJAN22.RData")
library(xtable)
xtable(final_table,digits = 4)

load("EXP3_r5000_outputs_gam.RData")
xtable(final_table,digits = 4)
#######################
load("EXP3_r5000_outputs_gam_power.RData")
xtable(final_table,digits = 4)
################gam power Jan 31
power_curve_final=final_table

conditions <- c(6, 7,8, 9,10)

power_by_time=power_curve_final[(power_curve_final$fl_choice %in% conditions & power_curve_final$num_subjects!=300),]

replacement_values <- c(0, 5,10,15,20)

# Use replace() to replace the names in the 'Names' column
power_by_time$fl_choice <- replace_all(power_by_time$fl_choice, 
                                       conditions, 
                                       replacement_values)

power_by_time_new=power_by_time[,c(1,3:5,7)]


library(data.table)
power_by_time_new_long <- melt(setDT(power_by_time_new), id.vars = c("fl_choice", "num_subjects", "num_timepoints"),
                               variable.name = "power_both")
#power_by_time_new_long
#######
power_by_time_new_long$num_timepoints_new <- c('"Timepoints = 90"', '"Timepoints = 180"')[as.factor(power_by_time_new_long$num_timepoints)]
power_by_time_new_long$power_both_new <- c('alpha*" = 0.05"', 'alpha*" = 0.10"')[power_by_time_new_long$power_both]

power_by_time_plot_new_gam=ggplot(power_by_time_new_long,
                              aes(x=fl_choice,y=unlist(value),color=as.factor(num_subjects),shape=as.factor(num_subjects)))+
    geom_line() +
    facet_grid(as.factor(num_timepoints_new)~power_both_new,
               labeller = label_parsed)+
    ylab("Power")+
    xlab(expression(~delta))+
    guides(color = guide_legend(title = "Subjects")) +
    theme(text = element_text(size = 20))  +
    theme(
        legend.position = c(.02, .98),
        legend.justification = c("left", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)
    )#+
#theme(strip.text.x = element_blank())
power_by_time_plot_new_gam
save(power_by_time_new_long,file="power_by_time_new_long.RData")
#################
conditions_more <- c(6, 7,8, 9,10,14,15)
power_by_test_type= power_curve_final[!power_curve_final$fl_choice %in% conditions_more &power_curve_final$num_timepoints==90,]
#power_by_test_type= power_curve_final[!power_curve_final$fl_choice %in% conditions_more &power_curve_final$num_timepoints==180,]
power_by_test_type

# Define the conditions and replacement values
conditions_rep=c(21,22,23,24,25)
replacement_values_fl <- c(0,10, 20,30,40)

# Use replace() to replace the names in the 'Names' column
power_by_test_type$fl_choice <- replace_all(power_by_test_type$fl_choice, 
                                            conditions_rep,
                                            replacement_values_fl)
power_by_test_type_subset=power_by_test_type[,c(1:3,5,7)]
power_by_test_type_subset
power_by_test_type_subset_long <- melt(setDT(power_by_test_type_subset), 
                                       id.vars = c("fl_choice", "num_subjects", "test_type"),
                                       variable.name = "power_both")
power_by_test_type_subset_long

power_by_test_type_subset_long$power_both_new <- c('alpha*" = 0.05"', 'alpha*" = 0.10"')[power_by_test_type_subset_long$power_both]
power_by_test_type_subset_long$test_type_new <- c('tilde(F[l])(t)*" = 0"', 'tilde(F[l])(t)*" = "*c')[power_by_test_type_subset_long$test_type]

save(power_by_test_type_subset_long,file="power_by_test_type_subset_long.RData")
power_by_test_type_subset_plot_new_gam=ggplot(power_by_test_type_subset_long[power_by_test_type_subset_long$num_subjects!=300,],
                                          aes(x=fl_choice,y=unlist(value),color=as.factor(num_subjects),shape=as.factor(num_subjects)))+
    geom_line() +
    facet_grid(  test_type_new~power_both_new,
                 labeller = label_parsed)+
    ylab("Power")+
    #xlab(expression(~delta~": Deivation from 0"))+
    xlab(expression(~1~"+"~delta~"t"))+
    guides(color = guide_legend(title = "Subjects")) +
    theme(text = element_text(size = 20))  +
    theme(
        legend.position = c(.02, .95),
        legend.justification = c("left", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)
    )#+
#theme(strip.text.x = element_blank())
power_by_test_type_subset_plot_new_gam

#######################################
#L2 penalty start
########################################
#########no gam l2 penalty type I error
load("EXP3_r5000_1febnogaml2penalty.RData")
# xtable(final_table,digits = 4)
# 
# febnogaml2penalty=final_table
# save(final_table,file="./simulation_data/febnogaml2penalty.RData")

table_list_column=final_table[,1:8]
table_part1_unlist=data.frame(matrix(unlist(table_list_column[,5:8]),nrow=30))
colnames(table_part1_unlist)=c("power","se","power_01","se_01")
table_part2=final_table[,1:4]
nogam_l2penalty=cbind(table_part1_unlist,table_part2)
nogam_l2penalty
#xtable(final_table,digits=4)
#typeI90and180RData=final_table
nogam_l2penaltytypeI90=nogam_l2penalty[nogam_l2penalty$num_timepoints==90,]
nogam_l2penaltytypeI180=nogam_l2penalty[nogam_l2penalty$num_timepoints==180,]

save(nogam_l2penaltytypeI90,nogam_l2penaltytypeI180,file="./simulation_data/nogaml2penaltytypeI.RData")
###############gam l2 penalty type I error
load("EXP3_r5000_outputs_gam_feb2.RData")
#xtable(final_table,digits = 4)

# febgaml2penalty=final_table
# save(final_table,file="./simulation_data/febgaml2penalty.RData")

table_list_column=final_table[,1:8]
table_part1_unlist=data.frame(matrix(unlist(table_list_column[,5:8]),nrow=30))
colnames(table_part1_unlist)=c("power","se","power_01","se_01")
table_part2=final_table[,1:4]
gam_l2penalty=cbind(table_part1_unlist,table_part2)
gam_l2penalty
#xtable(final_table,digits=4)
#typeI90and180RData=final_table
gam_l2penaltytypeI90=gam_l2penalty[gam_l2penalty$num_timepoints==90,]
gam_l2penaltytypeI180=gam_l2penalty[gam_l2penalty$num_timepoints==180,]

save(gam_l2penaltytypeI90,gam_l2penaltytypeI180,file="./simulation_data/gaml2penaltytypeI.RData")
#############no gam l2 penalty power
#power by time points
load("EXP3_r5000_cfdanongampower.RData")
power_curve_final=final_table

conditions <- c(6, 7,8, 9,10)

power_by_time=power_curve_final[power_curve_final$fl_choice %in% conditions,]

replacement_values <- c(0, 5,10,15,20)

# Use replace() to replace the names in the 'Names' column
power_by_time$fl_choice <- replace_all(power_by_time$fl_choice, 
                                       conditions, 
                                       replacement_values)

power_by_time_new=power_by_time[,c(1,3:5,7)]


library(data.table)
power_by_time_new_long <- melt(setDT(power_by_time_new), id.vars = c("fl_choice", "num_subjects", "num_timepoints"),
                               variable.name = "power_both")
#power_by_time_new_long
#######
power_by_time_new_long$num_timepoints_new <- c('"Timepoints = 90"', '"Timepoints = 180"')[as.factor(power_by_time_new_long$num_timepoints)]
power_by_time_new_long$power_both_new <- c('alpha*" = 0.05"', 'alpha*" = 0.10"')[power_by_time_new_long$power_both]

power_by_time_plot_new_gam=ggplot(power_by_time_new_long,
                                  aes(x=fl_choice,y=unlist(value),color=as.factor(num_subjects),shape=as.factor(num_subjects)))+
    geom_line() +
    facet_grid(as.factor(num_timepoints_new)~power_both_new,
               labeller = label_parsed)+
    ylab("Power")+
    xlab(expression(~delta))+
    guides(color = guide_legend(title = "Subjects")) +
    theme(text = element_text(size = 20))  +
    theme(
        legend.position = c(.02, .98),
        legend.justification = c("left", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)
    )#+
#theme(strip.text.x = element_blank())
power_by_time_plot_new_gam
#power by test type
conditions_more <- c(6, 7,8, 9,10,14,15)
power_by_test_type= power_curve_final[!power_curve_final$fl_choice %in% conditions_more &power_curve_final$num_timepoints==90,]
#power_by_test_type= power_curve_final[!power_curve_final$fl_choice %in% conditions_more &power_curve_final$num_timepoints==180,]
power_by_test_type

# Define the conditions and replacement values
conditions_rep=c(21,22,23,24,25)
replacement_values_fl <- c(0,10, 20,30,40)

# Use replace() to replace the names in the 'Names' column
power_by_test_type$fl_choice <- replace_all(power_by_test_type$fl_choice, 
                                            conditions_rep,
                                            replacement_values_fl)
power_by_test_type_subset=power_by_test_type[,c(1:3,5,7)]
power_by_test_type_subset
power_by_test_type_subset_long <- melt(setDT(power_by_test_type_subset), 
                                       id.vars = c("fl_choice", "num_subjects", "test_type"),
                                       variable.name = "power_both")
power_by_test_type_subset_long

power_by_test_type_subset_long$power_both_new <- c('alpha*" = 0.05"', 'alpha*" = 0.10"')[power_by_test_type_subset_long$power_both]
power_by_test_type_subset_long$test_type_new <- c('tilde(F[l])(t)*" = 0"', 'tilde(F[l])(t)*" = "*c')[power_by_test_type_subset_long$test_type]

#save(power_by_test_type_subset_long,file="power_by_test_type_subset_long.RData")
power_by_test_type_subset_plot_new_gam=ggplot(power_by_test_type_subset_long,
                                              aes(x=fl_choice,y=unlist(value),color=as.factor(num_subjects),shape=as.factor(num_subjects)))+
    geom_line() +
    facet_grid(  test_type_new~power_both_new,
                 labeller = label_parsed)+
    ylab("Power")+
    #xlab(expression(~delta~": Deivation from 0"))+
    xlab(expression(~1~"+"~delta~"t"))+
    guides(color = guide_legend(title = "Subjects")) +
    theme(text = element_text(size = 20))  +
    theme(
        legend.position = c(.02, .95),
        legend.justification = c("left", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)
    )#+
#theme(strip.text.x = element_blank())
power_by_test_type_subset_plot_new_gam
############# gam l2 penalty power
load("EXP3_r5000_cfdagampower.RData")
##power by timepoints
power_curve_final=final_table

conditions <- c(6, 7,8, 9,10)

power_by_time=power_curve_final[power_curve_final$fl_choice %in% conditions,]

replacement_values <- c(0, 5,10,15,20)

# Use replace() to replace the names in the 'Names' column
power_by_time$fl_choice <- replace_all(power_by_time$fl_choice, 
                                       conditions, 
                                       replacement_values)

power_by_time_new=power_by_time[,c(1,3:5,7)]


library(data.table)
power_by_time_new_long <- melt(setDT(power_by_time_new), id.vars = c("fl_choice", "num_subjects", "num_timepoints"),
                               variable.name = "power_both")
#power_by_time_new_long
#######
power_by_time_new_long$num_timepoints_new <- c('"Timepoints = 90"', '"Timepoints = 180"')[as.factor(power_by_time_new_long$num_timepoints)]
power_by_time_new_long$power_both_new <- c('alpha*" = 0.05"', 'alpha*" = 0.10"')[power_by_time_new_long$power_both]

power_by_time_plot_new_gam=ggplot(power_by_time_new_long,
                                  aes(x=fl_choice,y=unlist(value),color=as.factor(num_subjects),shape=as.factor(num_subjects)))+
    geom_line() +
    facet_grid(as.factor(num_timepoints_new)~power_both_new,
               labeller = label_parsed)+
    ylab("Power")+
    xlab(expression(~delta))+
    guides(color = guide_legend(title = "Subjects")) +
    theme(text = element_text(size = 20))  +
    theme(
        legend.position = c(.02, .98),
        legend.justification = c("left", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)
    )#+
#theme(strip.text.x = element_blank())
power_by_time_plot_new_gam
#power by test type
conditions_more <- c(6, 7,8, 9,10,14,15)
power_by_test_type= power_curve_final[!power_curve_final$fl_choice %in% conditions_more &power_curve_final$num_timepoints==90,]
#power_by_test_type= power_curve_final[!power_curve_final$fl_choice %in% conditions_more &power_curve_final$num_timepoints==180,]
power_by_test_type

# Define the conditions and replacement values
conditions_rep=c(21,22,23,24,25)
replacement_values_fl <- c(0,10, 20,30,40)

# Use replace() to replace the names in the 'Names' column
power_by_test_type$fl_choice <- replace_all(power_by_test_type$fl_choice, 
                                            conditions_rep,
                                            replacement_values_fl)
power_by_test_type_subset=power_by_test_type[,c(1:3,5,7)]
power_by_test_type_subset
power_by_test_type_subset_long <- melt(setDT(power_by_test_type_subset), 
                                       id.vars = c("fl_choice", "num_subjects", "test_type"),
                                       variable.name = "power_both")
power_by_test_type_subset_long

power_by_test_type_subset_long$power_both_new <- c('alpha*" = 0.05"', 'alpha*" = 0.10"')[power_by_test_type_subset_long$power_both]
power_by_test_type_subset_long$test_type_new <- c('tilde(F[l])(t)*" = 0"', 'tilde(F[l])(t)*" = "*c')[power_by_test_type_subset_long$test_type]

#save(power_by_test_type_subset_long,file="power_by_test_type_subset_long.RData")
power_by_test_type_subset_plot_new_gam=ggplot(power_by_test_type_subset_long,
                                              aes(x=fl_choice,y=unlist(value),color=as.factor(num_subjects),shape=as.factor(num_subjects)))+
    geom_line() +
    facet_grid(  test_type_new~power_both_new,
                 labeller = label_parsed)+
    ylab("Power")+
    #xlab(expression(~delta~": Deivation from 0"))+
    xlab(expression(~1~"+"~delta~"t"))+
    guides(color = guide_legend(title = "Subjects")) +
    theme(text = element_text(size = 20))  +
    theme(
        legend.position = c(.02, .95),
        legend.justification = c("left", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)
    )#+
#theme(strip.text.x = element_blank())
power_by_test_type_subset_plot_new_gam
#######################################
#L2 penalty end
########################################

#######################################
#Feb6, re do aRLRT, difference_penalty start, type I error
########################################
load("EXP3_aRLRTFeb5.RData")
table_list_column=final_table[,1:8]
table_part1_unlist=data.frame(matrix(unlist(table_list_column[,5:8]),nrow=30))
colnames(table_part1_unlist)=c("power","se","power_01","se_01")
table_part2=final_table[,1:4]
table_1_final=cbind(table_part1_unlist,table_part2)
#xtable(final_table,digits=4)
#typeI90and180RData=final_table
typeI90=table_1_final[table_1_final$num_timepoints==90,]
typeI180=table_1_final[table_1_final$num_timepoints==180,]
save(typeI90,typeI180,
     file="./simulation_data/aRLRT_typeI.RData")
#######################################
#aRLRT, difference_penalty end, type I error
########################################

#######################################
#Feb7, re do gam, difference_penalty start, type I error, bernoulli 1-y instead of n-y
########################################
load("./simulation_data/EXP3_r5000_differencepenalty_gam_bernoulli1.RData")
gam_differece_penalty=final_table[,1:8]
gam_differece_penalty90_typeI=gam_differece_penalty[gam_differece_penalty$num_timepoints==90,]
gam_differece_penalty180_typeI=gam_differece_penalty[gam_differece_penalty$num_timepoints==180,]
save(gam_differece_penalty90_typeI,gam_differece_penalty180_typeI,
     file="./simulation_data/gam_differece_penalty_typeI_ber1notn.RData")
load("./simulation_data/gam_differece_penalty_typeI_ber1notn.RData")
#######################################
#Feb7, re do gam, difference_penalty start, type I error, bernoulli 1-y instead of n-y
########################################


#######################################
#Feb7, re do gam, l2_penalty start, type I error, bernoulli 1-y instead of n-y
########################################
load("./simulation_data/EXP3_r5000_linear_gam_real.RData")
gam_l2_penalty=final_table[,1:8]
gam_l2_penalty90_typeI=gam_l2_penalty[gam_l2_penalty$num_timepoints==90,]
gam_l2_penalty180_typeI=gam_l2_penalty[gam_l2_penalty$num_timepoints==180,]
save(gam_l2_penalty90_typeI,gam_l2_penalty180_typeI,
     file="./simulation_data/gam_l2_penalty_typeI_ber1notn.RData")

load("./simulation_data/gam_l2_penalty_typeI_ber1notn.RData")
#######################################
#Feb7, re do gam, l2_penalty end, type I error, bernoulli 1-y instead of n-y
########################################


#######################################
#Feb8, re do gam, l2_penalty start, power, bernoulli 1-y instead of n-y
########################################
load("EXP3_r5000_linear_gam_real_power.RData")

#march 26th, after the type I error
#######################
load("EXP3_r5000_gampower90180inclusion.RData")
###############################
power_curve_gam_l2=final_table

conditions <- c(6, 7,8, 9,10)

power_by_time_gaml2=power_curve_gam_l2[power_curve_gam_l2$fl_choice %in% conditions,]

replacement_values <- c(0, 5,10,15,20)

# Use replace() to replace the names in the 'Names' column
power_by_time_gaml2$fl_choice <- replace_all(power_by_time_gaml2$fl_choice, 
                                       conditions, 
                                       replacement_values)

power_by_time_gaml2data=power_by_time_gaml2[,c(1,3:5,7)]


library(data.table)
power_by_time_gaml2long <- melt(setDT(power_by_time_gaml2data), id.vars = c("fl_choice", "num_subjects", "num_timepoints"),
                               variable.name = "power_both")
#power_by_time_new_long
#######
power_by_time_gaml2long$num_timepoints_new <- c('"Timepoints = 90"', '"Timepoints = 180"')[as.factor(power_by_time_gaml2long$num_timepoints)]
power_by_time_gaml2long$power_both_new <- c('alpha*" = 0.05"', 'alpha*" = 0.10"')[power_by_time_gaml2long$power_both]

save(power_by_time_gaml2long,file="./simulation_data/power_by_time_gaml2long.RData")
power_by_time_plot_gaml2=ggplot(power_by_time_gaml2long,
                                  aes(x=fl_choice,y=unlist(value),color=as.factor(num_subjects),shape=as.factor(num_subjects)))+
    geom_line() +
    facet_grid(as.factor(num_timepoints_new)~power_both_new,
               labeller = label_parsed)+
    ylab("Power")+
    xlab(expression(~delta))+
    guides(color = guide_legend(title = "Subjects")) +
    theme(text = element_text(size = 20))  +
    theme(
        legend.position = c(.01, .99),
        legend.justification = c("left", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)
    )#+
#theme(strip.text.x = element_blank())
power_by_time_plot_gaml2
ggsave("./simulation_data/power_by_time_plot_gaml2.png")
######not include 300
power_by_time_plot_gaml2_no300=ggplot(power_by_time_gaml2long[power_by_time_gaml2long$num_subjects!=300],
                                aes(x=fl_choice,y=unlist(value),color=as.factor(num_subjects),shape=as.factor(num_subjects)))+
    geom_line() +
    facet_grid(as.factor(num_timepoints_new)~power_both_new,
               labeller = label_parsed)+
    ylab("Power")+
    xlab(expression(~delta))+
    guides(color = guide_legend(title = "Subjects")) +
    theme(text = element_text(size = 20))  +
    theme(
        legend.position = c(.01, .99),
        legend.justification = c("left", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)
    )#+
#theme(strip.text.x = element_blank())
power_by_time_plot_gaml2_no300

######
#power by test type
conditions_more <- c(6, 7,8, 9,10,14,15)
power_by_test_type_gaml2= power_curve_gam_l2[!power_curve_gam_l2$fl_choice %in% conditions_more &power_curve_gam_l2$num_timepoints==90,]
#power_by_test_type= power_curve_final[!power_curve_final$fl_choice %in% conditions_more &power_curve_final$num_timepoints==180,]


# Define the conditions and replacement values
####current run for only 90, fl 11-15
# load("EXP3_r5000_gam_power_l2_bigger_fl.RData")
# power_by_test_type_gaml2=final_table


#####################
#powerbytesttype march 26
#load("EXP3_r5000_gampowertesttypemarch26.RData")
#######################
conditions_rep=c(21,22,23,24,25)
replacement_values_fl <- c(0,10, 20,30,40)

# conditions_rep=c(11,12,13,14,15)
# replacement_values_fl <- c(0,20, 40,60,80)

# Use replace() to replace the names in the 'Names' column
power_by_test_type_gaml2$fl_choice <- replace_all(power_by_test_type_gaml2$fl_choice, 
                                            conditions_rep,
                                            replacement_values_fl)
power_by_test_type_subset_gaml2=power_by_test_type_gaml2[,c(1:3,5,7)]


power_by_test_type_subset_long_gaml2 <- melt(setDT(power_by_test_type_subset_gaml2), 
                                       id.vars = c("fl_choice", "num_subjects", "test_type"),
                                       variable.name = "power_both")

power_by_test_type_subset_long_gaml2$power_both_new <- c('alpha*" = 0.05"', 'alpha*" = 0.10"')[power_by_test_type_subset_long_gaml2 $power_both]
power_by_test_type_subset_long_gaml2$test_type_new <- c('tilde(F[l])(t)*" = 0"', 'tilde(F[l])(t)*" = "*c')[power_by_test_type_subset_long_gaml2 $test_type]

save(power_by_test_type_subset_long_gaml2,file="./simulation_data/power_by_test_type_subset_long_gaml2.RData")
power_by_test_type_subset_plot_gaml2=ggplot(power_by_test_type_subset_long_gaml2,
                                              aes(x=fl_choice,y=unlist(value),color=as.factor(num_subjects),shape=as.factor(num_subjects)))+
    geom_line() +
    facet_grid(  test_type_new~power_both_new,
                 labeller = label_parsed)+
    ylab("Power")+
    #xlab(expression(~delta~": Deivation from 0"))+
    xlab(expression(~1~"+"~delta~"t"))+
    guides(color = guide_legend(title = "Subjects")) +
    theme(text = element_text(size = 20))  +
    theme(
        legend.position = c(.02, .95),
        legend.justification = c("left", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)
    )#+
#theme(strip.text.x = element_blank())
power_by_test_type_subset_plot_gaml2

########################################


##fl11-15
load("EXP3_r5000_gam_power_l2_bigger_fl.RData")
power_by_test_type_gaml2_bigfl=final_table

conditions_rep_newfl=c(11,12,13,14,15)
replacement_values_flnewfl <- c(0,20, 40,60,80)

# Use replace() to replace the names in the 'Names' column
power_by_test_type_gaml2_bigfl$fl_choice <- replace_all(power_by_test_type_gaml2_bigfl$fl_choice, 
                                                        conditions_rep_newfl,
                                                        replacement_values_flnewfl)
power_by_test_type_subset_gaml2bigfl=power_by_test_type_gaml2_bigfl[,c(1:3,5,7)]


power_by_test_type_subset_long_gaml2bigfl <- melt(setDT(power_by_test_type_subset_gaml2bigfl), 
                                             id.vars = c("fl_choice", "num_subjects", "test_type"),
                                             variable.name = "power_both")

power_by_test_type_subset_long_gaml2bigfl$power_both_new <- c('alpha*" = 0.05"', 'alpha*" = 0.10"')[power_by_test_type_subset_long_gaml2bigfl$power_both]
power_by_test_type_subset_long_gaml2bigfl$test_type_new <- c('tilde(F[l])(t)*" = 0"', 'tilde(F[l])(t)*" = "*c')[power_by_test_type_subset_long_gaml2bigfl$test_type]

save(power_by_test_type_subset_long_gaml2bigfl,file="./simulation_data/power_by_test_type_subset_long_gaml2bigfl.RData")

power_by_test_type_subset_plot_gaml2bigfl=ggplot(power_by_test_type_subset_long_gaml2bigfl,
                                            aes(x=fl_choice,y=unlist(value),color=as.factor(num_subjects),shape=as.factor(num_subjects)))+
    geom_line() +
    facet_grid(  test_type_new~power_both_new,
                 labeller = label_parsed)+
    ylab("Power")+
    #xlab(expression(~delta~": Deivation from 0"))+
    xlab(expression(~1~"+"~delta~"t"))+
    guides(color = guide_legend(title = "Subjects")) +
    theme(text = element_text(size = 20))  +
    theme(
        legend.position = c(.01, .98),
        legend.justification = c("left", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)
    )#+
#theme(strip.text.x = element_blank())
power_by_test_type_subset_plot_gaml2bigfl

###########use fl 0,10,20,30,40,60
load("./simulation_data/power_by_time_gaml2long.RData")
#power_by_test_type_subset_long_gaml2
power_by_test_type_subset_long_gaml2_sub=power_by_test_type_subset_long_gaml2[power_by_test_type_subset_long_gaml2$fl_choice %in% c(10,30) & power_by_test_type_subset_long_gaml2$num_subjects %in% c(1000,500,100),]
power_by_test_type_subset_long_gaml2bigfl_sub=power_by_test_type_subset_long_gaml2bigfl[power_by_test_type_subset_long_gaml2bigfl$fl_choice %in% c(0,20,40,60),]
power_test_type_finaldata=rbind(power_by_test_type_subset_long_gaml2_sub,power_by_test_type_subset_long_gaml2bigfl_sub)

save(power_test_type_finaldata,file="./simulation_data/power_test_type_finaldata.RData")
power_by_test_type_subset_plot_gaml2mixfl=ggplot(power_test_type_finaldata,
                                                 aes(x=fl_choice,y=unlist(value),color=as.factor(num_subjects),shape=as.factor(num_subjects)))+
    geom_line() +
    facet_grid(  test_type_new~power_both_new,
                 labeller = label_parsed)+
    ylab("Power")+
    #xlab(expression(~delta~": Deivation from 0"))+
    xlab(expression(~1~"+"~delta~"t"))+
    guides(color = guide_legend(title = "Subjects")) +
    theme(text = element_text(size = 20))  +
    theme(
        legend.position = c(.01, .98),
        legend.justification = c("left", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)
    )#+
#theme(strip.text.x = element_blank())
power_by_test_type_subset_plot_gaml2mixfl

#######################################
#Feb7, re do gam, l2_penalty end, power, bernoulli 1-y instead of n-y
########################################


#######################################
#Feb12, re do gam, l2_penalty start, power, bernoulli 1-y instead of n-y
########################################
load("EXP3_r5000_nonlinear_gam_power.RData")

#######################################
#March17, re do gam, l2_penalty end, power, bernoulli 1-y instead of n-y, type I
########################################
load("EXP3_r5000_cfdagamTypeImarch.RData")
march7_gam_typeI_90=final_table[final_table$num_timepoints==90,1:8]

march7_gam_typeI_90=march7_gam_typeI_90[!rownames(march7_gam_typeI_90)%in% c(
    "experiment_output.6","experiment_output.7","experiment_output.8",
    "experiment_output.18","experiment_output.19","experiment_output.20",
    "experiment_output.24","experiment_output.25","experiment_output.26",
    "experiment_output.30","experiment_output.31","experiment_output.32",
    "experiment_output.36","experiment_output.37","experiment_output.38"),]
march7_gam_typeI_90=march7_gam_typeI_90[,-4]
save(march7_gam_typeI_90,file="march7_gam_typeI_90.RData")

load("march7_gam_typeI_90.RData")
write.csv(march7_gam_typeI_90, "gam_90.csv",row.names = FALSE)
#march7_gam_typeI_180=final_table[final_table$num_timepoints==180,1:8]

####################################
#Final run for type I error gam L2 penalty March 23
###################################
load("EXP3_r5000_cfdagamtyieIL2March.RData")
RLRT_final_gam_type_I=final_table[,1:8]
RLRT_final_gam_type_I$power=unlist(RLRT_final_gam_type_I$power)
RLRT_final_gam_type_I$se=unlist(RLRT_final_gam_type_I$se)
RLRT_final_gam_type_I$power_01=unlist(RLRT_final_gam_type_I$power_01)
RLRT_final_gam_type_I$se01=unlist(RLRT_final_gam_type_I$se01)
write.csv(RLRT_final_gam_type_I, file = "RLRT_final_gam_type_I.csv")

#write a function to create csv file
load("EXP3_r5000_cfdagamtyieIL2March.RData")
gam_typeI_csv=function(final_table,number_time_points){
    RLRT_final_gam_type_I=final_table[,1:8]
    RLRT_final_gam_type_I$power=unlist(RLRT_final_gam_type_I$power)
    RLRT_final_gam_type_I$se=unlist(RLRT_final_gam_type_I$se)
    RLRT_final_gam_type_I$power_01=unlist(RLRT_final_gam_type_I$power_01)
    RLRT_final_gam_type_I$se01=unlist(RLRT_final_gam_type_I$se01)
    write.csv(RLRT_final_gam_type_I, file = paste0("RLRT_final_gam_type_I_",number_time_points,".csv") )
}
gam_typeI_csv(final_table,90)
###################################time points 180
load("EXP3_r5000_cfdagam180.RData")
gam_typeI_csv(final_table,180)
#user python code to write the table
######################################
#function graph power
replace_all <- function(old_vector, old_values, new_values){
    c(new_values, old_vector)[match(old_vector, c(old_values, old_vector))]
}
library(data.table)
library(ggplot2)

power_by_time_points_plot=function(final_table){
    power_curve_final=final_table
    conditions <- c(6, 7,8, 9,10)
    
    power_by_time=power_curve_final[power_curve_final$fl_choice %in% conditions,]
    
    replacement_values <- c(0, 5,10,15,20)
    
    # Use replace() to replace the names in the 'Names' column
    power_by_time$fl_choice <- replace_all(power_by_time$fl_choice, 
                                           conditions, 
                                           replacement_values)
    
    power_by_time_new=power_by_time[,c(1,3:5,7)]

    power_by_time_new_long <- melt(setDT(power_by_time_new), id.vars = c("fl_choice", "num_subjects", "num_timepoints"),
                                   variable.name = "power_both")
    #power_by_time_new_long
    #######
    power_by_time_new_long$num_timepoints_new <- c('"Timepoints = 90"', '"Timepoints = 180"')[as.factor(power_by_time_new_long$num_timepoints)]
    power_by_time_new_long$power_both_new <- c('alpha*" = 0.05"', 'alpha*" = 0.10"')[power_by_time_new_long$power_both]
    
    power_by_time_plot_new_gam=ggplot(power_by_time_new_long,
                                      aes(x=fl_choice,y=unlist(value),color=as.factor(num_subjects),shape=as.factor(num_subjects)))+
        geom_line() +
        facet_grid(as.factor(num_timepoints_new)~power_both_new,
                   labeller = label_parsed)+
        ylab("Power")+
        xlab(expression(~delta))+
        guides(color = guide_legend(title = "Subjects")) +
        theme(text = element_text(size = 20))  +
        theme(
            legend.position = c(.02, .98),
            legend.justification = c("left", "top"),
            legend.box.just = "right",
            legend.margin = margin(6, 6, 6, 6)
        )#+
    #theme(strip.text.x = element_blank())
    print(power_by_time_plot_new_gam)
}
load("EXP3_r5000_gampower90180inclusion.RData")
power_by_time_points_plot(final_table)

power_by_test_type_plot=function(final_table){
    power_by_test_type_gaml2=final_table
    conditions_rep=c(21,22,23,24,25,14)
    replacement_values_fl <- c(0,10, 20,30,40,60)
    
    # conditions_rep=c(11,12,13,14,15)
    # replacement_values_fl <- c(0,20, 40,60,80)
    
    # Use replace() to replace the names in the 'Names' column
    power_by_test_type_gaml2$fl_choice <- replace_all(power_by_test_type_gaml2$fl_choice, 
                                                      conditions_rep,
                                                      replacement_values_fl)
    power_by_test_type_subset_gaml2=power_by_test_type_gaml2[,c(1:3,5,7)]
    
    
    power_by_test_type_subset_long_gaml2 <- melt(setDT(power_by_test_type_subset_gaml2), 
                                                 id.vars = c("fl_choice", "num_subjects", "test_type"),
                                                 variable.name = "power_both")
    
    power_by_test_type_subset_long_gaml2$power_both_new <- c('alpha*" = 0.05"', 'alpha*" = 0.10"')[power_by_test_type_subset_long_gaml2 $power_both]
    power_by_test_type_subset_long_gaml2$test_type_new <- c('tilde(F[l])(t)*" = 0"', 'tilde(F[l])(t)*" = "*c')[power_by_test_type_subset_long_gaml2 $test_type]
    
    #save(power_by_test_type_subset_long_gaml2,file="./simulation_data/power_by_test_type_subset_long_gaml2.RData")
    power_by_test_type_subset_plot_gaml2=ggplot(power_by_test_type_subset_long_gaml2,
                                                aes(x=fl_choice,y=unlist(value),color=as.factor(num_subjects),shape=as.factor(num_subjects)))+
        geom_line() +
        facet_grid(  test_type_new~power_both_new,
                     labeller = label_parsed)+
        ylab("Power")+
        #xlab(expression(~delta~": Deivation from 0"))+
        xlab(expression(~1~"+"~delta~"t"))+
        guides(color = guide_legend(title = "Subjects")) +
        theme(text = element_text(size = 20))  +
        theme(
            legend.position = c(.02, .98),
            legend.justification = c("left", "top"),
            legend.box.just = "right",
            legend.margin = margin(6, 6, 6, 6)
        )#+
    #theme(strip.text.x = element_blank())
    print(power_by_test_type_subset_plot_gaml2)
}

load("EXP3_r5000_gampowertesttypemarch26.RData")
power_by_test_type_plot(final_table)

##3/25. T testing
#write a function to get all 5 dataset
load("./hazel_final_table_output/hazel16_100_1000_1000/Hazel_outputsTbootstrap_1_16_100_1000_1000_.RData")
n100_1=final_table$power

load("./hazel_final_table_output/hazel16_100_1000_1000/Hazel_outputsTbootstrap_2_16_100_1000_1000_.RData")
n100_2=final_table$power

load("./hazel_final_table_output/hazel16_100_1000_1000/Hazel_outputsTbootstrap_3_16_100_1000_1000_.RData")
n100_3=final_table$power

load("./hazel_final_table_output/hazel16_100_1000_1000/Hazel_outputsTbootstrap_4_16_100_1000_1000_.RData")
n100_4=final_table$power

load("./hazel_final_table_output/hazel16_100_1000_1000/Hazel_outputsTbootstrap_5_16_100_1000_1000_.RData")
n100_5=final_table$power

n_100_p=mean(unlist(c(n100_1,n100_2,n100_3,n100_4,n100_5)))
n_100_p #0.049
#[1] 0.048

n_100_se=sqrt(n_100_p*(1-n_100_p)/5000)
n_100_se #0.003052835 [1] 0.003023111

#######
# load("/Users/xzhao17/Documents/GitHub/xcj_hypothesis_test_cfd/hazel_final_table_output/testrun/Hazel_outputsTbootstrap_2_16300100100.RData")
# mean(unlist(final_table$power))

load("./hazel_final_table_output/Hazel_outputsTbootstrap_4_1630010001000.RData")
n300_4=final_table$power
mean(unlist(n300_4))
#[1] 0.01
########
# load("/Users/xzhao17/Documents/GitHub/xcj_hypothesis_test_cfd/hazel_final_table_output/testrun/Hazel_outputsTbootstrap_2_32500100100.RData")
# mean(unlist(final_table$power))
# [1] 0.05

load("./hazel_final_table_output/Hazel_outputsTbootstrap_1_3250010001000.RData")
n500_1=final_table$power

load("./hazel_final_table_output/Hazel_outputsTbootstrap_2_3250010001000.RData")
n500_2=final_table$power

load("./hazel_final_table_output/Hazel_outputsTbootstrap_3_3250010001000.RData")
n500_3=final_table$power

load("./hazel_final_table_output/Hazel_outputsTbootstrap_4_3250010001000.RData")
n500_4=final_table$power

load("./hazel_final_table_output/Hazel_outputsTbootstrap_5_3250010001000.RData")
n500_5=final_table$power

n_500_p=mean(unlist(c(n500_1,n500_2,n500_3,n500_4,n500_5)))
n_500_p #0.015

n_500_se=sqrt(n_500_p*(1-n_500_p)/5000)
n_500_se #0.003052835