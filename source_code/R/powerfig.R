
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

# Use replace() to replace the names in the 'Names' column
n100n500t90t300deltapower$fl_choice <- replace(n100n500t90t300deltapower$fl_choice, 
                                               n100n500t90t300deltapower$fl_choice %in% conditions, replacement_values)

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
    xlab(expression(~delta~": Slope deviation from 0"))+
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


n100n500t90fl21fl25_plot=ggplot(n100n500t90fl21tofl25,
                                    aes(x=fl_choice,y=unlist(power),
                                        color=as.factor(num_subjects),
                                        shape=as.factor(num_subjects)))+
    geom_line() +
    facet_grid(. ~ paste0("Test Type: ",test_type))+
    ylab("Power")+
    xlab(expression(~delta~": Slope deviation from 0"))+
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






