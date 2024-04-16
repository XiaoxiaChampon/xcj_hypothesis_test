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
# Date: 4/01/2024
#
##############################################################

calculate_power_stats <- function(power_values,alpha_level){
    power_mean <- mean(power_values<alpha_level)
    power_se <- sqrt(power_mean * (1 - power_mean) / length(power_values))
    return(list(mean = power_mean, standard_error = power_se))
}

calculate_stats <- function(directory_path) {
    power_values <- c()
    T_variable=c()
    
    # power_values2 <- c()
    # power_01_values2 <- c()
    
    
    # Get a list of all RData files in the specified directory by pattern
    files <- list.files(path = directory_path, pattern = "\\.RData$", full.names = TRUE)
    
    for (file_path in files) {
        load(file_path)
        #staicu
        power_values <- c(power_values, final_table$power[[1]])
        T_values <- c(  T_variable, final_table$T_rv[[1]])
        
        
        # power_values <- c(power_values, final_table$power[[1]])
        # power_01_values <- c(power_01_values, final_table$power_01[[1]])
        
        ###add one more power
        # power_values2 <- c(power_values2, final_table$power2[[1]])
        # power_01_values2 <- c(power_01_values2, final_table$power_012[[1]])
        ################################
    }
    
    cat("Total Length (power):", length(power_values))
    cat("\nTotal Length (power_01):", length(power_values))
    
    power_values_0.05 <- calculate_power_stats(power_values,0.05)
   # power_01_values <- calculate_power_stats(power_01_values)
    power_01_values <- calculate_power_stats(power_values,0.1)
    
    ######################
    # power_values2 <- calculate_power_stats(power_values2)
    # power_01_values2 <- calculate_power_stats(power_01_values2)
    #####################
    # return(list(power=power_values, power_01=power_01_values,
    #             power2=power_values2, power_012=power_01_values2))
    
    return(list(power=power_values_0.05, power_01=power_01_values
                ))
}


calculate_stats_percentile <- function(directory_path) {
    power_values <- c()
    T_variable=c()
    power_values_01 <- c()
    # power_values2 <- c()
    # power_01_values2 <- c()
    
    
    # Get a list of all RData files in the specified directory by pattern
    files <- list.files(path = directory_path, pattern = "\\.RData$", full.names = TRUE)
    
    for (file_path in files) {
        load(file_path)
        #staicu
        power_values <- c(power_values, final_table$power[[1]])
        power_values_01 <- c(power_values, final_table$power_01[[1]])
        T_values <- c(  T_variable, final_table$T_rv[[1]])
        
        
        # power_values <- c(power_values, final_table$power[[1]])
        # power_01_values <- c(power_01_values, final_table$power_01[[1]])
        
        ###add one more power
        # power_values2 <- c(power_values2, final_table$power2[[1]])
        # power_01_values2 <- c(power_01_values2, final_table$power_012[[1]])
        ################################
    }
    
    cat("Total Length (power):", length(power_values))
    cat("\nTotal Length (power_01):", length(power_values))
    
    power_values_0.05 <- calculate_power_stats_percentile(power_values)
    # power_01_values <- calculate_power_stats(power_01_values)
    power_01_values <- calculate_power_stats_percentile(power_values_01)
    
    ######################
    # power_values2 <- calculate_power_stats(power_values2)
    # power_01_values2 <- calculate_power_stats(power_01_values2)
    #####################
    # return(list(power=power_values, power_01=power_01_values,
    #             power2=power_values2, power_012=power_01_values2))
    
    return(list(power=power_values_0.05, power_01=power_01_values
    ))
}
calculate_power_stats_percentile <- function(power_values){
    power_mean <- mean(power_values)
    power_se <- sqrt(power_mean * (1 - power_mean) / length(power_values))
    return(list(mean = power_mean, standard_error = power_se))
}

# Example usage
directory_path <- "./hazel_final_table_output/hazel16_500_250_1000/"
stats <- calculate_stats(directory_path)
print(stats)
# Total Length (power): 5000
# Total Length (power_01): 5000> print(stats)
# $power
# $power$mean
# [1] 0.02
# 
# $power$standard_error
# [1] 0.001979899
# 
# 
# $power_01
# $power_01$mean
# [1] 0.04
# 
# $power_01$standard_error
# [1] 0.002771281

directory_path <- "./hazel_final_table_output/hazel16_100_1000_1000/"
stats <- calculate_stats(directory_path)
print(stats)

# Total Length (power): 5000
# Total Length (power_01): 5000> print(stats)
# $power
# $power$mean
# [1] 0.045
# 
# $power$standard_error
# [1] 0.002931723
# 
# 
# $power_01
# $power_01$mean
# [1] 0.073
# 
# $power_01$standard_error
# [1] 0.003678886

directory_path <- "./hazel_final_table_output/hazel16_300_1000_1000/"
stats <- calculate_stats(directory_path)
print(stats)

# Total Length (power): 5000
# Total Length (power_01): 5000> print(stats)
# $power
# $power$mean
# [1] 0.01
# 
# $power$standard_error
# [1] 0.001407125
# 
# 
# $power_01
# $power_01$mean
# [1] 0.03
# 
# $power_01$standard_error
# [1] 0.002412468


directory_path <- "./hazel_final_table_output/hazel_power/power_100/fl7/"
stats <- calculate_stats(directory_path)
 print(stats)
# 
# Total Length (power): 1000
# Total Length (power_01): 1000> print(stats)
# $power
# $power$mean
# [1] 0.053
# 
# $power$standard_error
# [1] 0.007084561
# 
# 
# $power_01
# $power_01$mean
# [1] 0.085
# 
# $power_01$standard_error
# [1] 0.008819014

directory_path <- "./hazel_final_table_output/hazel_power/power_100/fl8/"
stats <- calculate_stats(directory_path)
print(stats)
# $power
# $power$mean
# [1] 0.107
# 
# $power$standard_error
# [1] 0.009775019
# 
# 
# $power_01
# $power_01$mean
# [1] 0.146
# 
# $power_01$standard_error
# [1] 0.0111662

directory_path <- "./hazel_final_table_output/hazel_power/power_100/fl9/"
stats <- calculate_stats(directory_path)
print(stats)

# $power
# $power$mean
# [1] 0.154
# 
# $power$standard_error
# [1] 0.0114142
# 
# 
# $power_01
# $power_01$mean
# [1] 0.243
# 
# $power_01$standard_error
# [1] 0.01356285


directory_path <- "./hazel_final_table_output/hazel_power/power_100/fl10/"
stats <- calculate_stats(directory_path)
print(stats)


directory_path <- "./hazel_final_table_output/hazel_power/power_300/fl7/"
stats <- calculate_stats(directory_path)
print(stats)

directory_path <- "./hazel_final_table_output/hazel_power/power_300/fl8/"
stats <- calculate_stats(directory_path)
print(stats)

directory_path <- "./hazel_final_table_output/hazel_power/power_300/fl9/"
stats <- calculate_stats(directory_path)
print(stats)

directory_path <- "./hazel_final_table_output/hazel_power/power_300/fl10/"
stats <- calculate_stats(directory_path)
print(stats)

###########
directory_path <- "./hazel_final_table_output/hazel16_100_1000_1000/"
stats <- calculate_stats(directory_path)
print(stats)

#######partial 1000 users results
directory_path <- "./hazel_final_table_output/"
stats <- calculate_stats(directory_path)
print(stats)

############no shuffle
directory_path <- "./hazel_final_table_output/hazel_100_no_shuffle"
stats <- calculate_stats(directory_path)
print(stats)

directory_path <- "./hazel_final_table_output/hazel_300_no_shuffle"
stats <- calculate_stats(directory_path)
print(stats)

directory_path <- "./hazel_final_table_output/hazel_500_no_shuffle"
stats <- calculate_stats(directory_path)
print(stats)

# directory_path <- "/Users/xzhao17/Desktop/"
# stats <- calculate_stats(directory_path)
# print(stats)
######April 7, 2024
directory_path <- "./hazel_final_table_output/hazel_300_200_1000/"
stats <- calculate_stats(directory_path)
print(stats)

# Total Length (power): 4800
# Total Length (power_01): 4800> print(stats)
# $power
# $power$mean
# [1] 0.01458333
# 
# $power$standard_error
# [1] 0.001730285
# 
# 
# $power_01
# $power_01$mean
# [1] 0.03354167
# 
# $power_01$standard_error
# [1] 0.002598743
# 
# 
# $power2
# $power2$mean
# [1] 0.014375
# 
# $power2$standard_error
# [1] 0.001718063
# 
# 
# $power_012
# $power_012$mean
# [1] 0.03291667
# 
# $power_012$standard_error
# [1] 0.002575249

directory_path <- "./hazel_final_table_output/hazel_500_100_1000/"
stats <- calculate_stats(directory_path)
print(stats)

# Total Length (power): 5000
# Total Length (power_01): 5000> print(stats)
# $power
# $power$mean
# [1] 0.0176
# 
# $power$standard_error
# [1] 0.001859583
# 
# 
# $power_01
# $power_01$mean
# [1] 0.044
# 
# $power_01$standard_error
# [1] 0.002900483
# 
# 
# $power2
# $power2$mean
# [1] 0.0178
# 
# $power2$standard_error
# [1] 0.001869928
# 
# 
# $power_012
# $power_012$mean
# [1] 0.0444
# 
# $power_012$standard_error
# [1] 0.002913027


directory_path <- "./hazel_final_table_output/hazel_100_250_1000/"
stats <- calculate_stats(directory_path)
print(stats)
# $power
# $power$mean
# [1] 0.0564
# 
# $power$standard_error
# [1] 0.003262485
# 
# 
# $power_01
# $power_01$mean
# [1] 0.0832
# 
# $power_01$standard_error
# [1] 0.003905836

directory_path <- "./hazel_final_table_output/hazel_1000_100_1000/"
stats <- calculate_stats(directory_path)
print(stats)
# $power
# $power$mean
# [1] 0.0242
# 
# $power$standard_error
# [1] 0.002173217
# 
# 
# $power_01
# $power_01$mean
# [1] 0.0628
# 
# $power_01$standard_error
# [1] 0.003430923

##staicu
directory_path <- "./hazel_final_table_output/final_table_output_p_staicu_n100/"
stats <- calculate_stats(directory_path)
print(stats)

# $power
# $power$mean
# [1] 0.055
# 
# $power$standard_error
# [1] 0.003
# 
# 
# $power_01
# $power_01$mean
# [1] 0.0802
# 
# $power_01$standard_error
# [1] 0.0038


load("./hazel_final_table_output/final_table_output_p_staicu_n100/Hazel_outputsTbootstrap__100_250_1000_16_1.RData")
directory_path <- "./hazel_final_table_output/final_table_output_p_staicu_n300/"
stats <- calculate_stats(directory_path)
print(stats)
# $power
# $power$mean
# [1] 0.0118
# 
# $power$standard_error
# [1] 0.021597
# 
# 
# $power_01
# $power_01$mean
# [1] 1
# 
# $power_01$standard_error
# [1] 0
directory_path <- "./hazel_final_table_output/final_table_output_p_staicu_n500/"
stats <- calculate_stats(directory_path)
print(stats)
# $power
# $power$mean
# [1] 0.0182
# 
# $power$standard_error
# [1] 0.01890437
# 
# 

# $power_01
# $power_01$mean
# [1] 0.0432
# 
# $power_01$standard_error
# [1] 0.002875196

directory_path <- "./hazel_final_table_output/final_table_output_p_staicu_n1000/"
stats <- calculate_stats(directory_path)
print(stats)
# $power
# $power$mean
# [1] 0.02265306
# 
# $power$standard_error
# [1] 0.00212564
# 
# 
# $power_01
# $power_01$mean
# [1] 0.05979592
# 
# $power_01$standard_error
# [1] 0.003387262

#load("/Users/xzhao17/Documents/GitHub/xcj_hypothesis_test_cfd/hazel_final_table_output/final_table_output_p_staicu_n2000/Hazel_outputsTbootstrap__2000_40_1000_16_1.RData")
##
directory_path <- "./hazel_final_table_output/final_table_output_p_staicu_n2000/"
stats <- calculate_stats(directory_path)
print(stats)


########staicu no shuffle
directory_path <- "./hazel_final_table_output/final_table_output_p_staicu_noshffule_n100/"
stats <- calculate_stats(directory_path)
print(stats)
# $power
# $power$mean
# [1] 0.0552
# 
# $power$standard_error
# [1] 0.003229643
# 
# 
# $power_01
# $power_01$mean
# [1] 0.0824
# 
# $power_01$standard_error
# [1] 0.003888708

directory_path <- "./hazel_final_table_output/final_table_output_p_staicu_noshffule_n300/"
stats <- calculate_stats(directory_path)
print(stats)

# $power
# $power$mean
# [1] 0.0136
# 
# $power$standard_error
# [1] 0.001637989
# 
# 
# $power_01
# $power_01$mean
# [1] 0.0318
# 
# $power_01$standard_error
# [1] 0.002481482

directory_path <- "./hazel_final_table_output/final_table_output_p_staicu_noshffule_n500/"
stats <- calculate_stats(directory_path)
print(stats)
# $power
# $power$mean
# [1] 0.0172
# 
# $power$standard_error
# [1] 0.001838704
# 
# 
# $power_01
# $power_01$mean
# [1] 0.0434
# 
# $power_01$standard_error
# [1] 0.002881543

directory_path <- "./hazel_final_table_output/final_table_output_p_staicu_noshffule_n1000/"
stats <- calculate_stats(directory_path)
print(stats)

# $power
# $power$mean
# [1] 0.0234
# 
# $power$standard_error
# [1] 0.00213787
# 
# 
# $power_01
# $power_01$mean
# [1] 0.0622
# 
# $power_01$standard_error
# [1] 0.003415587


directory_path <- "./hazel_final_table_output/final_table_output_p_staicu_noshffule_n2000/"
stats <- calculate_stats(directory_path)
print(stats)


###################################
##power 4/10/2024
directory_path <- "./hazel_final_table_output/final_table_output_power_shuffle_n100_f7/"
stats <- calculate_stats_percentile(directory_path)
print(stats)

directory_path <- "./hazel_final_table_output/final_table_output_power_shuffle_n100_f8/"
stats <- calculate_stats_percentile(directory_path)
print(stats)

directory_path <- "./hazel_final_table_output/final_table_output_power_shuffle_n100_f9/"
stats <- calculate_stats_percentile(directory_path)
print(stats)

directory_path <- "./hazel_final_table_output/final_table_output_power_shuffle_n100_f10/"
stats <- calculate_stats_percentile(directory_path)
print(stats)
####
directory_path <- "./hazel_final_table_output/final_table_output_power_shuffle_n300_f7/"
stats <- calculate_stats_percentile(directory_path)
print(stats)

directory_path <- "./hazel_final_table_output/final_table_output_power_shuffle_n300_f8/"
stats <- calculate_stats_percentile(directory_path)
print(stats)

directory_path <- "./hazel_final_table_output/final_table_output_power_shuffle_n300_f9/"
stats <- calculate_stats_percentile(directory_path)
print(stats)
directory_path <- "./hazel_final_table_output/final_table_output_power_shuffle_n300_f10/"
stats <- calculate_stats_percentile(directory_path)
print(stats)
#####
directory_path <- "./hazel_final_table_output/final_table_output_power_shuffle_n500_f7/"
stats <- calculate_stats_percentile(directory_path)
print(stats)

directory_path <- "./hazel_final_table_output/final_table_output_power_shuffle_n500_f8/"
stats <- calculate_stats_percentile(directory_path)
print(stats)

directory_path <- "./hazel_final_table_output/final_table_output_power_shuffle_n500_f9/"
stats <- calculate_stats_percentile(directory_path)
print(stats)
directory_path <- "./hazel_final_table_output/final_table_output_power_shuffle_n500_f10/"
stats <- calculate_stats_percentile(directory_path)
print(stats)
########
directory_path <- "./hazel_final_table_output/final_table_output_power_shuffle_n1000_f7/"
stats <- calculate_stats_percentile(directory_path)
print(stats)

directory_path <- "./hazel_final_table_output/final_table_output_power_shuffle_n1000_f8/"
stats <- calculate_stats_percentile(directory_path)
print(stats)

directory_path <- "./hazel_final_table_output/final_table_output_power_shuffle_n1000_f9/"
stats <- calculate_stats_percentile(directory_path)
print(stats)

directory_path <- "./hazel_final_table_output/final_table_output_power_shuffle_n1000_f10/"
stats <- calculate_stats_percentile(directory_path)
print(stats)
###################################

#####boos
calculate_stats_boos <- function(directory_path) {
    power_values <- c()
    
    # Get a list of all RData files in the specified directory by pattern
    files <- list.files(path = directory_path, pattern = "\\.RData$", full.names = TRUE)
    
    for (file_path in files) {
        load(file_path)
        #staicu
        power_values <- c(power_values, unlist(final_table)[5:length(unlist(final_table))])
    }
    
    cat("Total Length (power):", length(power_values))
    
    power_values_0.05 <- calculate_power_stats(power_values,0.05)
    # power_01_values <- calculate_power_stats(power_01_values)
    power_01_values <- calculate_power_stats(power_values,0.1)
    
    ######################
    # power_values2 <- calculate_power_stats(power_values2)
    # power_01_values2 <- calculate_power_stats(power_01_values2)
    #####################
    # return(list(power=power_values, power_01=power_01_values,
    #             power2=power_values2, power_012=power_01_values2))
    
    return(list(power=power_values_0.05, power_01=power_01_values
    ))
}


#load("./hazel_final_table_output/final_table_output_p_boos_n100/Hazel_outputsTbootstrap__100_16_100_99_16_1.RData")

for(num_indvs in c(100, 300, 500))
{
    print("----------------")
    directory_path <- paste0("./test/final_table_output_p_boos_n", num_indvs, "/")
    print(directory_path)
    stats <- calculate_stats_boos(directory_path)
    print(paste0("", round(stats$power$mean, 3),
                 "(", round(stats$power$standard_error, 3), ")",
                 "    ", round(stats$power_01$mean, 3),
                 "(", round(stats$power_01$standard_error, 3), ")"))
    
    for(flchoice in c(7,8,9,10))
    {
        directory_path <- paste0("./test/final_table_output_p_boos_n", num_indvs, "_f", flchoice, "/")
        print(directory_path)
        stats <- calculate_stats_boos(directory_path)
        print(paste0("", round(stats$power$mean, 3),
                    "(", round(stats$power$standard_error, 3), ")",
                    "    ", round(stats$power_01$mean, 3),
                    "(", round(stats$power_01$standard_error, 3), ")"))
    }
    print("================")
}


#####
for(num_indvs in c(100, 300, 500,1000))
{
    print("----------------")
    directory_path <- paste0("./hazel_final_table_output/final_table_output_power_shuffle_n", num_indvs, "/")
    print(directory_path)
    stats <- calculate_stats_percentile(directory_path)
    print(paste0("", round(stats$power$mean, 3),
                 "(", round(stats$power$standard_error, 3), ")",
                 "    ", round(stats$power_01$mean, 3),
                 "(", round(stats$power_01$standard_error, 3), ")"))
    
    for(flchoice in c(7,8,9,10))
    {
        directory_path <- paste0("./hazel_final_table_output/final_table_output_power_shuffle_n", num_indvs, "_f", flchoice, "/")
        print(directory_path)
        stats <- calculate_stats_percentile(directory_path)
        print(paste0("", round(stats$power$mean, 3),
                     "(", round(stats$power$standard_error, 3), ")",
                     "    ", round(stats$power_01$mean, 3),
                     "(", round(stats$power_01$standard_error, 3), ")"))
    }
    print("================")
}
####table 2 first column
# 1] "----------------"
# [1] "./hazel_final_table_output/final_table_output_power_shuffle_n100/"
# Total Length (power): 0
# Total Length (power_01): 0[1] "NA(NA)    NA(NA)"
# [1] "./hazel_final_table_output/final_table_output_power_shuffle_n100_f7/"
# Total Length (power): 1000
# Total Length (power_01): 1000[1] "0.058(0.007)    0.059(0.007)"
# [1] "./hazel_final_table_output/final_table_output_power_shuffle_n100_f8/"
# Total Length (power): 1000
# Total Length (power_01): 1000[1] "0.144(0.011)    0.145(0.011)"
# [1] "./hazel_final_table_output/final_table_output_power_shuffle_n100_f9/"
# Total Length (power): 1000
# Total Length (power_01): 1000[1] "0.24(0.014)    0.242(0.013)"
# [1] "./hazel_final_table_output/final_table_output_power_shuffle_n100_f10/"
# Total Length (power): 1000
# Total Length (power_01): 1000[1] "0.35(0.015)    0.357(0.015)"
# [1] "================"
# [1] "----------------"
# [1] "./hazel_final_table_output/final_table_output_power_shuffle_n300/"
# Total Length (power): 0
# Total Length (power_01): 0[1] "NA(NA)    NA(NA)"
# [1] "./hazel_final_table_output/final_table_output_power_shuffle_n300_f7/"
# Total Length (power): 950
# Total Length (power_01): 950[1] "0.102(0.01)    0.11(0.01)"
# [1] "./hazel_final_table_output/final_table_output_power_shuffle_n300_f8/"
# Total Length (power): 1000
# Total Length (power_01): 1000[1] "0.467(0.016)    0.477(0.015)"
# [1] "./hazel_final_table_output/final_table_output_power_shuffle_n300_f9/"
# Total Length (power): 1000
# Total Length (power_01): 1000[1] "0.697(0.015)    0.7(0.014)"
# [1] "./hazel_final_table_output/final_table_output_power_shuffle_n300_f10/"
# Total Length (power): 1000
# Total Length (power_01): 1000[1] "0.78(0.013)    0.785(0.013)"
# [1] "================"
# [1] "----------------"
# [1] "./hazel_final_table_output/final_table_output_power_shuffle_n500/"
# Total Length (power): 0
# Total Length (power_01): 0[1] "NA(NA)    NA(NA)"
# [1] "./hazel_final_table_output/final_table_output_power_shuffle_n500_f7/"
# Total Length (power): 1000
# Total Length (power_01): 1000[1] "0.259(0.014)    0.263(0.014)"
# [1] "./hazel_final_table_output/final_table_output_power_shuffle_n500_f8/"
# Total Length (power): 1000
# Total Length (power_01): 1000[1] "0.745(0.014)    0.75(0.013)"
# [1] "./hazel_final_table_output/final_table_output_power_shuffle_n500_f9/"
# Total Length (power): 1000
# Total Length (power_01): 1000[1] "0.917(0.009)    0.917(0.009)"
# [1] "./hazel_final_table_output/final_table_output_power_shuffle_n500_f10/"
# Total Length (power): 1000
# Total Length (power_01): 1000[1] "0.93(0.008)    0.93(0.008)"
# [1] "================"
# [1] "----------------"
# [1] "./hazel_final_table_output/final_table_output_power_shuffle_n1000/"
# Total Length (power): 0
# Total Length (power_01): 0[1] "NA(NA)    NA(NA)"
# [1] "./hazel_final_table_output/final_table_output_power_shuffle_n1000_f7/"
# Total Length (power): 1000
# Total Length (power_01): 1000[1] "0.573(0.016)    0.578(0.015)"
# [1] "./hazel_final_table_output/final_table_output_power_shuffle_n1000_f8/"
# Total Length (power): 1000
# Total Length (power_01): 1000[1] "0.964(0.006)    0.966(0.006)"
# [1] "./hazel_final_table_output/final_table_output_power_shuffle_n1000_f9/"
# Total Length (power): 1000
# Total Length (power_01): 1000[1] "0.996(0.002)    0.996(0.002)"
# [1] "./hazel_final_table_output/final_table_output_power_shuffle_n1000_f10/"
# Total Length (power): 1000
# Total Length (power_01): 1000[1] "0.999(0.001)    0.999(0.001)"
# [1] "================"

###RLRT table 2 the 3rd column (only look at time points 90)
##march 26th, after the type I error
#######################
#load("EXP3_r5000_gampower90180inclusion.RData")
# final_table[,1:8]
#                       fl_choice test_type num_subjects num_timepoints  power           se power_01         se01
# experiment_output            6 Inclusion         1000             90 0.0588  0.003326937   0.1122  0.004463433
# experiment_output.1          7 Inclusion         1000             90 0.7934  0.005725669   0.8698  0.004759159
# experiment_output.2          8 Inclusion         1000             90 0.9988 0.0004896039   0.9994 0.0003463062
# experiment_output.3          9 Inclusion         1000             90      1            0        1            0
# experiment_output.4         10 Inclusion         1000             90      1            0        1            0
# experiment_output.5          6 Inclusion          500             90 0.0606  0.003374245   0.1116  0.004452986
# experiment_output.6          7 Inclusion          500             90 0.5038  0.007070864   0.6302  0.006827122
# experiment_output.7          8 Inclusion          500             90 0.9504    0.0030705   0.9766   0.00213787
# experiment_output.8          9 Inclusion          500             90 0.9986 0.0005287797   0.9994 0.0003463062
# experiment_output.9         10 Inclusion          500             90      1            0        1            0
# experiment_output.10         6 Inclusion          300             90 0.0512  0.003117004   0.1028  0.004294931
# experiment_output.11         7 Inclusion          300             90 0.3506  0.006748031   0.4696  0.007057986
# experiment_output.12         8 Inclusion          300             90 0.7954  0.005705065   0.8736  0.004699426
# experiment_output.13         9 Inclusion          300             90 0.9596  0.002784523   0.9792  0.002018284
# experiment_output.14        10 Inclusion          300             90  0.994  0.001092154    0.998 0.0006318228
# experiment_output.15         6 Inclusion          100             90 0.0576  0.003294912   0.1098  0.004421402
# experiment_output.16         7 Inclusion          100             90  0.137  0.004862736   0.2302  0.005953284
# experiment_output.17         8 Inclusion          100             90  0.368  0.006820205   0.4854  0.007068053
# experiment_output.18         9 Inclusion          100             90  0.564  0.007012902   0.6926  0.006525416
# experiment_output.19        10 Inclusion          100             90 0.7408  0.006197021   0.8298  0.005314733
# experiment_output.20         6 Inclusion         1000            180 0.0606  0.003374245   0.1076  0.004382288
# experiment_output.21         7 Inclusion         1000            180  0.567  0.007007296   0.6784  0.006605656
# experiment_output.22         8 Inclusion         1000            180 0.9722  0.002324958   0.9866  0.001626065
# experiment_output.23         9 Inclusion         1000            180 0.9998   0.00019998   0.9998   0.00019998
# experiment_output.24        10 Inclusion         1000            180      1            0        1            0
# experiment_output.25         6 Inclusion          500            180 0.0574  0.003289536   0.1124  0.004466906
# experiment_output.26         7 Inclusion          500            180 0.3346  0.006672973    0.451  0.007037031
# experiment_output.27         8 Inclusion          500            180  0.782   0.00583911   0.8602  0.004904201
# experiment_output.28         9 Inclusion          500            180 0.9538  0.002968689   0.9776  0.002092761
# experiment_output.29        10 Inclusion          500            180   0.99  0.001407125   0.9956 0.0009360171
# experiment_output.30         6 Inclusion          300            180 0.0578  0.003300278   0.1092   0.00441079
# experiment_output.31         7 Inclusion          300            180 0.2224  0.005881126   0.3256  0.006626985
# experiment_output.32         8 Inclusion          300            180  0.564  0.007012902   0.6858  0.006564729
# experiment_output.33         9 Inclusion          300            180   0.81  0.005547973     0.88   0.00459565
# experiment_output.34        10 Inclusion          300            180 0.9192  0.003854124   0.9564  0.002887873
# experiment_output.35         6 Inclusion          100            180 0.0582  0.003310974   0.1074  0.004378704
# experiment_output.36         7 Inclusion          100            180 0.1106  0.004435485   0.1852   0.00549365
# experiment_output.37         8 Inclusion          100            180 0.2386   0.00602777   0.3498  0.006744479
# experiment_output.38         9 Inclusion          100            180 0.3806  0.006866493    0.506  0.007070559
# experiment_output.39        10 Inclusion          100            180 0.5106  0.007069479   0.6436  0.006773168




