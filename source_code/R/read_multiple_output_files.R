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



calculate_power_stats <- function(power_values){
    power_mean <- mean(power_values)
    power_se <- sqrt(power_mean * (1 - power_mean) / length(power_values))
    return(list(mean = power_mean, standard_error = power_se))
}

calculate_stats <- function(directory_path) {
    power_values <- c()
    power_01_values <- c()
    
    power_values2 <- c()
    power_01_values2 <- c()
    
    
    # Get a list of all RData files in the specified directory by pattern
    files <- list.files(path = directory_path, pattern = "\\.RData$", full.names = TRUE)
    
    for (file_path in files) {
        load(file_path)
        power_values <- c(power_values, final_table$power[[1]])
        power_01_values <- c(power_01_values, final_table$power_01[[1]])
        
        ###add one more power
        power_values2 <- c(power_values2, final_table$power2[[1]])
        power_01_values2 <- c(power_01_values2, final_table$power_012[[1]])
        ################################
    }
    
    cat("Total Length (power):", length(power_values))
    cat("\nTotal Length (power_01):", length(power_01_values))
    
    power_values <- calculate_power_stats(power_values)
    power_01_values <- calculate_power_stats(power_01_values)
    
    ######################
    power_values2 <- calculate_power_stats(power_values2)
    power_01_values2 <- calculate_power_stats(power_01_values2)
    #####################
    return(list(power=power_values, power_01=power_01_values,
                power2=power_values2, power_012=power_01_values2))
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
