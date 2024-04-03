
calculate_power_stats <- function(power_values){
    power_mean <- mean(power_values)
    power_se <- sqrt(power_mean * (1 - power_mean) / length(power_values))
    return(list(mean = power_mean, standard_error = power_se))
}

calculate_stats <- function(directory_path) {
    power_values <- c()
    power_01_values <- c()
    
    # Get a list of all RData files in the specified directory by pattern
    files <- list.files(path = directory_path, pattern = "\\.RData$", full.names = TRUE)
    
    for (file_path in files) {
        load(file_path)
        power_values <- c(power_values, final_table$power[[1]])
        power_01_values <- c(power_01_values, final_table$power_01[[1]])
    }
    
    cat("Total Length (power):", length(power_values))
    cat("\nTotal Length (power_01):", length(power_01_values))
    
    power_values <- calculate_power_stats(power_values)
    power_01_values <- calculate_power_stats(power_01_values)
    
    return(list(power=power_values, power_01=power_01_values))
}

# Example usage
directory_path <- "./hazel_outputs"
stats <- calculate_stats(directory_path)
print(stats)