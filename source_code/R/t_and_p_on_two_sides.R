library(ggplot2)
library(tidyverse)
library(gridExtra)

files_nleq <- list.files(path = "./t_and_p_value_save_n100/", pattern = "pnleq.*\\.RData$", full.names = TRUE)
length(files_nleq)

files_leq <- list.files(path = "./t_and_p_value_save_n100/", pattern = "pleq.*\\.RData$", full.names = TRUE)
length(files_leq)

leq_T_values <- NULL
for (file_path in files_leq) {
  load(file_path)
  leq_T_values <- c(leq_T_values, T_value)
}

nleq_T_values <- NULL
for (file_path in files_nleq) {
  load(file_path)
  nleq_T_values <- c(nleq_T_values, T_value)
}

df <- tibble(
  T_val = c(leq_T_values, nleq_T_values),    # Combine vectors into one column
  Rejected = c(rep(TRUE, length(leq_T_values)), rep(FALSE, length(nleq_T_values)))  # Repeat TRUE for length of vector1, and FALSE for vector2
)


ggplot(df, aes(x = T_val, fill=Rejected) ) +
  geom_histogram(alpha=0.6) +
  labs(x = "Duration", y = "Frequency", title = "Histogram of T Value") +
  theme(plot.title = element_text(hjust = 0.5))+ 
  theme(text = element_text(size = 14) ) # Set the desired font size)


ggplot(df, aes(x = T_val) ) +
  geom_histogram(alpha=0.6) +
  facet_grid(. ~ Rejected) +
  labs(x = "Duration", y = "Frequency", title = "Histogram of T Value") +
  theme(plot.title = element_text(hjust = 0.5))+ 
  theme(text = element_text(size = 14) ) # Set the desired font size)
  


for (sam_file in sample(files_leq, size = 5)) {
  load(sam_file)
  
}

ggplot(data.frame(T_star_values), aes(x = T_star_values) ) +
  geom_histogram(alpha=0.6) +
  # labs(x = "Duration", y = "Frequency", title = "Histogram of T Star Values") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 14) ) + 
  geom_vline(xintercept = T_value, linetype = "dashed", color = "red")

# Set up plot area
par(mfrow = c(5, 2))  # 5 rows, 2 columns

file_leq_sam <- sample(files_leq, size = 5)
file_nleq_sam <- sample(files_nleq, size = 5)

for (idx in 1:5) {
  file <- file_leq_sam[idx]
  data <- load(file)  # assuming CSV format; adjust as needed
  
  # Create histogram and add vertical line for T_value
  hist(T_star_values, main = paste("Histogram for", basename(file)), xlab = "T_star values")
  abline(v = T_value, col = "red", lwd = 2)
  
  file <- file_nleq_sam[idx]
  data <- load(file)  # assuming CSV format; adjust as needed
  
  # Create histogram and add vertical line for T_value
  hist(T_star_values, main = paste("Histogram for", basename(file)), xlab = "T_star values")
  abline(v = T_value, col = "blue", lwd = 2)
}

# Reset plot settings
par(mfrow = c(1, 1))



# Function to create a plot for a given file
create_plot <- function(file) {
  data <- load(file)  # Load data, adjust if format differs
  
  # Create a ggplot histogram with a vertical line for T_value
  p <- ggplot(data.frame(T_star_values), aes(x = T_star_values)) +
    geom_histogram(binwidth = 1, fill = "grey", color = "black") +  # Adjust binwidth as necessary
    geom_vline(aes(xintercept = T_value), color = "red", linetype = "dashed", size = 1.5) +
    labs(title = paste("Histogram for", basename(file)), x = "T_star values") +
    theme_minimal()
  
  return(p)
}

# Generate plots for each file
plots_folderA <- lapply(sample(files_leq, size = 5), create_plot)
plots_folderB <- lapply(sample(files_nleq, size = 5), create_plot)

# Combine plots into a grid
grid_plots <- mapply(gridExtra::arrangeGrob, plots_folderA, plots_folderB, SIMPLIFY = FALSE)
combined_grid <- do.call(gridExtra::grid.arrange, c(grid_plots, ncol = 2))

