############################################################
# Copyright 2023 Xiaoxia Champon

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
# Purpose: Functions to perform Categorical Functional Data Hypothesis Testing
#           File that contains all the functions necessary to generate data 
# Author:  Xiaoxia Champon
# Date: 10/26/2023
#
##############################################################

library(matrixStats)

#' Function to find integral of penalty
#' #integral_penalty(time_interval,X_matrix[this_row,]*bspline[,this_col])$value
integral_penalty <- function(time_interval, function_value) { 
  #test
  # x_app <- 1:10
  # y_app <- rnorm(10)
  # time_interval=x_app
  # function_value=y_app
  number_points = max(256, length(time_interval))
  timeseries_length = length(time_interval)
  if (timeseries_length != length(function_value)) {stop("Unequal input vector lengths")}
  approx_value <- approx(time_interval, function_value, n = 2 * number_points + 1)
  delta_h = diff(approx_value$x)[1]
  # print(approx_value$y)
  integral = delta_h * (approx_value$y[2 * (1:number_points) - 1] + 4 * approx_value$y[2 * (1:number_points)] + 
                          approx_value$y[2 * (1:number_points) + 1])/3
  # print(integral)
  results = list(value = sum(integral), cdf = list(
    x = approx_value$x[2 * (1:number_points)], 
    y = cumsum(integral))
  )
  return(results)
}
#function_value_matrix is a matrxi t*R where t is the number of time points, R is the number of basis
#pl is a vecotr with length t, 
integral_penalty_matrix <- function(time_interval, function_value_matrix) { 
    #test
    # x_app <- 1:10
    # y_app <- rnorm(10)
    # time_interval=x_app
    # function_value=y_app
    number_points = max(256, length(time_interval))
    timeseries_length = length(time_interval)
    if (timeseries_length != dim(function_value_matrix)[1]) {stop("Unequal input vector lengths")}
    approx_value_y <- matrix(0, ncol=2 * number_points + 1 ,nrow= dim(function_value_matrix)[2])
    approx_value_x <- matrix(0, ncol=2 * number_points + 1 ,nrow= dim(function_value_matrix)[2])
    for(thing in 1:dim(function_value_matrix)[2]){
        new_value <- approx(time_interval, function_value_matrix[,thing], n = 2 * number_points + 1)
        approx_value_y[thing,] <- new_value$y
        approx_value_x[thing,] <-new_value$x
    }
    #print(dim(approx_value_y))
    delta_h = approx_value_x[,2] - approx_value_x[,1]
    integral = delta_h * (approx_value_y[,2 * (1:number_points) - 1] 
                          + 4 * approx_value_y[,2 * (1:number_points)] 
                          + approx_value_y[,2 * (1:number_points) + 1]) / 3
    results = list(value = rowSums(integral), cdf = list(
        x = approx_value_x[,2 * (1:number_points)], 
        y = rowCumsums(integral))
    )
    return(results)
}

integral_function <- function(time_interval, function_value){
        dt = diff(time_interval)/2.0
        m=length(time_interval)
        z = c(0, cumsum(dt*(function_value[1:(m-1)] + function_value[2:m])))
    return(z)
}