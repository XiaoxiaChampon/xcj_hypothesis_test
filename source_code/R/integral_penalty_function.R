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
  integral = delta_h * (approx_value$y[2 * (1:number_points) - 1] + 4 * approx_value$y[2 * (1:number_points)] + 
                          approx_value$y[2 * (1:number_points) + 1])/3
  results = list(value = sum(integral), cdf = list(
    x = approx_value$x[2 * (1:number_points)], 
    y = cumsum(integral))
  )
  return(results)
}