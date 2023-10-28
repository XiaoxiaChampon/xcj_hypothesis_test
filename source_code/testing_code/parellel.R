library(foreach)

set.seed(123)
foo <- function(){ return(runif(1)) }

J_matrix<-matrix(0,nrow=5,ncol=7) #empty
number_row <- nrow(J_matrix)
number_col <- ncol(J_matrix)

for(row in 1:number_row){
  for(col in 1:number_col){
    J_matrix[row,col] <- foo()
  }
}

set.seed(123)
J_list <- foreach(row = 1:number_row) %do%
  {
    a <- array(-123,col)
    for(col in 1:number_col){
      a[col] <- foo()
    }
    return(a)
  }
other <- do.call(rbind, J_list)

print(J_list)
print(J_matrix)
print(other)
