kron.prod <- function(y,matrices){
# INPUT
# matrices: list(matrix1, matrix2,..., matrixd)
# where dim(matrixj) = c(nj,nj) [SQUARE!]
# y: value vector
# with length(y) = n1*...*nd
#
# OUTPUT
# y1: kron(matrixd,kron(...,(kron(matrix2,matrix1))...))*y
#
# check user input
stopifnot(is.list(matrices))
# stop if not all of them are square
allsquare <- lapply(matrices,function(x){dim(x)[1]==dim(x)[2]})
stopifnot(all(unlist(allsquare)))
# get sizes
nmatrices <- length(matrices)
nall <- length(y)
nullvec <- rep(0,nall)
# compute product for first matrix
y0 <- y
stemp <- matrices[[1]]
n <- dim(stemp)[1]
m <- nall/n
y1 <- nullvec
for (i in 1:m){
y1[m*(0:(n-1)) + i] <- stemp %*% y0[(n*(i-1)) + (1:n)]
}
if (nmatrices > 1){
# for all other matrices
for(imat in 2:nmatrices){
y0 <- y1
stemp <- matrices[[imat]]
n <- dim(stemp)[1]
m <- nall/n
y1 <- nullvec
for (i in 1:m){
y1[m*(0:(n-1)) + i] <- stemp %*% y0[(n*(i-1)) + (1:n)]
}
}
}
return(y1)
}