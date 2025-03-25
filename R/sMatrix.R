#' @title Generate a difference matrix
#' @description This function generates a matrix for computing differences of a certain order, useful in numerical methods and for creating specific matrix patterns.
#' @param n The number of data points
#' @param q The order of the difference
#' @return  A matrix with dimensions \eqn{n-q-1} by \eqn{n}, whose elements correspond to the combinatorial values of \eqn{q}.
#' @examples
#' Mat1 <- DiffMat(n = 10, q = 0)
#' print(Mat1)
#'
#' Mat2 <- DiffMat(n = 15, q = 1)
#' print(Mat2)
#'
#' Mat3 <- DiffMat(n = 15, q = 2)
#' print(Mat3)
#'
#' @seealso \code{\link{XMat}}
#' @export
DiffMat <- function(n, q) {
  if( n <= q+1 ){
    stop("The number n should be greater than the order of the difference")
  }
  X <- matrix(0, (n-q-1), n)
  for (i in 1:(n-q-1)) {
    for (j in 1:n) {
      if (j >= i && j <= i + q + 1) {
        X[i, j] <- (-1)^(j - i + q - 1) * choose(q + 1, j - i)
      }
    }
  }
  return(X)
}



#' @title Generate an artificial design matrix
#' @description This matrix corresponds to the difference matrix, transforming the L0 trend filtering model into an inverse statistical problem.
#' @param n The number of data points
#' @param q The order of the difference
#' @return  A matrix with dimensions \eqn{n} by \eqn{n}, whose elements correspond to the difference matrix.
#' @examples
#' mat1 <- XMat(n = 10, q = 0)
#' print(mat1)
#'
#' mat2 <- XMat(n = 15, q = 1)
#' print(mat2)
#'
#' mat3 <- XMat(n = 15, q = 2)
#' print(mat3)
#'
#' Mat1 <- DiffMat(n = 10, q = 0)
#' Mat2 <- DiffMat(n = 15, q = 1)
#' Mat3 <- DiffMat(n = 15, q = 2)
#' print(Mat1%*%mat1)
#' print(Mat2%*%mat2)
#' print(Mat3%*%mat3)
#' @details Noticing the correspondence between \eqn{\boldsymbol{D}^{(q+1)}} and \eqn{\boldsymbol{X}^{(q+1)}}, the result of their matrix multiplication is a combination of a zero matrix and an identity matrix. Expressed as \eqn{\boldsymbol{D}^{(q+1)} \boldsymbol{X}^{(q+1)}=(\boldsymbol{O}_{(n-q-1)\times(q+1)},\quad \boldsymbol{I}_{(n-q-1)\times(n-q-1)})}. The result is advantageous for the invertible processing of the original L0 trend filtering problem.
#' @export
XMat <- function(n, q){
  X <- matrix(0,n,n)
  if(q == 0){
    for(i in 1:n){
      for(j in 1:i){
        X[i,j] <- 1
      }
    }
    return(X)
  }else{
    return(apply(XMat(n, q-1),2,cumsum))
  }
}






