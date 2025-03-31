#' @title Generate the inverse of the crossprod matrix
#' @description Generate the inverse matrix of \eqn{(\boldsymbol{X}^{(q+1)}_A)^T \boldsymbol{X}^{(q+1)}_A} for the cases where \eqn{q=0} or \eqn{q=1}, commonly employed in splicing algorithms. Note that an explicit solution exists for the inverse when \eqn{q=0}, but not when \eqn{q=1}.
#' @param n The number of data points
#' @param q The order of the difference, 0 or 1
#' @param A The set of indicators, a subset of \eqn{\{1,2,3,\dots,n\}}
#' @return The inverse matrix of \eqn{(\boldsymbol{X}^{(q+1)}_A)^T \boldsymbol{X}^{(q+1)}_A} for the cases where \eqn{q=} 0 or 1.
#' @examples
#' Mat1 <- XMat(n = 10, q = 0)
#' A1 = c(1,2,5,8)
#' mat1 = as.matrix(Mat1[,A1])
#' S1 <- solMat(n = 10, q = 0, A = A1)
#' print(S1)
#' print(round(S1%*%t(mat1)%*%mat1,10))
#'
#' Mat2 <- XMat(n = 15, q = 1)
#' A2 = c(1,3,8,10,15)
#' mat2 = as.matrix(Mat2[,A2])
#' S2 <- solMat(n = 15, q = 1, A = A2)
#' print(S2)
#' print(round(S2%*%t(mat2)%*%mat2,10))
#' @importFrom Matrix solve
#' @export
solMat <- function(n, q, A){
  if ( !(q %in% c(0,1)) ){
    stop("The specified order is not supported")
  }
  if( min(A)<1 | max(A)>n ){
    stop("
The maximum possible range for the set indicators should be within {1,2,...,n}")
  }
  k = length(A)
  m = as.numeric(k)
  m = n + 1 - A
  phi = matrix(0, k, k)
  if(q == 0){
    if(k == 1){
      phi[1,1] = 1/m[1]
      return(phi)
    }
    if(k == 2){
      phi[1,1] = m[2]
      phi[1,2] = -m[2]
      phi[2,1] = -m[2]
      phi[2,2] = m[1]
      phi = phi/(m[2]*(m[1]-m[2]))
      return(phi)
    }
    for(i in 1:k){
      if(i == 1){
        phi[1,1] = 1/(m[1]-m[2])
        phi[1,2] = -1/(m[1]-m[2])
        next
      }
      if(i == k){
        phi[k,(k-1)] = -1/(m[(k-1)]-m[k])
        phi[k,k] = 1/(m[(k-1)]-m[k])+1/m[k]
        break
      }
      phi[i,(i-1)] = -1/(m[(i-1)]-m[i])
      phi[i,i] = 1/(m[(i-1)]-m[i])+1/(m[i]-m[(i+1)])
      phi[i,(i+1)] = -1/(m[i]-m[(i+1)])
    }
    return(phi)
  }
  if (q == 1) {
    for (j in 1:k) {
      for (i in j:k) {
        if (j == i) {
          phi[i, j] <- m[j] * (m[j] + 1) * (2 * m[j] + 1) / 6
        } else {
          phi[i, j] <- m[i] * (m[i] + 1) * (3 * m[j] - m[i] + 1) / 6
        }
      }
    }
    phi <- phi + t(phi) - diag(diag(phi))
    phi <- solve(phi)
    return(phi)
  }
}



