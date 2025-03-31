#' @title Simulate Blocks Data
#' @description This function generates data points of piecewise constant trends.
#' @param n Number of data points
#' @param sigma Standard deviation of the noise added to the signal
#' @param seed An optional seed for random number generation to make results reproducible
#' @param tau The locations of change points in the underlying trend
#' @param h The constant values of the \eqn{length(tau)+1} segments of the underlying trend
#' @return
#' A list containing the piecewise constant simulated data and the underlying trend:
#' \item{x}{The set \{ \eqn{\frac{1}{n}}, \eqn{\frac{2}{n}}, \eqn{\frac{3}{n}},\dots, \eqn{1}\}}
#' \item{y}{The piecewise constant simulated data of length \eqn{n}}
#' \item{y0}{The underlying trend of length \eqn{n}}
#' \item{setA}{The set of position indicators of change points in the simulated data}
#' \item{tau}{The locations of change points in the underlying trend}
#'
#' @details
#' \itemize{
#' \item{}{To simplify the analysis, normalize the change point positions to a range between 0 and 1. Require that all elements of the input \eqn{tau} are within this range. Consequently, the change point positions in simulated data forms a subset of the set \{ \eqn{\frac{1}{n}}, \eqn{\frac{2}{n}}, \eqn{\frac{3}{n}},\dots, 1\}.}
#' \item{}{In fact, \eqn{length(tau)} change points can divide the interval into \eqn{length(tau)+1} segments of constant function values. Therefore, ensure that the length of vector \eqn{h} is \eqn{length(tau)+1}.}
#' }
#' @examples
#' tau = c(0.1, 0.3, 0.4, 0.7, 0.85)
#' h = c(-1, 5, 3, 0, -1, 2)
#' BlocksData <- SimuBlocksInv(n = 350, sigma = 0.1, seed = 50, tau = tau ,h = h)
#' plot(BlocksData$x, BlocksData$y, xlab="", ylab="")
#' lines(BlocksData$x, BlocksData$y0, col = "red")
#' print(BlocksData$setA)
#' print(BlocksData$tau)
#' @importFrom stats rnorm
#' @export
SimuBlocksInv <- function (n, sigma, seed = NA, tau, h ){
  if (!is.na(seed)) set.seed(seed)
  if( min(tau)<=0 | max(tau)>=1 ){
    stop("
The maximum possible range for the change points should be within [0,1]")
  }
  if (n < length(tau)+1){
    stop("The number of data points should be greater than the number of change points")
  }
  if (length(h)!=length(tau)+1){
    stop("The length of the vector does not meet the condition")
  }
  x = seq(1/n, 1,length.out = n)
  A = sapply(tau, function(z) which(x>=z)[1])
  beta = rep(0,n)
  beta[1] = h[1]
  beta[A+1] = diff(h)
  y0 = cumsum(beta)
  y = y0 + sigma*rnorm(n)
  return(list(x = x, y = y, y0 = y0, setA = A,  tau = tau))
}



#' @title Simulate Wave Data
#' @description This function generates data points of piecewise linear trends.
#' @param n Number of data points
#' @param sigma Standard deviation of the noise added to the signal
#' @param seed An optional seed for random number generation to make results reproducible
#' @param tau The locations of change points in the underlying trend
#' @param h The slope of the \eqn{length(tau)+1} segments of the underlying trend
#' @param a0 The initial point value
#' @return
#' A list containing the piecewise linear simulated data and the underlying trend:
#' \item{x}{The set \{ \eqn{\frac{1}{n}}, \eqn{\frac{2}{n}}, \eqn{\frac{3}{n}},\dots, \eqn{1} \}}
#' \item{y}{The piecewise linear simulated data of length \eqn{n}}
#' \item{y0}{The underlying trend of length \eqn{n}}
#' \item{setA}{The set of position indicators of change points in the simulated data}
#' \item{tau}{The locations of change points in the underlying trend}
#'
#' @examples
#' tau = c(0.1, 0.3, 0.4, 0.7, 0.85)
#' h = c(-1, 5, 3, 0, -1, 2)
#' a0 = -10
#' WaveData <- SimuWaveInv(n = 650, sigma = 0.1, seed = 50, tau = tau, h = h, a0 = a0)
#' plot(WaveData$x, WaveData$y, xlab="", ylab="")
#' lines(WaveData$x, WaveData$y0, col = "red")
#' print(WaveData$setA)
#' print(WaveData$tau)
#' @seealso \code{\link{SimuBlocksInv}}
#' @importFrom stats rnorm
#' @export
SimuWaveInv <- function (n, sigma, seed = NA, tau, h, a0 = 0 ){
  if (!is.na(seed)) set.seed(seed)
  if( min(tau)<=0 | max(tau)>=1 ){
    stop("
The maximum possible range for the change points should be within [0,1]")
  }
  if (n < length(tau)+2){
    stop("The number of data points should be greater than the number of change points")
  }
  if (length(h)!=length(tau)+1){
    stop("The length of the vector does not meet the condition")
  }
  x = seq(1/n, 1,length.out = n)
  A = sapply(tau, function(z) which(x>=z)[1])
  beta = rep(0,n)
  beta[1] = a0
  beta[2] = h[1]/n-beta[1]
  beta[A+1] = (diff(h))/n
  y0 = cumsum(cumsum(beta))
  y = y0 + sigma*rnorm(n)
  return(list(x = x, y = y, y0 = y0, setA = A,  tau = tau))
}

             
