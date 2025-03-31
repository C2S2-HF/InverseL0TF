
#' @title Print four metrics about change point detection results
#' @description Prints four metrics to compare the quality of change point detection results.
#' @param y0 The underlying trend
#' @param tau The locations of change points in the underlying trend
#' @param yhat The fitted trend
#' @param cpts The positions of the fitted change points
#' @return
#' \item{MSE}{The mean square error between the fitted trend and the underlying trend}
#' \item{MAD}{The median absolute deviation between the fitted trend and the underlying trend}
#' \item{dH}{Hausdorff Distance (dH) measures the accuracy of the estimated change points}
#' \item{nknot}{The number of detected change points}
#' @details
#' \eqn{\hat{\boldsymbol{\tau}}} represents the estimated change point positions, while \eqn{\boldsymbol{\tau}} denotes the locations of change points in the underlying trend.
#' \deqn{d_H=\frac{1}{n} \max \{\max_k \min_j |\tau_j-\hat{\tau}_k|,\max_j \min_k |\tau_j-\hat{\tau}_k|\}.}
#' Note that the number of \eqn{\hat{\boldsymbol{\tau}}} and \eqn{\boldsymbol{\tau}} does not need to be the same.
#' @examples
#'
#' tau = c(0.1, 0.3, 0.4, 0.7, 0.85)
#' h = c(-1, 5, 3, 0, -1, 2)
#' n = 500
#' BlocksData <- SimuBlocksInv(n = n, sigma = 0.2, seed = 50, tau = tau ,h = h)
#' res <- L0TFinv.opt(y=BlocksData$y, kmax=10, q=0, first=0.01, last=1, penalty="bic")
#' metrics <- TFmetrics(BlocksData$y0,BlocksData$tau,res$yopt,res$Aopt/n)
#' print(metrics)
#'
#' tau1 = c(0.1, 0.3, 0.4, 0.7, 0.85)
#' h1 = c(-1, 5, 3, 0, -1, 2)
#' a0 = -10
#' n1 = 2000
#' WaveData <- SimuWaveInv(n = n1, sigma = 0.1, seed = 50, tau = tau1, h = h1, a0 = a0)
#' res1 <- L0TFinv.fix(y=WaveData$y, k=20, q=1, first=0, last=0.99)
#' metrics1 <- TFmetrics(WaveData$y0,WaveData$tau,res1$y.all[,5],res1$A.all[[5]]/n1)
#' print(metrics1)
#'
#' @export
TFmetrics <- function(y0, tau=NULL, yhat, cpts = NULL){
  mse = mean((yhat-y0)^2)
  mad = mean(abs(yhat-y0))
  if(is.null(cpts)){
    tab = data.frame(MSE=mse, MAD=mad)
  }else{
    n = length(yhat)
    n.cpts = length(cpts)
    segments.endpoints.true = sort(unique(tau))
    segments.endpoints.est = sort(unique(cpts))
    distm = abs(matrix(rep(segments.endpoints.est, length(segments.endpoints.true)),
        nrow=length(segments.endpoints.est))-matrix(rep(segments.endpoints.true,
        length(segments.endpoints.est)), nrow=length(segments.endpoints.est), byrow=TRUE))
    screening.dist = max(apply(distm, 2, min))/n
    precision.dist = max(apply(distm, 1, min))/n
    haus.dist = max(screening.dist, precision.dist)
  }
  tab = data.frame(MSE=mse, MAD=mad, dH=haus.dist, nknot=n.cpts)
  return(tab)
}








