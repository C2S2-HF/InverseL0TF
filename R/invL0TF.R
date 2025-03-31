

ybeta <- function(beta=beta,q=q){
  if(q == 0){
    return(cumsum(beta))
  }
  if(q == 1){
    return(cumsum(cumsum(beta)))
  }
}

Dy <- function(y=y,q=q,U=U){
  if(q == 0){
    D = U-rev(cumsum(rev(y)))
    return(D)
  }
  if(q == 1){
    D = U-rev(cumsum(cumsum(rev(y))))
    return(D)
  }
}

Splicing <- function(A=A,I=I,k=k,y=y,q=q,H=H,U=U,first=first,last=last){
  n = length(y)
  All = 1:n
  S = 1:(q+1)
  AS = (q+2):n
  low = ceiling(first*n)
  high = floor(last*n)
  AS0 = intersect(low:high,AS)

  eps1 = rep(0,n)
  eps2 = rep(0,n)
  Al = sort(union(S,A),decreasing = FALSE)
  beta = rep(0,n)
  beta[Al] = solMat(n=n,q=q,A=Al)%*%U[Al]
  yhat = ybeta(beta=beta,q=q)
  D = Dy(y=yhat,q=q,U=U)/n
  eps1 = ((beta)^2)*H
  eps2 = D^2/H
  L = sum((y-yhat)^2)
  for(j in 1:(k%/%4)){
    A0 = as.vector(A[order(eps1[A],decreasing = F)[1:j]])
    I00 = intersect(I,AS0)
    I0 = as.vector(I00[order(eps2[I00],decreasing = T)[1:j]])
    Anew = c(setdiff(A,A0),I0)
    Inew = c(setdiff(I,I0),A0)
    Alnew = sort(union(Anew,S),decreasing = FALSE)
    betanew = rep(0,n)
    betanew[Alnew] = solMat(n=n,q=q,A=Alnew)%*%U[Alnew]
    yhatnew = ybeta(beta=betanew,q=q)
    Lnew = sum((y-yhatnew)^2)
    if(L-Lnew>0){
      A = Anew
      I = Inew
      beta = betanew
      L = Lnew
      yhat = ybeta(beta=beta,q=q)
      D = Dy(y=yhat,q=q,U=U)/n
      eps1 = ((beta)^2)*H
      eps2 = D^2/H
    }
  }
  return(list(A=A,I=I,beta=beta))
}

InvL0TFk <- function(A0=A0,y=y,q=q,k=k,H=H,U=U,first=first,last=last,max.step=50){
  n = length(y)
  All = 1:n
  S = 1:(q+1)
  AS = (q+2):n
  I0 = setdiff(AS,A0)
  for(j in 1:max.step){
    m = Splicing(A=A0,I=I0,k=k,y=y,q=q,H=H,U=U,first=first,last=last)
    A = m$A
    I = m$I
    beta = m$beta
    if(identical(A,A0) & identical(I,I0)){
      break
    }else{
      A0 = A
      I0 = I
    }
  }
  yhat = ybeta(beta=beta,q=q)
  Ahat = sort(A,decreasing = FALSE)
  return(list(betak=beta,yk=yhat,Ak=Ahat))
}

InverseL0TF <- function(y=y,kmax=kmax,q=q,first=0,last=1){
  n = length(y)
  All = 1:n
  S = 1:(q+1)
  AS = (q+2):n
  low = ceiling(first*n)
  high = floor(last*n)
  AS0 = intersect(low:high,AS)
  A0 = NULL
  I0 = setdiff(AS,A0)

  beta.all = NULL
  y.all = NULL
  A.all = list()
  mse = as.numeric(kmax)
  sic = as.numeric(kmax)
  eps = as.numeric(n)
  if(q == 0){
    H = n:1/n
    U = rev(cumsum(rev(y)))
  }
  if(q == 1){
    H = sapply(n:1, function(x) x*(x+1)*(2*x+1)/6)/n
    U = rev(cumsum(cumsum(rev(y))))
  }
  Al = sort(union(S,A0),decreasing = FALSE)
  beta = rep(0,n)
  beta[Al] = solMat(n=n,q=q,A=Al)%*%U[Al]
  D = Dy(y=ybeta(beta=beta,q=q),q=q,U=U)/n
  eps = D^2/H
  I00 = intersect(I0,AS0)
  A0 = union(A0,I00[which.max(eps[I00])])
  for(j in 1:kmax){
    result = InvL0TFk(A0=A0,y=y,q=q,k=j,H=H,U=U,first=first,last=last)
    beta.all = cbind(beta.all,result$betak)
    y.all = cbind(y.all,result$yk)
    A0 = result$Ak
    I0 = setdiff(AS,A0)
    A.all[[j]] = A0 - 1
    D = Dy(y=as.vector(y.all[,j]),q=q,U=U)/n
    eps = D^2/H
    I00 = intersect(I0,AS0)
    A0 = union(A0,I00[which.max(eps[I00])])
  }
  mse = colMeans((y-y.all)^2)
  df = 1:kmax + q + 1
  sic = n*log(mse) + 2*log(log(n))*log(n)*df
  bic = n*log(mse) + 2*log(n)*df
  return(list(beta.all=beta.all,y.all=y.all,A.all=A.all,
              sic=sic,bic=bic,mse=mse))
}


#' @title The inverse L0 trend filtering with fixed change points
#' @description Fit the input data points to a piecewise constant or piecewise linear trend with a given number of change points.
#' @param y The input data points
#' @param k The given number of change points
#' @param q 0 or 1. Correspond to a piecewise constant or piecewise linear trend
#' @param first The value ranges from 0 to 1. Represent the minimum percentile point where a change point may occur. If 'first' = 0.01, it means that change points cannot appear in the first 1\% of the data points. If 'first' = 0, there is no constraint on the position of the change point.
#' @param last The value ranges from 0 to 1. Represent the maximum percentile point where a change point may occur. If 'last' = 0.99, it means that change points cannot appear in the last 1\% of the data points. If 'last' = 1, there is no constraint on the position of the change point.
#' @return
#' An S3 object of type "L0TFinvfix". A list containing the fitted trend results:
#' \item{sic}{Information criterion value with a penalty term of \eqn{2\log(\log(n)) \times \log(n)}}
#' \item{bic}{Information criterion value with a penalty term of \eqn{2 \times \log(n)}}
#' \item{mse}{The mean square error between the fitted trend and the input data}
#' \item{y}{The input data points}
#' \item{betak}{The fitted \eqn{\hat{\boldsymbol{\beta}}} coefficients with the number of change points being \eqn{k} }
#' \item{yk}{The fitted trend with the number of change points being \eqn{k}}
#' \item{Ak}{The set of position indicators of the fitted change points with the number of change points being \eqn{k}}
#' \item{beta.all}{A data frame with dimensions \eqn{n \times k}, where each column represents the fitted \eqn{\hat{\boldsymbol{\beta}}} coefficients corresponding to a given number of change points}
#' \item{y.all}{A data frame with dimensions \eqn{n \times k}, where each column represents the fitted estimated trend corresponding to a given number of change points}
#' \item{A.all}{A list of length \eqn{k}, where each element corresponds to the set of position indicators of change points under a given number }
#'
#' @examples
#' tau = c(0.1, 0.3, 0.4, 0.7, 0.85)
#' h = c(-1, 5, 3, 0, -1, 2)
#' BlocksData <- SimuBlocksInv(n = 350, sigma = 0.2, seed = 50, tau = tau ,h = h)
#' res <- L0TFinv.fix(y=BlocksData$y, k=5, q=0, first=0.01, last=1)
#' print(res$Ak)
#' print(BlocksData$setA)
#' plot(BlocksData$x, BlocksData$y, xlab="", ylab="")
#' lines(BlocksData$x, BlocksData$y0, col = "red")
#' lines(BlocksData$x, res$yk, col = "lightgreen")
#'
#' tau1 = c(0.1, 0.3, 0.4, 0.7, 0.85)
#' h1 = c(-1, 5, 3, 0, -1, 2)
#' a0 = -10
#' WaveData <- SimuWaveInv(n = 2000, sigma = 0.1, seed = 50, tau = tau1, h = h1, a0 = a0)
#' res1 <- L0TFinv.fix(y=WaveData$y, k=5, q=1, first=0, last=0.99)
#' print(res1$Ak)
#' print(WaveData$setA)
#' plot(WaveData$x, WaveData$y, xlab="", ylab="")
#' lines(WaveData$x, WaveData$y0, col = "red")
#' lines(WaveData$x, res1$yk, col = "lightgreen")
#'
#' @seealso \code{\link{L0TFinv.opt}}
#' @export
L0TFinv.fix <- function(y=y, k=k, q=q, first=0, last=1){
  if ( !(q %in% c(0,1)) ){
    stop("The order is not supported")
  }
  if( first<0 | last>1 ){
    stop("
The maximum possible range for the change points should be within [0,1]")
  }
  if( length(y) <= k+q+1 ){
    stop("The number of data points should be greater than the number of change points")
  }
  res <- InverseL0TF(y=y, kmax=k, q=q, first=first, last=last)
  betak = res$beta.all[,k]
  Ak = sort(res$A.all[[k]])
  yk = res$y.all[,k]
  G = list(sic=res$sic,bic=res$bic,mse=res$mse,
           y=y,
           betak=betak,yk=yk,Ak=Ak,
           beta.all=res$beta.all,y.all=res$y.all,A.all=res$A.all)
  class(G) <- "L0TFinvfix"
  return(G)
}



#' @title The inverse L0 trend filtering with optimal change points
#' @description Fit the input data points to a piecewise constant or piecewise linear trend with optimal change points.
#' @param y The input data points
#' @param kmax The maximum number of change points
#' @param q 0 or 1. Correspond to a piecewise constant or piecewise linear trend
#' @param first The value ranges from 0 to 1. Represent the minimum percentile point where a change point may occur. If 'first' = 0.01, it means that change points cannot appear in the first 1\% of the data points. If 'first' = 0, there is no constraint on the position of the change point.
#' @param last The value ranges from 0 to 1. Represent the maximum percentile point where a change point may occur. If 'last' = 0.99, it means that change points cannot appear in the last 1\% of the data points. If 'last' = 1, there is no constraint on the position of the change point.
#' @param penalty 'sic' or 'bic' penalty
#' @return
#' An S3 object of type "L0TFinvopt". A list containing the fitted trend results:
#' \item{sic}{Information criterion value with a penalty term of \eqn{2\log(\log(n)) \times \log(n)}}
#' \item{bic}{Information criterion value with a penalty term of \eqn{2 \times \log(n)}}
#' \item{mse}{The mean square error between the fitted trend and the input data}
#' \item{y}{The input data points}
#' \item{betaopt}{The fitted \eqn{\hat{\boldsymbol{\beta}}} coefficients with optimal change points }
#' \item{yopt}{The fitted trend with optimal change points}
#' \item{Aopt}{The set of position indicators of the fitted change points with optimal change points}
#' \item{kopt}{Optimal number of change points}
#' \item{beta.all}{A data frame with dimensions \eqn{n \times k_{\text{max}}}, where each column represents the fitted \eqn{\hat{\boldsymbol{\beta}}} coefficients corresponding to a given number of change points}
#' \item{y.all}{A data frame with dimensions \eqn{n \times k_{\text{max}}}, where each column represents the fitted estimated trend corresponding to a given number of change points}
#' \item{A.all}{A list of length \eqn{k_{\text{max}}}, where each element corresponds to the set of position indicators of change points under a given number }
#'
#' @details
#' Let the fitted trend be denoted as \eqn{\hat{\boldsymbol{y}}}, then \deqn{\text{sic} = n \times \log(\frac{1}{n}\|\boldsymbol{y}-\hat{\boldsymbol{y}}\|_2^2)+2\log(\log(n)) \times \log(n) \times \text{df}(\hat{\boldsymbol{y}})} and \deqn{\text{bic} = n \times \log(\frac{1}{n}\|\boldsymbol{y}-\hat{\boldsymbol{y}}\|_2^2)+2\times \log(n) \times \text{df}(\hat{\boldsymbol{y}}).}
#' The term \eqn{\text{df}(\hat{\boldsymbol{y}})} represents the degrees of freedom for the estimated trend, where \eqn{\text{df}(\hat{\boldsymbol{y}})=k+q+1}. Here, \eqn{k} refers to the number of change points in the estimated trend.
#'
#' @examples
#' tau = c(0.1, 0.3, 0.4, 0.7, 0.85)
#' h = c(-1, 5, 3, 0, -1, 2)
#' BlocksData <- SimuBlocksInv(n = 500, sigma = 0.2, seed = 50, tau = tau ,h = h)
#' res <- L0TFinv.opt(y=BlocksData$y, kmax=20, q=0, first=0.01, last=1, penalty="bic")
#' print(res$Aopt)
#' print(BlocksData$setA)
#' plot(BlocksData$x, BlocksData$y, xlab="", ylab="")
#' lines(BlocksData$x, BlocksData$y0, col = "red")
#' lines(BlocksData$x, res$yopt, col = "lightgreen")
#'
#' tau1 = c(0.4, 0.6, 0.7)
#' h1 = c(-3, 5, -4, 6)
#' a0 = -10
#' WaveData <- SimuWaveInv(n = 500, sigma = 0.1, seed = 50, tau = tau1, h = h1, a0 = a0)
#' res1 <- L0TFinv.opt(y=WaveData$y, kmax=10, q=1, first=0, last=0.99, penalty="sic")
#' print(res1$Aopt)
#' print(WaveData$setA)
#' plot(WaveData$x, WaveData$y, xlab="", ylab="")
#' lines(WaveData$x, WaveData$y0, col = "red")
#' lines(WaveData$x, res1$yopt, col = "lightgreen")
#'
#' @export
L0TFinv.opt <- function(y=y, kmax=kmax, q=q, first=0, last=1, penalty="bic"){
  if ( !(q %in% c(0,1)) ){
    stop("The specified order is not supported")
  }
  if ( !(penalty %in% c("bic","sic")) ){
    stop("The specified penalty is not supported")
  }
  if( first<0 | last>1 ){
    stop("
The maximum possible range for the change points should be within [0,1]")
  }
  if( length(y) <= kmax+q+1 ){
    stop("The number of data points should be greater than the number of change points")
  }
  res <- InverseL0TF(y=y, kmax=kmax, q=q, first=first, last=last)
  if(penalty == "bic"){
    kopt = which.min(res$bic)
  }
  if(penalty == "sic"){
    kopt = which.min(res$sic)
  }
  betaopt = res$beta.all[,kopt]
  Aopt = sort(res$A.all[[kopt]])
  yopt = res$y.all[,kopt]
  G = list(sic=res$sic,bic=res$bic,mse=res$mse,
           y=y,
           betaopt=betaopt,yopt=yopt,Aopt=Aopt,kopt=kopt,
           beta.all=res$beta.all,y.all=res$y.all,A.all=res$A.all)
  class(G) <- "L0TFinvopt"
  return(G)
}

               

