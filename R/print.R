
#' @title Print L0TFinvfix or L0TFinvopt object
#' @description Prints a summary of L0TFinvfix or L0TFinvopt
#' @param x The output of L0TFinvfix or L0TFinvopt
#' @param ... ignore
#' @method print L0TFinvfix
#' @examples
#' library(ggplot2)
#'
#' tau = c(0.1, 0.3, 0.4, 0.7, 0.85)
#' h = c(-1, 5, 3, 0, -1, 2)
#' BlocksData <- SimuBlocksInv(n = 500, sigma = 0.2, seed = 50, tau = tau ,h = h)
#' res <- L0TFinv.opt(y=BlocksData$y, kmax=10, q=0, first=0.01, last=1, penalty="bic")
#' print(res)
#' coef(res,k=res$kopt)
#' plot(res,type="yhat")
#' plot(res,type="bic")
#'
#' tau1 = c(0.1, 0.3, 0.4, 0.7, 0.85)
#' h1 = c(-1, 5, 3, 0, -1, 2)
#' a0 = -10
#' WaveData <- SimuWaveInv(n = 2000, sigma = 0.1, seed = 50, tau = tau1, h = h1, a0 = a0)
#' res1 <- L0TFinv.fix(y=WaveData$y, k=20, q=1, first=0, last=0.99)
#' print(res1)
#' coef(res1,k=5)
#' plot(res1,type="yhat",k=5)
#' plot(res1,type="mse")
#'
#' @export
print.L0TFinvfix <- function(x, ...){
  return(list(sic=x$sic,bic=x$bic,mse=x$mse,
              A.all=x$A.all,beta.all=x$beta.all,y.all=x$y.all))
}

#' @rdname print.L0TFinvfix
#' @method print L0TFinvopt
#' @export
print.L0TFinvopt <- function(x, ...){
  return(list(sic=x$sic,bic=x$bic,mse=x$mse,kopt=x$kopt,
              A.all=x$A.all,beta.all=x$beta.all,y.all=x$y.all))
}


