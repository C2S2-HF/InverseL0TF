
#' @title Extract estimated trends
#' @description Extract the coefficients of the estimated trends under the constraint of a given number of change points.
#' @param object The output of L0TFinvfix or L0TFinvopt
#' @param k The given number of change points
#' @method coef L0TFinvfix
#' @seealso \code{\link{print.L0TFinvfix}}
#' @export
coef.L0TFinvfix <- function(object, k=NULL) {
  if(is.null(k)){
    return(list(A.all=object$A.all,beta.all=object$beta.all,y.all=object$y.all))
  }else{
    if(k > ncol(object$beta.all)){
      stop("The number of given change points exceeds the maximum range")
    }else{
      return(list(A=object$A.all[[k]],beta=object$beta.all[,k],yhat=object$y.all[,k]))
    }
  }
}

#' @rdname coef.L0TFinvfix
#' @method coef L0TFinvopt
#' @export
coef.L0TFinvopt <- function(object, k=NULL) {
  if(is.null(k)){
    return(list(A.all=object$A.all,beta.all=object$beta.all,y.all=object$y.all))
  }else{
    if(k > ncol(object$beta.all)){
      stop("The number of given change points exceeds the maximum range")
    }else{
      return(list(A=object$A.all[[k]],beta=object$beta.all[,k],yhat=object$y.all[,k]))
    }
  }
}
