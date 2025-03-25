

#' @title Plot L0TFinvfix or L0TFinvopt object
#'
#' @description Plots a summary of L0TFinvfix or L0TFinvopt
#' @param x The output of L0TFinvfix or L0TFinvopt
#' @param type The values are taken as c("\eqn{mse}", "\eqn{sic}", "\eqn{bic}", "\eqn{yhat}"). If \eqn{type} is "\eqn{mse}", plot the mse as it changes with change points.
#' The same applies to "\eqn{sic}" and "\eqn{bic}". If \eqn{type} is "\eqn{yhat}", plot the trend of the estimated values against the input data.
#' @param k Only used for \eqn{type} = "\eqn{yhat}". The given number of change points.
#' By default, the L0TFinvfix object outputs the estimated trend that corresponds to the fixed number of change points within the model. Conversely, the L0TFinvopt object provides the estimated trend based on the optimal change points.
#' @param ... ignore
#' @method plot L0TFinvfix
#' @import ggplot2
#' @seealso \code{\link{print.L0TFinvfix}}
#' @export
plot.L0TFinvfix <- function(x, type = NULL, k = NULL, ...){
  if ( !(type %in% c("mse","sic","bic","yhat")) ){
    stop("The specified type is not supported")
  }
  if(type == "mse"){
    num = 1:(length(x$mse))
    val = x$mse
    plotObject <- ggplot(data.frame(num, val), aes(x = num, y = val)) +
      geom_line() +
      geom_point() +
      labs(x = "The number of change points", y = "MSE") +
      theme_minimal()
  }
  if(type == "sic"){
    num = 1:(length(x$sic))
    val = x$sic
    plotObject <- ggplot(data.frame(num, val), aes(x = num, y = val)) +
      geom_line() +
      geom_point() +
      labs(x = "The number of change points", y = "SIC") +
      theme_minimal()
  }
  if(type == "bic"){
    num = 1:(length(x$bic))
    val = x$bic
    plotObject <- ggplot(data.frame(num, val), aes(x = num, y = val)) +
      geom_line() +
      geom_point() +
      labs(x = "The number of change points", y = "BIC") +
      theme_minimal()
  }
  if(type == "yhat"){
    n = length(x$y)
    num = 1:n
    if(is.null(k)){
      k = length(x$mse)
      val = x$y.all[,k]
    }else{
      if(k > ncol(x$beta.all)){
        stop("The number of given change points exceeds the maximum range")
      }
      val = x$y.all[,k]
    }
    input = x$y
    df <- data.frame(num = num, input = input, val = val)
    plotObject <- ggplot(data = df, aes(x = num)) +
      geom_point(aes(y = input), color = "black") +
      geom_line(aes(y = val), color = "lightgreen") +
      labs(x = "Position indicators", y = "Value") +
      theme_minimal()
  }
  return(plotObject)
}


#' @rdname plot.L0TFinvfix
#' @method plot L0TFinvopt
#' @export
plot.L0TFinvopt <- function(x, type = "mse", k = NULL, ...){
  if ( !(type %in% c("mse","sic","bic","yhat")) ){
    stop("The specified type is not supported")
  }
  if(type == "mse"){
    num = 1:(length(x$mse))
    val = x$mse
    plotObject <- ggplot(data.frame(num, val), aes(x = num, y = val)) +
      geom_line() +
      geom_point() +
      labs(x = "The number of change points", y = "MSE") +
      theme_minimal()
  }
  if(type == "sic"){
    num = 1:(length(x$sic))
    val = x$sic
    plotObject <- ggplot(data.frame(num, val), aes(x = num, y = val)) +
      geom_line() +
      geom_point() +
      labs(x = "The number of change points", y = "SIC") +
      theme_minimal()
  }
  if(type == "bic"){
    num = 1:(length(x$bic))
    val = x$bic
    plotObject <- ggplot(data.frame(num, val), aes(x = num, y = val)) +
      geom_line() +
      geom_point() +
      labs(x = "The number of change points", y = "BIC") +
      theme_minimal()
  }
  if(type == "yhat"){
    n = length(x$y)
    num = 1:n
    if(is.null(k)){
      k = x$kopt
      val = x$y.all[,k]
    }else{
      if(k > ncol(x$beta.all)){
        stop("The number of given change points exceeds the maximum range")
      }
      val = x$y.all[,k]
    }
    input = x$y
    df <- data.frame(num = num, input = input, val = val)
    plotObject <- ggplot(data = df, aes(x = num)) +
      geom_point(aes(y = input), color = "black") +
      geom_line(aes(y = val), color = "lightgreen") +
      labs(x = "Position indicators", y = "Value") +
      theme_minimal()
  }
  return(plotObject)
}
