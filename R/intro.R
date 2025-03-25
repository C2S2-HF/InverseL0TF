#' @_PACKAGE
#'
#' @name L0TFinv-package
#' @title A package for L0-regularized sparse approximation
#' @description
#' Trend filtering is a typical method for nonparametric regression.
#' The commonly used trend filtering models is the L1 trend filtering model \eqn{(a)} based on the difference matrix \eqn{\boldsymbol{D}^{(q+1)}}, as illustrated below.
#' \deqn{\min _{\boldsymbol{\alpha} \in \mathbb{R}^n} \frac{1}{2}\|\boldsymbol{y}-\boldsymbol{\alpha}\|_2^2 + \lambda\|\boldsymbol{D}^{(q+1)} \boldsymbol{\alpha}\|_{\ell_1}, \quad q=0,1,2, \ldots. \quad (a) }
#' L0 trend filtering \eqn{(b)} has a advantage over other trend filtering methods, especially in the detection of change points.
#' The expression for L0 trend filtering is as follows:
#' \deqn{\min _{\boldsymbol{\alpha} \in \mathbb{R}^n} \frac{1}{2}\|\boldsymbol{y}-\boldsymbol{\alpha}\|_2^2 + \lambda\|\boldsymbol{D}^{(q+1)} \boldsymbol{\alpha}\|_{\ell_0}. \quad (b)  }
#' We explore transforming the problem \eqn{(b)} into a L0-regularized sparse format \eqn{(c)} by introducing an artificial design matrix \eqn{\boldsymbol{X}^{(q+1)}} that corresponds to the difference matrix, thereby reformulating the L0 trend filtering problem into the following format.
#' \deqn{\min _{\boldsymbol{\beta} \in \mathbb{R}^n} \frac{1}{2}\|\boldsymbol{y}-\boldsymbol{X}^{(q+1)}\boldsymbol{\beta}\|_2^2 + \lambda \sum_{i=q+2}^n |\boldsymbol{\beta}_i|_{\ell_0}. \quad (c) }
#' In our practical approach, we consider the maximum number of change points \eqn{k_{\text{max}}} as a constraint, transforming the aforementioned L0 penalty problem \eqn{(c)} into the following L0 constraint problem.
#' \deqn{\text{ minimize }\frac{1}{2}\|\boldsymbol{y}-\boldsymbol{X}^{(q+1)}\boldsymbol{\beta}\|_2^2,\quad \text{ subject to } \sum_{i=q+2}^n |\boldsymbol{\beta}_i|_{\ell_0} \leq k_{\text{max}}. \quad (d)}
#' For such L0 constraint problems \eqn{(d)}, we employ a splicing-based approach to design algorithms for processing.
#' This package has the following seven main methods:
#' \itemize{
#' \item{\strong{matrix with special structure }}{\eqn{\quad}Generate \eqn{\boldsymbol{X}^{(q+1)}} or \eqn{\boldsymbol{D}^{(q+1)}} matrix.}
#' \item{\strong{inverse of the crossprod matrix }}{\eqn{\quad}Simplify the calculation of the inverse matrix of \eqn{(\boldsymbol{X}^{(q+1)}_A)^T \boldsymbol{X}^{(q+1)}_A} for the cases where \eqn{q=0} or \eqn{q=1}, which is frequently used in splicing algorithms.}
#' \item{\strong{inverse L0 trend filtering with fixed change points }}{\eqn{\quad}Fit a piecewise constant or piecewise linear estimated trend with a given number of change points.}
#' \item{\strong{inverse L0 trend filtering with optimal change points }}{\eqn{\quad}Fit a piecewise constant or piecewise linear estimated trend with a maximum number of change points, and select the optimal estimated trend using appropriate information criteria.}
#' \item{\strong{simulated data }}{\eqn{\quad}Generate piecewise constant or piecewise linear data.}
#' \item{\strong{print/coef}}{\eqn{\quad}Print a summary of the trend estimation results.}
#' \item{\strong{plot }}{\eqn{\quad}Plot a summary of the trend estimation results.}
#' }
#' @details
#' \itemize{
#' \item{}{In previous studies, algorithms solving trend filtering problems \eqn{(a)} necessitate the computation of \eqn{((\boldsymbol{D}^{(q+1)})^T \boldsymbol{D}^{(q+1)})^{-1}}.
#' When \eqn{n} is large, just fitting the matrix into memory becomes an issue.}
#' \item{}{In L0 trend filtering \eqn{(b)}, the positions of non-zero elements in the L0 norm correspond with the locations of change points.
#' We consider two subsets: the active set \eqn{A} for non-zero elements and the inactive set \eqn{I} for zero elements.
#' Despite this, computing \eqn{((\boldsymbol{D}^{(q+1)}_I)^T \boldsymbol{D}^{(q+1)}_I)^{-1}} remains a task involving a substantial matrix.}
#' \item{}{Due to the connection between L0 constraint problems and L0 penalty problems, and considering that the sparsity of \eqn{\boldsymbol{\beta}} is is more meaningful in practical applications than the selection of the hyperparameter \eqn{\lambda}.
#' We focus on the constraint that reflects our aim to achieve an estimated trend with a given number of change points.
#' So we transform the L0 penalty problem \eqn{(c)} into the L0 constraint problem \eqn{(d)}.}
#' }
#' @references
#' Kim SJ, Koh K, Boyd SP and Gorinevsky DM. L1 Trend Filtering. Society for Industrial and Applied Mathematics (2009).
#'
#' Wen C, Wang X and Zhang A. L0 Trend Filtering. INFORMS Journal on Computing (2023).
NULL
