#' Concentration curve of \emph{y} with respect to \emph{x}
#'
#' \code{Lorenz.curve} computes the concentration curve index of a vector \emph{y} with respect to another vector \emph{x}.
#' If \emph{y} and \emph{x} are identical, the obtained concentration curve boils down to the Lorenz curve.
#'
#' @param y variable of interest.
#' @param x variable to use for the ranking. By default \eqn{x=y}, and the obtained concentration curve is the Lorenz curve of \emph{y}.
#' @param na.rm should missing values be deleted. Default value is \code{TRUE}. If \code{FALSE} is selected, missing values generate an error message
#' @param ties.method What method should be used to break the ties in the rank index. Possible values are "mean" (default value) or "random". If "random" is selected, the ties are broken by further ranking in terms of a uniformly distributed random variable. If "mean" is selected, the average rank method is used.
#' @param seed seed imposed for the generation of the vector of uniform random variables used to break the ties. Default is NULL, in which case no seed is imposed.
#' @param weights vector of sample weights. By default, each observation is given the same weight.
#'
#' @return A function corresponding to the estimated Lorenz or concentration curve.
#'
#' @details The parameter \code{seed} allows for local seed setting to control randomness in the generation of the uniform random variables.
#' The specified seed is applied to the respective part of the computation, and the seed is reverted to its previous state after the operation.
#' This ensures that the seed settings do not interfere with the global random state or other parts of the code.
#'
#' @seealso \code{\link{Lorenz.graphs}}, \code{\link{Gini.coef}}
#'
#' @examples
#' data(Data.Incomes)
#' # We first compute the Lorenz curve of Income
#' Y <- Data.Incomes$Income
#' Lorenz.curve(y = Y)
#' # Then we compute the concentration curve of Income with respect to Age
#' X <- Data.Incomes$Age
#' Lorenz.curve(y = Y, x = X)
#'
#' @importFrom ggplot2 ggplot aes stat_function labs
#'
#' @export

Lorenz.curve <- function(y, x=y, na.rm=TRUE, ties.method=c("mean","random"), seed=NULL, weights=NULL){

  # 0. Preliminaries ----

  ties.method <- match.arg(ties.method)
  x <- as.vector(x); y <- as.vector(y)

  if(sum(is.na(c(x,y)))>0){
    if(na.rm){
      x.tmp <- x[!(is.na(x) | is.na(y))]
      y.tmp <- y[!(is.na(x) | is.na(y))]
      if (!is.null(weights)) weights <- weights[!(is.na(x) | is.na(y))]
      x <- x.tmp ; y <- y.tmp
    }else{
      stop("There are missing values in either x or y and na.rm is FALSE")
    }
  }

  n <- length(x)
  if (n < 1) stop("'x' must have 1 or more non-missing values")

  if(any(weights<0)) stop("Weights must be nonnegative")

  if(is.null(weights)){
    weights <- rep(1,n)
  }
  pi <- weights/sum(weights)

  # 1. Function ----

  if (ties.method == "random"){

    V <- runif_seed(n,seed = seed)

    y <- y[order(x,V)]
    pi <- pi[order(x,V)]
    y_mean <- as.numeric(pi%*%y)

    Fun <- stats::approxfun(cumsum(pi),cumsum(pi*y/y_mean) ,method = "linear", yleft = 0, yright = 1,ties = "ordered")

  }

  if (ties.method == "mean"){

    x_k <- sort(unique(x))
    pi_k <- sapply(1:length(x_k),function(k)sum(pi[x==x_k[k]]))
    y_k <- sapply(1:length(x_k),function(k)as.numeric(pi[x==x_k[k]]%*%y[x==x_k[k]]))

    y_mean <- as.numeric(pi%*%y)

    Fun <- stats::approxfun(c(0,cumsum(pi_k)),c(0,cumsum(y_k/y_mean)) ,method = "linear", yleft = 0, yright = 1,ties = "ordered")


  }

  return(Fun)

}

