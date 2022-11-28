#' Concentration index of \emph{y} wrt \emph{x}
#'
#' \code{Gini.coef} computes the concentration index of a vector \emph{y} with respect to another vector \emph{x}.
#' If \emph{y} and \emph{x} are identical, the obtained concentration index boils down to the Gini coefficient.
#'
#' @param y variable of interest.
#' @param x variable to use for the ranking. By default \eqn{x=y}, and the obtained concentration index is the Gini coefficient of \emph{y}.
#' @param na.rm should missing values be deleted. Default value is \code{TRUE}. If \code{FALSE} is selected, missing values generate an error message
#' @param ties.method What method should be used to break the ties in the rank index. Possible values are "mean" (default value) or "random". If "random" is selected, the ties are broken by further ranking in terms of a uniformly distributed random variable. If "mean" is selected, the average rank method is used.
#' @param seed seed imposed for the generation of the vector of uniform random variables used to break the ties. Default is NULL, in which case no seed is imposed.
#' @param weights vector of sample weights. By default, each observation is given the same weight.
#'
#'
#' @return The value of the concentration index (or Gini coefficient)
#'
#' @seealso \code{\link{Lorenz.curve}}, \code{\link{Lorenz.graphs}}
#'
#' @examples
#' data(Data.Incomes)
#' # We first compute the Gini coefficient of Income
#' Y <- Data.Incomes$Income
#' Gini.coef(y = Y)
#' # Then we compute the concentration index of Income with respect to Age
#' X <- Data.Incomes$Age
#' Gini.coef(y = Y, x = X)
#'
#' @export

Gini.coef <- function(y, x=y, na.rm=T, ties.method=c("mean","random"), seed=NULL, weights=NULL){

  tol <- floor(sqrt(.Machine$double.digits))
  x <- round(x,tol)

  ties.method <- match.arg(ties.method)

  n <- length(y)

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

  if(any(weights<0)) stop("Weights must be nonnegative")

  if(is.null(weights)){
    weights <- rep(1,n)
  }
  pi <- weights/sum(weights)

  if (ties.method == "random"){

    if(!is.null(seed)) set.seed(seed)
    V<-stats::runif(n)

    y <- y[order(x,V)]
    pi <- pi[order(x,V)]
    F_i <- cumsum(pi) - 0.5*pi # Ensures that sum(F_i*pi) = 0.5
    y_mean <- pi%*%y

  }

  if (ties.method == "mean"){

    x_k <- sort(unique(x))
    pi_k <- sapply(1:length(x_k),function(k)sum(pi[x==x_k[k]]))
    F_k <- cumsum(pi_k) - 0.5*pi_k
    F_i <- sapply(1:length(x),function(i)sum(F_k[x_k==x[i]])) # Ensures that sum(F_i*pi) = 0.5
    y_mean <- pi%*%y

  }

  Gini <- 2*((pi*y)%*%F_i)/y_mean - 1

  as.numeric(Gini)
}


