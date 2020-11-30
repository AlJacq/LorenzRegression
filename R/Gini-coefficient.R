#' Concentration index of \emph{y} wrt \emph{x}
#'
#' \code{Gini.coef} computes the concentration index of a vector \emph{y} with respect to another vector \emph{x}.
#' If \emph{y} and \emph{x} are identical, the obtained concentration index boils down to the Gini coefficient.
#'
#' @param y variable of interest.
#' @param x variable to use for the ranking. By default \eqn{x=y}, and the obtained concentration index is the Gini coefficient of \emph{y}.
#' @param na.rm should missing values be deleted. Default value is \code{TRUE}. If \code{FALSE} is selected, missing values generate an error message
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

Gini.coef <- function(y, x=y, na.rm=T){

  if(sum(is.na(c(x,y)))>0){
    if(na.rm){
     x.tmp <- x[!(is.na(x) | is.na(y))]
     y.tmp <- y[!(is.na(x) | is.na(y))]
     x <- x.tmp ; y <- y.tmp
    }else{
      stop("There are missing values in either x or y and na.rm is FALSE")
    }
  }

  n <- length(x)
  o.x <- order(x)
  Gini <- 2*(1:n%*%y[o.x])/n^2/mean(y) - (n+1)/n

  as.numeric(Gini)
}


