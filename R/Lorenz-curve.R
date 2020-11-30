#' Concentration curve of \emph{y} with respect to \emph{x}
#'
#' \code{Lorenz.curve} computes the concentration curve index of a vector \emph{y} with respect to another vector \emph{x}.
#' If \emph{y} and \emph{x} are identical, the obtained concentration curve boils down to the Lorenz curve.
#'
#' @param y variable of interest.
#' @param x variable to use for the ranking. By default \eqn{x=y}, and the obtained concentration curve is the Lorenz curve of \emph{y}.
#' @param graph whether a graph of the obtained concentration curve should be traced. Default value is FALSE.
#'
#' @return A function corresponding to the estimated Lorenz or concentration curve. If \code{graph} is TRUE, the curve is also plotted.
#'
#' @seealso \code{\link{Lorenz.graphs}}, \code{\link{Gini.coef}}
#'
#' @examples
#' data(Data.Incomes)
#' # We first compute the Lorenz curve of Income
#' Y <- Data.Incomes$Income
#' Lorenz.curve(y = Y, graph = TRUE)
#' # Then we compute the concentration curve of Income with respect to Age
#' X <- Data.Incomes$Age
#' Lorenz.curve(y = Y, x = X, graph = TRUE)
#'
#' @import ggplot2
#'
#' @export

Lorenz.curve <- function(y, x=y, graph=F){

  p <- NULL
  x <- as.vector(x); y <- as.vector(y)
  o <- order(x)
  x <- x[o] ; y <- y[o]
  n <- length(x)
  if (n < 1)
    stop("'x' must have 1 or more non-missing values")
  x_unique <- as.numeric(names(table(x)))
  tolerance <- .Machine$double.eps^0.5
  y_star <- c(0,sapply(1:length(x_unique),function(i)mean(y[abs(x-x_unique[i])<tolerance])))
  pi<-c(0,as.vector(table(x)/n))
  rval <- stats::approxfun(cumsum(pi),cumsum(pi*y_star/mean(y)) ,method = "linear", yleft = 0, yright = 1,ties = "ordered")

    if(graph){
    if(all.equal(x,y)==T){
      txt.title <- "Lorenz curve of y"
    }else{
      txt.title <- "Concentration curve of y wrt x"
    }
    print(ggplot2::ggplot(data.frame(p=c(0,1)), aes(p)) +
      stat_function(fun=function(p)rval(p), geom="line") +
      stat_function(fun=function(p)p, geom="line") +
      labs(x = "Cumulative share of the population",y = "Cumulative share of y",title = txt.title))
    }

  return(rval)

}

