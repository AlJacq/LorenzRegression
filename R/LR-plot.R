#' Plots for the Lorenz Regression
#'
#' \code{plot.LR} provides plots for an object of class \code{"LR"}.
#'
#' @param x An object of class \code{"LR"}.
#' @param ... Additional arguments
#'
#' @return The Lorenz curve of the response and concentration curve of the response with respect to the estimated index. The graph is as an object of class \code{"ggplot"}.
#'
#' @seealso \code{\link{Lorenz.Reg}}
#'
#' @examples
#' ## For examples see example(Lorenz.Reg)
#'
#' @import ggplot2
#'
#' @method plot LR
#' @export

plot.LR <- function(x, ...){

  if (!inherits(x, "LR")) stop("x must be of class 'LR'")
  if (is.null(x$theta)) stop("No plots are available for an empty model.")

  formula <- update.formula(x, . ~ index)
  data <- data.frame(x$y,x$index)
  names(data) <- all.vars(formula)

  g <- Lorenz.graphs(formula, data, weights = x$weights)
  g <- g + ggtitle("Observed and explained inequality")

  g

}
