#' Plots for the Unpenalized Lorenz Regression
#'
#' \code{plot.LR} provides plots for an object of class \code{LR}.
#'
#' @param x Output of a call to \code{\link{Lorenz.Reg}}, where \code{penalty=="none"}.
#' @param ... Additional arguments
#'
#' @return The Lorenz curve of the response and concentration curve of the response with respect to the estimated index
#'
#' @seealso \code{\link{Lorenz.Reg}}
#'
#' @examples
#' data(Data.Incomes)
#' NPLR <- Lorenz.Reg(Income ~ ., data = Data.Incomes, penalty = "none")
#' plot(NPLR)
#'
#' @import ggplot2
#'
#' @method plot LR
#' @export

plot.LR <- function(x, ...){

  p0 <- Lorenz.graphs(Response ~ ., x$Fit, weights = x$weights)
  p0 <- p0 + ggtitle("Observed and explained inequality")

  p0

}
