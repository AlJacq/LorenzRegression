#' Plots for the Unpenalized Lorenz Regression
#'
#' \code{plot.LR} provides plots for an object of class \code{LR}.
#'
#' @param LR Output of a call to \code{\link{Lorenz.Reg}}, where \code{penalty=="none"}.
#'
#' @return The Lorenz curve of the response and concentration curve of the response with respect to the estimated index
#'
#' @seealso \code{\link{Lorenz.Reg}}
#'
#' @examples
#' data(Data.Incomes)
#' NPLR <- Lorenz.Reg(Income ~ ., data = Data.Incomes, penalty = "none")
#' plot.PLR(NPLR)
#'
#' @import ggplot2
#'
#' @method plot LR
#' @export

plot.LR <- function(LR){

  p0 <- Lorenz.graphs(Response ~ ., LR$Fit, weights = LR$weights)
  p0 <- p0 + ggtitle("Observed and explained inequality")

  p0

}
