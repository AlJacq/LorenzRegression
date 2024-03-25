#' Printing method for the Lorenz Regression
#'
#' \code{print.LR} prints the arguments and estimated coefficients of an object of class \code{LR}.
#'
#' @param x Output of a call to \code{\link{Lorenz.Reg}}, where \code{penalty="none"}.
#' @param ... Additional arguments.
#'
#' @return No return value, called for printing an object of class \code{LR} to the console
#'
#' @seealso \code{\link{Lorenz.Reg}}
#'
#' @examples
#' data(Data.Incomes)
#' NPLR <- Lorenz.Reg(Income ~ ., data = Data.Incomes, penalty = "none")
#' print(NPLR)
#'
#' @import knitr
#'
#' @method print LR
#' @export

print.LR <- function(x, ...){

  cat("Call",
      x$call,
      sep="\n",
      "",
      "Coefficients",
      knitr::kable(x$theta))

}
