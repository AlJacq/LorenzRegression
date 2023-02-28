#' Printing method for the Penalized Lorenz Regression
#'
#' \code{print.PLR} prints the arguments and estimated coefficients of an object of class \code{PLR}.
#'
#' @param x Output of a call to \code{\link{Lorenz.Reg}}, where \code{penalty!="none"}.
#' @param ... Additional arguments.
#'
#' @return No return value, called for printing an object of class \code{PLR} to the console
#'
#' @seealso \code{\link{Lorenz.Reg}}
#'
#' @examples
#' data(Data.Incomes)
#' PLR <- Lorenz.Reg(Income ~ ., data = Data.Incomes, penalty = "SCAD",
#'                   sel.choice = c("BIC","CV"), h.grid = nrow(Data.Incomes)^(-1/5.5),
#'                   eps = 0.01, seed.CV = 123, nfolds = 5)
#' print(PLR)
#'
#' @import knitr
#'
#' @method print PLR
#' @export

print.PLR <- function(x, ...){

  cat("Call",
      x$call,
      sep="\n",
      "",
      "Coefficients",
      knitr::kable(t(x$theta)))

}
