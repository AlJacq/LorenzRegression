#' Printing method for the Lorenz Regression
#'
#' \code{print.LR} prints the arguments and estimated coefficients of an object of class \code{LR}.
#'
#' @param LR Output of a call to \code{\link{Lorenz.Reg}}, where \code{penalty="none"}.
#'
#' @seealso \code{\link{Lorenz.Reg}}
#'
#' @examples
#' data(Data.Incomes)
#' NPLR <- Lorenz.Reg(Income ~ ., data = Data.Incomes, penalty = "none")
#' print.LR(NPLR)
#'
#' @import knitr
#'
#' @method print LR
#' @export

print.LR <- function(LR){

  cat("Call",
      LR$call,
      sep="\n",
      "",
      "Coefficients",
      knitr::kable(LR$theta))

}
