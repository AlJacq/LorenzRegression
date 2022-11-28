#' Printing method for the Penalized Lorenz Regression
#'
#' \code{print.PLR} prints the arguments and estimated coefficients of an object of class \code{PLR}.
#'
#' @param PLR Output of a call to \code{\link{Lorenz.Reg}}, where \code{penalty!="none"}.
#'
#' @seealso \code{\link{Lorenz.Reg}}
#'
#' @examples
#' data(Data.Incomes)
#' PLR <- Lorenz.Reg(Income ~ ., data = Data.Incomes, penalty = "SCAD", lambda.choice = c("BIC","CV"), eps = 0.005, seed.CV = 123, nfolds = 5)
#' print.PLR(PLR)
#'
#' @import knitr
#'
#' @method print PLR
#' @export

print.PLR <- function(PLR){

  cat("Call",
      PLR$call,
      sep="\n",
      "",
      "Coefficients",
      knitr::kable(t(PLR$theta)))

}
