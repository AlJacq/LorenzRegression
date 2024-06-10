#' Printing method for the summary of a Penalized Lorenz Regression
#'
#' \code{print.summary.PLR} provides a printing method for an object of class \code{"summary.PLR"}.
#'
#' @param object An object of class \code{"summary.PLR"}.
#' @param ... Additional arguments
#'
#' @return No return value, called for printing an object of class \code{"PLR"} to the console.
#'
#' @seealso \code{\link{summary.PLR}}
#'
#' @examples
#' data(Data.Incomes)
#' NPLR <- Lorenz.Reg(Income ~ ., data = Data.Incomes, penalty = "none")
#' summary(NPLR)
#'
#' @method print summary.PLR
#' @export

print.summary.PLR <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

  if (!inherits(x, "summary.PLR")) stop("The object must be of class 'summary.PLR'")

  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")

  cat("Explained Inequality Table\n\n")

  printCoefmat(x$ineq, digits = digits, na.print = "NA", ...)

  cat("\n (Explained inequality is measured by the explained Gini coefficient. Total inequality is the Gini coefficient of the response. The Lorenz-R2 is the proportion of inequality explained by the model.")

  cat("\n\nEstimated Coefficients Table for the Single-Index Model\n")

  printCoefmat(x$coefficients, digits = digits, na.print = "NA", ...)

}
