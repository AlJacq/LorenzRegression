#' Printing method for the summary of a Lorenz Regression
#'
#' \code{print.summary.LR} provides a printing method for an object of class \code{"summary.LR"}.
#'
#' @param object An object of class \code{"summary.LR"}.
#' @param ... Additional arguments
#'
#' @return No return value, called for printing an object of class \code{"LR"} to the console.
#'
#' @seealso \code{\link{summary.LR}}
#'
#' @examples
#' data(Data.Incomes)
#' NPLR <- Lorenz.Reg(Income ~ ., data = Data.Incomes, penalty = "none")
#' summary(NPLR)
#'
#' @method print summary.LR
#' @export

print.summary.LR <- function(x, digits = max(3L, getOption("digits") - 3L), signif.stars = getOption("show.signif.stars"), ...){

  if (!inherits(x, "summary.LR")) stop("x must be of class 'summary.LR'")

  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")

  cat("Explained Inequality Table\n\n")

  printCoefmat(x$ineq, digits = digits, na.print = "NA", ...)

  cat("\n (Explained inequality is measured by the explained Gini coefficient. Total inequality is the Gini coefficient of the response. The Lorenz-R2 is the proportion of inequality explained by the model.")

  cat("\n\nCoefficients Table for the Single-Index Model\n")

  printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars,
               na.print = "NA", ...)

}
