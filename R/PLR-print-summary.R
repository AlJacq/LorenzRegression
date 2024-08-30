#' Printing method for the summary of a Penalized Lorenz Regression
#'
#' \code{print.summary.PLR} provides a printing method for an object of class \code{"summary.PLR"}.
#'
#' @aliases print.summary.PLR_boot print.summary.PLR_cv
#' @param x An object of class \code{"summary.PLR"}. The object might also have S3 class \code{"summary.PLR_boot"} and/or \code{"summary.PLR_cv"} (both inherit from class \code{"summary.LR"})
#' @param digits Number of significant digits to be passed.
#' @param ... Additional arguments passed to the function \code{\link{print}}.
#'
#' @return No return value, called for printing an object of class \code{"summary.PLR"} to the console.
#'
#' @seealso \code{\link{summary.PLR}}
#'
#' @examples
#' ## For examples see example(Lorenz.Reg), example(Lorenz.boot) and example(PLR.CV)
#'
#' @method print summary.PLR
#' @export

print.summary.PLR <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")

  cat("Explained Inequality Table\n\n")

  print(x$ineq, digits = digits)

  cat("\n (Explained inequality is measured by the explained Gini coefficient. Total inequality is the Gini coefficient of the response. The Lorenz-R2 is the proportion of inequality explained by the model.)")

  cat("\n\nEstimated Coefficients Table for the Single-Index Model\n")

  print(x$coefficients, digits = digits)

}

#' @method print summary.PLR_boot
#' @export

print.summary.PLR_boot <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  NextMethod("print")
}

#' @method print summary.PLR_cv
#' @export

print.summary.PLR_cv <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  NextMethod("print")
}
