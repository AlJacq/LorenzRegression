#' Printing method for the Lorenz Regression
#'
#' Prints the arguments, explained Gini coefficient and estimated coefficients of an object of class \code{"LR"}.
#'
#' @aliases print.LR_boot
#' @param x An object of class \code{"LR"}.
#' @param digits The number of significant digits to be passed.
#' @param ... Additional arguments.
#'
#' @return No return value, called for printing an object of class \code{"LR"} to the console.
#'
#' @seealso \code{\link{Lorenz.Reg}}
#'
#' @examples
#' ## For examples see example(Lorenz.Reg)
#'
#' @method print LR
#' @export

print.LR <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Explained Gini coefficient:", sprintf(paste0("%.", digits, "f"), ineqExplained.LR(x)), "\n")
  cat("\nCoefficients:\n")
  print.default(format(coef.LR(x), digits = digits), print.gap = 2L,
                quote = FALSE, ...)

}

#' @method print LR_boot
#' @export

print.LR_boot <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  NextMethod("print")
}
