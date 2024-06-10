#' Printing method for the Lorenz Regression
#'
#' \code{print.LR} prints the arguments, explained Gini coefficient and estimated coefficients of an object of class \code{"LR"}.
#'
#' @param x An object of class \code{"LR"}.
#' @param digits The number of significant digits to be passed.
#' @param ... Additional arguments.
#'
#' @return No return value, called for printing an object of class \code{"LR"} to the console.
#'
#' @seealso \code{\link{Lorenz.Reg}}
#'
#' @examples
#' data(Data.Incomes)
#' NPLR <- Lorenz.Reg(Income ~ ., data = Data.Incomes, penalty = "none")
#' print(NPLR)
#'
#' @method print LR
#' @export

print.LR <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

  if (!inherits(x, "LR")) stop("The object must be of class 'LR'")

  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Explained Gini coefficient:", sprintf(paste0("%.", digits, "f"), x$Gi.expl), "\n")
  cat("\nCoefficients:\n")
  print.default(format(coef(x), digits = digits), print.gap = 2L,
                quote = FALSE)

}
