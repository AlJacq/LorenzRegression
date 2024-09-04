#' Printing method for the Penalized Lorenz Regression
#'
#' Prints the arguments, explained Gini coefficient and estimated coefficients of an object of class \code{"PLR"}.
#'
#' @aliases print.PLR_boot print.PLR_cv
#' @param x An object of S3 class \code{"PLR"}. The object might also have S3 classes \code{"PLR_boot"} and/or \code{"PLR_cv"} (both inherit from class \code{"PLR"})
#' @param digits The number of significant digits to be passed.
#' @param ... Additional arguments.
#'
#' @return No return value, called for printing an object of class \code{"PLR"} to the console.
#'
#' @seealso \code{\link{Lorenz.Reg}}
#'
#' @details The explained Gini coefficient and estimated coefficients are returned for each available selection method, depending on the class of \code{x}.
#'
#' @examples
#' ## For examples see example(Lorenz.Reg), example(Lorenz.boot) and example(PLR.CV)
#'
#' @method print PLR
#' @export

print.PLR <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Explained Gini coefficient (BIC selection method):", sprintf(paste0("%.", digits, "f"), ineqExplained.PLR(x)), "\n")
  cat("\nCoefficients (BIC selection method):\n")
  print.default(format(coef.PLR(x), digits = digits), print.gap = 2L,
                quote = FALSE, ...)

}

#' @method print PLR_boot
#' @export

print.PLR_boot <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

  NextMethod("print")

  pars.boot <- c(x$grid.idx["Boot"],x$lambda.idx["Boot"])
  cat("\nExplained Gini coefficient (Bootstrap selection method):", sprintf(paste0("%.", digits, "f"), ineqExplained_PLR(x,type="Gini.explained", pars.idx=pars.boot)), "\n")
  cat("\nCoefficients (Bootstrap selection method):\n")
  print.default(format(coef_PLR(x, renormalize=TRUE, pars.idx=pars.boot), digits = digits), print.gap = 2L,
                quote = FALSE, ...)

}

#' @method print PLR_cv
#' @export

print.PLR_cv <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

  NextMethod("print")

  pars.cv <- c(x$grid.idx["CV"],x$lambda.idx["CV"])
  cat("\nExplained Gini coefficient (Cross-validation selection method):", sprintf(paste0("%.", digits, "f"), ineqExplained_PLR(x, type="Gini.explained", pars.idx=pars.cv)), "\n")
  cat("\nCoefficients (Cross-validation selection method):\n")
  print.default(format(coef_PLR(x, renormalize = TRUE, pars.idx=pars.cv), digits = digits), print.gap = 2L,
                quote = FALSE, ...)

}
