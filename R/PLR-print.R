#' Printing method for the Penalized Lorenz Regression
#'
#' \code{print.PLR} prints the arguments, explained Gini coefficient and estimated coefficients of an object of class \code{"PLR"}.
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
                quote = FALSE)

}

#' @method print PLR_boot
#' @rdname print.PLR
#' @export

print.PLR_boot <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

  NextMethod("print")


  cat("\nExplained Gini coefficient (Bootstrap selection method):", sprintf(paste0("%.", digits, "f"), ineqExplained.PLR_boot(x,pars.idx="Boot")), "\n")
  cat("\nCoefficients (Bootstrap selection method):\n")
  print.default(format(coef.PLR_boot(x,pars.idx="Boot"), digits = digits), print.gap = 2L,
                quote = FALSE)

}

#' @method print PLR_cv
#' @rdname print.PLR
#' @export

print.PLR_cv <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

  NextMethod("print")

  cat("\nExplained Gini coefficient (Cross-validation selection method):", sprintf(paste0("%.", digits, "f"), ineqExplained.PLR_cv(x,pars.idx="CV")), "\n")
  cat("\nCoefficients (Cross-validation selection method):\n")
  print.default(format(coef.PLR_cv(x,pars.idx="CV"), digits = digits), print.gap = 2L,
                quote = FALSE)

}
