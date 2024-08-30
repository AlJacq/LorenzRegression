#' Printing method for the summary of a Lorenz Regression
#'
#' \code{print.summary.LR} provides a printing method for an object of class \code{"summary.LR"}.
#'
#' @aliases print.summary.LR_boot
#' @param x An object of class \code{"summary.LR"}. The object might also have S3 class \code{"summary.LR_boot"} (which inherits from class \code{"summary.LR"})
#' @param digits Number of significant digits to be passed.
#' @param signif.stars Logical determining whether p-values should be also encoded visually. See the help of the function \code{\link{printCoefmat}} for more information.
#' This is only relevant if \code{x} inherits from \code{"summary.LR_boot"}.
#' @param ... Additional arguments passed to the function \code{\link{print}}.
#'
#' @return No return value, called for printing an object of class \code{"LR"} to the console.
#'
#' @seealso \code{\link{summary.LR}}
#'
#' @examples
#' ## For examples see example(Lorenz.Reg) and example(Lorenz.boot)
#'
#' @importFrom stats printCoefmat
#'
#' @method print summary.LR
#' @export

print.summary.LR <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")

  if(!is.null(x$ineq)){

    cat("Explained Inequality Table\n\n")

    print(x$ineq, digits = digits,...)

    cat("\n (Explained inequality is measured by the explained Gini coefficient. Total inequality is the Gini coefficient of the response. The Lorenz-R2 is the proportion of inequality explained by the model.)")

    cat("\n\nCoefficients Table for the Single-Index Model\n")

    printCoefmat(x$coefficients, digits = digits, na.print = "NA", ...)

  }else(

    cat("The model is empty. Therefore, explained inequality measured by the explained Gini coefficient is 0 and no single-index model is estimated.")

  )

}

#' @method print summary.LR_boot
#' @rdname print.summary.LR
#' @export

print.summary.LR_boot <- function(x, digits = max(3L, getOption("digits") - 3L), signif.stars = getOption("show.signif.stars"), ...){

  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")

  if(!is.null(x$ineq)){

    cat("Explained Inequality Table\n\n")

    print(x$ineq, digits = digits,...)

    cat("\n (Explained inequality is measured by the explained Gini coefficient. Total inequality is the Gini coefficient of the response. The Lorenz-R2 is the proportion of inequality explained by the model.)")

    cat("\n\nCoefficients Table for the Single-Index Model\n")

    printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars,
                 na.print = "NA", ...)

  }else(

    cat("The model is empty. Therefore, explained inequality measured by the explained Gini coefficient is 0 and no single-index model is estimated.")

  )

}

