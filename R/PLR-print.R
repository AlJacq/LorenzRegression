#' Printing method for the Penalized Lorenz Regression
#'
#' \code{print.PLR} prints the arguments, explained Gini coefficient and estimated coefficients of an object of class \code{"PLR"}.
#'
#' @param x An object of class \code{"PLR"}.
#' @param digits The number of significant digits to be passed.
#' @param ... Additional arguments.
#'
#' @return No return value, called for printing an object of class \code{"PLR"} to the console.
#'
#' @seealso \code{\link{Lorenz.Reg}}
#'
#' @examples
#' data(Data.Incomes)
#' PLR <- Lorenz.Reg(Income ~ ., data = Data.Incomes, penalty = "SCAD",
#'                   sel.choice = c("BIC","CV"), h.grid = nrow(Data.Incomes)^(-1/5.5),
#'                   eps = 0.01, seed.CV = 123, nfolds = 5)
#' print(PLR)
#'
#' @method print PLR
#' @export

print.PLR <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

  if (!inherits(x, "PLR")) stop("x must be of class 'PLR'")

  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")

  coef_x <- coef(x)

  if(inherits(x,c("PLR_cv","PLR_boot"))){

    cat("\nSelection methods for tuning and penalty parameter:\n",paste0(names(coef_x),collapse=","),"\n")
    cat("\nExplained Gini coefficient:\n")
    print.default(format(x$Gi.expl, digits = digits), print.gap = 2L,
                  quote = FALSE)
    for (s in names(coef_x)){
      cat("\nCoefficients obtained by",s,":\n")
      print.default(format(coef_x[[s]], digits = digits), print.gap = 2L,
                    quote = FALSE)
    }
  }else{

    cat("Explained Gini coefficient:", sprintf(paste0("%.", digits, "f"), x$Gi.expl), "\n")
    cat("\nCoefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2L,
                  quote = FALSE)

  }

}
