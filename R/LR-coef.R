#' Estimated coefficients for the Lorenz Regression
#'
#' \code{coef.LR} provides the estimated coefficients for an object of class \code{LR}.
#'
#' @param object Output of a call to \code{\link{Lorenz.Reg}}, where \code{penalty="none"}.
#' @param ... Additional arguments.
#'
#' @return a vector gathering the estimated coefficients
#'
#' @seealso \code{\link{Lorenz.Reg}}
#'
#' @examples
#' data(Data.Incomes)
#' NPLR <- Lorenz.Reg(Income ~ ., data = Data.Incomes, penalty = "none")
#' coef(NPLR)
#'
#' @method coef LR
#' @export

coef.LR <- function(object, ...){

  object$theta

}
