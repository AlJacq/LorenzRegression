#' Estimated coefficients for the Lorenz Regression
#'
#' \code{coef.LR} provides the estimated coefficients for an object of class \code{"LR"}.
#'
#' @param object An object of S3 class \code{"LR"}.
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
