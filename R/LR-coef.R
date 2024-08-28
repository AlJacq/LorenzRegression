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
#' ## For examples see example(Lorenz.Reg)
#'
#' @method coef LR
#' @export

coef.LR <- function(object, ...){

  object$theta

}

#' @method coef LR_boot
#' @export

coef.LR_boot <- function(object, ...){
  NextMethod("coef")
}

