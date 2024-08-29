#' Retrieve a measure of explained inequality from a model
#'
#' This generic function extracts a measure of explained inequality, such as the explained Gini coefficient or the Lorenz-R2, from a fitted model object.
#'
#' @param object An object for which the inequality metrics should be extracted.
#' @param type Character string specifying the type of inequality metric to retrieve. Options are \code{"Gini.explained"} for the explained Gini coefficient or \code{"Lorenz-R2"} for the Lorenz-\eqn{R^2}.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return The requested inequality metric.
#'
#' @seealso \code{\link{Lorenz.Reg}}
#'
#' @examples
#' ## For examples see example(Lorenz.Reg)
#'
#' @export

ineqExplained <- function(object, type = c("Gini.explained", "Lorenz-R2"), ...) {
  UseMethod("ineqExplained")
}

#' @describeIn ineqExplained Extract inequality metrics from an \code{"LR"} object
#'
#' \code{ineqExplained.LR} retrieves the explained Gini coefficient or the Lorenz-\eqn{R^2} from an object of class \code{"LR"}.
#'
#' @param object An object of S3 class \code{"LR"}.
#' @param type Character string specifying the type of inequality metric to retrieve. Options are \code{"Gini.explained"} and \code{"Lorenz-R2"}.
#' @param ... Additional arguments.
#'
#' @return A numeric value representing the requested inequality metric.
#'
#'
#' @method ineqExplained LR
#' @export

ineqExplained.LR <- function(object, type = c("Gini.explained","Lorenz-R2"), ...){

  type <- match.arg(type)

  if(type == "Gini.explained"){
    object$Gi.expl
  }else{
    object$LR2
  }

}

#' @method ineqExplained LR_boot
#' @export

ineqExplained.LR_boot <- function(object, ...){
  NextMethod("ineqExplained")
}

