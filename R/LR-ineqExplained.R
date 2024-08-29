#' Explained inequality metrics for the Lorenz regression
#'
#' \code{ineqExplained.LR} retrieves the explained Gini coefficient or the Lorenz-\eqn{R^2} from an object of class \code{"LR"}.
#'
#' @param object An object of S3 class \code{"LR"}.
#' @param type Character string specifying the type of inequality metric to retrieve. Options are \code{"Gini.explained"} and \code{"Lorenz-R2"}.
#' @param ... Additional arguments.
#'
#' @return A numeric value representing the requested inequality metric.
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

