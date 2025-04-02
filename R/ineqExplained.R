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

#' Explained inequality metrics for (generalized) linear models
#'
#' Retrieves the explained Gini coefficient or the Lorenz-\eqn{R^2} from an object of class \code{"lm"}.
#'
#' @param object An object of S3 class \code{"lm"}.
#' @param type Character string specifying the type of inequality metric to retrieve. Options are \code{"Gini.explained"} (for the explained Gini coefficient) and \code{"Lorenz-R2"} (for the Lorenz-\eqn{R^2}).
#' @param ... Additional arguments passed to \code{\link{Gini.coef}}.
#'
#' @return A numeric value representing the requested inequality metric.
#'
#' @method ineqExplained lm
#' @export

ineqExplained.lm <- function(object, type = c("Gini.explained","Lorenz-R2"), ...){

  type <- match.arg(type)

  if(!is.null(object$model)){

    y <- model.response(object$model)

  }else if (!is.null(object$y)){

    y <- object$y

  }else{

    stop("The response must be stored in object.")

  }

  yhat <- fitted(object)
  w <- object$weights

  Gi.expl <- Gini.coef(y, x = yhat, weights = w, ...)

  if(type == "Gini.explained"){

    return(Gi.expl)

  }else if(type == "Lorenz-R2"){

    Gi.y <- Gini.coef(y, weights = w, ...)
    LR2 <- Gi.expl/Gi.y
    return(LR2)

  }

}
