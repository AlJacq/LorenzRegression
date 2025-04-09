#' Prediction and fitted values for the Lorenz regression
#'
#' \code{predict} provides predictions for an object of class \code{"LR"},
#' while \code{fitted} extracts the fitted values.
#'
#' @aliases predict.LR_boot fitted.LR fitted.LR_boot
#' @param object An object of class \code{"LR"}.
#' @param newdata An optional data frame in which to look for variables with which to predict. If omitted, the original data are used.
#' @param type A character string indicating the type of prediction or fitted values. Possible values are \code{"response"} and \code{"index"} (the default).
#' In the first case, the prediction estimates the conditional expectation of the response given the covariates.
#' In the second case, the prediction estimates only the index of the single-index model.
#' @param ... Additional arguments passed to the function \code{\link{Rearrangement.estimation}}.
#'
#' @return A vector of predictions for \code{predict}, or a vector of fitted values for \code{fitted}.
#'
#' @details
#' The \code{type} argument distinguishes between two types of prediction outputs, aligned with the goals of the Lorenz regression.
#' When \code{type = "index"}, the function returns the estimated index \eqn{X^\top \theta} of the single-index model. This index captures the full ordering structure of the conditional expectation and is sufficient for computing the explained Gini coefficient, which is the primary focus of the method. Crucially, this estimation does not require recovering the full nonparametric link function.
#' When \code{type = "response"}, the function estimates the full conditional expectation \eqn{\mathbb{E}[Y | X]} by performing a second-stage estimation of the link function via \code{\link{Rearrangement.estimation}}. This is useful if fitted or predicted response values are needed for other purposes.
#'
#' @seealso \code{\link{Lorenz.Reg}}, \code{\link{Rearrangement.estimation}}
#'
#' @examples
#' ## For examples see example(Lorenz.Reg) and example(Lorenz.boot)
#'
#' @importFrom stats terms delete.response model.frame model.matrix
#'
#' @method predict LR
#' @export

predict.LR <- function(object, newdata, type=c("index","response"), ...){

  tt <- terms(object)
  if (is.null(object$theta)) stop("No prediction is available for an empty model.")
  type <- match.arg(type)
  noData <- (missing(newdata) || is.null(newdata))
  if(noData){
    x <- object$x
  }else{
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, xlev = object$xlevels)
    x <- model.matrix(Terms, m)[,-1,drop=FALSE]
  }
  object$index <- as.vector(object$theta%*%t(object$x)) # Necessarily on original data
  index <- as.vector(object$theta%*%t(x)) # Not necessarily on original data
  if(type=="index"){
    predictor <- index
  }else{
    predictor <- Rearrangement.estimation(object$y, object$index, t=index, weights=object$weights, ...)$H
    names(predictor) <- NULL
  }

  return(predictor)

}

#' @method predict LR_boot
#' @export

predict.LR_boot <- function(object, newdata, type=c("index","response"), ...){
  NextMethod("predict")
}

#' @method fitted LR
#' @rdname predict.LR
#' @export

fitted.LR <- function(object, type=c("index","response"), ...){
  type <- match.arg(type)
  predict.LR(object, type = type, ...)
}

#' @method fitted LR_boot
#' @export

fitted.LR_boot <- function(object, type=c("index","response"), ...){
  NextMethod("fitted")
}

#' Residuals for the Lorenz regression
#'
#' \code{residuals} provides residuals for an object of class \code{"LR"}.
#'
#' @aliases residuals.LR_boot
#' @param object An object of class \code{"LR"}.
#' @param ... Additional arguments passed to the function \code{\link{Rearrangement.estimation}}.
#'
#' @return A vector of residuals.
#'
#' @details Computing residuals entail to estimate the link function of the single-index model. This is done via the function \code{\link{Rearrangement.estimation}}.
#'
#' @seealso \code{\link{Lorenz.Reg}}, \code{\link{Rearrangement.estimation}}
#'
#' @examples
#' ## For examples see example(Lorenz.Reg) and example(Lorenz.boot)
#'
#' @method residuals LR
#' @export

residuals.LR <- function(object, ...){

  yhat <- fitted.LR(object, type = "response", ...)
  y <- object$y
  r <- y - yhat
  return(r)

}

#' @method residuals LR_boot
#' @export

residuals.LR_boot <- function(object, ...){
  NextMethod("residuals")
}
