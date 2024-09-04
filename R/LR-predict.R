#' Prediction for the Lorenz regression
#'
#' Provides predictions for an object of class \code{"LR"}.
#'
#' @aliases predict.LR_boot
#' @param object An object of class \code{"LR"}.
#' @param newdata An optional data frame in which to look for variables with which to predict. If omitted, the original data are used.
#' @param type A character string indicating the type of prediction. Possible values are \code{"response"} and \code{"index"} (the default).
#' In the first case, the prediction estimates the conditional expectation of the response given the covariates.
#' In the second case, the prediction estimates only the index of the single-index model.
#' @param ... Additional arguments passed to the function \code{\link{Rearrangement.estimation}}.
#'
#' @return a vector gathering the predictions
#'
#' @details If \code{type="Response"}, the link function of the single-index model must be estimated. This is done via the function \code{\link{Rearrangement.estimation}}.
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
