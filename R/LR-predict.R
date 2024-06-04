#' Prediction for the Lorenz regression
#'
#' \code{predict.LR} provides predictions for an object of class \code{LR}.
#'
#' @param object An object of class \code{"LR"}.
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
#' data(Data.Incomes)
#' NPLR <- Lorenz.Reg(Income ~ ., data = Data.Incomes, penalty = "none")
#' predict(NPLR)
#'
#' @method predict LR
#' @export

predict.LR <- function(object, newdata, type=c("index","response"), ...){

  tt <- terms(object)
  if (!inherits(object, "LR")) stop("The object must be of class 'LR'")
  type <- match.arg(type)
  noData <- (missing(newdata) || is.null(newdata))
  if(noData){
    x <- object$x
  }else{
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, xlev = object$xlevels)
    x <- model.matrix(Terms, m)[,-1,drop=FALSE]
  }
  index <- as.vector(object$theta%*%t(x))
  if(type=="index"){
    predictor <- index
  }else{
    predictor <- Rearrangement.estimation(object$y, object$index, t=index, weights=object$weights, ...)$H
  }

  return(predictor)

}
