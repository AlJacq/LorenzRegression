#' Prediction for the Penalized Lorenz regression
#'
#' \code{predict.PLR} provides predictions for an object of class \code{PLR}.
#'
#' @param object An object of class \code{"PLR"}.
#' @param type A character string indicating the type of prediction. Possible values are \code{"response"} and \code{"index"} (the default).
#' In the first case, the prediction estimates the conditional expectation of the response given the covariates.
#' In the second case, the prediction estimates only the index of the single-index model.
#' @param ... Additional arguments passed to the function \code{\link{Rearrangement.estimation}}.
#'
#' @return a vector gathering the predictions. If the object has also class \code{"PLR_boot"} and/or \code{"PLR_cv"}, the output is a matrix, where each column corresponds to a selection method.
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
#' @method predict PLR
#' @export

predict.PLR <- function(object, newdata, type=c("index","response"), ...){

  tt <- terms(object)
  if (!inherits(object, "PLR")) stop("The object must be of class 'PLR'")
  type <- match.arg(type)
  noData <- (missing(newdata) || is.null(newdata))
  if(noData){
    x <- object$x
  }else{
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, xlev = object$xlevels)
    x <- model_matrix_PLR(Terms,m)
  }
  index <- object$theta%*%t(x)
  if (inherits(object,c("PLR_boot","PLR_cv"))){
    if(type=="index"){
      predictor <- index
    }else{
      predictor <- t(sapply(1:nrow(index),function(i)Rearrangement.estimation(object$y, object$index[i,], t=index[i,], weights=object$weights, ...)$H))
      rownames(predictor) <- rownames(index)
    }
  }else{
    if(type=="index"){
      predictor <- as.vector(index)
    }else{
      predictor <- Rearrangement.estimation(object$y, object$index, t=as.vector(index), weights=object$weights, ...)$H
    }
  }

  return(predictor)

}
