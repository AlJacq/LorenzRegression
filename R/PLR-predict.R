#' Prediction for the Penalized Lorenz regression
#'
#' \code{predict.PLR} provides predictions for an object of class \code{"PLR"}.
#'
#' @param object An object of S3 class \code{"PLR"}.
#' @param newdata An optional data frame in which to look for variables with which to predict. If omitted, the original data are used.
#' @param type A character string indicating the type of prediction. Possible values are \code{"response"} and \code{"index"} (the default).
#' In the first case, the prediction estimates the conditional expectation of the response given the covariates.
#' In the second case, the prediction estimates only the index of the single-index model.
#' @param pars.idx A vector of size 2 specifying the index of the grid parameter (first element) and the index of the penalty parameter (second element) that should be selected.
#' Default is \code{NULL}, in which case the parameters are selected by the available methods : BIC (always), bootstrap (if \code{object} inherits from the \code{PLR_boot} class) and cross-validation (if \code{object} inherits from the \code{PLR_cv} class).
#' @param ... Additional arguments passed to the function \code{\link{Rearrangement.estimation}}.
#'
#' @return a vector gathering the predictions.
#' If the object has also class \code{"PLR_boot"} and/or \code{"PLR_cv"}, the output is a matrix, where each column corresponds to a selection method.
#'
#' @details If \code{type="response"}, the link function of the single-index model must be estimated. This is done via the function \code{\link{Rearrangement.estimation}}.
#'
#' @seealso \code{\link{Lorenz.Reg}}, \code{\link{Rearrangement.estimation}}
#'
#' @examples
#' ## For examples see example(Lorenz.Reg), example(Lorenz.boot) and example(PLR.CV)
#'
#' @method predict PLR
#' @export

predict.PLR <- function(object, newdata, type=c("index","response"), pars.idx = NULL, ...){

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
  if(!is.null(pars.idx)){
    l <- ncol(object$x)
    pth <- object$path[[pars.idx[1]]][,pars.idx[2]]
    object$theta <- pth[(length(pth)-l+1):length(pth)]
    object$index <- as.vector(object$theta%*%t(object$x))
    class(object) <- "PLR"
  }
  index <- object$theta%*%t(x)
  if (inherits(object,c("PLR_boot","PLR_cv"))){
    if(type=="index"){
      predictor <- t(index)
    }else{
      predictor <- sapply(1:nrow(index),function(i)Rearrangement.estimation(object$y, object$index[i,], t=index[i,], weights=object$weights, ...)$H)
      colnames(predictor) <- rownames(index)
      rownames(predictor) <- NULL
    }
  }else{
    if(type=="index"){
      predictor <- as.vector(index)
    }else{
      predictor <- Rearrangement.estimation(object$y, object$index, t=as.vector(index), weights=object$weights, ...)$H
      names(predictor) <- NULL
    }
  }

  return(predictor)

}
