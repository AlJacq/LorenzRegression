#' Prediction and fitted values for the penalized Lorenz regression
#'
#' \code{prediction} provides predictions for an object of class \code{"PLR"},
#' while \code{fitted} extracts the fitted values.
#'
#' @aliases predict.PLR_boot predict.PLR_cv fitted.PLR fitted.PLR_boot fitted.PLR_cv
#' @param object An object of S3 class \code{"PLR"}. The object might also have S3 classes \code{"PLR_boot"} and/or \code{"PLR_cv"} (both inherit from class \code{"PLR"})
#' @param newdata An optional data frame in which to look for variables with which to predict. If omitted, the original data are used.
#' @param type A character string indicating the type of prediction or fitted values. Possible values are \code{"response"} and \code{"index"} (the default).
#' In the first case, the conditional expectation of the response given the covariates is estimated.
#' In the second case, only the index of the single-index model is estimated.
#' @param pars.idx What grid and penalty parameters should be used for parameter selection. Either a character string specifying the selection method, where the possible values are:
#' \itemize{
#'    \item \code{"BIC"} (default) - Always available.
#'    \item \code{"Boot"} - Available if \code{object} inherits from \code{"PLR_boot"}.
#'    \item \code{"CV"} - Available if \code{object} inherits from \code{"PLR_cv"}.
#' }
#' Or a numeric vector of length 2, where the first element is the index of the grid parameter and the second is the index of the penalty parameter.
#' @param ... Additional arguments passed to the function \code{\link{Rearrangement.estimation}}.
#'
#' @return A vector of predictions for \code{predict}, or a vector of fitted values for \code{fitted}.
#'
#' @details If \code{type="response"}, the link function of the single-index model must be estimated. This is done via the function \code{\link{Rearrangement.estimation}}.
#'
#' @seealso \code{\link{Lorenz.Reg}}, \code{\link{Rearrangement.estimation}}
#'
#' @examples
#' ## For examples see example(Lorenz.Reg), example(Lorenz.boot) and example(PLR.CV)
#'
#' @importFrom stats terms delete.response model.frame
#'
#' @method predict PLR
#' @export

predict.PLR <- function(object, newdata, type=c("index","response"), pars.idx = "BIC", ...){

  type <- match.arg(type)

  if((is.numeric(pars.idx) & length(pars.idx)==2)){
    lth1 <- length(object$path)
    if(pars.idx[1] > lth1) stop("Index of grid parameter is out of bond.")
    lth2 <- ncol(object$path[[pars.idx[1]]])
    if(pars.idx[2] > lth2) stop("Index of lambda parameter is out of bond.")
  }else if(pars.idx == "BIC"){
    pars.idx <- c(object$grid.idx["BIC"],object$lambda.idx["BIC"])
  }else if(pars.idx == "Boot"){
    stop("object is not of class 'PLR_boot'. Therefore pars.idx cannot be set to 'Boot'.")
  }else if(pars.idx == "CV"){
    stop("object is not of class 'PLR_cv'. Therefore pars.idx cannot be set to 'CV'.")
  }else{
    stop("pars.idx does not have the correct format")
  }

  predict_PLR(object, newdata, type, pars.idx, ...)

}

#' @method predict PLR_boot
#' @export

predict.PLR_boot <- function(object, newdata, type=c("index","response"), pars.idx = "BIC", ...){

  type <- match.arg(type)

  if(all(pars.idx == "Boot")) pars.idx <- c(object$grid.idx["Boot"],object$lambda.idx["Boot"])
  NextMethod("predict")

}

#' @method predict PLR_cv
#' @export

predict.PLR_cv <- function(object, newdata, type=c("index","response"), pars.idx = "BIC", ...){

  type <- match.arg(type)

  if(all(pars.idx == "CV")) pars.idx <- c(object$grid.idx["CV"],object$lambda.idx["CV"])
  NextMethod("predict")

}

predict_PLR <- function(object, newdata, type, pars.idx, ...){

  # Data (re)-construction
  tt <- terms(object)
  noData <- (missing(newdata) || is.null(newdata))
  if(noData){
    x <- object$x
  }else{
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, xlev = object$xlevels)
    x <- model_matrix_PLR(Terms,m)
  }

  # Retrieving theta and index
  l <- ncol(object$x)
  pth <- object$path[[pars.idx[1]]][,pars.idx[2]]
  object$theta <- pth[(length(pth)-l+1):length(pth)]
  object$index <- as.vector(object$theta%*%t(object$x)) # Necessarily on the original x
  index <- object$theta%*%t(x) # on the x used for predictions

  # Defining the predictor
  if(type=="index"){
    predictor <- as.vector(index)
  }else{
    predictor <- Rearrangement.estimation(object$y, object$index, t=as.vector(index), weights=object$weights, ...)$H
    names(predictor) <- NULL
  }

  return(predictor)

}

#' @method fitted PLR
#' @rdname predict.PLR
#' @export

fitted.PLR <- function(object, type=c("index","response"), pars.idx = "BIC", ...){

  predict.PLR(object, type = type, pars.idx = pars.idx, ...)

}

#' @method fitted PLR_boot
#' @export

fitted.PLR_boot <- function(object, type=c("index","response"), pars.idx = "BIC", ...){

  if(all(pars.idx == "Boot")) pars.idx <- c(object$grid.idx["Boot"],object$lambda.idx["Boot"])
  NextMethod("fitted")

}

#' @method fitted PLR_cv
#' @export

fitted.PLR_cv <- function(object, type=c("index","response"), pars.idx = "BIC", ...){

  if(all(pars.idx == "CV")) pars.idx <- c(object$grid.idx["CV"],object$lambda.idx["CV"])
  NextMethod("fitted")

}
