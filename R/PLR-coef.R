#' Estimated coefficients for the Penalized Lorenz Regression
#'
#' \code{coef.PLR} provides the estimated coefficients for an object of class \code{"PLR"}.
#'
#' @aliases coef.PLR_boot coef.PLR_cv
#' @param object An object of S3 class \code{"PLR"}. The object might also have S3 classes \code{"PLR_boot"} and/or \code{"PLR_cv"} (both inherit from class \code{"PLR"})
#' @param renormalize A logical value determining whether the coefficient vector should be re-normalized to match the representation where the first category of each categorical variable is omitted. Default value is TRUE
#' @param pars.idx What grid and penalty parameters should be used for parameter selection. Either a character string specifying the selection method, where the possible values are:
#' \itemize{
#'    \item \code{"BIC"} (default) - Always available.
#'    \item \code{"Boot"} - Available if \code{object} inherits from \code{"PLR_boot"}.
#'    \item \code{"CV"} - Available if \code{object} inherits from \code{"PLR_cv"}.
#' }
#' Or a numeric vector of length 2, where the first element is the index of the grid parameter and the second is the index of the penalty parameter.
#' @param ... Additional arguments
#'
#' @return a vector gathering the estimated coefficients.
#'
#' @seealso \code{\link{Lorenz.Reg}}
#'
#' @examples
#' ## For examples see example(Lorenz.Reg), example(Lorenz.boot) and example(PLR.CV)
#'
#' @method coef PLR
#' @export

coef.PLR <- function(object, renormalize=TRUE, pars.idx="BIC", ...){

  if((is.numeric(pars.idx) & length(pars.idx)==2)){
    lth1 <- length(object$path)
    lth2 <- ncol(object$path[[lth1]])
    if(pars.idx[1] > lth1 | pars.idx[2] > lth2) stop("Indices in pars.idx are out of bounds.")
  }else if(pars.idx == "BIC"){
    pars.idx <- c(object$grid.idx["BIC"],object$lambda.idx["BIC"])
  }else if(pars.idx == "Boot"){
    stop("object is not of class 'PLR_boot'. Therefore pars.idx cannot be set to 'Boot'.")
  }else if(pars.idx == "CV"){
    stop("object is not of class 'PLR_cv'. Therefore pars.idx cannot be set to 'CV'.")
  }else{
    stop("pars.idx does not have the correct format")
  }

  coef_PLR(object, renormalize, pars.idx)

}

#' @method coef PLR_boot
#' @rdname coef.PLR
#' @export

coef.PLR_boot <- function(object, renormalize=TRUE, pars.idx="BIC", ...){

  if(pars.idx == "Boot"){
    pars.idx <- c(object$grid.idx["Boot"],object$lambda.idx["Boot"])
    coef_PLR(object, renormalize, pars.idx)
  }else{
    NextMethod("coef")
  }

}

#' @method coef PLR_cv
#' @rdname coef.PLR
#' @export

coef.PLR_cv <- function(object, renormalize=TRUE, pars.idx="BIC", ...){

  if(pars.idx == "CV"){
    pars.idx <- c(object$grid.idx["CV"],object$lambda.idx["CV"])
    coef_PLR(object, renormalize, pars.idx)
  }else{
    NextMethod("coef")
  }

}


coef_PLR <- function(object, renormalize, pars.idx){

  l <- ncol(object$x)
  pth <- object$path[[pars.idx[1]]][,pars.idx[2]]
  object$theta <- pth[(length(pth)-l+1):length(pth)]

  if(renormalize){
    m1 <- PLR.normalize(object)
  }else{
    m1 <- object$theta
  }

  m1

}
