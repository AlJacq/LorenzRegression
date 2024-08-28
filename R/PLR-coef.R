#' Estimated coefficients for the Penalized Lorenz Regression
#'
#' \code{coef.PLR} provides the estimated coefficients for an object of class \code{"PLR"}.
#'
#' @aliases coef.PLR_boot coef.PLR_cv
#' @param object An object of S3 class \code{"PLR"}.
#' @param renormalize A logical value determining whether the coefficient vector should be re-normalized to match the representation where the first category of each categorical variable is omitted. Default value is TRUE
#' @param pars.idx A vector of size 2 specifying the index of the grid parameter (first element) and the index of the penalty parameter (second element) that should be selected.
#' Default is \code{NULL}, in which case the parameters are selected by the available methods : BIC (always), bootstrap (if \code{object} inherits from the \code{PLR_boot} class) and cross-validation (if \code{object} inherits from the \code{PLR_cv} class).
#' @param ... Additional arguments
#'
#' @return a vector gathering the estimated coefficients.
#' If the object has also class \code{"PLR_boot"} and/or \code{"PLR_cv"}, the output is a matrix, where each column corresponds to a selection method.
#' If the argument \code{pars.idx} is specified, the output is a vector and the estimated coefficients correspond to the specified grid and penalty parameters.
#'
#' @seealso \code{\link{Lorenz.Reg}}
#'
#' @examples
#' ## For examples see example(Lorenz.Reg), example(Lorenz.boot) and example(PLR.CV)
#'
#' @method coef PLR
#' @export

coef.PLR <- function(object, renormalize=TRUE, pars.idx=NULL, ...){

  if(is.null(pars.idx)) pars.idx <- c(object$grid.idx["BIC"],object$lambda.idx["BIC"])

  coef_PLR(object, renormalize, pars.idx)

  # if(is.matrix(m1)){
  #
  #   m2 <- t(m1)
  #   l <- split(m2,rep(1:ncol(m2), each = nrow(m2)))
  #   names(l) <- colnames(m2)
  #   for (j in 1:length(l)) names(l[[j]]) <- rownames(m2)
  #
  #   return(l)
  #
  # }else{
  #
  #   return(m1)
  #
  # }

}

#' @method coef PLR_boot
#' @rdname coef.PLR
#' @export

coef.PLR_boot <- function(object, renormalize=TRUE, pars.idx=NULL, ...){

  m <- NextMethod("coef")

  if(is.null(pars.idx)){
    pars.idx <- c(object$grid.idx["Boot"],object$lambda.idx["Boot"])
    m1 <- coef_PLR(object, renormalize, pars.idx)
    if(is.list(m)){
      m$Boot <- m1
    }else{
      m <- list("BIC"=m,"Boot"=m1)
    }
  }

  return(m)

}

#' @method coef PLR_cv
#' @rdname coef.PLR
#' @export

coef.PLR_cv <- function(object, renormalize=TRUE, pars.idx=NULL, ...){

  m <- NextMethod("coef")

  if(is.null(pars.idx)){
    pars.idx <- c(object$grid.idx["CV"],object$lambda.idx["CV"])
    m1 <- coef_PLR(object, renormalize, pars.idx)
    if(is.list(m)){
      m$CV <- m1
    }else{
      m <- list("BIC"=m,"CV"=m1)
    }
  }

  return(m)

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
