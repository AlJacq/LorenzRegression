#' Explained inequality metrics for the penalized Lorenz regression
#'
#' Retrieves the explained Gini coefficient or the Lorenz-\eqn{R^2} from an object of class \code{"PLR"}.
#'
#' @aliases ineqExplained.PLR_boot ineqExplained.PLR_cv
#' @param object An object of S3 class \code{"PLR"}. The object might also have S3 classes \code{"PLR_boot"} and/or \code{"PLR_cv"} (both inherit from class \code{"PLR"})
#' @param type Character string specifying the type of inequality metric to retrieve. Options are \code{"Gini.explained"} (for the explained Gini coefficient) and \code{"Lorenz-R2"} (for the Lorenz-\eqn{R^2}).
#' @param pars.idx What grid and penalty parameters should be used for parameter selection. Either a character string specifying the selection method, where the possible values are:
#' \itemize{
#'    \item \code{"BIC"} (default) - Always available.
#'    \item \code{"Boot"} - Available if \code{object} inherits from \code{"PLR_boot"}.
#'    \item \code{"CV"} - Available if \code{object} inherits from \code{"PLR_cv"}.
#' }
#' Or a numeric vector of length 2, where the first element is the index of the grid parameter and the second is the index of the penalty parameter.
#' @param ... Additional arguments.
#'
#' @return A numeric value representing the requested inequality metric.
#'
#' @method ineqExplained PLR
#' @export

ineqExplained.PLR <- function(object, type = c("Gini.explained","Lorenz-R2"), pars.idx = "BIC", ...){

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

  ineqExplained_PLR(object, type, pars.idx)

}

#' @method ineqExplained PLR_boot
#' @export

ineqExplained.PLR_boot <- function(object, type = c("Gini.explained","Lorenz-R2"), pars.idx = "BIC", ...){

  type <- match.arg(type)

  if(all(pars.idx == "Boot")) pars.idx <- c(object$grid.idx["Boot"],object$lambda.idx["Boot"])
  NextMethod("ineqExplained")

}

#' @method ineqExplained PLR_cv
#' @export

ineqExplained.PLR_cv <- function(object, type = c("Gini.explained","Lorenz-R2"), pars.idx = "BIC", ...){

  type <- match.arg(type)

  if(all(pars.idx == "CV")) pars.idx <- c(object$grid.idx["CV"],object$lambda.idx["CV"])
  NextMethod("ineqExplained")

}

ineqExplained_PLR <- function(object, type, pars.idx){

  pth <- object$path[[pars.idx[1]]][,pars.idx[2]]

  if(type == "Gini.explained"){
    ie <- pth["Explained Gini"]
  }else{
    ie <- pth["Lorenz-R2"]
  }
  names(ie) <- NULL
  ie

}


