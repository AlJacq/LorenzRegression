#' Confidence intervals for the Penalized Lorenz Regression
#'
#' \code{confint.PLR_boot} provides bootstrap confidence intervals for the explained Gini coefficient and Lorenz-R2 for an object of class \code{"PLR_boot"}.
#'
#' @param object An object of class \code{"PLR_boot"}. The object might also have S3 class \code{"PLR_cv"}.
#' @param parm A character string determining whether the confidence interval is computed for the explained Gini coefficient or for the Lorenz-R2. Possible values are \code{"Gini"} (default, for the explained Gini) and \code{"LR2"} (for the Lorenz-R2)
#' @param level A numeric giving the level of the confidence interval. Default value is 0.95.
#' @param type A character string specifying the bootstrap method. Possible values are \code{"norm"}, \code{"basic"} and \code{"perc"}. For more information, see the argument \code{type} of the function \code{\link{boot.ci}} from the \emph{boot} library.
#' @param pars.idx What grid and penalty parameters should be used for parameter selection. Either a character string specifying the selection method, where the possible values are:
#' \itemize{
#'    \item \code{"BIC"} (default).
#'    \item \code{"Boot"}.
#'    \item \code{"CV"} - Available if \code{object} inherits from \code{"PLR_cv"}.
#' }
#' Or a numeric vector of length 2, where the first element is the index of the grid parameter and the second is the index of the penalty parameter.
#' @param bias.corr A logical determining whether bias correction should be performed. Only used if \code{type="norm"}. Default is \code{TRUE}.
#' @param ... Additional arguments.
#'
#' @return A vector providing the desired confidence interval.
#'
#' @seealso \code{\link{Lorenz.boot}}, \code{\link[boot]{boot.ci}}
#'
#' @examples
#' ## For examples see example(Lorenz.boot)
#'
#' @method confint PLR_boot
#' @export

confint.PLR_boot <- function(object, parm=c("Gini","LR2"), level=0.95, type=c("norm","basic","perc"), pars.idx = "BIC", bias.corr = TRUE, ...){

  parm <- match.arg(parm)
  type <- match.arg(type)

  if((is.numeric(pars.idx) & length(pars.idx)==2)){
    lth1 <- length(object$path)
    lth2 <- ncol(object$path[[lth1]])
    if(pars.idx[1] > lth1 | pars.idx[2] > lth2) stop("Indices in pars.idx are out of bounds.")
  }else if(pars.idx == "BIC"){
    pars.idx <- c(object$grid.idx["BIC"],object$lambda.idx["BIC"])
  }else if(pars.idx == "Boot"){
    pars.idx <- c(object$grid.idx["Boot"],object$lambda.idx["Boot"])
  }else if(pars.idx == "CV"){
    stop("object is not of class 'PLR_cv'. Therefore pars.idx cannot be set to 'CV'.")
  }else{
    stop("pars.idx does not have the correct format")
  }

  ci_boot_PLR(object, pars.idx, parm, level, type, bias.corr)

}

#' @method confint PLR_cv
#' @rdname confint.PLR_boot
#' @export

confint.PLR_cv <- function(object, parm=c("Gini","LR2"), level=0.95, type=c("norm","basic","perc"), pars.idx = "BIC", bias.corr = TRUE, ...){

  if(all(pars.idx == "CV")) pars.idx <- c(object$grid.idx["CV"],object$lambda.idx["CV"])
  NextMethod("confint")

}

#' @method confint PLR
#' @export

confint.PLR <- function(object, parm=c("Gini","LR2"), level=0.95, pars.idx = "BIC", ...){

  stop("The 'confint' method requires the object to inherit from 'PLR_boot'. The current implementation of the penalized Lorenz regression uses bootstrap for confidence interval construction. Please ensure the object is generated using bootstrap (i.e., it should be of class 'PLR_boot').")

}

ci_boot_PLR <- function(object, pars.idx, parm, level, type, bias.corr){
  type2 <- switch(type, "basic" = "basic", "norm" = "normal", "perc" = "percent")
  path.sizes <- sapply(object$path,ncol)
  path.size <- sum(path.sizes)
  lth.path <- length(path.sizes)
  idx <- lapply(1:lth.path,function(i)(cumsum(path.sizes)-path.sizes+1)[i]:cumsum(path.sizes)[i])
  i <- idx[[pars.idx[1]]][pars.idx[2]]
  if(parm == "LR2") i <- i + path.size
  ci <- boot.ci(object$boot_out, conf = level, type = type, index = i)
  ci <- ci[[type2]]
  ci <- ci[length(ci)-c(1,0)]
  if(!bias.corr & type=="norm") ci <- ci - mean(ci) + object$boot_out$t0[i]
  names(ci) <- paste0((c(0,level)+(1-level)/2)*100," %")
  return(ci)
}
