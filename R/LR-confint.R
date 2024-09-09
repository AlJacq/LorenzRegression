#' Confidence intervals for the Lorenz regression
#'
#' Provides bootstrap confidence intervals for the explained Gini coefficient, Lorenz-R2 and theta vector for an object of class \code{"LR_boot"}.
#'
#' @aliases confint.LR
#' @param object An object of class \code{"LR_boot"}. The current implementation requires bootstrap to construct confidence intervals. Hence, it is not sufficient that \code{object} inherits from \code{"LR"}.
#' @param parm A logical value determining whether the confidence interval is computed for the explained Gini coefficient, for the Lorenz-\eqn{R^2} or for the vector of coefficients of the single-index model. Possible values are \code{"Gini"} (default, for the explained Gini),\code{"LR2"} (for the Lorenz-\eqn{R^2}) and \code{"theta"} (for the index coefficients).
#' @param level A numeric giving the level of the confidence interval. Default value is 0.95.
#' @param type A character string specifying the bootstrap method. Possible values are \code{"norm"}, \code{"basic"} and \code{"perc"}. For more information, see the argument \code{type} of the function \code{\link[boot]{boot.ci}} from the \emph{boot} library.
#' @param bias.corr A logical determining whether bias correction should be performed. Only used if \code{type="norm"}. Default is \code{TRUE}.
#' @param ... Additional arguments.
#'
#' @return The desired confidence interval.
#' If \code{parm="Gini"} or \code{parm="LR2"}, the output is a vector.
#' If \code{parm="theta"}, it is a matrix where each row corresponds to a different coefficient.
#'
#' @seealso \code{\link{Lorenz.boot}}, \code{\link[boot]{boot.ci}}
#'
#' @examples
#' ## For examples see example(Lorenz.boot)
#'
#' @importFrom boot boot.ci
#'
#' @method confint LR_boot
#' @export

confint.LR_boot <- function(object, parm=c("Gini","LR2","theta"), level=0.95, type=c("norm","basic","perc"), bias.corr=TRUE, ...){

  parm <- match.arg(parm)
  type <- match.arg(type)

  ci_boot_LR(object, parm, level, type, bias.corr)

}

#' @method confint LR
#' @export

confint.LR <- function(object, parm=c("Gini","LR2","theta"), level=0.95, ...){

  stop("The 'confint' method requires the object to inherit from 'LR_boot'. The current implementation of the Lorenz regression uses bootstrap for confidence interval construction. Please ensure the object is generated using bootstrap (i.e., it should be of class 'LR_boot').")

}

ci_boot_LR <- function(object, parm, level, type, bias.corr){

  type2 <- switch(type, "basic" = "basic", "norm" = "normal", "perc" = "percent")

  ci.i <- function(i){
    ci <- boot.ci(object$boot_out, conf = level, type = type, index = i)
    ci <- ci[[type2]]
    ci <- ci[length(ci)-c(1,0)]
    if(!bias.corr & type=="norm") ci <- ci - mean(ci) + object$boot_out$t0[i]
    names(ci) <- paste0((c(0,level)+(1-level)/2)*100," %")
    return(ci)
  }

  if(parm == "Gini"){
    ci <- ci.i(1)
  }else if(parm == "LR2"){
    ci <- ci.i(2)
  }else{
    ci <- t(sapply(3:ncol(object$boot_out$t),ci.i))
    rownames(ci) <- names(object$theta)
  }

  return(ci)

}

