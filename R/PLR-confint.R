#' Confidence intervals for the Penalized Lorenz Regression
#'
#' \code{confint.PLR_boot} provides bootstrap confidence intervals for the explained Gini coefficient and Lorenz-R2 for an object of class \code{"PLR_boot"}.
#'
#' @param object An object of class \code{"PLR_boot"}.
#' @param parm A character string determining whether the confidence interval is computed for the explained Gini coefficient or for the Lorenz-R2. Possible values are \code{"Gini"} (default, for the explained Gini) and \code{"LR2"} (for the Lorenz-R2)
#' @param level A numeric giving the level of the confidence interval. Default value is 0.95.
#' @param type A character string specifying the bootstrap method. Possible values are \code{"norm"}, \code{"basic"} and \code{"perc"}. For more information, see the argument \code{type} of the function \code{\link{boot.ci}} from the \emph{boot} library.
#' @param which.pars A vector of size 2 specifying the index of the tuning parameter (first element) and the index of the penalty parameter (second element) that should be selected.
#' Default is \code{NULL}, in which case the parameters are selected by the available methods : BIC (always), bootstrap (if \code{object} inherits from the \code{PLR_boot} class) and cross-validation (if \code{object} inherits from the \code{PLR_cv} class).
#' @param bias.corr A logical determining whether bias correction should be performed. Only used if \code{type="norm"}. Default is \code{TRUE}.
#'
#' @return The desired confidence interval.
#' If \code{which.pars=NULL}, the output is a matrix where each row corresponds to a selection method.
#' Otherwise, the tuning and penalty parameters are specified by the vector \code{which.pars} and the output is a vector.
#'
#' @seealso \code{\link{Lorenz.boot}}, \code{\link[boot]{boot.ci}}
#'
#' @examples
#' ## For examples see example(Lorenz.boot)
#'
#' @method confint PLR_boot
#' @export

confint.PLR_boot <- function(object, parm=c("Gini","LR2"), level=0.95, type=c("norm","basic","perc"), which.pars = NULL, bias.corr = TRUE, ...){

  if (!inherits(object, "PLR_boot")) stop("The object must be of class 'PLR_boot'")

  parm <- match.arg(parm)
  type <- match.arg(type)
  type2 <- switch(type, "basic" = "basic", "norm" = "normal", "perc" = "percent")

  ci.i <- function(which.pars,parm){
    i <- idx[[which.pars[1]]][which.pars[2]]
    if(parm == "LR2") i <- i + path.size
    ci <- boot.ci(object$boot_out, conf = level, type = type, index = i)
    ci <- ci[[type2]]
    ci <- ci[length(ci)-c(1,0)]
    if(!bias.corr & type=="norm") ci <- ci - mean(ci) + object$boot_out$t0[i]
    names(ci) <- paste0((c(0,level)+(1-level)/2)*100," %")
    return(ci)
  }

  path.sizes <- sapply(object$path,ncol)
  path.size <- sum(path.sizes)
  lth.path <- length(path.sizes)
  idx <- lapply(1:lth.path,function(i)(cumsum(path.sizes)-path.sizes+1)[i]:cumsum(path.sizes)[i])

  if(!is.null(which.pars)){
    ci <- ci.i(which.pars,parm)
  }else{
    which.pars.m <- cbind(object$which.tuning,object$which.lambda)
    ci <- t(apply(which.pars.m,1,ci.i,parm=parm))
  }

  return(ci)

}
