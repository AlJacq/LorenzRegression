#' Confidence intervals for the Penalized Lorenz Regression
#'
#' \code{confint.PLR_boot} provides bootstrap confidence intervals for the explained Gini coefficient and Lorenz-R2 for an object of class \code{"PLR_boot"}.
#'
#' @param object An object of class \code{"PLR_boot"}.
#' @param parm A character string determining whether the confidence interval is computed for the explained Gini coefficient or for the Lorenz-R2. Possible values are \code{"Gini"} (default, for the explained Gini) and \code{"LR2"} (for the Lorenz-R2)
#' @param level A numeric giving the level of the confidence interval. Default value is 0.95.
#' @param type A character string specifying the bootstrap method. Possible values are \code{"norm"}, \code{"basic"} and \code{"perc"}. For more information, see the argument \code{type} of the function \code{\link{boot.ci}} from the \emph{boot} library.
#' @param pars.idx A vector of size 2 specifying the index of the grid parameter (first element) and the index of the penalty parameter (second element) that should be selected.
#' Default is \code{NULL}, in which case the parameters are selected by the available methods : BIC (always), bootstrap (if \code{object} inherits from the \code{PLR_boot} class) and cross-validation (if \code{object} inherits from the \code{PLR_cv} class).
#' @param bias.corr A logical determining whether bias correction should be performed. Only used if \code{type="norm"}. Default is \code{TRUE}.
#' @param ... Additional arguments.
#'
#' @return The desired confidence interval.
#' If \code{pars.idx=NULL}, the output is a matrix where each row corresponds to a selection method.
#' Otherwise, the grid and penalty parameters are specified by the vector \code{pars.idx} and the output is a vector.
#'
#' @seealso \code{\link{Lorenz.boot}}, \code{\link[boot]{boot.ci}}
#'
#' @examples
#' ## For examples see example(Lorenz.boot)
#'
#' @method confint PLR_boot
#' @export

confint.PLR_boot <- function(object, parm=c("Gini","LR2"), level=0.95, type=c("norm","basic","perc"), pars.idx = NULL, bias.corr = TRUE, ...){

  if (!inherits(object, "PLR_boot")) stop("The object must be of class 'PLR_boot'")

  parm <- match.arg(parm)
  type <- match.arg(type)
  type2 <- switch(type, "basic" = "basic", "norm" = "normal", "perc" = "percent")

  ci.i <- function(pars.idx,parm){
    i <- idx[[pars.idx[1]]][pars.idx[2]]
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

  if(!is.null(pars.idx)){
    ci <- ci.i(pars.idx,parm)
  }else{
    pars.idx.m <- cbind(object$grid.idx,object$lambda.idx)
    ci <- t(apply(pars.idx.m,1,ci.i,parm=parm))
  }

  return(ci)

}
