#' Confidence intervals for the Lorenz Regression
#'
#' \code{confint.LR} provides bootstrap confidence intervals for the explained Gini coefficient, Lorenz-R2 and theta vector for an object of class \code{"LR_boot"}.
#'
#' @param object An object of class \code{"LR_boot"}.
#' @param parm A logical value determining whether the confidence interval is computed for the explained Gini coefficient, for the Lorenz-R2 or for the vector of theta coefficients. Possible values are \code{"Gini"} (default, for the explained Gini),\code{"LR2"} (for the Lorenz-R2) and \code{"theta"} (for the vector theta).
#' @param level A numeric giving the level of the confidence interval. Default value is 0.95.
#' @param type A character string specifying the bootstrap method. Possible values are \code{"norm"}, \code{"basic"} and \code{"perc"}. For more information, see the argument \code{type} of the function \code{\link{boot.ci}} from the \code{\link{boot}} library.
#'
#' @param ... Additional arguments.
#'
#' @return The desired confidence interval.
#' If \code{parm="Gini"} or \code{parm="LR2"}, the output is a vector.
#' If \code{parm="theta"}, it is a matrix where each row corresponds to a different coefficient.
#'
#' @seealso \code{\link{Lorenz.boot}}, \code{\link[boot]{boot.ci}}
#'
#' @examples
#' \donttest{
#' # The following piece of code might take several minutes
#' data(Data.Incomes)
#' set.seed(123)
#' Data <- Data.Incomes[sample(1:nrow(Data.Incomes),50),]
#' NPLR <- Lorenz.Reg(Income ~ ., data = Data, penalty = "none",
#'                    seed.boot = 123, B = 40, Boot.inference = TRUE)
#' confint(NPLR)
#' }
#'
#' @importFrom boot boot.ci
#'
#' @method confint LR
#' @export

confint.LR <- function(object, parm=c("Gini","LR2","theta"), level=0.95, type=c("norm","basic","perc"), ...){

  if (!inherits(object, "LR_boot")) stop("The object must be of class 'LR_boot'")

  parm <- match.arg(parm)
  type <- match.arg(type)
  type2 <- switch(type, "basic" = "basic", "norm" = "normal", "perc" = "percent")

  ci.i <- function(i){
    ci <- boot.ci(object$boot_out, conf = level, type = type, index = i)
    ci <- ci[[type2]]
    ci <- ci[length(ci)-c(1,0)]
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
