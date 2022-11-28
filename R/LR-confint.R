#' Confidence intervals for the Lorenz Regression
#'
#' \code{confint.LR} provides confidence intervals for the explained Gini coefficient, Lorenz-R2 and theta vector for an object of class \code{LR}.
#'
#' @param object Output of a call to \code{\link{Lorenz.Reg}}, where \code{penalty="none"} and \code{Boot.inference=TRUE}.
#' @param parm Determines whether the confidence interval is computed for the explained Gini coefficient, for the Lorenz-R2 or for the vector of theta coefficients. Possible values are "Gini" (default, for the explained Gini),"LR2" (for the Lorenz-R2) and "theta" (for the vector theta).
#' @param level level of the confidence interval
#' @param boot.method What bootstrap method is used to construct the confidence interval. Default value is "Param", which exploits the asymptotic normality and only bootstraps the variance.
#' Other possible values are "Perc" (percentile bootstrap) and "Basic" (basic bootstrap). Percentile bootstrap directly plugs the quantiles of the bootstrap distribution.
#' Basic bootstrap is based on bootstrapping the whole distribution of the estimator.
#' @param ... Additional arguments.
#'
#' @return The desired confidence interval. If parm is set to either "Gini" or "LR2", the output is a vector. If parm is set to "theta", it is a matrix where each row corresponds to a different coefficient.
#'
#' @details Use this function only if Boot.inference was set to TRUE in the call to \code{\link{Lorenz.Reg}}. Otherwise, bootstrap was not computed and the confidence intervals cannot be determined.
#'
#' @seealso \code{\link{Lorenz.Reg}}
#'
#' @examples
#' data(Data.Incomes)
#' NPLR <- Lorenz.Reg(Income ~ ., data = Data.Incomes, penalty = "none", seed.boot = 123, B = 40, Boot.inference = TRUE)
#' confint(NPLR)
#'
#' @method confint LR
#' @export

confint.LR <- function(object, parm=c("Gini","LR2","theta"), level = 0.95, boot.method=c("Param","Basic","Perc"), ...){

  LR <- object
  parm <- match.arg(parm)
  boot.method <- match.arg(boot.method)
  alpha <- 1-level

  if(length(grep(".star",names(LR))) > 0){

    if(parm == "Gini") CI <- boot.confint(LR$Gi.expl, LR$Gi.star, alpha, boot.method)
    if(parm == "LR2") CI <- boot.confint(LR$LR2, LR$LR2.star, alpha, boot.method)
    if(parm == "theta"){
      CI <-t(sapply(1:length(LR$theta),function(i)boot.confint(LR$theta[i], LR$theta.star[,i], alpha, boot.method)))
      rownames(CI) <- names(LR$theta)
    }

    return(CI)

  }else{
    stop("The input must contain a bootstrap estimation. Consider turning the Boot.inference argument in the Lorenz.Reg function to TRUE")
  }

}
