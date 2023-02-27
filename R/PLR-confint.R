#' Confidence intervals for the Penalized Lorenz Regression
#'
#' \code{confint.PLR} provides confidence intervals for the explained Gini coefficient and Lorenz-R2 for an parm of class \code{PLR}.
#'
#' @param object Output of a call to \code{\link{Lorenz.Reg}}, where \code{penalty!="none"} and \code{Boot.inference=TRUE}.
#' @param parm Determines whether the confidence interval is computed for the explained Gini coefficient or for the Lorenz-R2. Possible values are "Gini" (default, for the explained Gini) and "LR2" (for the Lorenz-R2).
#' @param level level of the confidence interval
#' @param boot.method What bootstrap method is used to construct the confidence interval. Default value is "Param", which exploits the asymptotic normality and only bootstraps the variance.
#' Other possible values are "Perc" (percentile bootstrap) and "Basic" (basic bootstrap). Percentile bootstrap directly plugs the quantiles of the bootstrap distribution.
#' Basic bootstrap is based on bootstrapping the whole distribution of the estimator.
#' @param which.pars Which values of the bandwidth h and the penalty parameter lambda should be used. Default is NULL, in which case the optimal values are used.
#' @param ... Additional arguments.
#'
#' @return A matrix gathering the desired confidence intervals. Each row corresponds to a different selection method for the pair (h,lambda).
#'
#' @details Use this function only if Boot.inference was set to TRUE in the call to \code{\link{Lorenz.Reg}}. Otherwise, bootstrap was not computed and the confidence intervals cannot be determined.
#'
#' @seealso \code{\link{Lorenz.Reg}}
#'
#' @examples
#' data(Data.Incomes)
#' set.seed(123)
#' Data <- Data.Incomes[sample(1:nrow(Data.Incomes),50),]
#' PLR <- Lorenz.Reg(Income ~ ., data = Data, h.grid = nrow(Data)^(-1/5.5),
#'                   penalty = "SCAD", eps = 0.02, seed.boot = 123, B = 40, Boot.inference = TRUE)
#' confint(PLR)
#'
#' @method confint PLR
#' @export

confint.PLR <- function(object, parm=c("Gini","LR2"), level = 0.95, boot.method=c("Param","Basic","Perc"), which.pars = NULL, ...){

  PLR <- object
  parm <- match.arg(parm)
  boot.method <- match.arg(boot.method)
  alpha <- 1-level

  if(length(grep(".star",names(PLR))) > 0){

    if(is.null(which.pars)){
      which.h <- PLR$which.h
      which.lambda <- PLR$which.lambda
    }else{
      if(length(which.pars)!=2) stop("which.pars must either be NULL or a vector of size 2 where the first element is the index of the bandwidth and the second the index of lambda in the path.")
      which.h <- which.pars[1]
      which.lambda <- which.pars[2]
    }

    CI <- matrix(nrow = length(which.h), ncol = 2)
    colnames(CI) <- c("Lower bound", "Upper bound")
    if(is.null(which.pars)) rownames(CI) <- names(PLR$which.h)

    for (k in 1:length(which.h)){

      which.h.k <- which.h[k]
      which.lambda.k <- which.lambda[k]

      if(parm == "Gini") CI.k <- boot.confint(PLR$path[[which.h.k]]["Explained Gini",which.lambda.k],
                                              PLR$Gi.star[[which.h.k]][[which.lambda.k]],
                                              alpha, boot.method)
      if(parm == "LR2") CI.k <- boot.confint(PLR$path[[which.h.k]]["Lorenz-R2",which.lambda.k],
                                             PLR$LR2.star[[which.h.k]][[which.lambda.k]],
                                             alpha, boot.method)
      CI[k,] <- CI.k

    }

    return(CI)

  }else{
    stop("The input must contain a bootstrap estimation. Consider turning the Boot.inference argument in the Lorenz.Reg function to TRUE")
  }

}
