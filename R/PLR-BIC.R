#' Determines the regularization parameter (lambda) in a PLR via optimization of an information criterion.
#'
#' \code{PLR.BIC} takes as input a matrix of estimated parameter vectors, where each row corresponds to a covariate and each column corresponds to a value of lambda,
#' and returns the index of the optimal column by optimizing an information criterion. By default the BIC is used.
#'
#' @param YX_mat A matrix with the first column corresponding to the response vector, the remaining ones being the explanatory variables.
#' @param theta matrix gathering the path of estimated parameter vectors. Each row corresponds to a given covariate. Each column corresponds to a given value of lambda
#' @param weights vector of sample weights. By default, each observation is given the same weight.
#' @param IC indicates which information criterion is used. Possibles values are "BIC" (default) or "AIC".
#'
#' @return A list with two components
#' \describe{
#'    \item{\code{val}}{vector indicating the value attained by the information criterion for each value of lambda.}
#'    \item{\code{best}}{index of the value of lambda where the optimum is attained.}
#' }
#'
#' @seealso \code{\link{Lorenz.Reg}}, \code{\link{PLR.wrap}}, \code{\link{Lorenz.FABS}}, \code{\link{Lorenz.SCADFABS}}
#'
#' @section References:
#' Jacquemain, A., C. Heuchenne, and E. Pircalabelu (2022). A penalised bootstrap estimation procedure for the explained Gini coefficient.
#'
#' @examples
#' data(Data.Incomes)
#' YX_mat <- Data.Incomes
#' PLR <- PLR.wrap(YX_mat, h = nrow(YX_mat)^(-1/5.5), eps = 0.005)
#' PLR.BIC(YX_mat, PLR$theta)
#'
#' @export

PLR.BIC <- function(YX_mat, theta, weights=NULL, IC=c("BIC","AIC")){

  IC <- match.arg(IC)

  y <- YX_mat[,1]
  n <- length(y)
  X <- YX_mat[,-1]
  Index <- as.matrix(X)%*%theta

  score <- log(apply(Index, 2, function(t) Gini.coef(y, x=t, na.rm=T, ties.method="mean", weights=weights)))

  if (IC=="BIC") pen <- apply(theta,2,function(x)sum(x!=0))*log(n)/(2*n)
  if (IC=="AIC") pen <- apply(theta,2,function(x)sum(x!=0))/n

  IC.val <- score - pen
  IC.best <- which.max(IC.val)

  return.list <- list()
  return.list$val <- IC.val
  return.list$best <- IC.best

  return(return.list)

}
