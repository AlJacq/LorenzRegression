#' Determines the regularization parameter (lambda) in a PLR via optimization of an information criterion.
#'
#' \code{PLR.BIC} takes as input a matrix of estimated parameter vectors, where each row corresponds to a covariate and each column corresponds to a value of lambda,
#' and returns the index of the optimal column by optimizing an information criterion. By default the BIC is used.
#'
#' @param y a vector of responses
#' @param x a matrix of explanatory variables
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
#' @seealso \code{\link{Lorenz.Reg}}, \code{\link{Lorenz.FABS}}, \code{\link{Lorenz.SCADFABS}}
#'
#' @section References:
#' Jacquemain, A., C. Heuchenne, and E. Pircalabelu (2024). A penalised bootstrap estimation procedure for the explained Gini coefficient. \emph{Electronic Journal of Statistics 18(1) 247-300}.
#'
#' @examples
#' data(Data.Incomes)
#' y <- Data.Incomes[,1]
#' x <- as.matrix(Data.Incomes[,-c(1,2)])
#' PLR <- Lorenz.FABS(y, x, h = nrow(x)^(-1/5.5), eps = 0.005)
#' PLR.BIC(y, x, PLR$theta)
#'
#' @export

PLR.BIC <- function(y, x, theta, weights=NULL, IC=c("BIC","AIC")){

  IC <- match.arg(IC)

  n <- length(y)
  Index <- x%*%theta

  score <- log(apply(Index, 2, function(t) Gini.coef(y, x=t, na.rm=TRUE, ties.method="mean", weights=weights)))

  if (IC=="BIC") pen <- apply(theta,2,function(x)sum(x!=0))*log(n)/(2*n)
  if (IC=="AIC") pen <- apply(theta,2,function(x)sum(x!=0))/n

  IC.val <- score - pen
  IC.best <- which.max(IC.val)

  return.list <- list()
  return.list$val <- IC.val
  return.list$best <- IC.best

  return(return.list)

}
