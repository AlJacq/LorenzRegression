#' Summary for the Penalized Lorenz Regression
#'
#' \code{summary.PLR} provides a summary for an object of class \code{"PLR"}.
#'
#' @param object An object of class \code{"PLR"}.
#' @param ... Additional arguments
#'
#' @return An object of class \code{"summary.PLR"}, containing the following elements:
#' \describe{
#'    \item{\code{call}}{The matched call.}
#'    \item{\code{Gi.expl}}{The estimated explained Gini coefficient.}
#'    \item{\code{LR2}}{The Lorenz-\eqn{R^2} of the regression.}
#'    \item{\code{coefficients}}{A matrix providing information on the estimated coefficients. The first column gives the estimates.
#'    If the class of \code{object} contains \code{LR_boot}, bootstrap inference was performed and the matrix contains further information. The second column is the boostrap standard error. The third column is the z-value. Finally, the last column is the p-value.}
#' }
#'
#' @details
#' The inference provided in the \code{coefficients} matrix is obtained by using the asymptotic normality and estimating the asymptotic variance via bootstrap.
#'
#' @seealso \code{\link{Lorenz.Reg}}, \code{\link{Lorenz.boot}}
#'
#' @examples
#' data(Data.Incomes)
#' NPLR <- Lorenz.Reg(Income ~ ., data = Data.Incomes, penalty = "none")
#' summary(NPLR)
#'
#' @method summary PLR
#' @export

summary.PLR <- function(object, renormalize=TRUE, ...){

  if (!inherits(object, "PLR")) stop("The object must be of class 'PLR'")

  if(renormalize){
    m1 <- PLR.normalize(object)
  }else{
    m1 <- object$theta
  }

  ans <- list()
  ans$call <- object$call
  if(length(object$Gi.expl)==1) names(object$Gi.expl) <- ""
  ans$ineq <- matrix(c(object$Gi.expl,object$Gi.expl/object$LR2,object$LR2),
                     nrow = length(object$Gi.expl), ncol = 3,
                     dimnames = list(names(object$Gi.expl), c("Explained","Total","Lorenz-R2")))


  if (inherits(object, c("PLR_boot","PLR_cv"))){
    ans$coefficients <- t(m1)
  }else{
    ans$coefficients <- as.matrix(m1)
    colnames(ans$coefficients) <- "Estimate"
  }

  class(ans) <- "summary.PLR"

  return(ans)

}
