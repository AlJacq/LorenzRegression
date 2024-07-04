#' Summary for the Lorenz Regression
#'
#' \code{summary.LR} provides a summary for an object of class \code{"LR"}.
#'
#' @param object An object of class \code{"LR"}.
#' @param ... Additional arguments.
#'
#' @return An object of class \code{"summary.LR"}, containing the following elements:
#' \describe{
#'    \item{\code{call}}{The matched call.}
#'    \item{\code{ineq}}{A matrix with one row and three columns providing information on explained inequality. The first column gives the explained Gini coefficient, the second column gives the Gini coefficient of the response. The third column gives the Lorenz-R2.}
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
#' ## For examples see example(Lorenz.Reg) and example(Lorenz.boot)
#'
#' @method summary LR
#' @export

summary.LR <- function(object, ...){

  if (!inherits(object, "LR")) stop("The object must be of class 'LR'")

  ans <- list()
  ans$call <- object$call

  if(!is.null(object$theta)){

    ans$ineq <- matrix(c(object$Gi.expl,object$Gi.expl/object$LR2,object$LR2),
                       nrow = 1, ncol = 3,
                       dimnames = list("", c("Explained","Total","Lorenz-R2")))
    ans$coefficients <- as.matrix(object$theta)
    colnames(ans$coefficients) <- "Estimate"
    if (inherits(object, "LR_boot")){
      n <- nrow(object$x)
      p <- ncol(object$x)
      theta.boot <- object$boot_out$t[,3:ncol(object$boot_out$t)]
      Sigma.star <- n*stats::var(theta.boot)
      c.std <- sqrt(diag(Sigma.star)/n)
      ans$coefficients <- cbind(ans$coefficients, c.std)
      c.z <- object$theta/c.std
      ans$coefficients <- cbind(ans$coefficients, c.z)
      c.p <- sapply(1:p,function(k)2*stats::pnorm(abs(c.z[k]),lower.tail=FALSE))
      ans$coefficients <- cbind(ans$coefficients, c.p)
      colnames(ans$coefficients) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    }

  }else{

    ans$ineq <- NULL
    ans$coefficients <- NULL

  }



  class(ans) <- "summary.LR"

  return(ans)

}
