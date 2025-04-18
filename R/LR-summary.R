#' Summary for the Lorenz regression
#'
#' Provides a summary for an object of class \code{"LR"}.
#'
#' @aliases summary.LR_boot
#' @param object An object of class \code{"LR"}. The object might also have S3 class \code{"LR_boot"} (which inherits from class \code{"PLR"}).
#' @param ... Additional arguments.
#'
#' @return An object of class \code{"summary.LR"}, containing the following elements:
#' \describe{
#'    \item{\code{call}}{The matched call.}
#'    \item{\code{ineq}}{A matrix with one row and three columns providing information on explained inequality. The first column gives the explained Gini coefficient, the second column gives the Gini coefficient of the response. The third column gives the Lorenz-\eqn{R^2}.}
#'    \item{\code{coefficients}}{A matrix providing information on the estimated coefficients. The first column gives the estimates.
#'    If \code{object} inherits from \code{"LR_boot"}, bootstrap inference was performed and the matrix contains further information. The second column is the boostrap standard error. The third column is the z-value. Finally, the last column is the p-value.
#'    In this case, the class \code{"summary.LR_boot"} is added to the output.}
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

  ans <- list()
  ans$call <- object$call

  if(length(all.vars(object$terms))>1){

    ans$ineq <- matrix(c(ineqExplained.LR(object,type="Gini.explained"),
                         ineqExplained.LR(object,type="Gini.explained")/ineqExplained.LR(object,type="Lorenz-R2"),
                         ineqExplained.LR(object,type="Lorenz-R2")),
                       nrow = 1, ncol = 3,
                       dimnames = list("", c("Explained","Total","Lorenz-R2")))
    ans$coefficients <- as.matrix(coef.LR(object))
    colnames(ans$coefficients) <- "Estimate"

  }else{

    ans$ineq <- NULL
    ans$coefficients <- NULL

  }

  class(ans) <- "summary.LR"

  return(ans)

}

#' @method summary LR_boot
#' @export

summary.LR_boot <- function(object, ...){

  ans <- NextMethod("summary")

  if(length(all.vars(object$terms))>1){
    n <- nrow(object$x)
    p <- ncol(object$x)
    idx.theta <- (1:length(object$theta))+2
    theta.boot <- object$boot_out$t[,idx.theta]
    Sigma.star <- n*stats::var(theta.boot)
    c.std <- sqrt(diag(Sigma.star)/n)
    ans$coefficients <- cbind(ans$coefficients, c.std)
    c.z <- coef.LR(object)/c.std
    ans$coefficients <- cbind(ans$coefficients, c.z)
    c.p <- sapply(1:p,function(k)2*stats::pnorm(abs(c.z[k]),lower.tail=FALSE))
    ans$coefficients <- cbind(ans$coefficients, c.p)
    colnames(ans$coefficients) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")

    class(ans) <- c("summary.LR_boot",class(ans))

  }

  return(ans)

}
