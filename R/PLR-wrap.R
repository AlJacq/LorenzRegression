#' Wrapper for the \code{\link{Lorenz.SCADFABS}} and \code{\link{Lorenz.FABS}} functions
#'
#' \code{PLR.wrap} standardizes the covariates, run the penalized regression and spits out the path of parameter vectors.
#'
#' @param YX_mat a matrix with the first column corresponding to the response vector, the remaining ones being the explanatory variables.
#' @param standardize Should the variables be standardized before the estimation process? Default value is TRUE.
#' @param weights vector of sample weights. By default, each observation is given the same weight.
#' @param penalty penalty used in the Penalized Lorenz Regression. Possible values are "SCAD" (default) or "LASSO".
#' @param h bandwidth of the kernel, determining the smoothness of the approximation of the indicator function.
#' @param SCAD.nfwd optional tuning parameter used if penalty="SCAD". Default value is NULL. The larger the value of this parameter, the sooner the path produced by the SCAD will differ from the path produced by the LASSO.
#' @param eps step size in the FABS or SCADFABS algorithm. Default value is 0.005.
#' @param gamma value of the Lagrange multiplier in the loss function
#' @param ... Additional parameters corresponding to arguments passed in \code{\link{Lorenz.SCADFABS}} or \code{\link{Lorenz.FABS}} depending on the argument chosen in penalty.
#'
#' @return A list with several components:
#' \describe{
#'    \item{\code{lambda}}{vector gathering the different values of the regularization parameter}
#'    \item{\code{theta}}{matrix where column i provides the normalized estimated parameter vector corresponding to value lambda[i] of the regularization parameter.}
#'    \item{\code{LR2}}{vector where element i provides the Lorenz-\eqn{R^2} of the regression related to value lambda[i] of the regularization parameter.}
#'    \item{\code{Gi.expl}}{vector where element i provides the estimated explained Gini coefficient related to value lambda[i] of the regularization parameter.}
#' }
#'
#' @seealso \code{\link{Lorenz.SCADFABS}}, \code{\link{Lorenz.FABS}}
#'
#' @examples
#' data(Data.Incomes)
#' YX_mat <- Data.Incomes[,-2]
#' PLR.wrap(YX_mat, h = nrow(Data.Incomes)^(-1/5.5), eps = 0.005)
#'
#' @export

PLR.wrap <- function(YX_mat, standardize=TRUE, weights=NULL, penalty=c("SCAD","LASSO"), h, SCAD.nfwd = NULL, eps = 0.005, gamma = 0.05, kernel = c("Epan","Biweight"), ...){

  penalty <- match.arg(penalty)
  kernel <- match.arg(kernel)
  kernel <- switch(kernel, "Epan" = 1, "Biweight" = 2)

  n <- length(YX_mat[,1])
  p <- length(YX_mat[1,])-1

  # PRE-PLR ----

  if (standardize){

    X <- YX_mat[,-1]
    X.center <- colMeans(X)
    X <- X - rep(X.center, rep.int(n,p))
    # X.scale <- sqrt(colSums(X^2)/(n-1))
    X.scale <- sqrt(colSums(X^2)/(n)) # Changé le 25-04-2022 pour assurer l'équivalence au niveau des catégorielles
    X <- X / rep(X.scale, rep.int(n,p))

    YX_mat[,-1] <- X

  }

  # PLR ----

  if(penalty == "SCAD"){
    if(is.null(SCAD.nfwd)){
      PLR <- LorenzRegression::Lorenz.SCADFABS(YX_mat, weights=weights, eps=eps, h=h, gamma = gamma, kernel = kernel, ...)
    }else{
      if(is.null(weights)){
        weights <- rep(1,n)
      }
      b1 <- b0 <- rep(0,ncol(YX_mat)-1)
      Grad0 <- -.PLR_derivative_cpp(as.vector(YX_mat[,1]),
                                    as.matrix(YX_mat[,-1]),
                                    as.vector(weights/sum(weights)),
                                    as.vector(b0),
                                    as.double(h),
                                    as.double(gamma),
                                    as.integer(kernel))
      k0 <- which.max(abs(Grad0))
      b1[k0] <- - sign(Grad0[k0])*eps
      loss0 = .PLR_loss_cpp(as.matrix(YX_mat[,-1]),
                            as.vector(YX_mat[,1]),
                            as.vector(weights/sum(weights)),
                            as.vector(b0),
                            as.double(h),
                            as.double(gamma),
                            as.integer(kernel))
      loss1  = .PLR_loss_cpp(as.matrix(YX_mat[,-1]),
                             as.vector(YX_mat[,1]),
                             as.vector(weights/sum(weights)),
                             as.vector(b1),
                             as.double(h),
                             as.double(gamma),
                             as.integer(kernel))
      diff.loss.sqrt <- sqrt(loss0-loss1)
      eps.old <- eps
      eps <- diff.loss.sqrt/sqrt(SCAD.nfwd) + sqrt(.Machine$double.eps)
      h <- h*eps/eps.old
      PLR <- LorenzRegression::Lorenz.SCADFABS(YX_mat, weights=weights, eps=eps, h=h, gamma = gamma, kernel = kernel, ...)
    }
  }else if(penalty == "LASSO"){
    PLR <- LorenzRegression::Lorenz.FABS(YX_mat, weights=weights, eps=eps, h=h, gamma = gamma, kernel = kernel, ...)
  }

  # POST-PLR ----

  # We need to retrieve one parameter vector for each value of lambda

  iter.unique <- c(which(diff(PLR$lambda)<0),PLR$iter)

  theta <- PLR$theta[,iter.unique] # Only one value for each value of lambda
  if (standardize) theta <- theta/X.scale # Need to put back on the original scale
  theta <- apply(theta,2,function(x)x/sqrt(sum(x^2))) # Need to normalize
  # theta[abs(theta) < 10^(-10) ] <- 0

  lambda <- PLR$lambda[iter.unique]
  LR2 <- PLR$LR2[iter.unique]
  Gi.expl <- PLR$Gi.expl[iter.unique]

  return.list <- list(
    lambda = lambda,
    theta = theta,
    LR2=LR2,
    Gi.expl=Gi.expl
  )

  return(return.list)
}
