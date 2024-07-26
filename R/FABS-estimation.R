#' Solves the Penalized Lorenz Regression with Lasso penalty
#'
#' \code{Lorenz.FABS} solves the penalized Lorenz regression with (adaptive) Lasso penalty on a grid of lambda values.
#' For each value of lambda, the function returns estimates for the vector of parameters and for the estimated explained Gini coefficient, as well as the Lorenz-\eqn{R^2} of the regression.
#'
#' The regression is solved using the FABS algorithm developed by Shi et al (2018) and adapted to our case.
#' For a comprehensive explanation of the Penalized Lorenz Regression, see Jacquemain et al.
#' In order to ensure identifiability, theta is forced to have a L2-norm equal to one.
#'
#' @param y a vector of responses
#' @param x a matrix of explanatory variables
#' @param standardize Should the variables be standardized before the estimation process? Default value is TRUE.
#' @param weights vector of sample weights. By default, each observation is given the same weight.
#' @param kernel integer indicating what kernel function to use. The value 1 is the default and implies the use of an Epanechnikov kernel while the value of 2 implies the use of a biweight kernel.
#' @param h bandwidth of the kernel, determining the smoothness of the approximation of the indicator function. Default value is n^(-1/5.5) where n is the sample size.
#' @param gamma value of the Lagrange multiplier in the loss function
#' @param lambda this parameter relates to the regularization parameter. Several options are available.
#' \describe{
#'     \item{\code{grid}}{If \code{lambda="grid"}, lambda is defined on a grid, equidistant in the logarithmic scale.}
#'     \item{\code{Shi}}{If \code{lambda="Shi"}, lambda, is defined within the algorithm, as in Shi et al (2018).}
#'     \item{\code{supplied}}{If the user wants to supply the lambda vector himself}
#' }
#' @param w.adaptive vector of size equal to the number of covariates where each entry indicates the weight in the adaptive Lasso. By default, each covariate is given the same weight (Lasso).
#' @param eps step size in the FABS algorithm. Default value is 0.005.
#' @param iter maximum number of iterations. Default value is 10^4.
#' @param lambda.min lower bound of the penalty parameter. Only used if \code{lambda="Shi"}.
#'
#' @return A list with several components:
#' \describe{
#'    \item{\code{lambda}}{vector gathering the different values of the regularization parameter}
#'    \item{\code{theta}}{matrix where column i provides the vector of estimated coefficients corresponding to the value \code{lambda[i]} of the regularization parameter.}
#'    \item{\code{LR2}}{vector where element i provides the Lorenz-\eqn{R^2} attached to the value \code{lambda[i]} of the regularization parameter.}
#'    \item{\code{Gi.expl}}{vector where element i provides the estimated explained Gini coefficient related to the value \code{lambda[i]} of the regularization parameter.}
#' }
#'
#' @seealso \code{\link{Lorenz.Reg}}, \code{\link{Lorenz.SCADFABS}}
#'
#' @section References:
#' Jacquemain, A., C. Heuchenne, and E. Pircalabelu (2024). A penalised bootstrap estimation procedure for the explained Gini coefficient. \emph{Electronic Journal of Statistics 18(1) 247-300}.
#'
#' Shi, X., Y. Huang, J. Huang, and S. Ma (2018). A Forward and Backward Stagewise Algorithm for Nonconvex Loss Function with Adaptive Lasso, \emph{Computational Statistics & Data Analysis 124}, 235-251.
#'
#' @examples
#' data(Data.Incomes)
#' y <- Data.Incomes[,1]
#' x <- as.matrix(Data.Incomes[,-c(1,2)])
#' Lorenz.FABS(y, x)
#'
#' @import MASS
#'
#' @export

# Largely based on the code proposed by Xingjie Shi on github
Lorenz.FABS <- function(y, x, standardize = TRUE, weights=NULL,
                        kernel = 1, h=length(y)^(-1/5.5), gamma = 0.05,
                        lambda="Shi", w.adaptive=NULL,
                        eps=0.005, iter=10^4, lambda.min = 1e-7){

  n <- length(y)
  p <- ncol(x)

  # Standardization

  if (standardize){

    x.center <- colMeans(x)
    x <- x - rep(x.center, rep.int(n,p))
    # x.scale <- sqrt(colSums(x^2)/(n-1))
    x.scale <- sqrt(colSums(x^2)/(n)) # Changé le 25-04-2022 pour assurer l'équivalence au niveau des catégorielles
    x <- x / rep(x.scale, rep.int(n,p))

  }

  # Observation weights

  if(any(weights<0)) stop("Weights must be nonnegative")

  if(is.null(weights)){
    weights <- rep(1,n)
  }
  pi <- weights/sum(weights)

  # Adaptive Lasso weights

  if(is.null(w.adaptive)) w <- rep(1,p) else w <- w.adaptive

  # We are going to record the lambda and the direction (backward or forward) for each iteration
  lambda.out <- direction <- numeric(iter)

  # FABS > INITIALIZATION ----

  b0 <- rep(0,p)
  b <- matrix(0, ncol=iter, nrow=p)
  b[,1] <- b0

  # Computing k
  Grad0 <- -.PLR_derivative_cpp(as.vector(y),as.matrix(x),as.vector(pi),as.vector(b0),as.double(h),as.double(gamma),as.integer(kernel))
  k0 <- which.max(abs(Grad0)/w)
  A.set <- k0

  # Computing beta
  b[k0,1] <- - sign(Grad0[k0])/w[k0]*eps

  # Computing lambda and the direction
  loss0 = .PLR_loss_cpp(as.matrix(x), as.vector(y), as.vector(pi), as.vector(b0), as.double(h),as.double(gamma),as.integer(kernel))
  loss  = .PLR_loss_cpp(as.matrix(x), as.vector(y), as.vector(pi), as.vector(b[,1]), as.double(h),as.double(gamma),as.integer(kernel))

  if(length(lambda)==1){
    # Either lambda="grid" or lambda="Shi". in both cases, the starting lambda is the same
    lambda.out[1] <- (loss0-loss)/eps
  }else{
    # The full lambda vector is supplied by the user
    lambda.out[1] <- lambda[1]
  }
  direction[1] <- 1

  # If lambda is supplied by the user, that vector is used as a grid.
  if(length(lambda)>1){
    lambda.grid <- lambda
    lambda.pointer <- 1
  }
  # If lambda="grid", lambda is given by a grid
  if(all(lambda=="grid")){
    lambda.upper <- lambda.out[1]
    lambda.lower <- lambda.upper*eps*0.001
    lambda.grid <- exp(seq(to=log(lambda.lower),from=log(lambda.upper),length.out=100))
    lambda.pointer <- 1
  }

  # FABS > BACKWARD AND FORWARD STEPS ----

  loss.i <- loss
  for (i in 1:(iter-1))
  {
    b[,i+1] <- b[,i]
    Grad.i <- -.PLR_derivative_cpp(as.vector(y),as.matrix(x),as.vector(pi),as.vector(b[,i]),as.double(h),as.double(gamma),as.integer(kernel))

    # Backward direction
    k <- A.set[which.min(-Grad.i[A.set]*sign(b[A.set,i])/w[A.set])]
    Delta.k <- -sign(b[k,i])/w[k]
    b[k,i+1] <- b[k,i] + Delta.k*eps
    loss.back <- .PLR_loss_cpp(as.matrix(x), as.vector(y), as.vector(pi), as.vector(b[,i+1]), as.double(h),as.double(gamma),as.integer(kernel))
    back <- loss.back - loss.i - lambda.out[i]*eps*w[k] < -.Machine$double.eps^0.5
    if(back & (length(A.set)>1)){
      # Backward step
      lambda.out[i+1] <- lambda.out[i]
      direction[i+1] <- -1
      loss.i <- loss.back
      if(abs(b[k,i+1]) < .Machine$double.eps^0.5){
        b[k,i+1] <- 0
        A.set <- setdiff(A.set,k)
      }
    }else{
      # Forward step
      b[k,i+1] <- b[k,i] # We must take out what we did in the backward step
      k <- which.max(abs(Grad.i)/w)
      A.set <- union(A.set,k)
      b[k,i+1] <- b[k,i] - sign(Grad.i[k])/w[k]*eps
      loss.forward <- .PLR_loss_cpp(as.matrix(x), as.vector(y), as.vector(pi), as.vector(b[,i+1]), as.double(h),as.double(gamma),as.integer(kernel))
      if (lambda.out[i] > (loss.i-loss.forward)/eps){
        # It means that with this lambda, I can no longer improve the score function. Hence, I have to update lambda
        if(!all(lambda=="Shi")){
          lambda.pointer <- lambda.pointer + 1
          if (lambda.pointer > length(lambda.grid)) break
          lambda.out[i+1] <- lambda.grid[lambda.pointer]
        }else{
          lambda.out[i+1] <- (loss.i-loss.forward)/eps
        }
      }else{
        lambda.out[i+1] <- lambda.out[i]
      }
      direction[i+1] <- 1
      loss.i <- loss.forward
    }
    if(all(lambda=="Shi")){
      if (lambda.out[i+1] <= lambda.out[1]*lambda.min) break
    }
    if (i==(iter-1))
      warning("Solution path unfinished, more iterations are needed.")
  }

  # We retrieve the different values along the path until algo stops
  iter <- i
  lambda <- lambda.out[1:iter]
  theta <- b[,1:iter]
  Index.sol <- x%*%theta
  LR2.num <- apply(Index.sol, 2, function(t) Gini.coef(y, x=t, na.rm=TRUE, ties.method="mean", weights=weights))
  LR2.denom <- Gini.coef(y, na.rm=TRUE, ties.method="mean", weights=weights)
  LR2<-as.numeric(LR2.num/LR2.denom)
  Gi.expl<-as.numeric(LR2.num)

  # At this stage, there are several iterations for each value of lambda. We need to retrieve only the last one.
  iter.unique <- c(which(diff(lambda)<0),iter)
  theta <- theta[,iter.unique]
  if (standardize) theta <- theta/x.scale # Need to put back on the original scale
  theta <- apply(theta,2,function(x)x/sqrt(sum(x^2))) # Need to normalize
  lambda <- lambda[iter.unique]
  LR2 <- LR2[iter.unique]
  Gi.expl <- Gi.expl[iter.unique]

  # OUTPUT ----
  return.list <- list(
    lambda = lambda,
    theta = theta,
    LR2=LR2,
    Gi.expl=Gi.expl
  )

  return(return.list)
}
