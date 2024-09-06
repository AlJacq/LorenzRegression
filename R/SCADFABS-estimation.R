#' Estimates the parameter vector in a penalized Lorenz regression with SCAD penalty
#'
#' \code{Lorenz.SCADFABS} solves the penalized Lorenz regression with SCAD penalty on a grid of lambda values.
#' For each value of lambda, the function returns estimates for the vector of parameters and for the estimated explained Gini coefficient, as well as the Lorenz-\eqn{R^2} of the regression.
#'
#' The regression is solved using the SCAD-FABS algorithm developed by Jacquemain et al and adapted to our case.
#' For a comprehensive explanation of the Penalized Lorenz Regression, see Heuchenne et al.
#' In order to ensure identifiability, theta is forced to have a L2-norm equal to one.
#'
#' @param y a vector of responses
#' @param x a matrix of explanatory variables
#' @param standardize Should the variables be standardized before the estimation process? Default value is TRUE.
#' @param weights vector of sample weights. By default, each observation is given the same weight.
#' @param kernel integer indicating what kernel function to use. The value 1 is the default and implies the use of an Epanechnikov kernel while the value of 2 implies the use of a biweight kernel.
#' @param h bandwidth of the kernel, determining the smoothness of the approximation of the indicator function. Default value is n^(-1/5.5) where n is the sample size.
#' @param gamma value of the Lagrange multiplier in the loss function
#' @param a parameter of the SCAD penalty. Default value is 3.7.
#' @param lambda this parameter relates to the regularization parameter. Several options are available.
#' \describe{
#'     \item{\code{grid}}{If lambda="grid", lambda is defined on a grid, equidistant in the logarithmic scale.}
#'     \item{\code{Shi}}{If lambda="Shi", lambda, is defined within the algorithm, as in Shi et al (2018).}
#'     \item{\code{supplied}}{If the user wants to supply the lambda vector himself}
#' }
#' @param eps step size in the FABS algorithm. Default value is 0.005.
#' @param SCAD.nfwd optional tuning parameter used if penalty="SCAD". Default value is NULL. The larger the value of this parameter, the sooner the path produced by the SCAD will differ from the path produced by the LASSO.
#' @param iter maximum number of iterations. Default value is 10^4.
#' @param lambda.min lower bound of the penalty parameter. Only used if lambda="Shi".
#'
#' @return A list with several components:
#' \describe{
#'    \item{\code{lambda}}{vector gathering the different values of the regularization parameter}
#'    \item{\code{theta}}{matrix where column i provides the vector of estimated coefficients corresponding to the value \code{lambda[i]} of the regularization parameter.}
#'    \item{\code{LR2}}{vector where element i provides the Lorenz-\eqn{R^2} attached to the value \code{lambda[i]} of the regularization parameter.}
#'    \item{\code{Gi.expl}}{vector where element i provides the estimated explained Gini coefficient related to the value \code{lambda[i]} of the regularization parameter.}
#' }
#'
#' @seealso \code{\link{Lorenz.Reg}}, \code{\link{Lorenz.FABS}}
#'
#' @section References:
#' Jacquemain, A., C. Heuchenne, and E. Pircalabelu (2024). A penalised bootstrap estimation procedure for the explained Gini coefficient. \emph{Electronic Journal of Statistics 18(1) 247-300}.
#'
#' @examples
#' data(Data.Incomes)
#' y <- Data.Incomes[,1]
#' x <- as.matrix(Data.Incomes[,-c(1,2)])
#' Lorenz.SCADFABS(y, x)
#'
#' @import MASS
#'
#' @export

Lorenz.SCADFABS <- function(y, x, standardize = TRUE, weights=NULL,
                            kernel = 1, h=length(y)^(-1/5.5), gamma = 0.05,
                            a = 3.7, lambda="Shi",
                            eps = 0.005, SCAD.nfwd = NULL, iter=10^4, lambda.min = 1e-7){

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

  # We are going to record the lambda and the direction (backward or forward) for each iteration
  lambda.out <- direction <- numeric(iter)

  # SCAD-FABS > n_fwd argument

  if(!is.null(SCAD.nfwd)){

    b1 <- b0 <- rep(0,p)
    Grad0 <- -.PLR_derivative_cpp(as.vector(y),
                                  as.matrix(x),
                                  as.vector(weights/sum(weights)),
                                  as.vector(b0),
                                  as.double(h),
                                  as.double(gamma),
                                  as.integer(kernel))
    k0 <- which.max(abs(Grad0))
    b1[k0] <- - sign(Grad0[k0])*eps
    loss0 = .PLR_loss_cpp(as.matrix(x),
                          as.vector(y),
                          as.vector(weights/sum(weights)),
                          as.vector(b0),
                          as.double(h),
                          as.double(gamma),
                          as.integer(kernel))
    loss1  = .PLR_loss_cpp(as.matrix(x),
                           as.vector(y),
                           as.vector(weights/sum(weights)),
                           as.vector(b1),
                           as.double(h),
                           as.double(gamma),
                           as.integer(kernel))
    diff.loss.sqrt <- sqrt(loss0-loss1)
    eps.old <- eps
    eps <- diff.loss.sqrt/sqrt(SCAD.nfwd) + sqrt(.Machine$double.eps)
    h <- h*eps/eps.old

  }

  # SCAD-FABS > INITIALIZATION ----

  b0 <- rep(0,p)
  b <- matrix(0, ncol=iter, nrow=p)
  b[,1] <- b0

  # Computing k
  Grad0 <- -.PLR_derivative_cpp(as.vector(y),as.matrix(x),as.vector(pi),as.vector(b0),as.double(h),as.double(gamma),as.integer(kernel))
  k0 <- which.max(abs(Grad0))
  A.set <- k0
  B.set <- 1:p

  # Computing beta
  b[k0,1] <- - sign(Grad0[k0])*eps

  # Computing lambda and the direction
  loss0 = .PLR_loss_cpp(as.matrix(x), as.vector(y),as.vector(pi), as.vector(b0), as.double(h),as.double(gamma),as.integer(kernel))
  loss  = .PLR_loss_cpp(as.matrix(x), as.vector(y),as.vector(pi), as.vector(b[,1]), as.double(h),as.double(gamma),as.integer(kernel))
  direction[1] <- 1

  if(length(lambda)==1){
    # Either lambda="grid" or lambda="Shi". in both cases, the starting lambda is the same
    lambda.out[1] <- (loss0-loss)/eps
  }else{
    # The full lambda vector is supplied by the user
    lambda.out[1] <- lambda[1]
  }
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
    Grad.loss.i <- -.PLR_derivative_cpp(as.vector(y),as.matrix(x),as.vector(pi),as.vector(b[,i]),as.double(h),as.double(gamma),as.integer(kernel))
    Grad.Pen.i <- .SCAD_derivative_cpp(as.vector(abs(b[,i])), as.double(lambda.out[i]), as.double(a))
    # Backward direction
    Back.Obj <- -Grad.loss.i[A.set]*sign(b[A.set,i]) - Grad.Pen.i[A.set]
    k <- A.set[which.min(Back.Obj)]
    Delta.k <- -sign(b[k,i])
    b[k,i+1] <- b[k,i] + Delta.k*eps
    Back.Obj.opt <- -Grad.loss.i[k]*sign(b[k,i]) - Grad.Pen.i[k]
    loss.back <- .PLR_loss_cpp(as.matrix(x), as.vector(y),as.vector(pi), as.vector(b[,i+1]), as.double(h),as.double(gamma),as.integer(kernel))
    back <- loss.back - loss.i - Grad.Pen.i[k]*eps < -.Machine$double.eps^0.5
    if(back & (length(A.set)>1)){
      # Backward step
      lambda.out[i+1] <- lambda.out[i]
      direction[i+1] <- -1
      loss.i <- .PLR_loss_cpp(as.matrix(x), as.vector(y),as.vector(pi), as.vector(b[,i+1]), as.double(h),as.double(gamma),as.integer(kernel))
      if(abs(b[k,i+1]) < .Machine$double.eps^0.5){
        b[k,i+1] <- 0
        A.set <- setdiff(A.set,k)
      }
    }else{
      if(gamma > 0) B.set <- 1:p # We reset it at each iteration (except when gamma = 0 because the algo tends to go crazy)
      L_eps.check <- FALSE
      # Forward step
      while(!L_eps.check){
        b[k,i+1] <- b[k,i] # We must take out what we did in the backward step
        Fwd.Obj <- abs(Grad.loss.i[B.set]) - Grad.Pen.i[B.set]
        k <- B.set[which.max(Fwd.Obj)]
        A.set <- union(A.set,k)
        b[k,i+1] <- b[k,i] - sign(Grad.loss.i[k])*eps
        loss.forward <- .PLR_loss_cpp(as.matrix(x), as.vector(y),as.vector(pi), as.vector(b[,i+1]), as.double(h),as.double(gamma),as.integer(kernel))
        L_eps <- (loss.i-loss.forward)/eps
        L_eps.check <- L_eps > 0
        if(L_eps.check){
          B.set <- union(B.set,k)
        }else{
          B.set <- setdiff(B.set,k)
        }
        if (length(B.set)==0) L_eps.check <- TRUE
      }

      lambda_A <- L_eps
      lambda_B <- ((a-1)*L_eps+abs(b[k,i]))/a
      lambda.update <- (abs(b[k,i]) < a*lambda.out[i]) & (lambda.out[i] > max(lambda_A,lambda_B))
      if(lambda.update){
        if(!all(lambda=="Shi")){
          lambda.pointer <- lambda.pointer + 1
          if (lambda.pointer > length(lambda.grid)) break
          lambda.out[i+1] <- lambda.grid[lambda.pointer]
        }else{
          lambda.out[i+1] <- max(lambda_A,lambda_B)
        }
      }else{
        lambda.out[i+1] <- lambda.out[i]
      }

      direction[i+1] <- 1
      loss.i <- loss.forward
    }
    if(all(lambda=="Shi")){
      if (lambda.out[i+1] <= lambda.min) break
    }
    if (length(B.set)==0) break
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
