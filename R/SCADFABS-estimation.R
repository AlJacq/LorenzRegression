#' Solves the Penalized Lorenz Regression with SCAD penalty
#'
#' \code{Lorenz.SCADFABS} solves the penalized Lorenz regression with SCAD penalty on a grid of lambda values.
#' For each value of lambda, the function returns estimates for the vector of parameters and for the estimated explained Gini coefficient, as well as the Lorenz-\eqn{R^2} of the regression.
#'
#' The regression is solved using the SCAD-FABS algorithm developed by Jacquemain et al and adapted to our case.
#' For a comprehensive explanation of the Penalized Lorenz Regression, see Heuchenne et al.
#' In order to ensure identifiability, theta is forced to have a L2-norm equal to one.
#'
#' @param YX_mat a matrix with the first column corresponding to the response vector, the remaining ones being the explanatory variables.
#' @param weights vector of sample weights. By default, each observation is given the same weight.
#' @param h bandwidth of the kernel, determining the smoothness of the approximation of the indicator function.
#' @param eps step size in the FABS algorithm.
#' @param a parameter of the SCAD penalty. Default value is 3.7.
#' @param iter maximum number of iterations. Default value is 10^4.
#' @param lambda this parameter relates to the regularization parameter. Several options are available.
#' \describe{
#'     \item{\code{grid}}{If lambda="grid", lambda is defined on a grid, equidistant in the logarithmic scale.}
#'     \item{\code{Shi}}{If lambda="Shi", lambda, is defined within the algorithm, as in Shi et al (2018).}
#'     \item{\code{supplied}}{If the user wants to supply the lambda vector himself}
#' }
#' @param lambda.min lower bound of the penalty parameter. Only used if lambda="Shi".
#' @param gamma value of the Lagrange multiplier in the loss function
#'
#' @return A list with several components:
#' \describe{
#'    \item{\code{iter}}{number of iterations attained by the algorithm.}
#'    \item{\code{direction}}{vector providing the direction (-1 = backward step, 1 = forward step) for each iteration.}
#'    \item{\code{lambda}}{value of the regularization parameter for each iteration.}
#'    \item{\code{h}}{value of the bandwidth.}
#'    \item{\code{theta}}{matrix where column i provides the non-normalized estimated parameter vector for iteration i.}
#'    \item{\code{LR2}}{vector where element i provides the Lorenz-\eqn{R^2} of the regression for iteration i.}
#'    \item{\code{Gi.expl}}{vector where element i provides the estimated explained Gini coefficient for iteration i.}
#' }
#'
#' @seealso \code{\link{Lorenz.Reg}}, \code{\link{PLR.wrap}}, \code{\link{Lorenz.FABS}}
#'
#' @section References:
#' Jacquemain, A., C. Heuchenne, and E. Pircalabelu (2022). A penalised bootstrap estimation procedure for the explained Gini coefficient.
#'
#' @examples
#' data(Data.Incomes)
#' YX_mat <- Data.Incomes
#' Lorenz.SCADFABS(YX_mat, h = nrow(Data.Incomes)^(-1/5.5), eps = 0.005)
#'
#' @import MASS
#'
#' @export

Lorenz.SCADFABS <- function(YX_mat, weights=NULL, h, eps, a = 3.7,
                iter=10^4, lambda="Shi", lambda.min = 1e-7, gamma = 0.05){

  X <- YX_mat[,-1]
  y <- YX_mat[,1]

  n <- length(y)
  p <- ncol(X)

  # Observation weights

  if(any(weights<0)) stop("Weights must be nonnegative")

  if(is.null(weights)){
    weights <- rep(1,n)
  }
  pi <- weights/sum(weights)

  # We are going to record the lambda and the direction (backward or forward) for each iteration
  lambda.out <- direction <- numeric(iter)

  # SCAD-FABS > INITIALIZATION ----

  b0 <- rep(0,ncol(X))
  b <- matrix(0, ncol=iter, nrow=p)
  b[,1] <- b0

  # Computing k
  Grad0 <- -PLR_derivative_cpp(as.vector(y),as.matrix(X),as.vector(pi),as.vector(b0),as.double(h),as.double(gamma))
  k0 <- which.max(abs(Grad0))
  A.set <- k0
  B.set <- 1:p

  # Computing beta
  b[k0,1] <- - sign(Grad0[k0])*eps

  # Computing lambda and the direction
  loss0 = PLR_loss_cpp(as.matrix(X), as.vector(y),as.vector(pi), as.vector(b0), as.double(h),as.double(gamma))
  loss  = PLR_loss_cpp(as.matrix(X), as.vector(y),as.vector(pi), as.vector(b[,1]), as.double(h),as.double(gamma))
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
    Grad.loss.i <- -PLR_derivative_cpp(as.vector(y),as.matrix(X),as.vector(pi),as.vector(b[,i]),as.double(h),as.double(gamma))
    Grad.Pen.i <- SCAD_derivative_cpp(as.vector(abs(b[,i])), as.double(lambda.out[i]), as.double(a))
    # Backward direction
    Back.Obj <- -Grad.loss.i[A.set]*sign(b[A.set,i]) - Grad.Pen.i[A.set]
    k <- A.set[which.min(Back.Obj)]
    Delta.k <- -sign(b[k,i])
    b[k,i+1] <- b[k,i] + Delta.k*eps
    Back.Obj.opt <- -Grad.loss.i[k]*sign(b[k,i]) - Grad.Pen.i[k]
    loss.back <- PLR_loss_cpp(as.matrix(X), as.vector(y),as.vector(pi), as.vector(b[,i+1]), as.double(h),as.double(gamma))
    back <- loss.back - loss.i - Grad.Pen.i[k]*eps < -.Machine$double.eps^0.5
    if(back & (length(A.set)>1)){
      # Backward step
      lambda.out[i+1] <- lambda.out[i]
      direction[i+1] <- -1
      loss.i <- PLR_loss_cpp(as.matrix(X), as.vector(y),as.vector(pi), as.vector(b[,i+1]), as.double(h),as.double(gamma))
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
        loss.forward <- PLR_loss_cpp(as.matrix(X), as.vector(y),as.vector(pi), as.vector(b[,i+1]), as.double(h),as.double(gamma))
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

  # We compute the Lorenz-Rsquared and explained Gini coef along the path
  theta <- b[,1:i]
  Index.sol <- as.matrix(X)%*%theta

  LR2.num <- apply(Index.sol, 2, function(t) Gini.coef(y, x=t, na.rm=T, ties.method="mean", weights=weights))
  LR2.denom <- Gini.coef(y, na.rm=T, ties.method="mean", weights=weights)
  LR2<-as.numeric(LR2.num/LR2.denom)
  Gi.expl<-as.numeric(LR2.num)

  # WARNING ----

  # If eps is too large, the path may be really rough
  # if (length(unique(lambda.out[1:i]))<5) warning("The algorithm generated less than 5 different values for lambda. We suggest you to consider decreasing eps to have a finer grid")

  # OUTPUT ----

  return.list <- list(
    iter = i,
    direction = direction[1:i],
    lambda = lambda.out[1:i],
    h = h,
    theta = theta,
    LR2=LR2,
    Gi.expl=Gi.expl
  )

  return(return.list)
}
