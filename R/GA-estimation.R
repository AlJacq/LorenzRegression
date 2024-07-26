#' Estimates the parameter vector in Lorenz regression using a genetic algorithm
#'
#' \code{Lorenz.GA} estimates the coefficient vector of the single-index model.
#' It also returns the Lorenz-\eqn{R^2} of the regression as well as the estimated explained Gini coefficient.
#'
#' The genetic algorithm is solved using function \code{\link[GA]{ga}} from the \emph{GA} package. The fitness function is coded in Rcpp to speed up computation time.
#' When discrete covariates are introduced and ties occur in the index, the default option randomly breaks them, as advised in Section 3 of Heuchenne and Jacquemain (2022)
#'
#' @param y a vector of responses
#' @param x a matrix of explanatory variables
#' @param standardize Should the variables be standardized before the estimation process? Default value is TRUE.
#' @param weights vector of sample weights. By default, each observation is given the same weight.
#' @param popSize Size of the population of candidates in the genetic algorithm. Default value is 50.
#' @param maxiter Maximum number ot iterations in the genetic algorithm. Default value is 1500.
#' @param run Number of iterations without improvement in the best fitness necessary for the algorithm to stop. Default value is 150.
#' @param ties.method What method should be used to break the ties in optimization program. Possible values are "random" (default value) or "mean". If "random" is selected, the ties are broken by further ranking in terms of a uniformly distributed random variable. If "mean" is selected, the average rank method is used.
#' @param ties.Gini what method should be used to break the ties in the computation of the Gini coefficient at the end of the algorithm. Possible values and default choice are the same as above.
#' @param seed.random seed.random imposed for the generation of the vector of uniform random variables used to break the ties. Default is NULL, in which case no seed.random is imposed.
#' @param parallel.GA Whether parallel computing should be used to distribute the computations in the genetic algorithm. Either a logical value determining whether parallel computing is used (TRUE) or not (FALSE, the default value). Or a numerical value determining the number of cores to use.
#'
#' @return A list with several components:
#' \describe{
#'    \item{\code{theta}}{the estimated vector of parameters.}
#'    \item{\code{LR2}}{the Lorenz-\eqn{R^2} of the regression.}
#'    \item{\code{Gi.expl}}{the estimated explained Gini coefficient.}
#'    \item{\code{niter}}{number of iterations attained by the genetic algorithm.}
#'    \item{\code{fit}}{value attained by the fitness function at the optimum.}
#' }
#'
#' @seealso \code{\link{Lorenz.Reg}}, \code{\link[GA]{ga}}
#'
#' @section References:
#' Heuchenne, C. and A. Jacquemain (2022). Inference for monotone single-index conditional means: A Lorenz regression approach. \emph{Computational Statistics & Data Analysis 167(C)}.
#'
#' @examples
#' data(Data.Incomes)
#' y <- Data.Incomes$Income
#' x <- cbind(Data.Incomes$Age, Data.Incomes$Work.Hours)
#' Lorenz.GA(y, x, popSize = 40)
#'
#' @import GA
#'
#' @export

# unit-norm normalization ----
Lorenz.GA<-function(y, x, standardize=TRUE, weights=NULL, popSize=50, maxiter=1500, run=150, ties.method=c("random","mean"), ties.Gini=c("random","mean"), seed.random=NULL, parallel.GA = FALSE){

  # PRE-GA ----

  ties.method <- match.arg(ties.method)
  ties.Gini <- match.arg(ties.Gini)

  n <- length(y)
  p <- ncol(x)

  if(any(weights<0)) stop("Weights must be nonnegative")

  if(is.null(weights)){
    weights <- rep(1,n)
  }
  pi <- weights/sum(weights)

  if(p > 1){

    # PRE-GA > STANDARDIZE X ----

    if (standardize){

      x.center <- colMeans(x)
      x <- x - rep(x.center, rep.int(n,p))
      x.scale <- sqrt(colSums(x^2)/(n-1))
      x <- x / rep(x.scale, rep.int(n,p))

    }

    # GA ----

    if (ties.method == "random"){

      if(!is.null(seed.random)) set.seed(seed.random)
      V<-stats::runif(n)

      GA <- GA::ga(type = "real-valued",
                   population = Lorenz.Population,
                   fitness =  function(u).Fitness_cpp(u,as.vector(y),as.matrix(x),V,pi),
                   lower = rep(-1,p-1), upper = rep(1,p-1),
                   popSize = popSize, maxiter = maxiter, run = run, monitor = FALSE,
                   parallel = parallel.GA)

      # We need to take care of the fact that the first coefficient for theta may be positive or negative
      theta1<-c(GA@solution[1,],1-sum(abs(GA@solution[1,]))) #The theta solution if the last coeff is positive
      theta2<-c(GA@solution[1,],-(1-sum(abs(GA@solution[1,])))) #The theta solution if the last coeff is negative
      theta<-rbind(theta1,theta2)
      Index_1<-theta1%*%t(x)
      Y_1<-y[order(Index_1,V)]
      pi_1<-pi[order(Index_1,V)]
      rank_1<-cumsum(pi_1)-pi_1/2
      Index_2<-theta2%*%t(x)
      Y_2<-y[order(Index_2,V)]
      pi_2<-pi[order(Index_2,V)]
      rank_2<-cumsum(pi_2)-pi_2/2
      theta.argmax<-theta[which.max(c((Y_1*pi_1)%*%rank_1,(Y_2*pi_2)%*%rank_2)),]

      # We compute the Lorenz-Rsquared
      Index.sol<-x%*%theta.argmax

    }

    if (ties.method == "mean"){

      GA <- GA::ga(type = "real-valued",
                   population = Lorenz.Population,
                   fitness =  function(u).Fitness_meanrank(u,as.vector(y),as.matrix(x),pi),
                   lower = rep(-1,p-1), upper = rep(1,p-1),
                   popSize = popSize, maxiter = maxiter, run = run, monitor = FALSE,
                   parallel = parallel.GA)

      # We need to take care of the fact that the first coefficient for theta may be positive or negative
      theta1<-c(GA@solution[1,],1-sum(abs(GA@solution[1,]))) #The theta solution if the last coeff is positive
      theta2<-c(GA@solution[1,],-(1-sum(abs(GA@solution[1,])))) #The theta solution if the last coeff is negative
      theta<-rbind(theta1,theta2)
      index1<-theta1%*%t(x)
      index2<-theta2%*%t(x)
      index1_k <- sort(unique(index1))
      pi1_k <- sapply(1:length(index1_k),function(k)sum(pi[index1==index1_k[k]]))
      F1_k <- cumsum(pi1_k) - 0.5*pi1_k
      F1_i <- sapply(1:length(index1),function(i)sum(F1_k[index1_k==index1[i]])) # Ensures that sum(F_i*pi) = 0.5
      index2_k <- sort(unique(index2))
      pi2_k <- sapply(1:length(index2_k),function(k)sum(pi[index2==index2_k[k]]))
      F2_k <- cumsum(pi2_k) - 0.5*pi2_k
      F2_i <- sapply(1:length(index2),function(i)sum(F2_k[index2_k==index2[i]])) # Ensures that sum(F_i*pi) = 0.5
      theta.argmax<-theta[which.max(c((pi*y)%*%F1_i,(pi*y)%*%F2_i)),]

      # We compute the Lorenz-Rsquared
      Index.sol<-x%*%theta.argmax

    }

    # POST-LR ----

    if(ties.Gini == "random"){

      LR2.num<-Gini.coef(y, x=Index.sol, na.rm=TRUE, ties.method="random", seed=seed.random, weights=weights)
      LR2.denom<-Gini.coef(y, na.rm=TRUE, ties.method="random", seed=seed.random, weights=weights)
      LR2<-as.numeric(LR2.num/LR2.denom)
      Gi.expl<-as.numeric(LR2.num)

    }else{

      LR2.num<-Gini.coef(y, x=Index.sol, na.rm=TRUE, ties.method="mean", seed=seed.random, weights=weights)
      LR2.denom<-Gini.coef(y, na.rm=TRUE, ties.method="mean", seed=seed.random, weights=weights)
      LR2<-as.numeric(LR2.num/LR2.denom)
      Gi.expl<-as.numeric(LR2.num)

    }

    if (standardize) theta.argmax <- theta.argmax/x.scale # Need to put back on the original scale
    theta <- theta.argmax/sqrt(sum(theta.argmax^2))
    niter <- length(GA@summary[,1])
    fit <- GA@fitnessValue

  }else{

    # Empty model
    if(all(x==1)){

      theta <- niter <- fit <- NULL
      Gi.expl <- LR2 <- 0

    # Only one covariate
    }else{

      if(ties.Gini == "random"){

        CI<-Gini.coef(y, x=x[,1], na.rm=TRUE, ties.method="random", seed=seed.random, weights=weights)
        LR2.denom<-Gini.coef(y, na.rm=TRUE, ties.method="random", seed=seed.random, weights=weights)
        Gi.expl<-as.numeric(abs(CI))
        LR2<-as.numeric(Gi.expl/LR2.denom)


      }else{

        CI<-Gini.coef(y, x=x[,1], na.rm=TRUE, ties.method="mean", seed=seed.random, weights=weights)
        LR2.denom<-Gini.coef(y, na.rm=TRUE, ties.method="mean", seed=seed.random, weights=weights)
        Gi.expl<-as.numeric(abs(CI))
        LR2<-as.numeric(Gi.expl/LR2.denom)

      }

      niter <- fit <- NULL
      theta <- sign(CI)

    }

  }

  Return.list <- list()
  Return.list$theta <- theta
  Return.list$LR2 <- LR2
  Return.list$Gi.expl <- Gi.expl
  Return.list$niter <- niter
  Return.list$fit <- fit

  return(Return.list)

}
