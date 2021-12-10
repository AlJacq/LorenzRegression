#' Estimates the parameter vector in Lorenz regression using a genetic algorithm
#'
#' \code{Lorenz.GA.cpp} estimates the vector of parameters in Lorenz regression using the unit-norm normalization
#' It also returns the Lorenz-\eqn{R^2} of the regression as well as the estimated explained Gini coefficient.
#'
#' The genetic algorithm is solved using function \code{\link[GA]{ga}} from the \emph{GA} package. The fitness function is coded in Rcpp to speed up computation time.
#' When discrete covariates are introduced and ties occur in the index, the default option randomly breaks them, as advised in Section 3 of Heuchenne and Jacquemain (2020)
#'
#' @param YX_mat A matrix with the first column corresponding to the response vector, the remaining ones being the explanatory variables.
#' @param popSize Size of the population of candidates in the genetic algorithm. Default value is 50.
#' @param maxiter Maximum number ot iterations in the genetic algorithm. Default value is 1500.
#' @param run Number of iterations without improvement in the best fitness necessary for the algorithm to stop. Default value is 150.
#' @param ties.method What method should be used to break the ties in the rank index. Possible values are "random" (default value) or "mean". If "random" is selected, the ties are broken by further ranking in terms of a uniformly distributed random variable. If "mean" is selected, the average rank method is used.
#' @param seed seed imposed for the generation of the vector of uniform random variables used to break the ties. Default is NULL, in which case no seed is imposed.
#' @param weights vector of sample weights. By default, each observation is given the same weight.
#' @param parallel Whether parallel computing should be used in the genetic algorithm. Default value is FALSE.
#'
#' @return A list with several components:
#' \describe{
#'    \item{\code{sol}}{the estimated vector of parameters.}
#'    \item{\code{LR2}}{the Lorenz-\eqn{R^2} of the regression.}
#'    \item{\code{Gi.expl}}{the estimated explained Gini coefficient.}
#'    \item{\code{niter}}{number of iterations attained by the genetic algorithm.}
#'    \item{\code{fit}}{value attained by the fitness function at the optimum.}
#' }
#'
#' @seealso \code{\link{Lorenz.Reg}}, \code{\link[GA]{ga}}
#'
#' @section References:
#' Heuchenne, C. and A. Jacquemain. 2020. “Inference for monotone single-index conditional means: a Lorenz regression approach”
#'
#' @examples
#' data(Data.Incomes)
#' YX_mat <- cbind(Data.Incomes$Income, Data.Incomes$Age, Data.Incomes$Work.Hours)
#' Lorenz.GA.cpp(YX_mat, popSize = 40)
#'
#' @import GA
#'
#' @export

# unit-norm normalization ----
Lorenz.GA.cpp<-function(YX_mat,popSize=50,maxiter=1500,run=150, ties.method=c("random","mean"), seed=NULL, weights=NULL, parallel = F){

  ties.method <- match.arg(ties.method)

  n.param<-length(YX_mat[1,-1])
  n<-length(YX_mat[,1])

  if(any(weights<0)) stop("Weights must be nonnegative")

  if(is.null(weights)){
    weights <- rep(1,n)
  }
  pi <- weights/sum(weights)

  # The GA in itself
  if (ties.method == "random"){

    if(!is.null(seed)) set.seed(seed)
    V<-stats::runif(n)

    GA <- GA::ga(type = "real-valued",
                 population = Lorenz.Population,
                 fitness =  function(u)Fitness_cpp(u,as.vector(YX_mat[,1]),as.matrix(YX_mat[,-1]),V,pi),
                 lower = rep(-1,n.param-1), upper = rep(1,n.param-1),
                 popSize = popSize, maxiter = maxiter, run = run, monitor = FALSE,
                 parallel = parallel)

    # We need to take care of the fact that the first coefficient for theta may be positive or negative
    theta1<-c(GA@solution[1,],1-sum(abs(GA@solution[1,]))) #The theta solution if the last coeff is positive
    theta2<-c(GA@solution[1,],-(1-sum(abs(GA@solution[1,])))) #The theta solution if the last coeff is negative
    theta<-rbind(theta1,theta2)
    Index_1<-theta1%*%t(YX_mat[,-1])
    Y_1<-YX_mat[order(Index_1,V),1]
    pi_1<-pi[order(Index_1,V)]
    rank_1<-cumsum(pi_1)-pi_1/2
    Index_2<-theta2%*%t(YX_mat[,-1])
    Y_2<-YX_mat[order(Index_2,V),1]
    pi_2<-pi[order(Index_2,V)]
    rank_2<-cumsum(pi_2)-pi_2/2
    theta.argmax<-theta[which.max(c((Y_1*pi_1)%*%rank_1,(Y_2*pi_2)%*%rank_2)),]

    # We compute the Lorenz-Rsquared
    Index.sol<-theta.argmax%*%t(YX_mat[,-1])
    Y<-YX_mat[,1]

    LR2.num<-Gini.coef(Y, x=Index.sol, na.rm=T, ties.method="random", seed=seed, weights=weights)
    LR2.denom<-Gini.coef(Y, na.rm=T, ties.method="random", seed=seed, weights=weights)
    LR2<-as.numeric(LR2.num/LR2.denom)
    Gi.expl<-as.numeric(LR2.num)

  }

  if (ties.method == "mean"){

    GA <- GA::ga(type = "real-valued",
                 population = Lorenz.Population,
                 fitness =  function(u)Fitness_meanrank(u,as.vector(YX_mat[,1]),as.matrix(YX_mat[,-1]),pi),
                 lower = rep(-1,n.param-1), upper = rep(1,n.param-1),
                 popSize = popSize, maxiter = maxiter, run = run, monitor = FALSE,
                 parallel = parallel)

    # We need to take care of the fact that the first coefficient for theta may be positive or negative
    theta1<-c(GA@solution[1,],1-sum(abs(GA@solution[1,]))) #The theta solution if the last coeff is positive
    theta2<-c(GA@solution[1,],-(1-sum(abs(GA@solution[1,])))) #The theta solution if the last coeff is negative
    theta<-rbind(theta1,theta2)
    Y <- YX_mat[,1]
    index1<-theta1%*%t(YX_mat[,-1])
    index2<-theta2%*%t(YX_mat[,-1])
    index1_k <- sort(unique(index1))
    pi1_k <- sapply(1:length(index1_k),function(k)sum(pi[index1==index1_k[k]]))
    F1_k <- cumsum(pi1_k) - 0.5*pi1_k
    F1_i <- sapply(1:length(index1),function(i)sum(F1_k[index1_k==index1[i]])) # Ensures that sum(F_i*pi) = 0.5
    index2_k <- sort(unique(index2))
    pi2_k <- sapply(1:length(index2_k),function(k)sum(pi[index2==index2_k[k]]))
    F2_k <- cumsum(pi2_k) - 0.5*pi2_k
    F2_i <- sapply(1:length(index2),function(i)sum(F2_k[index2_k==index2[i]])) # Ensures that sum(F_i*pi) = 0.5
    theta.argmax<-theta[which.max(c((pi*Y)%*%F1_i,(pi*Y)%*%F2_i)),]

    # We compute the Lorenz-Rsquared
    Index.sol<-theta.argmax%*%t(YX_mat[,-1])
    Y<-YX_mat[,1]

    LR2.num<-Gini.coef(Y, x=Index.sol, na.rm=T, ties.method="mean", seed=seed, weights=weights)
    LR2.denom<-Gini.coef(Y, na.rm=T, ties.method="mean", seed=seed, weights=weights)
    LR2<-as.numeric(LR2.num/LR2.denom)
    Gi.expl<-as.numeric(LR2.num)

  }

  # We spit out the actual solution, with the right sign for the last coeff as well as the Lorenz-Rsquared
  result <- list(sol=theta.argmax,LR2=LR2,Gi.expl=Gi.expl,niter=length(GA@summary[,1]),fit=GA@fitnessValue)

  return(result)

}
