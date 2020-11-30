#' Estimates the parameter vector in Lorenz regression using a genetic algorithm
#'
#' \code{Lorenz.GA.cpp} estimates the vector of parameters in Lorenz regression using the unit-norm normalization
#' It also returns the Lorenz-\eqn{R^2} of the regression.
#'
#' The genetic algorithm is solved using function \code{\link[GA]{ga}} from the \emph{GA} package. The fitness function is coded in Rcpp to speed up computation time.
#' When discrete covariates are introduced and ties occur in the index, they are randomly broken, as advised in Section 3 of Heuchenne and Jacquemain (2020)
#'
#' @param YX_mat A matrix with the first column corresponding to the response vector, the remaining ones being the explanatory variables.
#' @param parallel Whether parallel computing should be used in the genetic algorithm. Default value is FALSE.
#' @param popSize Size of the population of candidates in the genetic algorithm. Default value is 50.
#' @param maxiter Maximum number ot iterations in the genetic algorithm. Default value is 1500.
#' @param run Number of iterations without improvement in the best fitness necessary for the algorithm to stop. Default value is 150.
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
Lorenz.GA.cpp<-function(YX_mat,popSize=50,maxiter=1500,run=150, parallel = F){

  YX_mat<-YX_mat[order(YX_mat[,1]),]

  n.param<-length(YX_mat[1,-1])
  n<-length(YX_mat[,1])
  V<-stats::runif(n)

  # The GA in itself
  GA <- GA::ga(type = "real-valued",
           population = Lorenz.Population,
           fitness =  function(u)Fitness_cpp(u,as.vector(YX_mat[,1]),as.matrix(YX_mat[,-1]),V),
           lower = rep(-1,n.param-1), upper = rep(1,n.param-1),
           popSize = popSize, maxiter = maxiter, run = run, monitor = FALSE,
           parallel = parallel)

  # We need to take care of the fact that the first coefficient for theta may be positive or negative
  theta1<-c(GA@solution[1,],1-sum(abs(GA@solution[1,]))) #The theta solution if the last coeff is positive
  theta2<-c(GA@solution[1,],-(1-sum(abs(GA@solution[1,])))) #The theta solution if the last coeff is negative
  theta<-rbind(theta1,theta2)
  Index_1<-theta1%*%t(YX_mat[,-1])
  Y_1<-YX_mat[order(Index_1,V),1]
  Index_2<-theta2%*%t(YX_mat[,-1])
  Y_2<-YX_mat[order(Index_2,V),1]
  theta.argmax<-theta[which.max(c(Y_1%*%seq(from=1,to=n),Y_2%*%seq(from=1,to=n))),]

  # We compute the Lorenz-Rsquared
  Index.sol<-theta.argmax%*%t(YX_mat[,-1])
  Y.sol<-YX_mat[order(Index.sol,V),1]
  Y<-YX_mat[,1]
  LR2.num<-2*(Y.sol%*%seq(from=1,to=n))/n^2/mean(Y) - (n+1)/n
  LR2.denom<-2*(Y%*%seq(from=1,to=n))/n^2/mean(Y) - (n+1)/n
  LR2<-as.numeric(LR2.num/LR2.denom)
  Gi.expl<-as.numeric(LR2.num)

  # We spit out the actual solution, with the right sign for the last coeff as well as the Lorenz-Rsquared
  result <- list(sol=theta.argmax,LR2=LR2,Gi.expl=Gi.expl,niter=length(GA@summary[,1]),fit=GA@fitnessValue)

  return(result)

}
