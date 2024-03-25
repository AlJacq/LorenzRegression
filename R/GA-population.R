#' Defines the population used in the genetic algorithm
#'
#' \code{Lorenz.Population} creates the initial population of the genetic algorithm used to solve the Lorenz regression.
#'
#' Note that this population produces an initial solution ensuring a unit norm.
#'
#' @param object An object of class "\code{ga}", resulting from a call to function \code{ga}.
#'
#' @seealso \code{\link{Lorenz.GA}}
#'
#' @return A matrix of dimension \code{object@popSize} times the number of explanatory variables minus one, gathering the initial population.

Lorenz.Population<-function (object){
  min <- object@lower
  max <- object@upper
  nvars <- length(min)+1
  population<-matrix(stats::runif(object@popSize*nvars,-1,1),ncol=nvars)
  population<-t(apply(population,1,function(x)x/sum(abs(x))))[,-nvars]
  if(nvars==2){
    population<-as.matrix(population)
  }
  return(population)
}
