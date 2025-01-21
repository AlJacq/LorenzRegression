#' Call to the genetic algorithm for the Lorenz regression
#'
#' \code{Lorenz.ga.call} encapsulates the call to ga for a local management of seed setting
#'
#' @param ties.method Either \code{"mean"} or \code{"random"}.
#' @param y vector of responses.
#' @param x matrix of covariates.
#' @param pi sample weights (normalized).
#' @param V vector of uniformly distributed rvs.
#' @param popSize passed to \code{ga}.
#' @param maxiter passed to \code{ga}.
#' @param run passed to \code{ga}.
#' @param parallel.GA passed to \code{ga}.
#' @param suggestions passed to \code{ga}.
#' @param seed An optional integer for setting the seed for random number generation. Default is \code{NULL}.
#'
#' @importFrom GA ga
#' @importFrom stats runif
#'
#' @return The fitted genetic algorithm
#' @keywords internal

Lorenz.ga.call <- function(ties.method, y, x, pi, V, popSize, maxiter, run, parallel.GA, suggestions, seed = NULL){

  p <- ncol(x)

  if (!is.null(seed)){
    if(exists(".Random.seed")){
      old <- .Random.seed
    }else{
      runif(1)
      old <- .Random.seed
    }
    on.exit( { .Random.seed <<- old } )
  }

  tolerance <- sqrt(.Machine$double.eps)
  if (ties.method == "random") fitness_function <- function(u) .Fitness_cpp(u,as.vector(y),as.matrix(x),V,pi,tolerance)
  if (ties.method == "mean") fitness_function <- function(u) .Fitness_meanrank(u,as.vector(y),as.matrix(x),pi,tolerance)

  GA <- GA::ga(type = "real-valued",
               population = Lorenz.Population,
               fitness =  fitness_function,
               lower = rep(-1,p-1), upper = rep(1,p-1),
               popSize = popSize, maxiter = maxiter, run = run, monitor = FALSE,
               parallel = parallel.GA, seed = seed, suggestions = suggestions)

  return(GA)
}
