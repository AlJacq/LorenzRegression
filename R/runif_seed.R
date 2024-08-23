#' Generates a sample of uniform random variables with a specific seed
#'
#' \code{runif_seed} generates a vector of uniform(0,1) random variables with a specific seed. The seed is only used locally.
#'
#' @param n the sample size
#' @param seed the seed to use
#'
#' @return A vector with the generated random variables
#' @keywords internal

runif_seed <- function(n, min = 0, max = 1, seed = NULL) {
  if (!is.null(seed)){
    if(exists(".Random.seed")){
      old <- .Random.seed
    }else{
      runif(1)
      old <- .Random.seed
    }
    on.exit( { .Random.seed <<- old } )
    set.seed(seed)
  }
  return(stats::runif(n, min = min, max = max))
}
