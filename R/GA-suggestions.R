#' Defines the suggestions used in the genetic algorithm
#'
#' \code{Lorenz.Suggestions} creates suggestions for the genetic algorithm used to solve the Lorenz regression.
#'
#' @param suggestions either a character string 'OLS' or a numeric matrix with at most \code{popSize} rows and \code{ncol(x)} columns.
#' @param popSize population size of the genetic algorithm
#' @param y vector of responses
#' @param x matrix of covariates (after standardization if applied)
#' @param pi vector of normalized weights
#' @param x.scale vector of standard deviations of the covariates
#' @param seed seed used in the generation of the suggestions
#'
#' @seealso \code{\link{Lorenz.GA}}
#'
#' @return A matrix with at most \code{popsize} rows and with a number of columns equal to the number of explanatory variables minus one.
#'
#' @importFrom stats lm median rnorm
#'
#' @keywords internal

Lorenz.Suggestions <- function(suggestions, popSize, y, x, pi, x.scale, seed){

  p <- ncol(x)
  sig <- 1
  prop <- 0.5

  # 1. Checks ----

  type_matrix <- type_OLS <- FALSE

  if (is.character(suggestions)) {
    if (!identical(suggestions, "OLS")) {
      stop("suggestions must be 'OLS' if it is a character string.")
    }else{
      type_OLS <- TRUE
    }
  } else if (is.matrix(suggestions)) {
    if (!is.numeric(suggestions)) {
      stop("If suggestions is a matrix, it must be numeric.")
    }
    if (ncol(suggestions) != p) {
      stop(paste0("The number of columns of the suggestion matrix must match the number of explanatory variables."))
    }
    if (nrow(suggestions) > popSize) {
      stop(paste0("The number of rows of the suggestion matrix cannot exceed the population size of the genetic algorithm. "))
    }
    type_matrix <- TRUE
  } else {
    stop("suggestions must be either a character string 'OLS' or a numeric matrix.")
  }

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

  # 2. OLS suggestions ----

  if (type_OLS){

    b <- lm(y ~ x, weights = pi)$coefficients
    b <- b[-1] # Delete intercept
    b_m <- median(abs(b))
    b <- matrix(b, nrow = floor(prop*popSize), ncol = length(b), byrow = TRUE)
    eps <- rnorm(n = floor(prop*popSize) - 1, mean = 0, sd = sig*b_m)
    mod_cols <- sample(1:ncol(b), size = floor(prop*popSize) - 1, replace = TRUE)
    for (i in 2:nrow(b)) {
      b[i, mod_cols[i]] <- eps[i]
    }
    sugg_matrix <- t(apply(b,1,function(x)x/sum(abs(x))))[,-ncol(b)]

  }

  # 3. User suggestions ----

  if (type_matrix){

    b <- t(apply(suggestions,1,function(x)x*x.scale))
    sugg_matrix <- t(apply(b,1,function(x)x/sum(abs(x))))[,-ncol(b)]

  }

  return(sugg_matrix)

}
