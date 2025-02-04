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
#' @param suggestions Initial guesses used in the genetic algorithm. The default value is \code{NULL}, meaning no suggestions are passed.
#' Other possible values are a numeric matrix with at most \code{popSize} rows and \code{ncol(x)} columns, or a character string "OLS".
#' In the latter case, \code{0.5*popSize} suggestions are created as random perturbations of the OLS solutions.
#' @param ties.method What method should be used to break the ties in optimization program. Possible values are "random" (default value) or "mean". If "random" is selected, the ties are broken by further ranking in terms of a uniformly distributed random variable. If "mean" is selected, the average rank method is used.
#' @param ties.Gini what method should be used to break the ties in the computation of the Gini coefficient at the end of the algorithm. Possible values and default choice are the same as above.
#' @param seed.random An optional seed for generating the vector of uniform random variables used to break ties in the genetic algorithm. Defaults to \code{NULL}, which means no specific seed is set.
#' @param seed.Gini An optional seed for generating the vector of uniform random variables used to break ties in the computation of the Gini coefficient. Defaults to \code{NULL}, meaning no specific seed is applied.
#' @param seed.GA An optional seed for \code{\link[GA]{ga}}, used during the fitting of the genetic algorithm. Defaults to \code{NULL}, implying that no specific seed is set.
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
#' @details The parameters \code{seed.random}, \code{seed.Gini}, and \code{seed.GA} allow for local seed setting to control randomness in specific parts of the function.
#' Each seed is applied to the respective part of the computation, and the seed is reverted to its previous state after the operation.
#' This ensures that the seed settings do not interfere with the global random state or other parts of the code.
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
#' @export

PLR.fit <- function(y, x, weights = NULL, penalty, grid.arg, grid.value, lambda.list, ...){

  # 1. Model fitting ----

  if(is.null(grid.value)){
    lth.path <- 1
  }else{
    lth.path <- length(grid.value)
  }
  fun <- switch(penalty,
                "LASSO" = Lorenz.FABS,
                "SCAD" = Lorenz.SCADFABS)
  arg.list <- lapply(1:lth.path,function(z)list(y = y, x = x, weights = w))
  for (i in 1:lth.path){
    if(!is.null(lambda.list)) arg.list[[i]]$lambda <- lambda.list[[i]]
    if(!is.null(grid.value)) arg.list[[i]][grid.arg] <- grid.value[i]
  }
  dots <- list(...)
  call.list <- lapply(1:lth.path,function(i)c(arg.list[[i]],dots))
  LR <- lapply(1:lth.path,function(i)do.call(fun,call.list[[i]]))

  # 2. Return ----

  return.list <- list()

  # Construction of the path > Number of selected vars
  n_selected <- lapply(1:lth.path,function(i)apply(LR[[i]]$theta,2,function(x)sum(abs(x) > 10^(-10))))
  # Construction of the path > Main objects
  Path <- lapply(1:lth.path,function(i)rbind(LR[[i]]$lambda, LR[[i]]$LR2, LR[[i]]$Gi.expl, n_selected[[i]]))
  for(i in 1:lth.path) rownames(Path[[i]]) <- c("lambda","Lorenz-R2","Explained Gini", "Number of nonzeroes")
  # Construction of the path > BIC score
  Path_BIC <- lapply(1:lth.path,function(i)PLR.BIC(y, x, LR[[i]]$theta, weights = w))
  best.BIC <- lapply(1:lth.path,function(i)Path_BIC[[i]]$best)
  val.BIC <- lapply(1:lth.path,function(i)Path_BIC[[i]]$val)
  for (i in 1:lth.path){
    Path[[i]] <- rbind(Path[[i]], val.BIC[[i]])
    rownames(Path[[i]])[nrow(Path[[i]])] <- "BIC score"
  }
  # Construction of the path > theta's
  for (i in 1:lth.path){
    lth <- nrow(Path[[i]])
    Path[[i]] <- rbind(Path[[i]], LR[[i]]$theta)
    rownames(Path[[i]])[(lth+1):nrow(Path[[i]])] <- colnames(x)
  }
  return.list$path <- Path
  # Optimum grid params for BIC
  # grid refers either to h or to SCAD.nfwd
  grid.idx <- which.max(sapply(1:lth.path,function(i)max(val.BIC[[i]])))
  lambda.idx <- best.BIC[[grid.idx]]
  names(grid.idx) <- names(lambda.idx) <- "BIC"
  return.list$grid.idx <- grid.idx
  return.list$lambda.idx <- lambda.idx
  return.list$grid.value <- grid.value

  return(return.list)
}
