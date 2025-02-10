#' Penalized Lorenz Regression Fit Function
#'
#' \code{PLR.fit} fits a penalized Lorenz regression model using either the LASSO or SCAD penalty.
#' It serves as an internal wrapper that applies the fit function over a grid of tuning parameter values.
#'
#' @param y A numeric vector representing the response variable.
#' @param x A numeric matrix of covariates.
#' @param weights An optional numeric vector of sample weights. Default is \code{NULL}.
#' @param penalty A character string specifying the penalty type. Possible values are \code{"LASSO"} and \code{"SCAD"}.
#' @param grid.arg A character string specifying the tuning parameter for which a grid is constructed.
#' @param grid.value A numeric vector specifying the grid values for \code{grid.arg}. If \code{NULL}, no grid is constructed.
#' @param lambda.list An optional list specifying penalty values (\eqn{\lambda}) to be used for each grid value.
#' @param ... Additional arguments passed to \code{\link{Lorenz.FABS}} or \code{\link{Lorenz.SCADFABS}}, depending on the penalty type.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{path}}{A list of matrices, where each element corresponds to a grid value. Each matrix contains lambda values, Lorenz-\eqn{R^2}, explained Gini coefficients, BIC scores, and estimated coefficients.}
#'   \item{\code{grid.idx}}{The index of the optimal grid parameter selected by the BIC criterion.}
#'   \item{\code{lambda.idx}}{The index of the optimal \eqn{\lambda} selected by the BIC criterion.}
#'   \item{\code{grid.value}}{The grid values used for \code{grid.arg}.}
#'   \item{\code{lambda.list}}{A list of \eqn{\lambda} values along the solution paths.}
#'   \item{\code{grid.arg}}{The tuning parameter for which the grid was constructed.}
#' }
#'
#' @details
#' The function applies either \code{\link{Lorenz.FABS}} (for LASSO) or \code{\link{Lorenz.SCADFABS}} (for SCAD) for each grid value.
#' The best model is selected based on the BIC score.
#'
#' @seealso \code{\link{Lorenz.FABS}}, \code{\link{Lorenz.SCADFABS}}, \code{\link{Lorenz.boot}}, \code{\link{Lorenz.Reg}}
#'
#' @examples
#' data(Data.Incomes)
#' y <- Data.Incomes$Income
#' x <- as.matrix(Data.Incomes[,-c(1,2)])
#' PLR.fit(y, x, penalty = "SCAD", grid.arg = "eps", grid.value = c(0.2,0.5), lambda.list = NULL)
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
  arg.list <- lapply(1:lth.path,function(z)list(y = y, x = x, weights = weights))
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
  Path_BIC <- lapply(1:lth.path,function(i)PLR.BIC(y, x, LR[[i]]$theta, weights = weights))
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
  return.list$lambda.list <- lapply(Path,function(x)x["lambda",])
  return.list$grid.arg <- grid.arg

  return(return.list)
}
