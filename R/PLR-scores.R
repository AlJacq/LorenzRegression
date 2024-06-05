#' Computes Gini scores for the Penalized Lorenz Regression
#'
#' \code{PLR.scores} computes the Gini scores (either OOB-scores or CV-scores) obtained for a specific validation sample and associated to a list of parameters obtained by the Penalized Lorenz Regression.
#'
#' @param y the vector of responses
#' @param x the design matrix (after data management steps, i.e. standardization and transformations of the categorical covariates into binaries)
#' @param weights vector of sample weights. By default, each observation is given the same weight.
#' @param theta.list list of matrices. Each element of the list correspond to a value of the tuning parameter. The columns of the matrices correspond to values of the penalty parameters. The rows correspond to the different covariates.
#'
#' @return A list of vectors gathering the Gini scores. Each element of the list corresponds to a value of the tuning parameter and each element of the vector corresponds to a value of the penalization parameter.

PLR.scores <- function(y, x, weights, theta.list){

  Gi.list <- lapply(1:length(theta.list),function(i)apply(theta.list[[i]],2,function(index)Gini.coef(y = y, x = x%*%index, na.rm=TRUE, ties.method = "mean", weights = weights)))

  return(Gi.list)

}
