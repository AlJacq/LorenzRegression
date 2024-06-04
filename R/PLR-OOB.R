#' Computes the OOB for the Penalized Lorenz Regression
#'
#' \code{PLR.OOB} computes the out-of-bag score obtained for a specific validation sample and associated to a list of parameters obtained by the Penalized Lorenz Regression.
#'
#' @param y the vector of responses
#' @param x the design matrix (after data management steps, i.e. standardization and transformations of the categorical covariates into binaries)
#' @param weights vector of sample weights related to the validation sample. By default, each observation is given the same weight.
#' @param theta.boot list of matrices. Each element of the list correspond to a value of the tuning parameter. The columns of the matrices correspond to values of the penalty parameters. The rows correspond to the different covariates.
#'
#' @return A list of vectors gathering the OOB-score. Each element of the list corresponds to a value of the tuning parameter and each element of the vector corresponds to a value of the penalization parameter.
#' @importFrom parsnip contr_one_hot

PLR.OOB <- function(y, x, weights, theta.boot){

  OOB.list <- lapply(1:length(theta.boot),function(i)apply(theta.boot[[i]],2,function(index)Gini.coef(y = y, x = x%*%index, na.rm=TRUE, ties.method = "mean", weights = weights)))

  return(OOB.list)

}
