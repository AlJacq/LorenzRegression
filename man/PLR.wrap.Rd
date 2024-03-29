% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PLR-wrap.R
\name{PLR.wrap}
\alias{PLR.wrap}
\title{Wrapper for the \code{\link{Lorenz.SCADFABS}} and \code{\link{Lorenz.FABS}} functions}
\usage{
PLR.wrap(
  YX_mat,
  standardize = TRUE,
  weights = NULL,
  penalty = c("SCAD", "LASSO"),
  h,
  SCAD.nfwd = NULL,
  eps = 0.005,
  gamma = 0.05,
  kernel = c("Epan", "Biweight"),
  ...
)
}
\arguments{
\item{YX_mat}{a matrix with the first column corresponding to the response vector, the remaining ones being the explanatory variables.}

\item{standardize}{Should the variables be standardized before the estimation process? Default value is TRUE.}

\item{weights}{vector of sample weights. By default, each observation is given the same weight.}

\item{penalty}{penalty used in the Penalized Lorenz Regression. Possible values are "SCAD" (default) or "LASSO".}

\item{h}{bandwidth of the kernel, determining the smoothness of the approximation of the indicator function.}

\item{SCAD.nfwd}{optional tuning parameter used if penalty="SCAD". Default value is NULL. The larger the value of this parameter, the sooner the path produced by the SCAD will differ from the path produced by the LASSO.}

\item{eps}{step size in the FABS or SCADFABS algorithm. Default value is 0.005.}

\item{gamma}{value of the Lagrange multiplier in the loss function}

\item{...}{Additional parameters corresponding to arguments passed in \code{\link{Lorenz.SCADFABS}} or \code{\link{Lorenz.FABS}} depending on the argument chosen in penalty.}
}
\value{
A list with several components:
\describe{
   \item{\code{lambda}}{vector gathering the different values of the regularization parameter}
   \item{\code{theta}}{matrix where column i provides the normalized estimated parameter vector corresponding to value lambda[i] of the regularization parameter.}
   \item{\code{LR2}}{vector where element i provides the Lorenz-\eqn{R^2} of the regression related to value lambda[i] of the regularization parameter.}
   \item{\code{Gi.expl}}{vector where element i provides the estimated explained Gini coefficient related to value lambda[i] of the regularization parameter.}
}
}
\description{
\code{PLR.wrap} standardizes the covariates, run the penalized regression and spits out the path of parameter vectors.
}
\examples{
data(Data.Incomes)
YX_mat <- Data.Incomes[,-2]
PLR.wrap(YX_mat, h = nrow(Data.Incomes)^(-1/5.5), eps = 0.005)

}
\seealso{
\code{\link{Lorenz.SCADFABS}}, \code{\link{Lorenz.FABS}}
}
