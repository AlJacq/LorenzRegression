% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PLR-CV.R
\name{PLR.CV}
\alias{PLR.CV}
\title{Determines the regularization parameter (lambda) in a PLR via cross-validation}
\usage{
PLR.CV(
  formula,
  data,
  penalty = "SCAD",
  h,
  SCAD.nfwd = NULL,
  PLR.est = NULL,
  standardize = TRUE,
  weights = NULL,
  eps,
  nfolds = 10,
  foldID = NULL,
  seed.CV = NULL,
  parallel = FALSE,
  ...
)
}
\arguments{
\item{formula}{A formula object of the form \emph{response} ~ \emph{other_variables}.}

\item{data}{A data frame containing the variables displayed in the formula.}

\item{penalty}{penalty used in the Penalized Lorenz Regression. Possible values are "SCAD" (default) or "LASSO".}

\item{h}{bandwidth of the kernel, determining the smoothness of the approximation of the indicator function.}

\item{SCAD.nfwd}{optional tuning parameter used if penalty="SCAD". Default value is NULL. The larger the value of this parameter, the sooner the path produced by the SCAD will differ from the path produced by the LASSO.}

\item{PLR.est}{Output of a call to \code{\link{PLR.wrap}} corresponding to the estimation of the Penalized Lorenz Regression on the full sample. Default value is NULL in which case the estimation on the full sample is run internally.}

\item{standardize}{Should the variables be standardized before the estimation process? Default value is TRUE.}

\item{weights}{vector of sample weights. By default, each observation is given the same weight.}

\item{eps}{Step size in the FABS or SCADFABS algorithm. Default value is 0.005.}

\item{nfolds}{Number of folds. Default value is 10.}

\item{foldID}{vector taking value from 1 to nfolds specifying the fold index of each observation. Default value is NULL in which case the folds are defined internally.}

\item{seed.CV}{Should a specific seed be used in the definition of the folds. Default value is NULL in which case no seed is imposed.}

\item{parallel}{Whether parallel computing should be used to distribute the \code{nfolds} computations on different CPUs. Either a logical value determining whether parallel computing is used (TRUE) or not (FALSE, the default value). Or a numerical value determining the number of cores to use.}

\item{...}{Additional parameters corresponding to arguments passed in \code{\link{Lorenz.SCADFABS}} or \code{\link{Lorenz.FABS}} depending on the argument chosen in penalty.}
}
\value{
A list with two components
\describe{
   \item{\code{val}}{vector indicating the CV-score for each value of lambda.}
   \item{\code{best}}{index where the optimum is attained.}
}
}
\description{
\code{PLR.CV} undertakes k-fold cross-validation for a Penalized Lorenz Regression. It returns the CV-score associated to each value of the regularization parameter and the index of the optimum.
}
\section{References}{

Jacquemain, A., C. Heuchenne, and E. Pircalabelu (2024). A penalised bootstrap estimation procedure for the explained Gini coefficient. \emph{Electronic Journal of Statistics 18(1) 247-300}.
}

\examples{
YX_mat <- Data.Incomes[,-2]
PLR <- PLR.wrap(YX_mat, h = nrow(YX_mat)^(-1/5.5), eps=0.01)
PLR.CV(Income ~ ., Data.Incomes, PLR.est = PLR,
       h = nrow(Data.Incomes)^(-1/5.5), eps = 0.01, nfolds = 5)

}
\seealso{
\code{\link{Lorenz.Reg}}, \code{\link{PLR.wrap}}, \code{\link{Lorenz.FABS}}, \code{\link{Lorenz.SCADFABS}}
}
