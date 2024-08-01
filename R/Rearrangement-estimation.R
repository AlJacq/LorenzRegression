#' Estimates a monotonic regression curve via Chernozhukov et al (2009)
#'
#' \code{Rearrangement.estimation} estimates the increasing link function of a single index model via the methodology proposed in Chernozhukov et al (2009).
#'
#' A first estimator of the link function, neglecting the assumption of monotonicity, is obtained with function \code{\link[locpol]{locpol}} from the \emph{locpol} package.
#' The final estimator is obtained through the rearrangement operation explained in Chernozhukov et al (2009). This operation is carried out with function \code{\link[Rearrangement]{rearrangement}} from package \emph{Rearrangement}.
#'
#' @param Y The response variable.
#' @param Index The estimated index. The user may obtain it using function \code{\link{Lorenz.Reg}}.
#' @param t A vector of points over which the link function \eqn{H(.)} should be estimated. Default is the estimated index.
#' @param weights vector of sample weights. By default, each observation is given the same weight.
#' @param degree.pol degree of the polynomial used in the local polynomial regression. Default value is 1.
#'
#' @return A list with the following components
#' \describe{
#'     \item{\code{t}}{the points over which the estimation has been undertaken.}
#'     \item{\code{H}}{the estimated link function evaluated at \emph{t}.}
#' }
#'
#' @seealso \code{\link{Lorenz.Reg}}, \code{\link[locpol]{locpol}}, \code{\link[Rearrangement]{rearrangement}}
#'
#' @section References:
#' Chernozhukov, V., I. Fernández-Val, and A. Galichon (2009). Improving Point and Interval Estimators of Monotone Functions by Rearrangement. \emph{Biometrika 96 (3)}. 559–75.
#'
#' @examples
#' data(Data.Incomes)
#' PLR <- Lorenz.Reg(Income ~ ., data = Data.Incomes,
#'                   penalty = "SCAD", eps = 0.01)
#' Y <- PLR$y
#' Index <- PLR$index
#' Rearrangement.estimation(Y = Y, Index = Index)
#'
#' @import Rearrangement
#' @import locpol
#'
#' @export

Rearrangement.estimation<-function(Y, Index, t=Index, weights=NULL, degree.pol=1){

  #################
  #ESTIMATION OF H#
  #################

  if(is.null(weights)){
    n <- length(Y)
    weights <- rep(1,n)
  }

  # Original estimator ----

  Index.std <- Index/max(Index) # It seems to be better for numerical stability
  t.std <- t/max(Index)

  Data.locpol <- data.frame(Y = Y, Index.std = Index.std)

  H.LocPolyn <- locpol::locpol(Y ~ Index.std, Data.locpol, xeval = t.std, deg = degree.pol, weig = weights)$lpFit$Y # Watch out ... the response is spitted out as if it was computed on sort(t.std).

  o <- rank(t.std,ties.method="first")
  H.LocPolyn <- H.LocPolyn[o]

  # Transformed estimator ----

  H.Rear <- Rearrangement::rearrangement(as.data.frame(t),H.LocPolyn)

  result <- list(t=t,H=H.Rear)

  return(result)
}
