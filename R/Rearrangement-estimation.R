#' Estimates a monotonic regression curve via Chernozhukov et al (2009)
#'
#' \code{Rearrangement.estimation} estimates the increasing link function of a single index model via the methodology proposed in Chernozhukov et al (2009).
#'
#' A first estimator of the link function, neglecting the assumption of monotonicity, is obtained with function \code{\link[NonpModelCheck]{localpoly.reg}} from the \emph{NonpModelCheck} package.
#' The final estimator is obtained through the rearrangement operation explained in Chernozhukov et al (2009). This operation is carried out with function \code{\link[Rearrangement]{rearrangement}} from package \code{Rearrangement}.
#'
#' @param Y The response variable.
#' @param Index The estimated index. The user may obtain it using function \code{\link{Lorenz.Reg}}.
#' @param t A vector of points over which the link function \eqn{H(.)} should be estimated. Default is the estimated index.
#' @param band the method used to estimate the bandwidth in the local polynomial regression
#' @param degree.pol degree of the polynomial used in the local polynomial regression
#'
#' @return A list with the following components
#' \describe{
#'     \item{\code{t}}{the points over which the estimation has been undertaken.}
#'     \item{\code{H}}{the estimated link function evaluated at \emph{t}.}
#' }
#'
#' @seealso \code{\link{Lorenz.Reg}}, \code{\link[NonpModelCheck]{localpoly.reg}}, \code{\link[Rearrangement]{rearrangement}}
#'
#' @section References:
#' Chernozhukov, V., I. Fernández-Val, and A. Galichon. 2009. “Improving Point and Interval Estimators of Monotone Functions by Rearrangement.” Biometrika 96 (3): 559–75.
#'
#' @examples
#' data(Data.Incomes)
#' LR <- Lorenz.Reg(Income ~ Age + Work.Hours, data = Data.Incomes)
#' Y <- LR$Fit[,1]
#' Index <- LR$Fit[,2]
#' Rearrangement.estimation(Y = Y, Index = Index)
#'
#' @import Rearrangement
#' @import NonpModelCheck
#'
#' @export

Rearrangement.estimation<-function(Y, Index, t=Index, band="CV", degree.pol=1){

  #################
  #ESTIMATION OF H#
  #################

  # Original estimator ----
  H.LocPolyn<-NonpModelCheck::localpoly.reg(X=Index,Y=Y,points=t,
                            bandwidth=band,degree.pol=degree.pol,kernel.type="gaussian")$predicted

  # Transformed estimator ----

  H.Rear<-Rearrangement::rearrangement(as.data.frame(t),H.LocPolyn)

  result<-list(t=t,H=H.Rear)

  return(result)
}
