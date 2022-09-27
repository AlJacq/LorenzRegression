#' Estimated coefficients for the Lorenz Regression
#'
#' \code{coef.LR} provides the estimated coefficients for an object of class \code{LR}.
#'
#' @param LR Output of a call to \code{\link{Lorenz.Reg}}, where \code{penalty="none"}.
#'
#' @return a vector gathering the estimated coefficients
#'
#' @seealso \code{\link{Lorenz.Reg}}
#'
#' @examples
#' data(Data.Incomes)
#' NPLR <- Lorenz.Reg(Income ~ ., data = Data.Incomes, penalty = "none")
#' coef.LR(NPLR)
#'
#' @method coef LR
#' @export

coef.LR <- function(LR){

  LR$theta

}
