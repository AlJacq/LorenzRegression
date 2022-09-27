#' Confidence intervals for the Penalized Lorenz Regression
#'
#' \code{confint.PLR} provides confidence intervals for the explained Gini coefficient and Lorenz-R2 for an object of class \code{PLR}.
#'
#' @param PLR Output of a call to \code{\link{Lorenz.Reg}}, where \code{penalty!="none"} and \code{Boot.inference=TRUE}.
#'
#' @return A list with two elements
#' \describe{
#'    \item{\code{Gi}}{a list of confidence intervals for the explained Gini coefficient. Each element of the list correspond to a bootstrap method.}
#'    \item{\code{LR2}}{a list of confidence intervals for the Lorenz-\eqn{R^2}. Each element of the list correspond to a bootstrap method.}
#' }
#'
#' @details This function is only meant to display the confidence intervals obtained from an already computed PLR.
#' The level and bootstrap methods must be chosen via the arguments \code{alpha} and \code{which.CI} in the \code{Lorenz.Reg} function.
#'
#' @seealso \code{\link{Lorenz.Reg}}
#'
#' @examples
#' data(Data.Incomes)
#' PLR <- Lorenz.Reg(Income ~ ., data = Data.Incomes, h.grid = nrow(Data.Incomes)^(-1/5.5), penalty = "SCAD", eps = 0.01, seed.boot = 123, B = 50, Boot.inference = TRUE)
#' confint.PLR(PLR)
#'
#' @method confint PLR
#' @export

confint.PLR <- function(PLR){

  if(length(grep("CI",names(PLR))) > 0){
    CI.list <- list()
    CI.list$Gi <- PLR$CI.Gi
    CI.list$LR2 <- PLR$CI.LR2
    return(CI.list)
  }else{
    stop("The input must contain a bootstrap estimation. Consider turning the Boot.inference argument in the Lorenz.Reg function to TRUE")
  }

}
