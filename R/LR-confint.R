#' Confidence intervals for the Lorenz Regression
#'
#' \code{confint.LR} provides confidence intervals for the explained Gini coefficient and Lorenz-R2 for an object of class \code{LR}.
#'
#' @param LR Output of a call to \code{\link{Lorenz.Reg}}, where \code{penalty="none"} and \code{Boot.inference=TRUE}.
#'
#' @return A list with two elements
#' \describe{
#'    \item{\code{Gi}}{a matrix where each row is a confidence interval for the explained Gini coefficient. The different rows correspond to different bootstrap methods.}
#'    \item{\code{LR2}}{a matrix where each row is a confidence interval for the Lorenz-\eqn{R^2}. The different rows correspond to different bootstrap methods.}
#' }
#'
#' @details This function is only meant to display the confidence intervals obtained from an already computed LR.
#' The level and bootstrap methods must be chosen via the arguments \code{alpha} and \code{which.CI} in the \code{Lorenz.Reg} function.
#'
#' @seealso \code{\link{Lorenz.Reg}}
#'
#' @examples
#' data(Data.Incomes)
#' NPLR <- Lorenz.Reg(Income ~ ., data = Data.Incomes, penalty = "none", seed.boot = 123, B = 40, Boot.inference = TRUE)
#' confint.LR(NPLR)
#'
#' @method confint LR
#' @export

confint.LR <- function(LR){

  if(length(grep("CI",names(LR))) > 0){
    CI.list <- list()
    CI.list$Gi <- LR$CI.Gi
    CI.list$LR2 <- LR$CI.LR2
    return(CI.list)
  }else{
    stop("The input must contain a bootstrap estimation. Consider turning the Boot.inference argument in the Lorenz.Reg function to TRUE")
  }

}
