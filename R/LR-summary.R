#' Summary for the Lorenz Regression
#'
#' \code{summary.LR} provides a summary for an object of class \code{LR}.
#'
#' @param LR Output of a call to \code{\link{Lorenz.Reg}}, where \code{penalty=="none"}.
#'
#' @return A summary displaying the explained Gini coefficient, Lorenz-\eqn{R^2} and a table gathering the estimated coefficients.
#'
#' @seealso \code{\link{Lorenz.Reg}}
#'
#' @examples
#' data(Data.Incomes)
#' NPLR <- Lorenz.Reg(Income ~ ., data = Data.Incomes, penalty = "none")
#' summary.LR(NPLR)
#'
#' @import knitr
#'
#' @method summary LR
#' @export

summary.LR <- function(LR){

  theta.mat <- as.matrix(LR$theta)
  colnames(theta.mat) <- c("estimate")
  theta.table <- knitr::kable(theta.mat)

  cat(paste0("The explained Gini coefficient is of ",round(LR$Gi.expl,5)),
            "",
            paste0("The Lorenz-R2 is of ",round(LR$LR2,5)),
            "",
            "Estimated coefficients",
            theta.table,
            sep = "\n")

}
