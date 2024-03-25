#' Summary for the Penalized Lorenz Regression
#'
#' \code{summary.PLR} provides a summary for an object of class \code{PLR}.
#'
#' @param object Output of a call to \code{\link{Lorenz.Reg}}, where \code{penalty!="none"}.
#' @param renormalize whether the coefficient vector should be re-normalized to match the representation where the first category of each categorical variable is omitted. Default value is TRUE
#' @param ... Additional arguments.
#'
#' @return A summary displaying two tables: a summary of the model and the estimated coefficients.
#'
#' @seealso \code{\link{Lorenz.Reg}}
#'
#' @examples
#' data(Data.Incomes)
#' PLR <- Lorenz.Reg(Income ~ ., data = Data.Incomes, penalty = "SCAD",
#'                   sel.choice = c("BIC","CV"), h.grid = nrow(Data.Incomes)^(-1/5.5),
#'                   eps = 0.01, seed.CV = 123, nfolds = 5)
#' summary(PLR)
#'
#' @import knitr
#'
#' @method summary PLR
#' @export

summary.PLR <- function(object, renormalize=TRUE, ...){

  PLR <- object
  sum.table <- knitr::kable(PLR$summary)

  if (renormalize){
    theta <- PLR.normalize(PLR)
  }else{
    theta <- PLR$theta
  }
  theta.table <- knitr::kable(t(theta))

  cat("Summary of the model fit",
      sum.table,
      "",
      "Estimated coefficients",
      theta.table,
      sep = "\n")

}
