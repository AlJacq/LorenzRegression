#' Summary for the Penalized Lorenz Regression
#'
#' \code{summary.PLR} provides a summary for an object of class \code{PLR}.
#'
#' @param PLR Output of a call to \code{\link{Lorenz.Reg}}, where \code{penalty!="none"}.
#' @param renormalize whether the coefficient vector should be re-normalized to match the representation where the first category of each categorical variable is omitted. Default value is TRUE
#'
#' @return A summary displaying two tables: a summary of the model and the estimated coefficients.
#'
#' @seealso \code{\link{Lorenz.Reg}}
#'
#' @examples
#' data(Data.Incomes)
#' PLR <- Lorenz.Reg(Income ~ ., data = Data.Incomes, penalty = "SCAD", lambda.choice = c("BIC","CV"), eps = 0.005, seed.CV = 123, nfolds = 5)
#' summary.PLR(PLR)
#'
#' @import knitr
#'
#' @method summary PLR
#' @export

summary.PLR <- function(PLR, renormalize=TRUE){

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
