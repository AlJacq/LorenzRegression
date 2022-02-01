#' Trace plot of a Penalized Lorenz Regression
#'
#' \code{PLR.traceplot} displays the trace plot of a Penalized Lorenz Regression, i.e. it displays the evolution of the estimated coefficient attached to each covariate as a function of the regularization parameter.
#'
#' @param PLR Output of a call to \code{\link{Lorenz.Reg}}, where \code{penalty!="none"}.
#'
#' @return A plot where the horizontal axis is -log(lambda), lambda being the value of the regularization parameter, and where the vertical axis gives the size of the coefficient attached to each covariate.
#'
#' @seealso \code{\link{Lorenz.Reg}}
#'
#' @examples
#' data(Data.Incomes)
#' PLR <- Lorenz.Reg(Income ~ ., data = Data.Incomes, penalty = "SCAD", lambda.choice = "BIC", eps = 0.005, nfolds = 5)
#' PLR.traceplot(PLR)
#'
#' @import ggplot2
#'
#' @export

PLR.traceplot <- function(PLR, type.var=NULL){

  Path <- PLR$path
  names.var <- rownames(Path)[-c(1:3)]
  n.iter <- ncol(Path)
  Plot.data <- data.frame(
    "Variable" = rep(names.var,n.iter),
    "theta" = as.vector(Path[-c(1:3),]),
    "minloglambda" = rep(-log(Path[1,]),each=length(names.var))
  )
  ggplot2::ggplot(Plot.data) +
    aes(x = minloglambda, y = theta, colour = Variable) +
    geom_line(size = 1L) +
    scale_color_hue() +
    labs(x = expression(paste("-log(", symbol(lambda), ")",sep="")),
         y = expression(symbol(theta)[k])) +
    theme_minimal()

}
