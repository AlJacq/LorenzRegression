#' Plots for the Lorenz Regression
#'
#' \code{autoplot.LR} generates a plot for an object of class \code{"LR"} and returns it as a \code{ggplot} object.
#' The \code{plot.LR} method is a wrapper around \code{autoplot.LR} that directly displays the plot,
#' providing a more familiar interface for users accustomed to base R plotting.
#'
#' @aliases plot.LR
#' @param x An object of class \code{"LR"}.
#' @param ... Additional arguments passed to \code{\link{Lorenz.graphs}}.
#'
#' @return \code{autoplot.LR} returns a \code{ggplot} object representing the Lorenz curve of the response and the concentration curve of the response with respect to the estimated index. \code{plot.LR} directly displays this plot.
#'
#' @seealso \code{\link{Lorenz.Reg}}
#'
#' @examples
#' ## For examples see example(Lorenz.Reg)
#'
#' @importFrom ggplot2 ggtitle autoplot
#' @importFrom stats update.formula
#'
#' @method autoplot LR
#' @export

autoplot.LR <- function(x, ...){

  if (is.null(x$theta)) stop("No plots are available for an empty model.")

  formula <- update.formula(x, . ~ index)
  data <- data.frame(x$y,x$index)
  names(data) <- all.vars(formula)

  g <- Lorenz.graphs(formula, data, weights = x$weights, ...)
  g <- g + ggtitle("Observed and explained inequality")

  g

}

#' @importFrom graphics plot
#' @method plot LR
#' @rdname autoplot.LR
#' @export
plot.LR <- function(x, ...) {
  print(autoplot(x, ...))
}

#' @method autoplot LR_boot
#' @export
autoplot.LR_boot <- function(x, ...){
  NextMethod("autoplot")
}

#' @method plot LR_boot
#' @export
plot.LR_boot <- function(x, ...) {
  print(autoplot(x, ...))
}

