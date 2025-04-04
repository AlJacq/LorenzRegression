#' Plots for the Lorenz regression
#'
#' \code{autoplot} generates a plot for an object of class \code{"LR"} and returns it as a \code{ggplot} object.
#' The \code{plot} method is a wrapper around \code{autoplot} that directly displays the plot,
#' providing a more familiar interface for users accustomed to base R plotting.
#'
#' @aliases plot.LR autoplot.LR_boot plot.LR_boot
#' @param x An object of class \code{"LR"}.
#' @param object An object of class \code{"LR"}.
#' @param band.level Confidence level for the bootstrap confidence intervals.
#' @param ... Additional arguments passed to \code{\link{Lorenz.graphs}}.
#'
#' @return \code{autoplot} returns a \code{ggplot} object representing the Lorenz curve of the response and the concentration curve of the response with respect to the estimated index.
#' If \code{object} inherits from \code{"LR_boot"} and \code{LC_store} was set to \code{TRUE} in \code{\link{Lorenz.boot}}, pointwise confidence intervals for the concentration curve are added. Their confidence level is set via the argument \code{band.level}.
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

autoplot.LR <- function(object, band.level = 0.95, ...){

  if (is.null(object$theta)) stop("No plots are available for an empty model.")

  formula <- update.formula(object, . ~ index)
  data <- data.frame(object$y,predict.LR(object))
  names(data) <- all.vars(formula)

  g <- Lorenz.graphs(formula, data, weights = object$weights, ...)
  g <- g + ggtitle("Observed and explained inequality")

  return(g)

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
autoplot.LR_boot <- function(object, band.level = 0.95, ...){

  g <- NextMethod("autoplot")

  if(object$store_LC){

    LC_start <- ncol(object$x)+2 # LC ordinates are stored after Gi.expl, LR2 and theta vector
    LC_lth <- 100
    LC_ordinates <- object$boot_out$t[,LC_start + 1:LC_lth]

    g <- Lorenz.bands(g, LC_ordinates, level = band.level, ...)

  }

  return(g)

}

#' @method plot LR_boot
#' @export
plot.LR_boot <- function(x, ...) {
  print(autoplot(x, ...))
}

