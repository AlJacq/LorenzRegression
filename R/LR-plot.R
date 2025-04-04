#' Plots for the Lorenz regression
#'
#' \code{autoplot} generates a plot for an object of class \code{"LR"} and returns it as a \code{ggplot} object.
#' The \code{plot} method is a wrapper around \code{autoplot} that directly displays the plot,
#' providing a more familiar interface for users accustomed to base R plotting.
#'
#' @aliases plot.LR autoplot.LR_boot plot.LR_boot
#' @param x An object of class \code{"LR"}.
#' @param object An object of class \code{"LR"}.
#' @param type A character string indicating the type of plot. Possible values are \code{"explained"} and \code{"residuals"}.
#' \itemize{
#' \item If \code{"explained"} is selected, the graph displays the Lorenz curve of the response and concentration curve of the response with respect to the estimated index.
#' If \code{object} inherits from \code{"PLR_boot"} and \code{LC_store} was set to \code{TRUE} in \code{\link{Lorenz.boot}}, pointwise confidence intervals for the concentration curve are added. Their confidence level is set via the argument \code{band.level}.
#' \item If \code{"residuals"} is selected, the graph displays a scatterplot of residuals with respect to the estimated index.
#' Obtaining residuals entail to estimate the link function of the single-index. This is performed via the function \code{\link{Rearrangement.estimation}}, as explained in \code{\link{predict.LR}}.
#' }
#' @param band.level Confidence level for the bootstrap confidence intervals.
#' @param ... Additional arguments passed to \code{\link{Lorenz.graphs}}.
#'
#' @return \code{autoplot} returns a \code{ggplot} object representing the desired graph. \code{plot} directly displays this plot.
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

autoplot.LR <- function(object, type = c("explained","residuals"), band.level = 0.95, ...){

  type <- match.arg(type)

  if (is.null(object$theta)) stop("No plots are available for an empty model.")

  if (type == "explained"){

    formula <- update.formula(object, . ~ index)
    data <- data.frame(object$y,fitted.LR(object))
    names(data) <- all.vars(formula)

    g <- Lorenz.graphs(formula, data, weights = object$weights, ...)
    g <- g + ggtitle("Observed and explained inequality")

  }else if (type == "residuals"){

    data <- data.frame(index = fitted.LR(object),
                       resid = residuals.LR(object))

    # Scatterplot
    # Should I include sample weights in the scatterplot ?

  }

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
autoplot.LR_boot <- function(object, type = c("explained","residuals"), band.level = 0.95, ...){

  g <- NextMethod("autoplot")

  if (type == "explained"){

    if(object$store_LC){

      LC_start <- ncol(object$x)+2 # LC ordinates are stored after Gi.expl, LR2 and theta vector
      LC_lth <- 100
      LC_ordinates <- object$boot_out$t[,LC_start + 1:LC_lth]

      g <- Lorenz.bands(g, LC_ordinates, level = band.level, ...)

    }

  }

  return(g)

}

#' @method plot LR_boot
#' @export
plot.LR_boot <- function(x, ...) {
  print(autoplot(x, ...))
}

