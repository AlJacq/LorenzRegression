#' Estimates a monotonic regression curve via Chernozhukov et al (2009)
#'
#' \code{Rearrangement.estimation} estimates the increasing link function of a single index model via the methodology proposed in Chernozhukov et al (2009).
#'
#' A first estimator of the link function, neglecting the assumption of monotonicity, is obtained using the procedure chosen via \code{method}.
#' The final estimator is obtained through the rearrangement operation explained in Chernozhukov et al (2009). This operation is carried out with function \code{\link[Rearrangement]{rearrangement}} from package \emph{Rearrangement}.
#'
#' @param y The response variable.
#' @param index The estimated index. The user may obtain it using function \code{\link{Lorenz.Reg}}.
#' @param t A vector of points over which the link function \eqn{H(.)} should be estimated. Default is the estimated index.
#' @param weights A vector of sample weights. By default, each observation is given the same weight.
#' @param method Either a character string specifying a smoothing method (e.g., \code{"loess"}, the default), or a list containing two elements:
#' \itemize{
#'   \item \code{fit_fun}: a function used to fit the model, typically taking arguments \code{y}, \code{x}, and optionally \code{weights}, and returning a fitted object.
#'   \item \code{predict_fun}: a function to generate predictions from the fitted object. It must accept the fitted object and an argument \code{newdata}.
#' }
#' The specification of a custom method is illustrated in the Examples section.
#' If a character string is provided, a \code{predict} method must exist for that function (e.g., \code{\link[stats]{predict.loess}}). If \code{weights} are provided but unsupported by the fit function, a warning is issued and the weights are ignored.
#'
#' @param ... Additional arguments passed to the fit function defined by \code{method}.
#'
#' @return A list with the following components
#' \describe{
#'     \item{\code{t}}{the points over which the estimation has been undertaken.}
#'     \item{\code{H}}{the estimated link function evaluated at \emph{t}.}
#' }
#'
#' @seealso \code{\link{Lorenz.Reg}}, \code{\link[Rearrangement]{rearrangement}}
#'
#' @section References:
#' Chernozhukov, V., I. Fernández-Val, and A. Galichon (2009). Improving Point and Interval Estimators of Monotone Functions by Rearrangement. \emph{Biometrika 96 (3)}. 559–75.
#'
#' @examples
#' data(Data.Incomes)
#' PLR <- Lorenz.Reg(Income ~ ., data = Data.Incomes,
#'                   penalty = "SCAD", eps = 0.01)
#' y <- PLR$y
#' index <- predict(PLR)
#' # Default method where the first step is obtained with loess()
#' Rearrangement.estimation(y = y, index = index, method = "loess")
#' # Custom method, where the first step is obtained with ksmooth()
#' # ksmooth() lacks from a separation between fitting and predicting interfaces
#' ksmooth_method <- list(
#'   fit_fun = function(y, x, ...) {
#'     list(y = y, x = x, args = list(...))
#'   },
#'   predict_fun = function(fit, newdata) {
#'     if(missing(newdata)){
#'       x.points <- fit$x
#'     } else {
#'       x.points <- newdata[,1]
#'     }
#'     o <- order(order(x.points))
#'     yhat <- do.call(ksmooth, c(list(x = fit$x, y = fit$y, x.points = x.points), fit$args))$y
#'     yhat[o]
#'   }
#' )
#' Rearrangement.estimation(y = y, index = index, method = ksmooth_method, bandwidth = 0.1)
#'
#' @importFrom Rearrangement rearrangement
#' @importFrom stats loess
#'
#' @export

Rearrangement.estimation<-function(y, index, t = index, weights = NULL, method = "loess", ...){

  # Handling method argument
  if (is.character(method)) {
    fit_fun <- get(method, mode = "function")
    predict_fun <- utils::getS3method("predict", method, optional = TRUE)
    if (is.null(predict_fun)) {
      stop(sprintf("No predict.%s method found. Either use another method or provide a custom predict function.", method))
    }
  }else if (is.list(method) && all(c("fit_fun", "predict_fun") %in% names(method))) {
    fit_fun <- method$fit_fun
    predict_fun <- method$predict_fun
  } else {
    stop("Argument 'method' must be either a character string or a list with elements 'fit_fun' and 'predict_fun'")
  }
  # Handling weights argument
  fit_args <- list(...)
  fit_formals <- names(formals(fit_fun))
  if (!is.null(weights)) {
    if("weights" %in% fit_formals){
      fit_args$weights <- weights
    }else{
      warning("Weights are supplied but the proposed method does not have a 'weights' argument. The weights are therefore unused")
    }
  }

  # Original estimator ----

  # Fitting
  if("formula" %in% fit_formals){
    fit <- do.call(fit_fun, c(list(formula = y ~ index), fit_args))
  }else if (all(c("y", "x") %in% fit_formals)){
    fit <- do.call(fit_fun, c(list(y = y, x = index), fit_args))
  }else{
    stop("The supplied method must support either a formula/data interface or y/x interface.")
  }
  # Predicting
  predict_formals <- names(formals(predict_fun))
  if("newdata" %in% predict_formals){
    H.orig <- do.call(predict_fun, list(fit, newdata = data.frame(index = t)))
  }else{
    stop("The supplied prediction method must have a 'newdata' argument.")
  }
  if (!is.numeric(H.orig)) {
    stop("The prediction method must return a numeric output.")
  }
  H.orig <- as.vector(H.orig)
  if(length(H.orig) != length(t)) stop("Prediction output length does not match newdata length.")

  # Transformed estimator ----

  H.Rear <- rearrangement(as.data.frame(t),H.orig)
  result <- list(t=t,H=H.Rear)

  return(result)
}
