#' Retrieve a measure of explained inequality from a model
#'
#' This generic function extracts a measure of explained inequality, such as the explained Gini coefficient or the Lorenz-R2, from a fitted model object.
#'
#' @param object An object for which the inequality metrics should be extracted.
#' @param type Character string specifying the type of inequality metric to retrieve. Options are \code{"Gini.explained"} for the explained Gini coefficient or \code{"Lorenz-R2"} for the Lorenz-\eqn{R^2}.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return The requested inequality metric.
#'
#' @seealso \code{\link{Lorenz.Reg}}
#'
#' @examples
#' ## For examples see example(Lorenz.Reg)
#'
#' @export

ineqExplained <- function(object, type = c("Gini.explained", "Lorenz-R2"), ...) {
  UseMethod("ineqExplained")
}
