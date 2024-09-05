#' Summary for the Penalized Lorenz Regression
#'
#' Provides a summary for an object of class \code{"PLR"}.
#'
#' @aliases summary.PLR_boot summary.PLR_cv
#' @param object An object of class \code{"PLR"}. The object might also have S3 classes \code{"PLR_boot"} and/or \code{"PLR_cv"} (both inherit from class \code{"PLR"})
#' @param renormalize A logical value determining whether the coefficient vector should be re-normalized to match the representation where the first category of each categorical variable is omitted. Default value is TRUE
#' @param ... Additional arguments
#'
#' @return An object of class \code{"summary.PLR"}, which contains:
#' \describe{
#'    \item{\code{call}}{The matched call.}
#'    \item{\code{ineq}}{A table of explained inequality metrics. The columns display the explained Gini coefficient, the Gini coefficient of the response, and the Lorenz-R2. The first row contains the results obtained by BIC.}
#'    \item{\code{coefficients}}{A matrix with estimated coefficients, each row corresponding to a specific coefficient. The first column contains the results obtained by BIC.}
#' }
#' If the object inherits from \code{"PLR_boot"}, \code{ineq} and \code{coefficients} also include results from bootstrap, and the class \code{"summary.PLR_boot"} is added to the output.
#' Similarly, if the object inherits from \code{"PLR_cv"}, \code{ineq} and \code{coefficients} also include results from cross-validation, and the class \code{"summary.PLR_cv"} is added to the output.
#'
#' @seealso \code{\link{Lorenz.Reg}}, \code{\link{Lorenz.boot}}, \code{\link{PLR.CV}}
#'
#' @examples
#' ## For examples see example(Lorenz.Reg), example(Lorenz.boot) and example(PLR.CV)
#'
#' @method summary PLR
#' @export

summary.PLR <- function(object, renormalize=TRUE, ...){

  ans <- list()
  ans$call <- object$call
  ans$ineq <- matrix(c(ineqExplained.PLR(object, type = "Gini.explained", pars.idx = "BIC"),
                       ineqExplained.PLR(object, type = "Gini.explained", pars.idx = "BIC")/ineqExplained.PLR(object, type = "Lorenz-R2", pars.idx = "BIC"),
                       ineqExplained.PLR(object, type = "Lorenz-R2", pars.idx = "BIC")),
                     nrow = 1, ncol = 3,
                     dimnames = list("BIC", c("Explained","Total","Lorenz-R2")))
  ans$coefficients <- as.matrix(coef.PLR(object, renormalize = renormalize, pars.idx = "BIC"))
  colnames(ans$coefficients) <- ""

  class(ans) <- "summary.PLR"

  return(ans)

}

#' @method summary PLR_boot
#' @export

summary.PLR_boot <- function(object, renormalize=TRUE, ...){

  pars.boot <- c(object$grid.idx["Boot"],object$lambda.idx["Boot"])
  ans <- NextMethod("summary")
  ans$ineq <- rbind(ans$ineq,
                    c(ineqExplained_PLR(object, type = "Gini.explained", pars.idx = pars.boot),
                      ineqExplained_PLR(object, type = "Gini.explained", pars.idx = pars.boot)/ineqExplained_PLR(object, type = "Lorenz-R2", pars.idx = pars.boot),
                      ineqExplained_PLR(object, type = "Lorenz-R2", pars.idx = pars.boot)))
  ans$coefficients <- cbind(ans$coefficients,
                            coef_PLR(object, renormalize = renormalize, pars.idx = pars.boot))

  if(ncol(ans$coefficients)==2){
    rownames(ans$ineq) <- c("BIC","Boot")
    colnames(ans$coefficients) <- c("BIC","Boot")
  }else{
    rownames(ans$ineq)[nrow(ans$ineq)] <- "Boot"
    colnames(ans$coefficients)[ncol(ans$coefficients)] <- "Boot"
  }

  class(ans) <- c("summary.PLR_boot",class(ans))

  return(ans)

}

#' @method summary PLR_cv
#' @export

summary.PLR_cv <- function(object, renormalize=TRUE, ...){

  pars.cv <- c(object$grid.idx["CV"],object$lambda.idx["CV"])
  ans <- NextMethod("summary")
  ans$ineq <- rbind(ans$ineq,
                    c(ineqExplained_PLR(object, type = "Gini.explained", pars.idx = pars.cv),
                      ineqExplained_PLR(object, type = "Gini.explained", pars.idx = pars.cv)/ineqExplained_PLR(object, type = "Lorenz-R2", pars.idx = pars.cv),
                      ineqExplained_PLR(object, type = "Lorenz-R2", pars.idx = pars.cv)))
  ans$coefficients <- cbind(ans$coefficients,
                            coef_PLR(object, renormalize = renormalize, pars.idx = pars.cv))

  if(ncol(ans$coefficients)==2){
    rownames(ans$ineq) <- c("BIC","CV")
    colnames(ans$coefficients) <- c("BIC","CV")
  }else{
    rownames(ans$ineq)[nrow(ans$ineq)] <- "CV"
    colnames(ans$coefficients)[ncol(ans$coefficients)] <- "CV"
  }

  class(ans) <- c("summary.PLR_cv",class(ans))

  return(ans)

}
