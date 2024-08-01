#' Summary for the Penalized Lorenz Regression
#'
#' \code{summary.PLR} provides a summary for an object of class \code{"PLR"}.
#'
#' @param object An object of class \code{"PLR"}.
#' @param renormalize A logical value determining whether the coefficient vector should be re-normalized to match the representation where the first category of each categorical variable is omitted. Default value is TRUE
#' @param pars.idx A vector of size 2 specifying the index of the grid parameter (first element) and the index of the penalty parameter (second element) that should be selected.
#' Default is \code{NULL}, in which case the parameters are selected by the available methods : BIC (always), bootstrap (if \code{object} inherits from the \code{PLR_boot} class) and cross-validation (if \code{object} inherits from the \code{PLR_cv} class).
#' @param ... Additional arguments
#'
#' @return An object of class \code{"summary.PLR"}, containing the following elements:
#' \describe{
#'    \item{\code{call}}{The matched call.}
#'    \item{\code{ineq}}{The table of explained inequality. The first column displays the explained Gini coefficient, the second displays the Gini coefficient of the response and the last column displays the Lorenz-R2. Each row corresponds to a different selection method.}
#'    \item{\code{coefficients}}{A matrix providing information on the estimated coefficients. Each column corresponds to a different selection method}
#' }
#'
#' @seealso \code{\link{Lorenz.Reg}}, \code{\link{Lorenz.boot}}
#'
#' @examples
#' ## For examples see example(Lorenz.Reg), example(Lorenz.boot) and example(PLR.CV)
#'
#' @method summary PLR
#' @export

summary.PLR <- function(object, renormalize=TRUE, pars.idx=NULL, ...){

  if (!inherits(object, "PLR")) stop("The object must be of class 'PLR'")

  if(!is.null(pars.idx)){
    l <- ncol(object$x)
    pth <- object$path[[pars.idx[1]]][,pars.idx[2]]
    object$theta <- pth[(length(pth)-l+1):length(pth)]
    object$Gi.expl <- object$path[[pars.idx[1]]]["Explained Gini",pars.idx[2]]
    object$LR2 <- object$path[[pars.idx[1]]]["Lorenz-R2",pars.idx[2]]
    class(object) <- "PLR"
  }

  if(renormalize){
    m1 <- PLR.normalize(object)
  }else{
    m1 <- object$theta
  }

  ans <- list()
  ans$call <- object$call
  if(length(object$Gi.expl)==1) names(object$Gi.expl) <- ""
  ans$ineq <- matrix(c(object$Gi.expl,object$Gi.expl/object$LR2,object$LR2),
                     nrow = length(object$Gi.expl), ncol = 3,
                     dimnames = list(names(object$Gi.expl), c("Explained","Total","Lorenz-R2")))


  if (inherits(object, c("PLR_boot","PLR_cv"))){
    ans$coefficients <- t(m1)
  }else{
    ans$coefficients <- as.matrix(m1)
    colnames(ans$coefficients) <- "Estimate"
  }

  class(ans) <- "summary.PLR"

  return(ans)

}
