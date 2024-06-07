#' Estimated coefficients for the Penalized Lorenz Regression
#'
#' \code{coef.PLR} provides the estimated coefficients for an object of class \code{"PLR"}.
#'
#' @param object An object of S3 class \code{"PLR"}.
#' @param renormalize A logical value determining whether the coefficient vector should be re-normalized to match the representation where the first category of each categorical variable is omitted. Default value is TRUE
#' @param ... Additional arguments
#'
#' @return a vector gathering the estimated coefficients.
#' If the object has also class \code{"PLR_boot"} and/or \code{"PLR_cv"}, the output is a matrix, where each column corresponds to a selection method.
#'
#' @seealso \code{\link{Lorenz.Reg}}
#'
#' @examples
#' data(Data.Incomes)
#' PLR <- Lorenz.Reg(Income ~ ., data = Data.Incomes, penalty = "SCAD",
#'                   h.grid = nrow(Data.Incomes)^(-1/5.5), sel.choice = c("BIC","CV"),
#'                   eps = 0.01, seed.CV = 123, nfolds = 5)
#' coef(PLR)
#'
#' @method coef PLR
#' @export

coef.PLR <- function(object, renormalize=TRUE, ...){

  if (!inherits(object, "PLR")) stop("The object must be of class 'PLR'")

  if(renormalize){
    m1 <- PLR.normalize(object)
  }else{
    m1 <- object$theta
  }

  if(is.matrix(m1)){

    m2 <- t(m1)
    l <- split(m2,rep(1:ncol(m2), each = nrow(m2)))
    names(l) <- colnames(m2)
    for (j in 1:length(l)) names(l[[j]]) <- rownames(m2)

    return(l)

  }else{

    return(m1)

  }

}
