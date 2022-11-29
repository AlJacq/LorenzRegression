#' Estimated coefficients for the Penalized Lorenz Regression
#'
#' \code{coef.PLR} provides the estimated coefficients for an object of class \code{PLR}.
#'
#' @param object Output of a call to \code{\link{Lorenz.Reg}}, where \code{penalty!="none"}.
#' @param renormalize whether the coefficient vector should be re-normalized to match the representation where the first category of each categorical variable is omitted. Default value is TRUE
#' @param ... Additional arguments
#'
#' @return If the PLR was fitted with only one selection method, the output is a vector gathering the estimated coefficients.
#' If several selection methods were selected, it outputs a list of vectors, where each element of the list corresponds to a different selection method.
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

  PLR <- object
  if(renormalize){
    m1 <- PLR.normalize(PLR)
  }else{
    m1 <- PLR$theta
  }

  if(nrow(m1) > 1){

    m2 <- t(m1)
    l <- split(m2,rep(1:ncol(m2), each = nrow(m2)))
    names(l) <- rownames(PLR$theta)
    for (j in 1:length(l)) names(l[[j]]) <- rownames(m2)

    return(l)

  }else{

    return(m1[1,])

  }

}
