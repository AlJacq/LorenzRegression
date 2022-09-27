#' Estimated coefficients for the Penalized Lorenz Regression
#'
#' \code{coef.PLR} provides the estimated coefficients for an object of class \code{PLR}.
#'
#' @param PLR Output of a call to \code{\link{Lorenz.Reg}}, where \code{penalty!="none"}.
#'
#' @return If the PLR was fitted with only one selection method, the output is a vector gathering the estimated coefficients.
#' If several selection methods were selected, it outputs a list of vectors, where each element of the list corresponds to a different selection method.
#'
#' @seealso \code{\link{Lorenz.Reg}}
#'
#' @examples
#' data(Data.Incomes)
#' PLR <- Lorenz.Reg(Income ~ ., data = Data.Incomes, penalty = "SCAD", lambda.choice = c("BIC","CV"), eps = 0.005, seed.CV = 123, nfolds = 5)
#' coef.PLR(PLR)
#'
#' @method coef PLR
#' @export

coef.PLR <- function(PLR){

  m1 <- PLR$theta

  if(nrow(m1) > 1){

    m2 <- t(PLR$theta)
    l <- split(m2,rep(1:ncol(m2), each = nrow(m2)))
    names(l) <- rownames(PLR$theta)

    return(l)

  }else{

    return(m1[1,])

  }

}
