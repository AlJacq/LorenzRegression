#' Bootstrap confidence intervals
#'
#' \code{boot.confint} computes boostrap confidence intervals given an estimation on the original sample and on the boostrap samples
#'
#' @param x.hat estimator on the original sample.
#' @param x.star vector gathering the estimation on the bootstrapped sample.
#' @param alpha 1-level of the confidence interval
#' @param boot.method bootstrap method.
#'
#' @return A vector of dimension two with the desired confidence interval

boot.confint <- function(x.hat, x.star, alpha, boot.method){

  if(boot.method %in% c("Basic","Perc")){

    B <- length(x.star)
    l1 <- floor(alpha/2*B+1)
    l2 <- ceiling((1-alpha/2)*B)

    x1 <- sort(x.star)[l1]
    x2 <- sort(x.star)[l2]

    if(boot.method=="Perc") CI <- c(x1,x2)
    if(boot.method=="Basic") CI <- 2*x.hat - c(x2,x1)

  }else{

    Sigma.star <- var(x.star,na.rm=T)
    CI <- x.hat+c(-1,1)*sqrt(Sigma.star)*qnorm(1-alpha/2)

  }

  names(CI) <- c("Lower bound", "Upper bound")
  return(CI)

}
