#' Undertakes a Lorenz regression
#'
#' \code{Lorenz.Reg} performs the Lorenz regression of a response with respect to several covariates.
#'
#' @param formula A formula object of the form \emph{response} ~ \emph{other_variables}.
#' @param data A data frame containing the variables displayed in the formula.
#' @param standardize Should the variables be standardized before the estimation process? Default value is TRUE.
#' @param ... Additional parameters corresponding to arguments passed in \code{\link{Lorenz.GA.cpp}} (i.e. parameters of the genetic algorithm).
#'
#' @return A list with several components:
#' \describe{
#'    \item{\code{theta}}{the estimated vector of parameters (on the original scale, even if \code{standardize} is TRUE).}
#'    \item{\code{expl.Gini}}{the estimated explained Gini coefficient}
#'    \item{\code{LorenzR2}}{the Lorenz-\eqn{R^2} of the regression.}
#'    \item{\code{MRS}}{the matrix of estimated marginal rates of substitution. More precisely, if we want the MRS of X1 (numerator) with respect to X2 (denominator),
#'    we should look for row corresponding to X1 and column corresponding to X2}.
#'    \item{\code{Fit}}{A data frame containing the response (first column) and the estimated index (second column).}
#' }
#'
#' @seealso \code{\link{Lorenz.GA.cpp}}, \code{\link{Lorenz.boot}}
#'
#' @examples
#' data(Data.Incomes)
#' Lorenz.Reg(Income ~ Age + Work.Hours, data = Data.Incomes)
#'
#'
#' @export

Lorenz.Reg <- function(formula, data, standardize=T, ...){

  return.list <- list()

  # 0. Obtain YX_mat ----

  # Transform the formula into dataframe
  Data.temp.X <- as.data.frame(stats::model.matrix(formula,data=data)[,-1])
  Data.temp.Y <- stats::model.frame(formula,data=data)[,1]
  Data.temp <- cbind(Data.temp.Y,Data.temp.X)
  colnames(Data.temp)[1] <- colnames(stats::model.frame(formula,data=data))[1]

  # Put the dataframe in the right format
  n.param<-length(Data.temp[1,])-1

  are.factor<-sapply(1:n.param,function(i)is.factor(Data.temp[,i+1])) #Before anything we need to treat the categorical variables

  if(sum(are.factor)!=0){

    length.factor<-sapply(which(are.factor),function(i)length(levels(Data.temp[,i+1])))

    YX_mat<-Data.temp[,-(which(are.factor)+1)]

    for (f in 1:sum(are.factor)){
      for (l in 2:length.factor[f]){
        tmp.var<-ifelse(Data.temp[,which(are.factor)[f]+1]==levels(Data.temp[,which(are.factor)[f]+1])[l],1,0)
        name.tmp.var<-paste(colnames(Data.temp)[which(are.factor)[f]+1],".",levels(Data.temp[,which(are.factor)[f]+1])[l],sep="")
        YX_mat<-cbind(YX_mat,tmp.var)
        colnames(YX_mat)[length(YX_mat[1,])]<-name.tmp.var
      }
    }
  }else{
    YX_mat <- Data.temp
  }

  n <- length(YX_mat[,1])
  p <- length(YX_mat[1,])-1

  # 1. Standardize X ----

  if (standardize){
    X <- YX_mat[,-1]
    X.temp <- stats::setNames(data.frame(matrix(ncol = p, nrow = n)), colnames(X))
    which.disc <- which(apply(X,2,function(x)length(unique(x)))<=2)
    n.disc <- length(which.disc)
    X.scale <- rep(NA,p)

    if(n.disc>0){
        X.temp[,which.disc] <- X[,which.disc]
        X.scale[which.disc] <- 1
    }
      if( (n.disc<p) & (n.disc>0) ){
        if (n.disc<(p-1)){
          X.center <- colMeans(X[,-which.disc])
        }else{
          X.center <- mean(X[,-which.disc])
        }
        X[,-which.disc] <- X[,-which.disc] - rep(X.center, rep.int(n,p-n.disc))
        if (n.disc<(p-1)){
          X.scale[-which.disc] <- sqrt(colSums(X[,-which.disc]^2)/(n-1))
        }else{
          X.scale[-which.disc] <- sqrt(sum(X[,-which.disc]^2)/(n-1))
        }
        X[,-which.disc] <- X[,-which.disc] / rep(X.scale[-which.disc], rep.int(n,p-n.disc))
    }
      if(n.disc==0){
        X.center <- colMeans(X)
        X <- X - rep(X.center, rep.int(n,p))
        X.scale <- sqrt(colSums(X^2)/(n-1))
        X <- X / rep(X.scale, rep.int(n,p))
    }
    YX_mat[,-1] <- X
  }

  # 2. Estimation of theta ----

  LR <- Lorenz.GA.cpp(YX_mat,...)

  # Estimation of the explained Gini coef
  Gi.expl <- LR$Gi.expl
  # Computation of Lorenz-R2
  LorenzR2 <- LR$LR2
  # Vector of estimated thetas
  theta <- LR$sol
  if (standardize){# Should take into account the standardization
    theta.standardize <- theta
    names(theta.standardize) <- colnames(YX_mat[,-1])
    MRS.standardize <- outer(theta.standardize,theta.standardize,"/")
    theta <- theta/X.scale
    theta <- theta/sum(abs(theta))
  }
  names(theta) <- colnames(YX_mat[,-1])
  # Matrix of MRS
  MRS <- outer(theta,theta,"/")
  # Estimated index
  Fit <- data.frame(Response = Data.temp[,1], Index = as.vector(theta%*%t(Data.temp[,-1])))

  return.list$theta <- theta
  return.list$expl.Gini <- Gi.expl
  return.list$LorenzR2 <- LorenzR2
  return.list$MRS <- MRS
  return.list$Fit <- Fit

  # 3. Estimation of regression curve ----

  # 4. Estimation of marginal impacts ----

  return(return.list)
}

