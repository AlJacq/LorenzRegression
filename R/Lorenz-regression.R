#' Undertakes a Lorenz regression
#'
#' \code{Lorenz.Reg} performs the Lorenz regression of a response with respect to several covariates.
#'
#' @param formula A formula object of the form \emph{response} ~ \emph{other_variables}.
#' @param data A data frame containing the variables displayed in the formula.
#' @param standardize Should the variables be standardized before the estimation process? Default value is TRUE.
#' @param weights vector of sample weights. By default, each observation is given the same weight.
#' @param parallel Whether parallel computing should be used. Default value is FALSE.
#' @param penalty should the regression include a penalty on the coefficients size.
#' If "none" is chosen, a non-penalized Lorenz regression is computed using function \code{\link{Lorenz.GA.cpp}}.
#' If "SCAD" is chosen, a penalized Lorenz regression with SCAD penalty is computed using function \code{\link{Lorenz.SCADFABS}}.
#' IF "LASSO" is chosen, a penalized Lorenz regression with LASSO penalty is computed using function \code{\link{Lorenz.FABS}}.
#' @param eps Only used if penalty="SCAD" or penalty="LASSO". Step size in the FABS or SCADFABS algorithm. Default value is 0.005.
#' @param lambda.choice Only used if penalty="SCAD" or penalty="LASSO". Determines what method is used to determine the optimal regularization parameter. Possibles values are "BIC" (Default), "CV" or c("BIC","CV"). In the last case, both methods are used.
#' @param nfolds Only used if lambda.choice contains "CV". Number of folds in the cross-validation.
#' @param seed.CV Only used if lambda.choice contains "CV". Should a specific seed be used in the definition of the folds. Default value is NULL in which case no seed is imposed.
#' @param foldID vector taking value from 1 to nfolds specifying the fold index of each observation. Default value is NULL in which case the folds are defined internally.
#' @param ... Additional parameters corresponding to arguments passed in \code{\link{Lorenz.GA.cpp}}, \code{\link{Lorenz.SCADFABS}} or \code{\link{Lorenz.FABS}} depending on the argument chosen in penalty.
#'
#' @return For the Non-penalized Lorenz Regression, a list with the following elements :
#' \describe{
#'    \item{\code{theta}}{the estimated vector of parameters.}
#'    \item{\code{summary}}{a vector including the estimated explained Gini coefficient and the Lorenz-\eqn{R^2}.}
#'    \item{\code{Gi.expl}}{the estimated explained Gini coefficient}
#'    \item{\code{LR2}}{the Lorenz-\eqn{R^2} of the regression.}
#'    \item{\code{MRS}}{the matrix of estimated marginal rates of substitution. More precisely, if we want the MRS of X1 (numerator) with respect to X2 (denominator),
#'    we should look for row corresponding to X1 and column corresponding to X2.}
#'    \item{\code{Fit}}{A data frame containing the response (first column) and the estimated index (second column).}
#' }
#' For the Penalized Lorenz Regression, a list with the following elements.
#' \describe{
#'    \item{\code{path}}{a matrix where the first line displays the path of regularization parameters. The second and third lines display the evolution of the Lorenz-\eqn{R^2} and explained Gini coefficient along that path. The remaining lines display the evolution of the estimated parameter vector.}
#'    \item{\code{theta}}{the estimated vector of parameters at the optimal value of the regularization parameter. If lambda.choice=c("BIC","CV"), it becomes a matrix where the first line corresponds to cross-validation, the second to BIC.}
#'    \item{\code{summary}}{a vector including the estimated explained Gini coefficient and the Lorenz-\eqn{R^2} at the optimal value of the regularization parameter. If lambda.choice=c("BIC","CV"), it becomes a matrix where the first line corresponds to cross-validation, the second to BIC.}
#'    \item{\code{Gi.expl}}{the estimated explained Gini coefficient at the optimal value of the regularization parameter. If lambda.choice=c("BIC","CV"), it becomes a vector where the first element corresponds to cross-validation, the second to BIC.}
#'    \item{\code{LR2}}{the Lorenz-\eqn{R^2} of the regression. If lambda.choice=c("BIC","CV"), it becomes a vector where the first element corresponds to cross-validation, the second to BIC.}
#'    \item{\code{MRS}}{the matrix of estimated marginal rates of substitution for non-zero coefficients at the optimal value of the regularization parameter. If lambda.choice=c("BIC","CV"), it becomes a list where the first element corresponds to cross-validation, the second to BIC.}
#'    \item{\code{Fit}}{A data frame containing the response (first column) and the estimated index (second column) at the optimal value of the regularization parameter. If lambda.choice=c("BIC","CV"), the dataframe contains three columns: the first corresponding to the response, the second to the estimated index obtained using CV. Finally, the third corresponds to the estimated index using BIC.}
#' }
#'
#' @seealso \code{\link{Lorenz.GA.cpp}}, \code{\link{Lorenz.SCADFABS}}, \code{\link{Lorenz.FABS}}, \code{\link{PLR.wrap}}
#'
#' @section References:
#' Heuchenne, C. and A. Jacquemain (2022). Inference for monotone single-index conditional means: A Lorenz regression approach. \emph{Computational Statistics & Data Analysis 167(C)}.
#'
#' @examples
#' data(Data.Incomes)
#' # 1. Non-penalized regression
#' NPLR <- Lorenz.Reg(Income ~ ., data = Data.Incomes, penalty = "none")
#' # 2. Penalized regression
#' PLR <- Lorenz.Reg(Income ~ ., data = Data.Incomes, penalty = "SCAD", lambda.choice = c("BIC","CV"), eps = 0.01, nfolds = 5)
#' # Comparison
#' NPLR$theta/sqrt(sum(NPLR$theta**2));PLR$theta
#' NPLR$summary;PLR$summary
#'
#' @export

Lorenz.Reg <- function(formula, data, standardize=T, weights=NULL, parallel=F, penalty=c("none","SCAD","LASSO"), eps=0.005, lambda.choice=c("BIC","CV")[1], nfolds=10, seed.CV=NULL, foldID=NULL, ...){

  # Check on penalty
  penalty <- match.arg(penalty)

  # Check on lambda.choice
  if( !all(lambda.choice%in%c("BIC","CV")) ) stop("lambda.choice should be a subvector of c(\"BIC\",\"CV\")")

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

  # 1. (Penalized) Lorenz Regression ----

  if(penalty == "none"){
    LR <- Lorenz.GA.cpp(YX_mat, standardize=standardize, weights=weights, parallel=parallel, ...)
  }else{
    LR <- PLR.wrap(YX_mat, standardize=standardize, weights=weights, penalty=penalty, eps=eps, ...)
  }

  # 2. Output of the PLR ----

  if(penalty == "none"){

    # Estimation of theta
    theta <- LR$theta # Vector of estimated coefficients
    names(theta) <- colnames(YX_mat[,-1])
    # Summary
    summary <- c()
    summary["Explained Gini"] <- LR$Gi.expl
    summary["Lorenz-R2"] <- LR$LR2
    # Matrix of MRS
    MRS <- outer(theta,theta,"/")
    # Estimated index
    Fit <- data.frame(Response = Data.temp[,1], Index = as.vector(theta%*%t(Data.temp[,-1])))
    # Return
    return.list$theta <- theta
    return.list$summary <- summary
    return.list$Gi.expl <- LR$Gi.expl
    return.list$LR2 <- LR$LR2
    return.list$MRS <- MRS
    return.list$Fit <- Fit

  }else{

    # Number of variables selected
    n_selected <- apply(LR$theta,2,function(x)sum(abs(x) > 10^(-10)))
    # Path
    Path <- rbind(LR$lambda, LR$LR2, LR$Gi.expl, n_selected)
    rownames(Path) <- c("lambda","Lorenz-R2","Explained Gini", "Number of nonzeroes")
    # BIC and/or CV
    if ("BIC" %in% lambda.choice){
      Path_BIC <- PLR.BIC(YX_mat, LR$theta, weights = weights)
      best.BIC <- Path_BIC$best
      val.BIC <- Path_BIC$val
      Path <- rbind(Path, val.BIC)
      rownames(Path)[nrow(Path)] <- "BIC score"
    }
    if ("CV" %in% lambda.choice){
      Path_CV <- PLR.CV(formula, data, penalty = penalty, PLR.est = LR, standardize = standardize, weights = weights, eps = eps, nfolds = nfolds, parallel = parallel, seed.CV = seed.CV, foldID = foldID, ...)
      best.CV <- Path_CV$best
      val.CV <- Path_CV$val
      Path <- rbind(Path, val.CV)
      rownames(Path)[nrow(Path)] <- "CV score"
    }
    lth <- nrow(Path)
    Path <- rbind(Path, LR$theta)
    rownames(Path)[(lth+1):nrow(Path)] <- colnames(YX_mat[,-1])

    # BIC only output
    if ( ("BIC" %in% lambda.choice) & !("CV" %in% lambda.choice) ){
      # Estimation of theta
      theta <- LR$theta[,best.BIC] # Vector of estimated coefficients
      names(theta) <- colnames(YX_mat[,-1])
      # Summary
      summary <- c()
      summary["Explained Gini"] <- LR$Gi.expl[best.BIC]
      summary["Lorenz-R2"] <- LR$LR2[best.BIC]
      summary["lambda"] <- LR$lambda[best.BIC]
      summary["Number of variables"] <- n_selected[best.BIC]
      summary["BIC score"] <- val.BIC[best.BIC]
      # Matrix of MRS
      theta.MRS <- theta[theta!=0]
      MRS <- outer(theta.MRS,theta.MRS,"/")
      # Estimated index
      Fit <- data.frame(Response = Data.temp[,1], Index = as.vector(theta%*%t(Data.temp[,-1])))
      # Return
      return.list$path <- Path
      return.list$theta <- theta
      return.list$summary <- summary
      return.list$Gi.expl <- LR$Gi.expl[best.BIC]
      return.list$LR2 <- LR$LR2[best.BIC]
      return.list$MRS <- MRS
      return.list$Fit <- Fit
    }
    # CV only output
    if ( ("CV" %in% lambda.choice) & !("BIC" %in% lambda.choice) ){
      # Estimation of theta
      theta <- LR$theta[,best.CV] # Vector of estimated coefficients
      names(theta) <- colnames(YX_mat[,-1])
      # Summary
      summary <- c()
      summary["Explained Gini"] <- LR$Gi.expl[best.CV]
      summary["Lorenz-R2"] <- LR$LR2[best.CV]
      summary["lambda"] <- LR$lambda[best.CV]
      summary["Number of variables"] <- n_selected[best.CV]
      summary["CV score"] <- val.CV[best.CV]
      # Matrix of MRS
      theta.MRS <- theta[theta!=0]
      MRS <- outer(theta.MRS,theta.MRS,"/")
      # Estimated index
      Fit <- data.frame(Response = Data.temp[,1], Index = as.vector(theta%*%t(Data.temp[,-1])))
      # Return
      return.list$path <- Path
      return.list$theta <- theta
      return.list$summary <- summary
      return.list$Gi.expl <- LR$Gi.expl[best.CV]
      return.list$LR2 <- LR$LR2[best.CV]
      return.list$MRS <- MRS
      return.list$Fit <- Fit
    }
    # CV and BIC output
    if ( ("CV" %in% lambda.choice) & ("BIC" %in% lambda.choice) ){
      # Estimation of theta
      theta <- t(LR$theta[,c(best.CV,best.BIC)]) # Vector of estimated coefficients
      colnames(theta) <- colnames(YX_mat[,-1])
      rownames(theta) <- c("CV","BIC")
      # Summary
      summary <- matrix(nrow=2,ncol=6)
      summary[,1] <- LR$Gi.expl[c(best.CV, best.BIC)]
      summary[,2] <- LR$LR2[c(best.CV, best.BIC)]
      summary[,3] <- LR$lambda[c(best.CV, best.BIC)]
      summary[,4] <- n_selected[c(best.CV, best.BIC)]
      summary[,5] <- val.BIC[c(best.CV, best.BIC)]
      summary[,6] <- val.CV[c(best.CV, best.BIC)]
      rownames(summary) <- c("CV","BIC")
      colnames(summary) <- c("Explained Gini", "Lorenz-R2", "lambda", "Number of variables", "BIC score", "CV score")
      # Matrix of MRS
      MRS <- list()
      theta.MRS.CV <- theta["CV",][theta["CV",]!=0]
      MRS$CV <- outer(theta.MRS.CV,theta.MRS.CV,"/")
      theta.MRS.BIC <- theta["BIC",][theta["BIC",]!=0]
      MRS$BIC <- outer(theta.MRS.BIC,theta.MRS.BIC,"/")
      # Estimated index
      Fit <- data.frame(Response = Data.temp[,1], Index.CV = as.vector(theta["CV",]%*%t(Data.temp[,-1])), Index.BIC = as.vector(theta["BIC",]%*%t(Data.temp[,-1])))
      # Return
      return.list$path <- Path
      return.list$theta <- theta
      return.list$summary <- summary
      return.list$Gi.expl <- summary[,1]
      return.list$LR2 <- summary[,2]
      return.list$MRS <- MRS
      return.list$Fit <- Fit
    }
  }

  return(return.list)
}

