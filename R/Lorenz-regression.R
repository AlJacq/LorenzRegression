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
#' @param lambda.choice Only used if penalty="SCAD" or penalty="LASSO". Determines what method is used to determine the optimal regularization parameter. Possibles values are any subvector of c("BIC","CV","Boot"). Default is "BIC". Notice that "Boot" is necessarily added if Boot.inference is set to TRUE.
#' @param nfolds Only used if lambda.choice contains "CV". Number of folds in the cross-validation.
#' @param seed.CV Only used if lambda.choice contains "CV". Should a specific seed be used in the definition of the folds. Default value is NULL in which case no seed is imposed.
#' @param foldID vector taking value from 1 to nfolds specifying the fold index of each observation. Default value is NULL in which case the folds are defined internally.
#' @param Boot.inference should bootstrap inference be produced ? Default is FALSE. It is automatically turned to TRUE if lambda.choice contains "Boot".
#' @param B Only used if Boot.inference is TRUE. Number of bootstrap resamples. Default is 400.
#' @param which.CI Only used if Boot.inference is TRUE. Determines what method is used to construct the bootstrap confidence intervals. Sub-vector of c("Basic","Perc","Param"), the default being the whole vector.
#' "Basic" corresponds to basic bootstrap, based on bootstrapping the whole distribution of the estimated explained Gini coefficient.
#' "Perc" corresponds to percentile bootstrap, where the quantiles of bootstrap distribution of the estimated explained Gini are directly plugged in the bounds.
#' "Param" corresponds to parametric bootstrap, which exploits the asymptotic normality and only bootstraps the asymptotic variance.
#' @param alpha Only used if Boot.inference is TRUE. significance level for the bootstrap confidence intervals. Default is 0.05.
#' @param bootID Only used if Boot.inference is TRUE. matrix where each row provides the ID of the observations selected in each bootstrap resample. Default is NULL, in which case these are defined internally.
#' @param seed.boot Only used if Boot.inference is TRUE. Should a specific seed be used in the definition of the folds. Default value is NULL in which case no seed is imposed.
#' @param
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
#'    \item{\code{CI.Gi}}{Only returned if Boot.inference is TRUE. A matrix where each line corresponds to a bootstrap method, as specified in which.CI, and the columns provide the bounds of the intervals for the explained Gini coefficient.}
#'    \item{\code{CI.LR2}}{Only returned if Boot.inference is TRUE. A matrix where each line corresponds to a bootstrap method, as specified in which.CI, and the columns provide the bounds of the intervals for the Lorenz-\eqn{R^2}.}
#' }
#' For the Penalized Lorenz Regression, a list with the following elements.
#' \describe{
#'    \item{\code{path}}{a matrix where the first line displays the path of regularization parameters. The second and third lines display the evolution of the Lorenz-\eqn{R^2} and explained Gini coefficient along that path. The next lines display the evolution of the scores of the methods chosen in lambda.choice. The remaining lines display the evolution of the estimated parameter vector.}
#'    \item{\code{theta}}{a matrix where the different lines correspond to the methods chosen in lambda.choice. Each line provides the estimated vector of parameters at the optimal value of the regularization parameter.
#'    \item{\code{summary}}{a matrix where the different lines correspond to the methods chosen in lambda.choice. Each line provides the estimated explained Gini coefficient, the Lorenz-\eqn{R^2} and the scores at the optimal value of the regularization parameter.}
#'    \item{\code{Gi.expl}}{a vector providing the estimated explained Gini coefficient at the optimal value of the regularization parameter for each method in lambda.choice.}
#'    \item{\code{LR2}}{a vector providing the Lorenz-\eqn{R^2} at the optimal value of the regularization parameter for each method in lambda.choice.}
#'    \item{\code{MRS}}{a list where the different elements correspond to a method in lambda.choice. Each element is a matrix of estimated marginal rates of substitution for non-zero coefficients at the optimal value of the regularization parameter.}
#'    \item{\code{Fit}}{A data frame containing the response (first column). The remaining columns give the estimated index at the optimal value of the regularization parameter, for each method chosen in lambda.choice.}
#'    \item{\code{CI.Gi}}{Only returned if Boot.inference is TRUE. A list of matrices where the different elements correspond to the bootstrap method chosen with which.CI. For each matrix, the rows correspond to the method chosen for the selection of the regularization parameter and the columns provide the bounds of the confidence interval for the explained Gini coefficient.}
#'    \item{\code{CI.LR2}}{Only returned if Boot.inference is TRUE. A list of matrices where the different elements correspond to the bootstrap method chosen with which.CI. For each matrix, the rows correspond to the method chosen for the selection of the regularization parameter and the columns provide the bounds of the confidence interval for the Lorenz-\eqn{R^2}.}
#'    \item{\code{CI.Gi.path}}{Only returned if Boot.inference is TRUE. A list of matrices where the different elements of the lists correspond to a different lambda.value. For each matrix, the rows correspond to the bootstrap methods chosen in which.CI and the columns provide the bounds of the confidence interval for the explained Gini coefficient.}
#'    \item{\code{CI.LR2.path}}{Only returned if Boot.inference is TRUE. A list of matrices where the different elements of the lists correspond to a different lambda.value. For each matrix, the rows correspond to the bootstrap methods chosen in which.CI and the columns provide the bounds of the confidence interval for the Lorenz-\eqn{R^2}.}
#' }
#'
#' @seealso \code{\link{Lorenz.GA.cpp}}, \code{\link{Lorenz.SCADFABS}}, \code{\link{Lorenz.FABS}}, \code{\link{PLR.wrap}}, \code{\link{PLR.boot}}
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

Lorenz.Reg <- function(formula,
                       data,
                       standardize=T,
                       weights=NULL,
                       parallel=F,
                       penalty=c("none","SCAD","LASSO"),
                       eps=0.005,
                       lambda.choice=c("BIC","CV","Boot")[1],
                       nfolds=10,
                       seed.CV=NULL,
                       foldID=NULL,
                       Boot.inference=FALSE,
                       B=400,
                       which.CI=c("Basic","Perc","Param"),
                       alpha=0.05,
                       bootID=NULL,
                       seed.boot=NULL,
                       ...){

  # Check on penalty
  penalty <- match.arg(penalty)

  # Check on lambda.choice
  if( !all(lambda.choice%in%c("BIC","CV","Boot")) ) stop("lambda.choice should be a subvector of c(\"BIC\",\"CV\",\"Boot\")")

  # Check on bootstrap
  # The idea is that if bootstrap is required, we provide all the results that are available anyway
  if( Boot.inference ) lambda.choice <- union(lambda.choice,"Boot")
  if ("Boot"%in%lambda.choice) Boot.inference <- T

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
    # Bootstrap
    LR_Boot <- Lorenz.boot(formula, data, standardize = standardize, weights = weights, LR.est = LR, penalty = penalty, which.CI = which.CI, alpha = alpha, B = B, bootID = bootID, seed.boot = seed.boot, parallel = parallel, ...)
    # Return
    return.list$theta <- theta
    return.list$summary <- summary
    return.list$Gi.expl <- LR$Gi.expl
    return.list$LR2 <- LR$LR2
    return.list$MRS <- MRS
    return.list$Fit <- Fit
    return.list$CI.Gi <- LR_Boot$CI.Gi
    return.list$CI.LR2 <- LR_Boot$CI.LR2

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
    if ("Boot" %in% lambda.choice){
      Path_Boot <- Lorenz.boot(formula, data, standardize = standardize, weights = weights, LR.est = LR, penalty = penalty, eps = eps, which.CI = which.CI, alpha = alpha, B = B, bootID = bootID, seed.boot = seed.boot, parallel = parallel, ...)
      best.Boot <- Path_Boot$OOB.best
      val.Boot <- Path_Boot$OOB.total
      CI.Gi <- Path_Boot$CI.Gi
      CI.LR2 <- Path_Boot$CI.LR2
      Path <- rbind(Path, val.Boot)
      rownames(Path)[nrow(Path)] <- "Boot score"
    }

    lth <- nrow(Path)
    Path <- rbind(Path, LR$theta)
    rownames(Path)[(lth+1):nrow(Path)] <- colnames(YX_mat[,-1])

    # Output without bootstrap
    theta <- matrix(nrow = length(lambda.choice), ncol = ncol(YX_mat)-1)
    colnames(theta) <- colnames(YX_mat[,-1])
    summary <- matrix(nrow = length(lambda.choice), ncol = 4 + length(lambda.choice))
    colnames(summary) <- c("Explained Gini", "Lorenz-R2", "lambda", "Number of variables", paste0(lambda.choice," score"))
    Gi.expl <- rep(NA,length(lambda.choice))
    LR2 <- rep(NA,length(lambda.choice))
    MRS <- lapply(1:length(lambda.choice),function(x)matrix(nrow = ncol(YX_mat)-1, ncol = ncol(YX_mat)-1))
    Fit <- data.frame(Response = Data.temp[,1])
    names(MRS) <- rownames(theta) <- rownames(summary) <- names(Gi.expl) <- names(LR2) <- lambda.choice

    if( "BIC" %in% lambda.choice ){

      i.BIC <- which(lambda.choice == "BIC")
      # theta
      theta[i.BIC,] <- LR$theta[,best.BIC]
      # summary
      summary[i.BIC,1] <- LR$Gi.expl[best.BIC]
      summary[i.BIC,2] <- LR$LR2[best.BIC]
      summary[i.BIC,3] <- LR$lambda[best.BIC]
      summary[i.BIC,4] <- n_selected[best.BIC]
      summary[i.BIC,"BIC score"] <- val.BIC[best.BIC]
      if ("Boot" %in% lambda.choice) summary[i.BIC,"Boot score"] <- val.Boot[best.BIC]
      if ("CV" %in% lambda.choice) summary[i.BIC,"CV score"] <- val.CV[best.BIC]
      # Gi.expl and LR2
      Gi.expl[i.BIC] <- summary[i.BIC,1]
      LR2[i.BIC] <- summary[i.BIC,2]
      # Matrix of MRS
      theta.MRS.BIC <- theta[i.BIC,][theta[i.BIC,]!=0]
      MRS$BIC <- outer(theta.MRS.BIC,theta.MRS.BIC,"/")
      # Fit
      Fit$Index.BIC <- as.vector(theta[i.BIC,]%*%t(Data.temp[,-1]))

    }

    if( "Boot" %in% lambda.choice ){

      i.Boot <- which(lambda.choice == "Boot")
      # theta
      theta[i.Boot,] <- LR$theta[,best.Boot]
      # summary
      summary[i.Boot,1] <- LR$Gi.expl[best.Boot]
      summary[i.Boot,2] <- LR$LR2[best.Boot]
      summary[i.Boot,3] <- LR$lambda[best.Boot]
      summary[i.Boot,4] <- n_selected[best.Boot]
      summary[i.Boot,"Boot score"] <- val.Boot[best.Boot]
      if ("BIC" %in% lambda.choice) summary[i.Boot,"BIC score"] <- val.BIC[best.Boot]
      if ("CV" %in% lambda.choice) summary[i.Boot,"CV score"] <- val.CV[best.Boot]
      # Gi.expl and LR2
      Gi.expl[i.Boot] <- summary[i.Boot,1]
      LR2[i.Boot] <- summary[i.Boot,2]
      # Matrix of MRS
      theta.MRS.Boot <- theta[i.Boot,][theta[i.Boot,]!=0]
      MRS$Boot <- outer(theta.MRS.Boot,theta.MRS.Boot,"/")
      # Fit
      Fit$Index.Boot <- as.vector(theta[i.Boot,]%*%t(Data.temp[,-1]))

    }

    if( "CV" %in% lambda.choice ){

      i.CV <- which(lambda.choice == "CV")
      # theta
      theta[i.CV,] <- LR$theta[,best.CV]
      # summary
      summary[i.CV,1] <- LR$Gi.expl[best.CV]
      summary[i.CV,2] <- LR$LR2[best.CV]
      summary[i.CV,3] <- LR$lambda[best.CV]
      summary[i.CV,4] <- n_selected[best.CV]
      summary[i.CV,"CV score"] <- val.CV[best.CV]
      if ("BIC" %in% lambda.choice) summary[i.CV,"BIC score"] <- val.BIC[best.CV]
      if ("Boot" %in% lambda.choice) summary[i.CV,"Boot score"] <- val.Boot[best.CV]
      # Gi.expl and LR2
      Gi.expl[i.CV] <- summary[i.CV,1]
      LR2[i.CV] <- summary[i.CV,2]
      # Matrix of MRS
      theta.MRS.CV <- theta[i.CV,][theta[i.CV,]!=0]
      MRS$CV <- outer(theta.MRS.CV,theta.MRS.CV,"/")
      # Fit
      Fit$Index.CV <- as.vector(theta[i.CV,]%*%t(Data.temp[,-1]))

    }

    # Return
    return.list$path <- Path
    return.list$theta <- theta
    return.list$summary <- summary
    return.list$Gi.expl <- Gi.expl
    return.list$LR2 <- LR2
    return.list$MRS <- MRS
    return.list$Fit <- Fit

    # Output of the bootstrap

    if (Boot.inference){

      CI.Gi.list <- lapply(1:length(which.CI),function(x)matrix(nrow = length(lambda.choice), ncol = 2))
      names(CI.Gi.list) <- which.CI
      for (l in 1:length(CI.Gi.list)){
        rownames(CI.Gi.list[[l]]) <- lambda.choice
        colnames(CI.Gi.list[[l]]) <- c("Lower bound","Upper bound")
      }
      CI.LR2.list <- CI.Gi.list

      if( "BIC" %in% lambda.choice ){
        for (l in 1:length(CI.Gi.list)){
          CI.Gi.list[[which.CI[l]]]["BIC",] <- CI.Gi[[best.BIC]][which.CI[l],]
          CI.LR2.list[[which.CI[l]]]["BIC",] <- CI.LR2[[best.BIC]][which.CI[l],]
        }
      }

      if( "Boot" %in% lambda.choice ){
        for (l in 1:length(CI.Gi.list)){
          CI.Gi.list[[which.CI[l]]]["Boot",] <- CI.Gi[[best.Boot]][which.CI[l],]
          CI.LR2.list[[which.CI[l]]]["Boot",] <- CI.LR2[[best.Boot]][which.CI[l],]
        }
      }

      if( "CV" %in% lambda.choice ){
        for (l in 1:length(CI.Gi.list)){
          CI.Gi.list[[which.CI[l]]]["CV",] <- CI.Gi[[best.CV]][which.CI[l],]
          CI.LR2.list[[which.CI[l]]]["CV",] <- CI.LR2[[best.CV]][which.CI[l],]
        }
      }

      return.list$CI.Gi <- CI.Gi.list
      return.list$CI.LR2 <- CI.LR2.list
      return.list$CI.Gi.path <- CI.Gi
      return.list$CI.LR2.path <- CI.LR2

    }

  }

  return(return.list)
}

