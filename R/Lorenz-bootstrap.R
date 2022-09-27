#' Produces bootstrap-based inference for (penalized) Lorenz regression
#'
#' \code{Lorenz.boot} constructs bootstrap confidence intervals for the explained Gini coefficient and Lorenz-\eqn{R^2} and, if applies, selects the regularization parameter.
#'
#' @param formula A formula object of the form \emph{response} ~ \emph{other_variables}.
#' @param data A data frame containing the variables displayed in the formula.
#' @param standardize Should the variables be standardized before the estimation process? Default value is TRUE.
#' @param weights vector of sample weights. By default, each observation is given the same weight.
#' @param LR.est Estimation on the original sample. Output of a call to \code{\link{Lorenz.GA}} or \code{\link{PLR.wrap}}.
#' @param penalty should the regression include a penalty on the coefficients size.
#' If "none" is chosen, a non-penalized Lorenz regression is computed using function \code{\link{Lorenz.GA}}.
#' If "SCAD" is chosen, a penalized Lorenz regression with SCAD penalty is computed using function \code{\link{Lorenz.SCADFABS}}.
#' IF "LASSO" is chosen, a penalized Lorenz regression with LASSO penalty is computed using function \code{\link{Lorenz.FABS}}.
#' @param h Only used if penalty="SCAD" or penalty="LASSO". Bandwidth of the kernel, determining the smoothness of the approximation of the indicator function. Default value is NULL (unpenalized case) but has to be specified if penalty="LASSO" or penalty="SCAD".
#' @param eps Only used if penalty="SCAD" or penalty="LASSO". Step size in the FABS or SCADFABS algorithm. Default value is 0.005.
#' @param which.CI Determines what method is used to construct the bootstrap confidence intervals. Sub-vector of c("Basic","Perc","Param"), the default being the whole vector.
#' "Basic" corresponds to basic bootstrap, based on bootstrapping the whole distribution of the estimated explained Gini coefficient.
#' "Perc" corresponds to percentile bootstrap, where the quantiles of bootstrap distribution of the estimated explained Gini are directly plugged in the bounds.
#' "Param" corresponds to parametric bootstrap, which exploits the asymptotic normality and only bootstraps the asymptotic variance.
#' @param alpha significance level for the bootstrap confidence intervals. Default is 0.05.
#' @param B Number of bootstrap resamples. Default is 400.
#' @param bootID matrix where each row provides the ID of the observations selected in each bootstrap resample. Default is NULL, in which case these are defined internally.
#' @param seed.boot Should a specific seed be used in the definition of the folds. Default value is NULL in which case no seed is imposed.
#' @param parallel Whether parallel computing should be used to distribute the \code{B} computations on different CPUs. Either a logical value determining whether parallel computing is used (TRUE) or not (FALSE, the default value). Or a numerical value determining the number of cores to use.
#' @param ... Additional parameters corresponding to arguments passed in \code{\link{Lorenz.GA}}, \code{\link{Lorenz.SCADFABS}} or \code{\link{Lorenz.FABS}} depending on the argument chosen in penalty.
#' @return A list with several components:
#' \describe{
#'    \item{\code{LR.est}}{Estimation on the original sample.}
#'    \item{\code{CI.Gi}}{In the unpenalized case, a matrix where the rows correspond to the bootstrap methods and the columns correspond to the bounds of the confidence interval for the explained Gini coefficient. In the penalized case, a list of those matrices, where each element of the list corresponds to a different value of the regularization parameter.}
#'    \item{\code{CI.LR2}}{In the unpenalized case, a matrix where the rows correspond to the bootstrap methods and the columns correspond to the bounds of the confidence interval for the Lorenz-\eqn{R^2}. In the penalized case, a list of those matrices, where each element of the list corresponds to a different value of the regularization parameter.}
#'    \item{\code{Gi.star}}{In the unpenalized case, a vector gathering the bootstrap estimators of the explained Gini coefficient. In the penalized case, it becomes a list of vectors.}
#'    \item{\code{LR2.star}}{In the unpenalized case, a vector gathering the bootstrap estimators of the Lorenz-\eqn{R^2}. In the penalized case, it becomes a list of vectors.}
#'    \item{\code{OOB.total}}{In the penalized case only. Vector gathering the OOB-score for each lambda value.}
#'    \item{\code{OOB.best}}{In the penalized case only. index of the lambda value attaining the highest OOB-score.}
#' }
#'
#' @seealso \code{\link{Lorenz.Reg}}, \code{\link{Lorenz.GA}}, \code{\link{Lorenz.SCADFABS}}, \code{\link{Lorenz.FABS}}, \code{\link{PLR.wrap}}
#'
#' @section References:
#' Heuchenne, C. and A. Jacquemain (2022). Inference for monotone single-index conditional means: A Lorenz regression approach. \emph{Computational Statistics & Data Analysis 167(C)}.
#' Jacquemain, A., C. Heuchenne, and E. Pircalabelu (2022). A penalised bootstrap estimation procedure for the explained Gini coefficient.
#'
#' @examples
#' data(Data.Incomes)
#' Lorenz.boot(Income ~ ., data = Data.Incomes, penalty = "SCAD", h = nrow(Data.Incomes)^(-1/5.5), eps = 0.01, B = 20, seed.boot = 123)
#'
#'
#' @export

Lorenz.boot<-function(formula,
                      data,
                      standardize=T,
                      weights=NULL,
                      LR.est=NULL,
                      penalty=c("none","SCAD","LASSO"),
                      h=NULL,
                      eps=0.005,
                      which.CI=c("Basic","Perc","Param"),
                      alpha=0.05,
                      B = 400,
                      bootID = NULL,
                      seed.boot = NULL,
                      parallel=F,
                      ...
){

  # Check on penalty
  penalty <- match.arg(penalty)

  # Check on which.CI
  if( !all(which.CI%in%c("Basic","Perc","Param")) ) stop("which.CI should be a subvector of c(\"Basic\",\"Perc\",\"Param\")")

  # Number of bootstrap resamples
  if( !is.null(bootID) ) B <- nrow(bootID)

  # Check on h
  if( is.null(h) & penalty %in% c("SCAD","LASSO") ) stop("h has to be specified in the penalized case")

  # PRE-BOOT ----

  Basic.T <- ("Basic"%in%which.CI)
  Perc.T <- ("Perc"%in%which.CI)
  Param.T <- ("Param"%in%which.CI)

  # PRE-BOOT > GETTING YX_MAT ----

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

  # PRE-BOOT > STANDARDIZE X ----
  # No need to do it since it will be dealt with later in

  # PRE-BOOT > (P)LR ----

  if(is.null(LR.est)){
    if(penalty == "none"){
      LR.est <- LorenzRegression::Lorenz.GA(YX_mat, standardize = standardize, weights = weights, parallel = parallel, ...)
    }else{
      LR.est <- LorenzRegression::PLR.wrap(YX_mat, standardize = standardize, weights = weights, h = h, penalty = penalty, eps = eps, ...)
    }
  }
  theta.hat <- LR.est$theta
  Gi.hat <- LR.est$Gi.expl
  LR2 <- LR.est$LR2

  # PRE-BOOT > Boot ID ----

  if (is.null(bootID)){
    if(!is.null(seed.boot)) set.seed(seed.boot)
    bootID <- t(sapply(1:B,function(b)sample(1:n, replace = TRUE)))
  }

  # BOOT > INNER ----

  Boot.inner <- function(b, ...){

    Return.list <- list()

    # Construct Test and Validation bootstrap samples
    idx.boot <- bootID[b,]
    YX_mat.test <- YX_mat[idx.boot,]
    weights.test <- weights[idx.boot]
    if(penalty != "none"){
      YX_mat.valid <- YX_mat[-unique(idx.boot),]
      weights.valid <- weights[-unique(idx.boot)]
    }

    # Perform the estimation
    if(penalty == "none"){
      LR.est.star <- LorenzRegression::Lorenz.GA(YX_mat.test, standardize = standardize, weights = weights.test, parallel=parallel, ...)
    }else{
      LR.est.star <- LorenzRegression::PLR.wrap(YX_mat.test, standardize = standardize, weights = weights.test, penalty = penalty, h = h, eps = eps, lambda = LR.est$lambda, ...)
      lambda.star <- LR.est.star$lambda
    }
    theta.star <- LR.est.star$theta
    Gi.star <- LR.est.star$Gi.expl
    LR2.star <- LR.est.star$LR2

    # Compute the OOB-score (if penalized LR)
    if(penalty != "none"){
      y.valid <- YX_mat.valid[,1]
      X.valid <- as.matrix(YX_mat.valid[,-1])
      n.valid <- length(y.valid)
      OOB.score <- apply(theta.star,2,function(x)Gini.coef(y = y.valid, x = X.valid%*%x, na.rm=T, ties.method = "mean", weights = weights.valid))
      # With SCAD, the algorithm may stop sooner than in original sample. Hence, the lambda vectors might be different (the bootstrapped might be shorter)
      if( length(lambda.star) != length(LR.est$lambda) ){
        diff.lengths <- length(LR.est$lambda)-length(lambda.star)
        theta.star <- cbind(theta.star, matrix(NA,p,diff.lengths))
        Gi.star <- c(Gi.star,rep(NA,diff.lengths))
        LR2.star <- c(LR2.star,rep(NA,diff.lengths))
        OOB.score <- c(OOB.score,rep(NA,diff.lengths))
      }
      Return.list$OOB.score <- OOB.score
    }

    Return.list$theta.star <- theta.star
    Return.list$Gi.star <- Gi.star
    Return.list$LR2.star <- LR2.star

    return(Return.list)

  }

  # BOOT > ITERATIONS ----
  if(parallel){
    if(is.numeric(parallel)){
      registerDoParallel(parallel)
    }else{
      numCores <- detectCores()
      registerDoParallel(numCores-1)
    }
    Boot.b <- foreach(b=1:B) %dopar% {
      Boot.inner(b)
    }
    stopImplicitCluster()
  }else{
    Boot.b <- foreach(b=1:B) %do% {
      Boot.inner(b)
    }
  }

  # BOOT > ITERATIONS > RETRIEVE ----
  if(penalty != "none"){

    OOB.matrix <- t(sapply(1:B,function(b)Boot.b[[b]]$OOB.score))
    L <- ncol(OOB.matrix)
    iter.ok <- which(apply(OOB.matrix,2,function(x)mean(is.na(x)))< 0.05)
    OOB.total <- c()
    for (i in 1:L){
      if (i %in% iter.ok){
        OOB.total[i] <- mean(OOB.matrix[,i], na.rm=T)
      }else{
        OOB.total[i] <- -Inf
      }
    }
    OOB.best <- which.max(OOB.total)
    lambda.OOB <- LR.est$lambda[OOB.best]
    Gi.star <- lapply(1:L,function(i)sapply(1:B,function(b)Boot.b[[b]]$Gi.star[i]))
    LR2.star <- lapply(1:L,function(i)sapply(1:B,function(b)Boot.b[[b]]$LR2.star[i]))

  }

  # BOOT > POST-ITER ----
  if(penalty != "none"){

    CI.construction <- function(i){

      if (i %in% iter.ok){

        Gi.star.i <- Gi.star[[i]]
        LR2.star.i <- LR2.star[[i]]
        B.i <- length(sort(Gi.star.i)) # Watch out, it might be different than B because of NA's !!
        l1.i <- floor(alpha/2*B.i+1)
        l2.i <- ceiling((1-alpha/2)*B.i)
        CI.Gi.i <- matrix(nrow=length(which.CI),ncol=2)
        colnames(CI.Gi.i) <- c("Lower.limit","Upper.limit")
        rownames(CI.Gi.i) <- which.CI
        CI.LR2.i <- CI.Gi.i

        G1.i <- sort(Gi.star.i)[l1.i]
        G2.i <- sort(Gi.star.i)[l2.i]
        LR2_1.i <- sort(LR2.star.i)[l1.i]
        LR2_2.i <- sort(LR2.star.i)[l2.i]
        if ("Perc" %in% which.CI){
          CI.Gi.i["Perc",] <- c(G1.i,G2.i)
          CI.LR2.i["Perc",] <- c(LR2_1.i,LR2_2.i)
        }
        if ("Basic" %in% which.CI){
          CI.Gi.i["Basic",] <- 2*Gi.hat[i] - c(G2.i,G1.i)
          CI.LR2.i["Basic",] <- 2*LR2[i] - c(LR2_2.i,LR2_1.i)
        }
        if ("Param" %in% which.CI){
          Sigma.Gi.star.i <- n*var(Gi.star.i,na.rm=T)
          CI.Gi.i["Param",] <- Gi.hat[i]+c(-1,1)*sqrt(Sigma.Gi.star.i)/sqrt(n)*qnorm(1-alpha/2)
          Sigma.LR2.star.i <- n*var(LR2.star.i,na.rm=T)
          CI.LR2.i["Param",] <- LR2[i]+c(-1,1)*sqrt(Sigma.LR2.star.i)/sqrt(n)*qnorm(1-alpha/2)
        }

      }else{

        CI.Gi.i <- NA
        CI.LR2.i <- NA

      }

      return.list <- list()
      return.list$CI.Gi.i <- CI.Gi.i
      return.list$CI.LR2.i <- CI.LR2.i

      return(return.list)

    }

    CIs <- lapply(1:L,CI.construction)
    CI.Gi <- lapply(1:L, function(i)CIs[[i]]$CI.Gi.i)
    CI.LR2 <- lapply(1:L, function(i)CIs[[i]]$CI.LR2.i)

  }else{

    Gi.star <- sapply(1:B,function(b)Boot.b[[b]]$Gi.star)
    LR2.star <- sapply(1:B,function(b)Boot.b[[b]]$LR2.star)
    l1 <- floor(alpha/2*B+1)
    l2 <- ceiling((1-alpha/2)*B)
    CI.Gi <- matrix(nrow=length(which.CI),ncol=2)
    colnames(CI.Gi) <- c("Lower.limit","Upper.limit")
    rownames(CI.Gi) <- which.CI
    CI.LR2 <- CI.Gi

    G1 <- sort(Gi.star)[l1]
    G2 <- sort(Gi.star)[l2]
    LR2_1 <- sort(LR2.star)[l1]
    LR2_2 <- sort(LR2.star)[l2]
    if ("Perc" %in% which.CI){
      CI.Gi["Perc",] <- c(G1,G2)
      CI.LR2["Perc",] <- c(LR2_1,LR2_2)
    }
    if ("Basic" %in% which.CI){
      CI.Gi["Basic",] <- 2*Gi.hat - c(G2,G1)
      CI.LR2["Basic",] <- 2*LR2 - c(LR2_2,LR2_1)
    }
    if ("Param" %in% which.CI){
      Sigma.Gi.star <- n*var(Gi.star,na.rm=T)
      CI.Gi["Param",] <- Gi.hat+c(-1,1)*sqrt(Sigma.Gi.star)/sqrt(n)*qnorm(1-alpha/2)
      Sigma.LR2.star <- n*var(LR2.star,na.rm=T)
      CI.LR2["Param",] <- LR2+c(-1,1)*sqrt(Sigma.LR2.star)/sqrt(n)*qnorm(1-alpha/2)
    }

  }

  # RETURN LIST ----

  Return.list <- list()
  Return.list$LR.est <- LR.est
  Return.list$CI.Gi <- CI.Gi
  Return.list$CI.LR2 <- CI.LR2
  Return.list$Gi.star <- Gi.star
  Return.list$LR2.star <- LR2.star

  if(penalty != "none"){
    Return.list$OOB.best <- OOB.best
    Return.list$OOB.total <- OOB.total
  }

  return(Return.list)

}
