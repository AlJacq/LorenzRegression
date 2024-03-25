#' Produces bootstrap-based inference for (penalized) Lorenz regression
#'
#' \code{Lorenz.boot} determines bootstrap estimators for the weight vector, explained Gini coefficient and Lorenz-\eqn{R^2} and, if applies, selects the regularization parameter.
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
#' @param SCAD.nfwd optional tuning parameter used if penalty="SCAD". Default value is NULL. The larger the value of this parameter, the sooner the path produced by the SCAD will differ from the path produced by the LASSO.
#' @param eps Only used if penalty="SCAD" or penalty="LASSO". Step size in the FABS or SCADFABS algorithm. Default value is 0.005.
#' @param B Number of bootstrap resamples. Default is 500.
#' @param bootID matrix where each row provides the ID of the observations selected in each bootstrap resample. Default is NULL, in which case these are defined internally.
#' @param seed.boot Should a specific seed be used in the definition of the folds. Default value is NULL in which case no seed is imposed.
#' @param parallel Whether parallel computing should be used to distribute the \code{B} computations on different CPUs. Either a logical value determining whether parallel computing is used (TRUE) or not (FALSE, the default value). Or a numerical value determining the number of cores to use.
#' @param ... Additional parameters corresponding to arguments passed in \code{\link{Lorenz.GA}}, \code{\link{Lorenz.SCADFABS}} or \code{\link{Lorenz.FABS}} depending on the argument chosen in penalty.
#'
#' @return A list with several components:
#' \describe{
#'    \item{\code{LR.est}}{Estimation on the original sample.}
#'    \item{\code{Gi.star}}{In the unpenalized case, a vector gathering the bootstrap estimators of the explained Gini coefficient. In the penalized case, it becomes a list of vectors. Each element of the list corresponds to a different value of the penalization parameter}
#'    \item{\code{LR2.star}}{In the unpenalized case, a vector gathering the bootstrap estimators of the Lorenz-\eqn{R^2}. In the penalized case, it becomes a list of vectors.}
#'    \item{\code{theta.star}}{In the unpenalized case, a matrix gathering the bootstrap estimators of theta (rows correspond to bootstrap iterations and columns refer to the different coefficients). In the penalized case, it becomes a list of matrices.}
#'    \item{\code{OOB.total}}{In the penalized case only. Vector gathering the OOB-score for each lambda value.}
#'    \item{\code{OOB.best}}{In the penalized case only. index of the lambda value attaining the highest OOB-score.}
#' }
#'
#' @seealso \code{\link{Lorenz.Reg}}, \code{\link{Lorenz.GA}}, \code{\link{Lorenz.SCADFABS}}, \code{\link{Lorenz.FABS}}, \code{\link{PLR.wrap}}
#'
#' @section References:
#' Heuchenne, C. and A. Jacquemain (2022). Inference for monotone single-index conditional means: A Lorenz regression approach. \emph{Computational Statistics & Data Analysis 167(C)}.
#' Jacquemain, A., C. Heuchenne, and E. Pircalabelu (2024). A penalised bootstrap estimation procedure for the explained Gini coefficient. \emph{Electronic Journal of Statistics 18(1) 247-300}.
#'
#' @examples
#' data(Data.Incomes)
#' set.seed(123)
#' Data <- Data.Incomes[sample(1:nrow(Data.Incomes),50),]
#' Lorenz.boot(Income ~ ., data = Data,
#'             penalty = "SCAD", h = nrow(Data)^(-1/5.5),
#'             eps = 0.02, B = 40, seed.boot = 123)
#'
#'
#' @export

Lorenz.boot<-function(formula,
                      data,
                      standardize=TRUE,
                      weights=NULL,
                      LR.est=NULL,
                      penalty=c("none","SCAD","LASSO"),
                      h=NULL,
                      SCAD.nfwd=NULL,
                      eps=0.005,
                      B = 500,
                      bootID = NULL,
                      seed.boot = NULL,
                      parallel=FALSE,
                      ...
){

  # Check on penalty
  penalty <- match.arg(penalty)

  # Number of bootstrap resamples
  if( !is.null(bootID) ) B <- nrow(bootID)

  # Check on h
  if( is.null(h) & penalty %in% c("SCAD","LASSO") ) stop("h has to be specified in the penalized case")

  # PRE-BOOT > GETTING YX_MAT ----

  # Transform the formula into dataframe
  if (penalty=="none"){
    YX_mat <- as.data.frame(stats::model.matrix(formula,data=data)[,-1])
  }else{

    YX_mat <- as.data.frame(stats::model.matrix(formula,data=data,
                                                contrasts.arg=lapply(data[,sapply(data,is.factor),drop=FALSE],stats::contrasts,contrasts=FALSE))[,-1])

    which.factor <- which(sapply(1:ncol(data),function(i)class(data[,i])=="factor" & length(levels(data[,i]))==2))
    binary.exclude <- sapply(which.factor,function(i)paste0(colnames(data)[i],levels(data[,i])[1]))

    YX_mat <- YX_mat[,!(names(YX_mat) %in% binary.exclude)]


  }
  YX_mat <- cbind(stats::model.frame(formula,data=data)[,1],YX_mat)
  colnames(YX_mat)[1] <- colnames(stats::model.frame(formula,data=data))[1]
  n <- length(YX_mat[,1])
  p <- length(YX_mat[1,])-1

  # PRE-BOOT > (P)LR ----

  if(is.null(LR.est)){
    if(penalty == "none"){
      LR.est <- LorenzRegression::Lorenz.GA(YX_mat, standardize = standardize, weights = weights, parallel = parallel, ...)
    }else{
      LR.est <- LorenzRegression::PLR.wrap(YX_mat, standardize = standardize, weights = weights, h = h, SCAD.nfwd = SCAD.nfwd, penalty = penalty, eps = eps, ...)
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

  args.ellipsis <- list(...)
  args.GA <- formalArgs(LorenzRegression::Lorenz.GA)
  args.wrap <- c(formalArgs(LorenzRegression::PLR.wrap),
                 formalArgs(LorenzRegression::Lorenz.SCADFABS),
                 formalArgs(LorenzRegression::Lorenz.FABS))
  args.list.GA <- args.ellipsis[names(args.ellipsis)%in%args.GA]
  args.list.wrap <- args.ellipsis[names(args.ellipsis)%in%args.wrap]

  Boot.inner <- function(b){

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
      LR.est.star <- do.call(LorenzRegression::Lorenz.GA, c(list(YX_mat = YX_mat.test, standardize = standardize, weights = weights.test, parallel=FALSE), args.list.GA)) # parallel turned to FALSE because we already use parallel computing to distribute the boot iterations
    }else{
      LR.est.star <- do.call(LorenzRegression::PLR.wrap, c(list(YX_mat = YX_mat.test, standardize = standardize, weights = weights.test, penalty = penalty, h = h, SCAD.nfwd = SCAD.nfwd, eps = eps, lambda = LR.est$lambda), args.list.wrap))
      lambda.star <- LR.est.star$lambda
    }
    theta.star <- LR.est.star$theta
    Gi.star <- LR.est.star$Gi.expl
    LR2.star <- LR.est.star$LR2

    # With SCAD, the algorithm may stop sooner than in original sample. Hence, the lambda vectors might be different (the bootstrapped might be shorter)
    if(penalty != "none"){
      diff.lengths <- length(LR.est$lambda)-length(lambda.star)
      if(diff.lengths > 0){
        l.path <- length(Gi.star)
        theta.star <- cbind(theta.star, matrix(replicate(diff.lengths,theta.star[,l.path]),nrow=nrow(theta.star)))
        Gi.star <- c(Gi.star,rep(Gi.star[l.path],diff.lengths))
        LR2.star <- c(LR2.star,rep(LR2.star[l.path],diff.lengths))
      }
    }

    # Compute the OOB-score (if penalized LR)
    if(penalty != "none"){
      y.valid <- YX_mat.valid[,1]
      X.valid <- as.matrix(YX_mat.valid[,-1])
      n.valid <- length(y.valid)
      OOB.score <- apply(theta.star,2,function(x)Gini.coef(y = y.valid, x = X.valid%*%x, na.rm=TRUE, ties.method = "mean", weights = weights.valid))
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
    OOB.total <- colMeans(OOB.matrix)
    OOB.best <- which.max(OOB.total)
    lambda.OOB <- LR.est$lambda[OOB.best]
    Gi.star <- lapply(1:L,function(i)sapply(1:B,function(b)Boot.b[[b]]$Gi.star[i]))
    LR2.star <- lapply(1:L,function(i)sapply(1:B,function(b)Boot.b[[b]]$LR2.star[i]))
    theta.star <- lapply(1:L,function(i)t(sapply(1:B,function(b)Boot.b[[b]]$theta.star[,i])))

  }else{

    Gi.star <- sapply(1:B,function(b)Boot.b[[b]]$Gi.star)
    LR2.star <- sapply(1:B,function(b)Boot.b[[b]]$LR2.star)
    theta.star <- t(sapply(1:B,function(b)Boot.b[[b]]$theta.star))

  }

  # RETURN LIST ----

  Return.list <- list()
  Return.list$LR.est <- LR.est
  Return.list$Gi.star <- Gi.star
  Return.list$LR2.star <- LR2.star
  Return.list$theta.star <- theta.star

  if(penalty != "none"){
    Return.list$OOB.best <- OOB.best
    Return.list$OOB.total <- OOB.total
  }

  return(Return.list)

}
