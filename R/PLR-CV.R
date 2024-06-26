#' Determines the regularization parameter (lambda) in a PLR via cross-validation
#'
#' \code{PLR.CV} undertakes k-fold cross-validation for a Penalized Lorenz Regression. It returns the CV-score associated to each value of the regularization parameter and the index of the optimum.
#'
#' @param formula A formula object of the form \emph{response} ~ \emph{other_variables}.
#' @param data A data frame containing the variables displayed in the formula.
#' @param penalty penalty used in the Penalized Lorenz Regression. Possible values are "SCAD" (default) or "LASSO".
#' @param h bandwidth of the kernel, determining the smoothness of the approximation of the indicator function.
#' @param SCAD.nfwd optional tuning parameter used if penalty="SCAD". Default value is NULL. The larger the value of this parameter, the sooner the path produced by the SCAD will differ from the path produced by the LASSO.
#' @param PLR.est Output of a call to \code{\link{PLR.wrap}} corresponding to the estimation of the Penalized Lorenz Regression on the full sample. Default value is NULL in which case the estimation on the full sample is run internally.
#' @param standardize Should the variables be standardized before the estimation process? Default value is TRUE.
#' @param weights vector of sample weights. By default, each observation is given the same weight.
#' @param eps Step size in the FABS or SCADFABS algorithm. Default value is 0.005.
#' @param nfolds Number of folds. Default value is 10.
#' @param foldID vector taking value from 1 to nfolds specifying the fold index of each observation. Default value is NULL in which case the folds are defined internally.
#' @param seed.CV Should a specific seed be used in the definition of the folds. Default value is NULL in which case no seed is imposed.
#' @param parallel Whether parallel computing should be used to distribute the \code{nfolds} computations on different CPUs. Either a logical value determining whether parallel computing is used (TRUE) or not (FALSE, the default value). Or a numerical value determining the number of cores to use.
#' @param ... Additional parameters corresponding to arguments passed in \code{\link{Lorenz.SCADFABS}} or \code{\link{Lorenz.FABS}} depending on the argument chosen in penalty.
#'
#' @return A list with two components
#' \describe{
#'    \item{\code{val}}{vector indicating the CV-score for each value of lambda.}
#'    \item{\code{best}}{index where the optimum is attained.}
#' }
#'
#' @seealso \code{\link{Lorenz.Reg}}, \code{\link{PLR.wrap}}, \code{\link{Lorenz.FABS}}, \code{\link{Lorenz.SCADFABS}}
#'
#' @section References:
#' Jacquemain, A., C. Heuchenne, and E. Pircalabelu (2024). A penalised bootstrap estimation procedure for the explained Gini coefficient. \emph{Electronic Journal of Statistics 18(1) 247-300}.
#'
#' @examples
#' YX_mat <- Data.Incomes[,-2]
#' PLR <- PLR.wrap(YX_mat, h = nrow(YX_mat)^(-1/5.5), eps=0.01)
#' PLR.CV(Income ~ ., Data.Incomes, PLR.est = PLR,
#'        h = nrow(Data.Incomes)^(-1/5.5), eps = 0.01, nfolds = 5)
#'
#' @import parallel
#' @import foreach
#' @import doParallel
#'
#' @export

PLR.CV<-function(formula,
                 data,
                 penalty="SCAD",
                 h,
                 SCAD.nfwd=NULL,
                 PLR.est=NULL,
                 standardize=TRUE,
                 weights=NULL,
                 eps,
                 nfolds=10,
                 foldID=NULL,
                 seed.CV=NULL,
                 parallel=FALSE,
                 ...
){

  if(!is.null(foldID)){
    nfolds <- length(unique(foldID))
  }

  # PRE-CV > GETTING YX_mat ----

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

  # PRE-CV > STANDARDIZE X ----
  # No need to do it since it will be dealt with in PLR.wrap

  # PRE-CV > INITIAL EST ----
  if(is.null(PLR.est)) PLR.est <- PLR.wrap(YX_mat, standardize = standardize, weights = weights, h = h, SCAD.nfwd = SCAD.nfwd, penalty = penalty, eps = eps, ...)

  # CV > INNER ----

  if (is.null(foldID)){
    if(!is.null(seed.CV)) set.seed(seed.CV)
    folds <- cut(sample(seq(1,n)),breaks=nfolds,labels=FALSE)
  }else{
    folds <- foldID
  }

  args.ellipsis <- list(...)
  args.wrap <- c(formalArgs(LorenzRegression::PLR.wrap),
                 formalArgs(LorenzRegression::Lorenz.SCADFABS),
                 formalArgs(LorenzRegression::Lorenz.FABS))
  args.list.wrap <- args.ellipsis[names(args.ellipsis)%in%args.wrap]

  CV.inner <- function(k){

    # Construct Test and Validation bootstrap samples
    fold.k <- folds==k
    YX_mat.train <- YX_mat[!fold.k,]
    weights.train <- weights[!fold.k]
    YX_mat.valid <- YX_mat[fold.k,]
    weights.valid <- weights[fold.k]
    # Perform the estimation
    PLR.est.k <- do.call(PLR.wrap,c(list(YX_mat = YX_mat.train, standardize = standardize, weights = weights.train, h=h, SCAD.nfwd=SCAD.nfwd, penalty = penalty, eps = eps, lambda = PLR.est$lambda), args.list.wrap))
    theta.k <- PLR.est.k$theta
    lambda.k <- PLR.est.k$lambda
    # Compute the CV-score
    y.valid <- YX_mat.valid[,1]
    X.valid <- as.matrix(YX_mat.valid[,-1])
    CV.score <- apply(theta.k,2,function(x)Gini.coef(y.valid, X.valid%*%x, na.rm=TRUE, ties.method="mean", weights=weights.valid))
    # With SCAD, the algorithm may stop sooner than in original sample. Hence, the lambda vectors might be different
    if( length(lambda.k) != length(PLR.est$lambda) ){
      diff.lengths <- length(PLR.est$lambda)-length(lambda.k)
      if(diff.lengths > 0) CV.score <- c(CV.score,rep(CV.score[length(CV.score)],diff.lengths))
    }

    return(CV.score)

  }

  # CV > ITERATIONS ----

  if(parallel){
    if(is.numeric(parallel)){
      registerDoParallel(parallel)
    }else{
      numCores <- detectCores()
      registerDoParallel(numCores-1)
    }
    CV.k <- foreach(k=1:nfolds) %dopar% {
      CV.inner(k)
    }
    stopImplicitCluster()
  }else{
    CV.k <- foreach(k=1:nfolds) %do% {
      CV.inner(k)
    }
  }

  # CV > POST-ITER > SCORE ----

  CV.matrix <- t(sapply(1:nfolds,function(k)CV.k[[k]]))
  CV.total <- colMeans(CV.matrix)
  CV.best <- which.max(CV.total)

  # CV > POST-ITER > RETURN ----

  Return.list <- list()
  Return.list$val <- CV.total
  Return.list$best <- CV.best

  return(Return.list)

}
