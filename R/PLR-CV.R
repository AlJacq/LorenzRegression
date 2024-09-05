#' Cross-validation for penalized Lorenz regression
#'
#' \code{PLR.CV} selects the grid and penalization parameters of the penalized Lorenz regression by cross-validation.
#'
#' @param object An object with S3 class \code{"PLR"}, i.e. the return of a call to the \code{\link{Lorenz.Reg}} function where \code{penalty=="SCAD"} or \code{penalty=="LASSO"}.
#' @param k An integer indicating the number of folds in the k-fold cross-validation
#' @param data.orig A data frame corresponding to the original dataset, used in the \code{\link{Lorenz.Reg}} call.
#' @param seed.CV An optional seed that is used internally for the creation of the folds. Default is \code{NULL}, in which case no seed is imposed.
#' @param parallel Whether parallel computing should be used to distribute the cross-validation computations. Either a logical value determining whether parallel computing is used (TRUE) or not (FALSE, the default value). Or a numerical value determining the number of cores to use.
#' @param ... Additional parameters corresponding to arguments passed to the function \code{\link{vfold_cv}} from the \emph{rsample} library.
#'
#' @return An object of class \code{c("PLR_cv", "PLR")}. The object is a list containing the following components:
#' \describe{
#'    \item{\code{path}}{See the \code{\link{Lorenz.Reg}} function for the documentation of the original path. To this path is added the CV-score.}
#'    \item{\code{lambda.idx}}{A vector indicating the index of the optimal lambda obtained by each selection method.}
#'    \item{\code{grid.idx}}{A vector indicating the index of the optimal grid parameter obtained by each selection method.}
#' }
#' Note: The returned object may have additional classes such as \code{"PLR_boot"} if bootstrap was performed.
#'
#' @details The parameter \code{seed.CV} allows for local seed setting to control randomness in the generation of the folds.
#' The specified seed is applied to the respective part of the computation, and the seed is reverted to its previous state after the operation.
#' This ensures that the seed settings do not interfere with the global random state or other parts of the code.
#'
#' @seealso \code{\link{Lorenz.Reg}}, \code{\link{Lorenz.SCADFABS}}, \code{\link{Lorenz.FABS}}, \code{\link{Lorenz.boot}}
#'
#' @section References:
#' Jacquemain, A., C. Heuchenne, and E. Pircalabelu (2024). A penalised bootstrap estimation procedure for the explained Gini coefficient. \emph{Electronic Journal of Statistics 18(1) 247-300}.
#'
#' @examples
#' \dontshow{
#' utils::example(Lorenz.Reg, echo = FALSE)
#' }
#' # Continuing the  Lorenz.Reg(.) example:
#' PLR_CV <- PLR.CV(PLR, k = 5, data.orig = data, seed.CV = 123)
#' # The object now inherits from the class "PLR_CV".
#' # Hence the methods (also) display the results obtained by cross-validation.
#' print(PLR_CV)
#' summary(PLR_CV)
#' coef(PLR_CV)
#' predict(PLR_CV)
#' plot(PLR_CV)
#' plot(PLR_CV, type = "diagnostic") # Plot of the scores depending on the grid and penalty parameters
#'
#' @importFrom rsample vfold_cv analysis
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom parallel detectCores
#' @importFrom foreach foreach '%do%' '%dopar%'
#' @importFrom stats setNames
#'
#' @export

PLR.CV<-function(object,
                 k,
                 data.orig,
                 seed.CV=NULL,
                 parallel=FALSE,
                 ...
){

  # 0. Checks ----
  if(!inherits(object,"PLR")) stop("object must be the output of a penalized Lorenz regression.")

  # 1. statistic for cv ----
  cv.f <- function(split) {
    train.sample <- analysis(split)
    indices <- split$in_id
    train.call <- object$call
    if(!is.null(object$weights)) train.call$weights <- train.sample$weights_CV
    train.call$data <- quote(train.sample)
    train.call$lambda.list <- lapply(object$path,function(x)x["lambda",])
    train.LR <- eval(train.call)
    # With penalized reg, the algorithm may stop sooner than in the original sample.
    # Therefore the paths would be shorter and the objects would not have the same size
    compare.paths <- function(path.long,path.short){
      lth.diff <- ncol(path.long) - ncol(path.short)
      if(lth.diff > 0) path.short <- cbind(path.short,replicate(lth.diff,path.short[,ncol(path.short)]))
      return(path.short)
    }
    train.LR$path <- lapply(1:length(object$path),function(i)compare.paths(object$path[[i]],train.LR$path[[i]]))
    # Computation of the CV score
    test.x <- object$x[-unique(indices),]
    test.y <- object$y[-unique(indices)]
    if(!is.null(object$weights)){
      test.weights <- object$weights[-unique(indices)]
    }else{
      test.weights <- NULL
    }
    theta.train <- lapply(train.LR$path,function(x)x[(nrow(x)-ncol(object$x)+1):nrow(x),])
    cv.score <- PLR.scores(test.y,test.x,test.weights,theta.train)
    cv.vec <- unlist(cv.score)

    return(cv.vec)

  }

  # 2. cv folds ----

  if (!is.null(object$weights)){
    data.orig$weights_CV <- object$weights
  }
  if (!is.null(seed.CV)) {
    if(exists(".Random.seed")){
      old <- .Random.seed
    }else{
      runif(1)
      old <- .Random.seed
    }
    on.exit( { .Random.seed <<- old } )
    set.seed(seed.CV)
  }
  cv_folds <- vfold_cv(data.orig, v = k, ...)

  # 3. cv computations ----

  j <- NULL

  if(parallel){
    if(is.numeric(parallel)){
      registerDoParallel(parallel)
    }else{
      numCores <- detectCores()
      registerDoParallel(numCores-1)
    }
    cv.j <- foreach(j=1:k) %dopar% {
      cv.f(cv_folds$splits[[j]])
    }
    stopImplicitCluster()
  }else{
    cv.j <- foreach(j=1:k) %do% {
      cv.f(cv_folds$splits[[j]])
    }
  }

  cv_out <- t(sapply(1:k,function(j)cv.j[[j]]))

  # 4. Adding to the PLR object ----
  path.sizes <- sapply(object$path,ncol)
  path.size <- sum(path.sizes)
  lth.path <- length(path.sizes)
  # the cv score is the mean of the cv scores across the folds
  cv_total <- colMeans(cv_out)
  # Adding cv score to the path
  idx <- lapply(1:lth.path,function(i)(cumsum(path.sizes)-path.sizes+1)[i]:cumsum(path.sizes)[i])
  val.cv <- lapply(idx,function(i)cv_total[i])
  lth.theta <- ncol(object$x)
  lth <- nrow(object$path[[1]]) # Same for all anyway (what changes is ncol)
  for (i in 1:lth.path){
    path.tmp <- rbind(object$path[[i]][1:(lth-lth.theta),],
                      "CV score" = val.cv[[i]])
    object$path[[i]] <- rbind(path.tmp,
                              object$path[[i]][(lth-lth.theta+1):lth,])
  }
  lth <- lth + 1
  # Best pair (grid,lambda) in terms of CV score
  path.wl <- unlist(sapply(path.sizes,function(x)1:x))
  path.wt <- rep(1:lth.path,times=path.sizes)
  wl <- path.wl[which.max(cv_total)]
  wt <- path.wt[which.max(cv_total)]
  object$grid.idx <- c(object$grid.idx,"CV"=wt)
  object$lambda.idx <- c(object$lambda.idx,"CV"=wl)
  object$splits <- cv_folds$splits

  class(object) <- c("PLR_cv",class(object))

  return(object)

}
