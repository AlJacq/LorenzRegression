#' Bootstrap for (penalized) Lorenz regression
#'
#' \code{Lorenz.boot} determines bootstrap estimators for the weight vector, explained Gini coefficient and Lorenz-\eqn{R^2} and, if applies, selects the regularization parameter.
#'
#' @param object An object with S3 class \code{"LR"} or \code{"PLR"}, i.e. the return of a call to the \code{\link{Lorenz.Reg}} function.
#' @param R An integer indicating the number of bootstrap replicates.
#' @param data.orig A data frame corresponding to the original dataset, used in the \code{\link{Lorenz.Reg}} call.
#' @param ... Additional parameters corresponding to arguments passed to the function \code{\link{boot}} from the \code{\link{boot}} library.
#'
#' @return An object of class \code{c("LR_boot", "LR")} or \code{c("PLR_boot", "PLR")}, depending on whether a non-penalized or penalized regression was fitted. The object contains:
#' \describe{
#'    \item{\code{theta}}{The estimated vector of parameters. In the penalized case, it is a matrix where each row corresponds to a different selection method (e.g., BIC, bootstrap, cross-validation).}
#'    \item{\code{Gi.expl}}{The estimated explained Gini coefficient. In the penalized case, it is a vector, where each element corresponds to a different selection method.}
#'    \item{\code{LR2}}{The Lorenz-\eqn{R^2} of the regression. In the penalized case, it is a vector, where each element corresponds to a different selection method.}
#'    \item{\code{MRS}}{The matrix of estimated marginal rates of substitution. In the penalized case, it is a list where each element corresponds to a different selection method.}
#'    \item{\code{index}}{The estimated index. In the penalized case, it is a matrix where each row corresponds to a different selection method.}
#'    \item{\code{boot_out}}{An object of class \code{"boot"} containing the output of the bootstrap calculation.}
#' }
#' For the Penalized Lorenz Regression, the list also contains the following elements.
#' \describe{
#'    \item{\code{path}}{See the \code{\link{Lorenz.Reg}} function for the original path. To this path is added the OOB-score.}
#'    \item{\code{which.lambda}}{A vector indicating the index of the optimal lambda obtained by each selection method.}
#'    \item{\code{which.tuning}}{A vector indicating the index of the optimal tuning parameter obtained by each selection method.}
#' }
#' Note: The returned object may have additional classes such as \code{"PLR_cv"} if cross-validation was performed and used as a selection method in the penalized case.
#'
#' @seealso \code{\link{Lorenz.Reg}}, \code{\link{Lorenz.GA}}, \code{\link{Lorenz.SCADFABS}}, \code{\link{Lorenz.FABS}}, \code{\link{PLR.wrap}}, \code{\link{PLR.CV}}, \code{\link[boot]{boot}}
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
#' @importFrom boot boot
#'
#' @export

Lorenz.boot <- function(object, R, data.orig, ...){

  # 0. Checks ----
  if(!inherits(object,c("LR","PLR"))) stop("object must be the output of a (penalized) Lorenz regression.")

  if(inherits(object,"PLR")){
    method <- "PLR"
  }else{
    method <- "LR"
  }

  # 1. statistic in boot() ----
  boot.f <- function(data, indices){

    # Construction similar to the "Boot" function in library "car".
    # We want to avoid recomputation on the original sample
    first <- all(indices == seq(length(indices)))
    if(first){
      result <- object
    }else{
      boot.sample <- data[indices, ]
      boot.call <- object$call
      boot.call$data <- quote(boot.sample)
      if(method == "LR") boot.call$parallel.GA <- quote(FALSE) # parallel will be used for bootstrap
      if(method == "PLR") boot.call$lambda.list <- lapply(object$path,function(x)x["lambda",])
      boot.LR <- eval(boot.call)
      if(method == "PLR"){
        # With penalized reg, the algorithm may stop sooner than in the original sample.
        # Therefore the paths would be shorter and the objects would not have the same size
        compare.paths <- function(path.long,path.short){
          lth.diff <- ncol(path.long) - ncol(path.short)
          if(lth.diff > 0) path.short <- cbind(path.short,replicate(lth.diff,path.short[,ncol(path.short)]))
          return(path.short)
        }
        boot.LR$path <- lapply(1:length(object$path),function(i)compare.paths(object$path[[i]],boot.LR$path[[i]]))
        # Computation of the OOB score
        OOB.x <- object$x[-unique(indices),]
        OOB.y <- object$y[-unique(indices)]
        if(!is.null(object$weights)){
          OOB.weights <- object$weights[-unique(indices)]
        }else{
          OOB.weights <- NULL
        }
        theta.boot <- lapply(boot.LR$path,function(x)x[(nrow(x)-length(boot.LR$theta)+1):nrow(x),])
        OOB.score <- PLR.scores(OOB.y,OOB.x,OOB.weights,theta.boot)
      }
      result <- boot.LR
    }

    # All objects that require bootstrapping are stacked in a vector
    if(method == "LR"){
      boot.vec <- c("Gi.expl"=result$Gi.expl,
                    "LR2"=result$LR2,
                    result$theta)
    }else{
      Gi.vec <- unlist(sapply(result$path,function(x)x["Explained Gini",]))
      LR2.vec <- unlist(sapply(result$path,function(x)x["Lorenz-R2",]))
      if(first){
        OOB.vec <- rep(0,length(Gi.vec))
      }else{
        OOB.vec <- unlist(OOB.score)
      }
      boot.vec <- c(Gi.vec,LR2.vec,OOB.vec)
    }

    return(boot.vec)

  }

  # 3. boot() ----
  boot_out <- boot(data = data.orig, statistic = boot.f, R = R, ...)
  object$boot_out <- boot_out

  # 4. PLR specifics ----
  if(method == "PLR"){

    path.sizes <- sapply(object$path,ncol)
    path.size <- sum(path.sizes)
    lth.path <- length(path.sizes)
    # the OOB score is the mean of the OOB scores across the bootstrap samples
    OOB_matrix <- boot_out$t[,(ncol(boot_out$t)-path.size+1):ncol(boot_out$t)]
    OOB_total <- colMeans(OOB_matrix)
    # Adding OOB score to the path
    idx <- lapply(1:lth.path,function(i)(cumsum(path.sizes)-path.sizes+1)[i]:cumsum(path.sizes)[i])
    val.OOB <- lapply(idx,function(i)OOB_total[i])
    lth.theta <- ifelse(is.vector(object$theta),length(object$theta),ncol(object$theta))
    lth <- nrow(object$path[[1]]) # Same for all anyway (what changes is ncol)
    for (i in 1:lth.path){
      path.tmp <- rbind(object$path[[i]][1:(lth-lth.theta),],
                        "OOB score" = val.OOB[[i]])
      object$path[[i]] <- rbind(path.tmp,
                                object$path[[i]][(lth-lth.theta+1):lth,])
    }
    lth <- lth + 1
    # Best pair (tuning,lambda) in terms of OOB score
    path.wl <- unlist(sapply(path.sizes,function(x)1:x))
    path.wt <- rep(1:lth.path,times=path.sizes)
    wl <- path.wl[which.max(OOB_total)]
    wt <- path.wt[which.max(OOB_total)]
    if(length(class(object))==1){
      names(object$Gi.expl) <-
        names(object$LR2) <-
        names(object$which.lambda) <-
        names(object$which.tuning) <-
        "BIC"
      object$theta <- t(as.matrix(object$theta))
      rownames(object$theta) <- "BIC"
      object$index <- t(as.matrix(object$index))
      rownames(object$index) <- "BIC"
      object$MRS <- list("BIC" = object$MRS)
    }
    object$which.tuning <- c(object$which.tuning,"Boot"=wt)
    object$which.lambda <- c(object$which.lambda,"Boot"=wl)
    object$Gi.expl <- setNames(c(object$Gi.expl, object$path[[wt]]["Explained Gini", wl]), c(names(object$Gi.expl), "Boot"))
    object$LR2 <- setNames(c(object$LR2, object$path[[wt]]["Lorenz-R2", wl]), c(names(object$LR2), "Boot"))
    theta.boot <- object$path[[wt]][(lth-lth.theta+1):lth, wl]
    object$theta <- rbind(object$theta, "Boot" = theta.boot)
    theta.boot.nz <- theta.boot[theta.boot!=0]
    object$MRS$Boot <- outer(theta.boot.nz,theta.boot.nz,"/")
    index.boot <- as.vector(theta.boot%*%t(object$x))
    object$index <- rbind(object$index, "Boot" = index.boot)
  }

  # 5. Return ----
  if(method == "LR"){
    class(object) <- c(class(object),"LR_boot")
  }else{
    class(object) <- c(class(object),"PLR_boot")
  }

  return(object)

}
