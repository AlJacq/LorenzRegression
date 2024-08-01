#' Combines bootstrap Lorenz regressions
#'
#' \code{Lorenz.boot.combine} combine outputs of different instances of the \code{\link{Lorenz.boot}} function.
#'
#' @param boot_list list of objects, each element being the output of a call to the function \code{\link{Lorenz.boot}}.
#'
#' @return An object of class \code{c("LR_boot", "LR")} or \code{c("PLR_boot", "PLR")}, depending on whether a non-penalized or penalized regression was fitted.
#'
#' The method \code{confint} is used on an object of class \code{"LR_boot"} or \code{"PLR_boot"} to obtain bootstrap inference on the model parameters.
#'
#' The object returned is a list that contains the following components:
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
#'    \item{\code{path}}{See the \code{\link{Lorenz.Reg}} function for the original path. To this path is added the out-of-bag (OOB) score.}
#'    \item{\code{lambda.idx}}{A vector indicating the index of the optimal lambda obtained by each selection method.}
#'    \item{\code{grid.idx}}{A vector indicating the index of the optimal grid parameter obtained by each selection method.}
#' }
#' Note: The returned object may have additional classes such as \code{"PLR_cv"} if cross-validation was performed and used as a selection method in the penalized case.
#'
#' @seealso \code{\link{Lorenz.boot}}
#'
#' @section References:
#' Heuchenne, C. and A. Jacquemain (2022). Inference for monotone single-index conditional means: A Lorenz regression approach. \emph{Computational Statistics & Data Analysis 167(C)}.
#'
#' Jacquemain, A., C. Heuchenne, and E. Pircalabelu (2024). A penalised bootstrap estimation procedure for the explained Gini coefficient. \emph{Electronic Journal of Statistics 18(1) 247-300}.
#'
#' @examples
#' \dontshow{
#' utils::example(Lorenz.Reg, echo = FALSE)
#' }
#' # Continuing the Lorenz.Reg(.) example for the penalized regression:
#' boot_list <- list()
#' set.seed(123)
#' boot_list[[1]] <- Lorenz.boot(PLR, R = 20, data.orig = data, boot_out_only = TRUE)
#' set.seed(456)
#' boot_list[[2]] <- Lorenz.boot(PLR, R = 20, data.orig = data, boot_out_only = TRUE)
#' PLR_boot <- Lorenz.boot.combine(boot_list)
#' summary(PLR_boot)
#'
#' @importFrom boot boot
#' @importFrom stats setNames
#'
#' @export

Lorenz.boot.combine <- function(boot_list){

  if(!length(unique(sapply(boot_list,function(x)x$call)))==1)
    stop("The elements of boot_list were not obtained from the same original call.")
  if(!length(unique(lapply(boot_list,function(x)class(x))))==1)
    stop("The elements of boot_list do not have the same class")
  object <- boot_list[[1]]
  if(!inherits(object,c("LR","PLR")))
    stop("object must be the output of a (penalized) Lorenz regression.")
  if(inherits(object,"PLR")){
    method <- "PLR"
  }else{
    method <- "LR"
  }

  boot_out <- do.call(c, lapply(boot_list, function(x) x$boot_out))
  object$boot_out <- boot_out

  # PLR specifics ----
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
    # Best pair (grid,lambda) in terms of OOB score
    path.wl <- unlist(sapply(path.sizes,function(x)1:x))
    path.wt <- rep(1:lth.path,times=path.sizes)
    wl <- path.wl[which.max(OOB_total)]
    wt <- path.wt[which.max(OOB_total)]
    if(length(class(object))==1){
      names(object$Gi.expl) <-
        names(object$LR2) <-
        names(object$lambda.idx) <-
        names(object$grid.idx) <-
        "BIC"
      object$theta <- t(as.matrix(object$theta))
      rownames(object$theta) <- "BIC"
      object$index <- t(as.matrix(object$index))
      rownames(object$index) <- "BIC"
      object$MRS <- list("BIC" = object$MRS)
    }
    object$grid.idx <- c(object$grid.idx,"Boot"=wt)
    object$lambda.idx <- c(object$lambda.idx,"Boot"=wl)
    object$Gi.expl <- setNames(c(object$Gi.expl, object$path[[wt]]["Explained Gini", wl]), c(names(object$Gi.expl), "Boot"))
    object$LR2 <- setNames(c(object$LR2, object$path[[wt]]["Lorenz-R2", wl]), c(names(object$LR2), "Boot"))
    theta.boot <- object$path[[wt]][(lth-lth.theta+1):lth, wl]
    object$theta <- rbind(object$theta, "Boot" = theta.boot)
    theta.boot.nz <- theta.boot[theta.boot!=0]
    object$MRS$Boot <- outer(theta.boot.nz,theta.boot.nz,"/")
    index.boot <- as.vector(theta.boot%*%t(object$x))
    object$index <- rbind(object$index, "Boot" = index.boot)
  }

  # Class ----
  if(method == "LR"){
    class(object) <- c(class(object),"LR_boot")
  }else{
    class(object) <- c(class(object),"PLR_boot")
  }

  return(object)

}
