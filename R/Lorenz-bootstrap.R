#' Bootstrap for the (penalized) Lorenz regression
#'
#' \code{Lorenz.boot} determines bootstrap estimators for the vector of coefficients of the single-index model, explained Gini coefficient and Lorenz-\eqn{R^2}. In the penalized case, it also provides a selection method.
#'
#' @param object An object with S3 class \code{"LR"} or \code{"PLR"}, i.e. the return of a call to the \code{\link{Lorenz.Reg}} function.
#' @param R An integer indicating the number of bootstrap replicates.
#' @param data.orig A data frame corresponding to the original dataset, used in the \code{\link{Lorenz.Reg}} call.
#' @param boot_out_only A logical determining whether the function should return raw bootstrap results. This is an advanced feature that helps save computation time in certain instances. See Details.
#' @param ... Additional parameters corresponding to arguments passed to the function \code{\link{boot}} from the \emph{boot} library.
#'
#' @return An object of class \code{c("LR_boot", "LR")} or \code{c("PLR_boot", "PLR")}, depending on whether a non-penalized or penalized regression was fitted.
#'
#' The method \code{confint} is used on an object of class \code{"LR_boot"} or \code{"PLR_boot"} to construct confidence intervals for the model parameters.
#'
#' For the non-penalized Lorenz regression, the returned object is a list containing the following components:
#' \describe{
#'    \item{\code{theta}}{The estimated vector of parameters. In the penalized case, it is a matrix where each row corresponds to a different selection method (e.g., BIC, bootstrap, cross-validation).}
#'    \item{\code{Gi.expl}}{The estimated explained Gini coefficient. In the penalized case, it is a vector, where each element corresponds to a different selection method.}
#'    \item{\code{LR2}}{The Lorenz-\eqn{R^2} of the regression. In the penalized case, it is a vector, where each element corresponds to a different selection method.}
#'    \item{\code{boot_out}}{An object of class \code{"boot"} containing the output of the bootstrap calculation.}
#' }
#' For the penalized Lorenz regression, the returned object is a list containing the following components:
#' \describe{
#'    \item{\code{path}}{See \code{\link{Lorenz.Reg}} for the original path. To this path is added the out-of-bag (OOB) score.}
#'    \item{\code{lambda.idx}}{A vector indicating the index of the optimal lambda obtained by each selection method.}
#'    \item{\code{grid.idx}}{A vector indicating the index of the optimal grid parameter obtained by each selection method.}
#' }
#' Note: in the penalized case, the returned object may have additional classes such as \code{"PLR_cv"} if cross-validation was performed and used as a selection method.
#'
#' @details
#' Users that want to perform parallel computing have two options. The first and most obvious option is to use the facilities provided by the function \code{\link[boot]{boot}}.
#' Indeed, arguments such as \code{parallel}, \code{ncpus} and \code{cl} can be passed through the \code{...}.
#' Alternatively, users might want to run different instances of the function, each taking care of a portion of the bootstrap samples.
#' The argument \code{boot_out_only} can be set to \code{TRUE} to avoid unnecessary computations. If so, the returned object does not inherit from the class \code{"LR_boot"} or \code{"PLR_boot"}. The function simply returns the original \code{object}, to which is added the \code{boot_out} object.
#' If this second option is chosen, the instances have to be combined using the function \code{\link{Lorenz.boot.combine}}.
#'
#' @seealso \code{\link{Lorenz.Reg}}, \code{\link{Lorenz.GA}}, \code{\link{Lorenz.SCADFABS}}, \code{\link{Lorenz.FABS}}, \code{\link{PLR.CV}}, \code{\link[boot]{boot}}
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
#' # Continuing the Lorenz.Reg(.) example for the non-penalized regression:
#' # This example is not run as it takes > 5 seconds to run.
#' \dontrun{
#' set.seed(123)
#' NPLR_boot <- Lorenz.boot(NPLR, R = 40, data.orig = data)
#' # The method confint() is available to objects of class "LR_boot".
#' confint(NPLR_boot)
#' summary(NPLR_boot)
#' }
#' # Continuing the Lorenz.Reg(.) example for the penalized regression:
#' set.seed(123)
#' PLR_boot <- Lorenz.boot(PLR, R = 40, data.orig = data)
#' # The object now inherits from the class "PLR_boot"
#' # Hence the methods (also) display the results obtained by bootstrap.
#' print(PLR_boot)
#' summary(PLR_boot)
#' coef(PLR_boot)
#' predict(PLR_boot)
#' plot(PLR_boot)
#' # Plot of the scores for each selection method depending on the grid and penalty parameters
#' plot(PLR_boot, type = "diagnostic")
#' # The method confint() is available to objects of class "PLR_boot".
#' confint(PLR_boot)
#'
#' @importFrom boot boot
#' @importFrom stats setNames
#'
#' @export

Lorenz.boot <- function(object, R, data.orig, boot_out_only = FALSE, ...){

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
      if(!is.null(object$weights)) boot.call$weights <- object$weights[indices]
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
        theta.boot <- lapply(boot.LR$path,function(x)x[(nrow(x)-ncol(object$x)+1):nrow(x),])
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

  if(!boot_out_only){

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
      lth.theta <- ncol(object$x)
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
      object$grid.idx <- c(object$grid.idx,"Boot"=wt)
      object$lambda.idx <- c(object$lambda.idx,"Boot"=wl)

    }

    # 5. Class ----
    if(method == "LR"){
      class(object) <- c("LR_boot",class(object))
    }else{
      lth.class <- length(class(object))
      if(lth.class==1){
        class(object) <- c("PLR_boot",class(object))
      }else{
        # PLR_boot must come right after "PLR"
        class(object) <- c(class(object)[-lth.class],"PLR_boot","PLR")
      }
    }

  }

  return(object)

}
