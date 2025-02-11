#' Bootstrap for the (penalized) Lorenz regression
#'
#' \code{Lorenz.boot} performs bootstrap estimation for the vector of coefficients of the single-index model, the explained Gini coefficient, and the Lorenz-\eqn{R^2}. In the penalized case, it also provides a selection method.
#'
#' @param object An object of class \code{"LR"} or \code{"PLR"}, i.e., the output of a call to \code{\link{Lorenz.Reg}}.
#' @param R An integer specifying the number of bootstrap replicates.
#' @param boot_out_only A logical value indicating whether the function should return only the raw bootstrap output. This advanced feature can help save computation time in specific use cases. See Details.
#' @param ... Additional arguments passed to either the bootstrap function \code{\link[boot]{boot}} from the \pkg{boot} package or the underlying fit functions (\code{\link{Lorenz.GA}}, \code{\link{Lorenz.FABS}}, or \code{\link{Lorenz.SCADFABS}}). By default, the fit function uses the same parameters as in the original call to \code{Lorenz.Reg}, but these can be overridden by explicitly passing them in \code{...}.
#'
#' @return An object of class \code{c("LR_boot", "LR")} or \code{c("PLR_boot", "PLR")}, depending on whether a non-penalized or penalized regression was fitted.
#'
#' The methods \code{\link{confint.LR}} and \code{\link{confint.PLR}} can be used on objects of class \code{"LR_boot"} or \code{"PLR_boot"} to construct confidence intervals for the model parameters.
#'
#' For the non-penalized Lorenz regression, the returned object is a list containing:
#' \describe{
#'    \item{\code{theta}}{The estimated vector of parameters. In the penalized case, this is a matrix where each row corresponds to a different selection method (e.g., BIC, bootstrap, cross-validation).}
#'    \item{\code{Gi.expl}}{The estimated explained Gini coefficient. In the penalized case, this is a vector, where each element corresponds to a different selection method.}
#'    \item{\code{LR2}}{The Lorenz-\eqn{R^2} of the regression. In the penalized case, this is a vector, where each element corresponds to a different selection method.}
#'    \item{\code{boot_out}}{An object of class \code{"boot"} containing the raw bootstrap output.}
#' }
#' For the penalized Lorenz regression, the returned object includes:
#' \describe{
#'    \item{\code{path}}{See \code{\link{Lorenz.Reg}} for the original path. The out-of-bag (OOB) score is added.}
#'    \item{\code{lambda.idx}}{A vector indicating the index of the optimal lambda obtained by each selection method.}
#'    \item{\code{grid.idx}}{A vector indicating the index of the optimal grid parameter obtained by each selection method.}
#' }
#' Note: In the penalized case, the returned object may have additional classes such as \code{"PLR_cv"} if cross-validation was performed and used for selection.
#'
#' @details
#' The function supports parallel computing in two ways:
#' \enumerate{
#'   \item Using the built-in parallelization options of \code{\link[boot]{boot}}, which can be controlled via the \code{...} arguments such as \code{parallel}, \code{ncpus}, and \code{cl}.
#'   \item Running multiple independent instances of \code{Lorenz.boot()}, each handling a subset of the bootstrap samples. In this case, setting \code{boot_out_only = TRUE} ensures that the function only returns the raw bootstrap results. These results can be merged using \code{\link{Lorenz.boot.combine}}.
#' }
#'
#' \strong{Handling of additional arguments (\code{...}):}
#' The function allows for two types of arguments through \code{...}:
#' \itemize{
#'   \item Arguments for \code{\link[boot]{boot}}, used to control the bootstrap procedure.
#'   \item Arguments for the underlying fit functions (\code{\link{Lorenz.GA}}, \code{\link{Lorenz.FABS}}, or \code{\link{Lorenz.SCADFABS}}). By default, the function retrieves these parameters from the original \code{Lorenz.Reg} call. However, users can override them by explicitly specifying new values in \code{...}.
#' }
#'
#' @seealso \code{\link{Lorenz.Reg}}, \code{\link{Lorenz.GA}}, \code{\link{Lorenz.SCADFABS}}, \code{\link{Lorenz.FABS}}, \code{\link{PLR.CV}}, \code{\link[boot]{boot}}
#'
#' @section References:
#' Heuchenne, C. and A. Jacquemain (2022). Inference for monotone single-index conditional means: A Lorenz regression approach. \emph{Computational Statistics & Data Analysis 167(C)}.
#'
#' Jacquemain, A., C. Heuchenne, and E. Pircalabelu (2024). A penalized bootstrap estimation procedure for the explained Gini coefficient. \emph{Electronic Journal of Statistics 18(1) 247-300}.
#'
#' @examples
#' \dontshow{
#' utils::example(Lorenz.Reg, echo = FALSE)
#' }
#' # Non-penalized regression example (not run due to execution time)
#' \dontrun{
#' set.seed(123)
#' NPLR_boot <- Lorenz.boot(NPLR, R = 30)
#' confint(NPLR_boot) # Confidence intervals
#' summary(NPLR_boot)
#' }
#'
#' # Penalized regression example:
#' set.seed(123)
#' PLR_boot <- Lorenz.boot(PLR, R = 20)
#' print(PLR_boot)
#' summary(PLR_boot)
#' coef(PLR_boot, pars.idx = "Boot")
#' predict(PLR_boot, pars.idx = "Boot")
#' plot(PLR_boot)
#' plot(PLR_boot, type = "diagnostic")
#'
#' # Confidence intervals for different selection methods:
#' confint(PLR_boot, pars.idx = "BIC")  # Using BIC-selected tuning parameters
#' confint(PLR_boot, pars.idx = "Boot") # Using bootstrap-selected tuning parameters
#'
#' @importFrom boot boot
#' @importFrom stats setNames
#' @importFrom utils modifyList
#'
#' @export

Lorenz.boot <- function(object, R, boot_out_only = FALSE, ...){

  # 0. Checks ----
  if(!inherits(object,c("LR","PLR"))) stop("object must be the output of a (penalized) Lorenz regression.")

  if(inherits(object,"PLR")){
    method <- "PLR"
  }else{
    method <- "LR"
  }

  args <- list(...)

  # 1. Arguments of the bootstrap ----
  data <- cbind(object$y, object$x)
  boot_args <- args[names(args) %in% names(formals(boot))]

  if(object$penalty == "none"){
    fit_fun <- Lorenz.GA
  }else{
    fit_fun <- PLR.fit
    PLR_args <- list("penalty"=object$penalty, "grid.arg"=object$grid.arg,
                     "grid.value"=object$grid.value, "lambda.list"=object$lambda.list)
  }
  fit_formals <- switch(object$penalty,
                        "none" = names(formals(fit_fun)),
                        "LASSO" = names(formals(Lorenz.FABS)),
                        "SCAD" = names(formals(Lorenz.SCADFABS)))
  fit_args <- modifyList(object$fit_args, args[names(args) %in% fit_formals])
  if(object$penalty != "none") fit_args <- c(fit_args, PLR_args)
  if(object$penalty == "none") fit_args$parallel.GA <- FALSE

  # 2. statistic in boot() ----
  boot.f <- function(data, indices){
    # Construction similar to the "Boot" function in library "car".
    # We want to avoid recomputation on the original sample
    first <- all(indices == seq(length(indices)))
    if(first){
      result <- object
    }else{
      x.boot <- data[indices,-1,drop=FALSE]
      y.boot <- data[indices,1]
      if(!is.null(object$weights)){
        w.boot <- object$weights[indices]
      }else{
        w.boot <- NULL
      }
      if(method == "PLR"){
        x.oob <- data[-unique(indices),-1,drop=FALSE]
        y.oob <- data[-unique(indices),1]
        if(!is.null(object$weights)){
          w.oob <- object$weights[-unique(indices)]
        }else{
          w.oob <- NULL
        }
      }
      boot.LR <- do.call(fit_fun, c(list(y = y.boot, x = x.boot, weights = w.boot), fit_args))
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
        theta.boot <- lapply(boot.LR$path,function(x)x[(nrow(x)-ncol(object$x)+1):nrow(x),])
        OOB.score <- PLR.scores(y.oob,x.oob,w.oob,theta.boot)
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
  boot_out <- do.call(boot, c(list(data = data, statistic = boot.f, R = R), boot_args))
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
