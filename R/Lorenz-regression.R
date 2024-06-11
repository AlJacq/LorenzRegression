#' Undertakes a Lorenz regression
#'
#' \code{Lorenz.Reg} performs the Lorenz regression of a response with respect to several covariates.
#'
#' @param formula A formula object of the form \emph{response} ~ \emph{other_variables}.
#' @param data A data frame containing the variables displayed in the formula.
#' @param standardize A logical determining whether the variables should be standardized before the estimation process. Default value is \code{TRUE}.
#' @param weights A vector of sample weights. By default, each observation is given the same weight.
#' @param penalty A character string specifying the type of penalty on the coefficients size.
#' If \code{"none"} is chosen, a non-penalized Lorenz regression is computed using function \code{\link{Lorenz.GA}}.
#' If \code{"SCAD"} is chosen, a penalized Lorenz regression with SCAD penalty is computed using function \code{\link{Lorenz.SCADFABS}}.
#' If \code{"LASSO"} is chosen, a penalized Lorenz regression with LASSO penalty is computed using function \code{\link{Lorenz.FABS}}.
#' @param h.grid Only used if \code{penalty="SCAD"} or \code{penalty="LASSO"}. A vector (grid) of values for the bandwidth of the kernel, determining the smoothness of the approximation of the indicator function. Default value is (0.1,0.2,1,2,5)*n^(-1/5.5), where n is the sample size.
#' @param SCAD.nfwd.grid Only used if \code{penalty="SCAD"}. A vector (grid) of values for the \code{SCAD.nfwd} argument used in the \code{\link{PLR.wrap}} function. Default value is \code{NULL}. If a vector is supplied, the argument \code{h.grid} is modified to only use the first value of the vector.
#' @param eps Only used if \code{penalty="SCAD"} or \code{penalty="LASSO"}. A scalar indicating the step size in the FABS or SCADFABS algorithm. Default value is 0.005.
#' @param ... Additional parameters corresponding to arguments passed in \code{\link{Lorenz.GA}}, \code{\link{Lorenz.SCADFABS}} or \code{\link{Lorenz.FABS}} depending on the argument chosen in penalty.
#'
#' @return An object of class \code{"LR"} for the non-penalized Lorenz regression or of class \code{"PLR"} for a penalized Lorenz regression.
#'
#' For both classes, several methods are available. The function \code{summary} is used to summarize the model fits.
#' Information on the coefficient of the single-index model is obtained via \code{coef}.
#' The method \code{predict} is used to predict either the response or the index of the model.
#' A visual representation of explained inequality through Lorenz curves is provided with the method \code{plot}.
#'
#' The object of class \code{"LR"} is a list containing the following components:
#' \describe{
#'    \item{\code{theta}}{The estimated vector of parameters.}
#'    \item{\code{Gi.expl}}{The estimated explained Gini coefficient.}
#'    \item{\code{LR2}}{The Lorenz-\eqn{R^2} of the regression.}
#'    \item{\code{MRS}}{The matrix of estimated marginal rates of substitution. More precisely, if we want the MRS of X1 (numerator) with respect to X2 (denominator),
#'    we should look for row corresponding to X1 and column corresponding to X2.}
#'    \item{\code{index}}{A vector gathering the estimated index.}
#' }
#' For the Penalized Lorenz Regression, the tuning parameter (i.e. the value on the grid \code{h.grid} or \code{SCAD.nfwd.grid}) and the penalization parameter are chosen optimally via the BIC method.
#' The object of class \code{"PLR"} is a list containing the same components as previously, and in addition :
#' \describe{
#'    \item{\code{path}}{A list where the different elements correspond to the values of the tuning parameter. Each element is a matrix where the first line displays the vector of lambda values. The second and third lines display the evolution of the Lorenz-\eqn{R^2} and explained Gini coefficient along that vector. The next lines display the evolution of the BIC score. The remaining lines display the evolution of the estimated parameter vector.}
#'    \item{\code{which.lambda}}{the index of the optimal lambda obtained by the BIC method}
#'    \item{\code{which.tuning}}{the index of the optimal tuning parameter obtained by the BIC method.}
#' }
#' In both cases, the list also provide technical information, such as the specified formula, weights and call, as well as the design matrix \code{x} and the response vector \code{y}.
#'
#'
#' @seealso \code{\link{Lorenz.GA}}, \code{\link{Lorenz.SCADFABS}}, \code{\link{Lorenz.FABS}}, \code{\link{PLR.wrap}}, \code{\link{Lorenz.boot}}
#'
#' @section References:
#' Heuchenne, C. and A. Jacquemain (2022). Inference for monotone single-index conditional means: A Lorenz regression approach. \emph{Computational Statistics & Data Analysis 167(C)}.
#' Jacquemain, A., C. Heuchenne, and E. Pircalabelu (2024). A penalised bootstrap estimation procedure for the explained Gini coefficient. \emph{Electronic Journal of Statistics 18(1) 247-300}.
#'
#' @examples
#' data(Data.Incomes)
#' set.seed(123)
#' Data <- Data.Incomes[sample(1:nrow(Data.Incomes),50),]
#' # 1. Non-penalized regression
#' NPLR <- Lorenz.Reg(Income ~ ., data = Data, penalty = "none",
#'                    popSize = 30)
#' # 2. Penalized regression
#' PLR <- Lorenz.Reg(Income ~ ., data = Data, penalty = "SCAD",
#'                   h.grid = nrow(Data.Incomes)^(-1/5.5),
#'                   sel.choice = c("BIC","CV"), eps = 0.01, nfolds = 5)
#' # Comparison
#' NPLR$theta;PLR$theta
#' NPLR$summary;PLR$summary
#'
#'
#' @export

Lorenz.Reg <- function(formula,
                       data,
                       weights,
                       na.action,
                       standardize = TRUE,
                       penalty=c("none","SCAD","LASSO"),
                       h.grid=c(0.1,0.2,1,2,5)*nrow(data)^(-1/5.5),
                       SCAD.nfwd.grid = NULL,
                       eps=0.005,
                       lambda.list=NULL,
                       ...){

  # 0 > Calls ----
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "weights", "na.action"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")

  # 0 > Checks ----
  # Check on penalty
  penalty <- match.arg(penalty)
  # Check on SCAD.fwd.grid and h.grid
  # Choose either a grid for the bandwidth (via h.grid) or a grid for n_fwd (via SCAD.nfwd.grid), but not both.
  if(length(h.grid)>1 & length(SCAD.nfwd.grid)>1){
    warning("To avoid enormous computation time, the code does not accept a grid for h and nfwd at the same time. As such, only the first value for h.grid is used, while the whole vector is used for SCAD.nfwd.grid")
    h.grid <- h.grid[1]
  }

  # 0 > Response and Design ----

  y <- model.response(mf, "numeric")
  w <- as.vector(model.weights(mf))
  if (!is.null(w) && !is.numeric(w))
    stop("'weights' must be a numeric vector")

  if (is.empty.model(mt)) {
    # In this case, we do not explain anything : explained Gini = Lorenz R2 = 0
  }
  else {
    # We need to distinguish between LR and PLR because specific treatment of categorical in PLR
    if (penalty == "none"){
      # In LR, only need to exclude the intercept
      x <- model.matrix(mt, mf)[,-1,drop=FALSE]
    }else{
      x <- model_matrix_PLR(mt, mf)
    }

  }

  n <- nrow(x)
  p <- ncol(x)

  # 0 > Return ----
  return.list <- list()
  return.list$na.action <- attr(mf, "na.action")
  return.list$call <- cl
  return.list$terms <- mt
  return.list$xlevels <- .getXlevels(mt, mf)
  return.list$weights <- w
  return.list$x <- x
  return.list$y <- y

  # 1. (Penalized) Lorenz Regression ----

  if(penalty == "none"){
    LR <- Lorenz.GA(y, x, standardize=standardize, weights=w, ...)
  }else{
    if(is.null(SCAD.nfwd.grid)|penalty != "SCAD"){
      n.h <- length(h.grid)
      if(is.null(lambda.list)){
        LR <- lapply(1:n.h,function(i)PLR.wrap(y, x, standardize=standardize, weights=w, penalty=penalty, h = h.grid[i], SCAD.nfwd = NULL, eps=eps, ...))
      }else{
        LR <- lapply(1:n.h,function(i)PLR.wrap(y, x, standardize=standardize, weights=w, penalty=penalty, h = h.grid[i], SCAD.nfwd = NULL, eps=eps, lambda = lambda.list[[i]], ...))
      }
    }else{
      n.c <- length(SCAD.nfwd.grid)
      if(is.null(lambda.list)){
        LR <- lapply(1:n.c,function(i)PLR.wrap(y, x, standardize=standardize, weights=w, penalty=penalty, h = h.grid[1], SCAD.nfwd = SCAD.nfwd.grid[i], eps=eps, ...))
      }else{
        LR <- lapply(1:n.c,function(i)PLR.wrap(y, x, standardize=standardize, weights=w, penalty=penalty, h = h.grid[1], SCAD.nfwd = SCAD.nfwd.grid[i], eps=eps, lambda = lambda.list[[i]], ...))
      }
    }
  }

  # 2. Output of the (P)LR ----

  if(penalty == "none"){
    theta <- LR$theta
    Gi.expl <- LR$Gi.expl
    LR2 <- LR$LR2
    class(return.list) <- "LR"
  }else{
    # Construction of the path
    lth.path <- ifelse(!is.null(SCAD.nfwd.grid) & penalty=="SCAD",n.c,n.h)
    # Construction of the path > Number of selected vars
    n_selected <- lapply(1:lth.path,function(i)apply(LR[[i]]$theta,2,function(x)sum(abs(x) > 10^(-10))))
    # Construction of the path > Main objects
    Path <- lapply(1:lth.path,function(i)rbind(LR[[i]]$lambda, LR[[i]]$LR2, LR[[i]]$Gi.expl, n_selected[[i]]))
    for(i in 1:lth.path) rownames(Path[[i]]) <- c("lambda","Lorenz-R2","Explained Gini", "Number of nonzeroes")
    # Construction of the path > BIC score
    Path_BIC <- lapply(1:lth.path,function(i)PLR.BIC(y, x, LR[[i]]$theta, weights = w))
    best.BIC <- lapply(1:lth.path,function(i)Path_BIC[[i]]$best)
    val.BIC <- lapply(1:lth.path,function(i)Path_BIC[[i]]$val)
    for (i in 1:lth.path){
      Path[[i]] <- rbind(Path[[i]], val.BIC[[i]])
      rownames(Path[[i]])[nrow(Path[[i]])] <- "BIC score"
    }
    # Construction of the path > theta's
    for (i in 1:lth.path){
      lth <- nrow(Path[[i]])
      Path[[i]] <- rbind(Path[[i]], LR[[i]]$theta)
      rownames(Path[[i]])[(lth+1):nrow(Path[[i]])] <- colnames(x)
    }
    return.list$path <- Path
    # Optimum tuning params for BIC
    # tuning refers either to h or to SCAD.nfwd
    which.tuning <- which.max(sapply(1:lth.path,function(i)max(val.BIC[[i]])))
    which.lambda <- best.BIC[[which.tuning]]
    return.list$which.tuning <- which.tuning
    return.list$which.lambda <- which.lambda

    theta <- LR[[which.tuning]]$theta[,which.lambda]
    Gi.expl <- LR[[which.tuning]]$Gi.expl[which.lambda]
    LR2 <- LR[[which.tuning]]$LR2[which.lambda]
    class(return.list) <- "PLR"
  }

  names(theta) <- colnames(x)
  # Matrix of MRS
  if(penalty == "none"){
    MRS <- outer(theta,theta,"/")
  }else{
    theta_nz <- theta[theta!=0]
    MRS <- outer(theta_nz,theta_nz,"/")
  }
  # Estimated index
  index <- as.vector(theta%*%t(x))
  # Return fitted objects
  return.list$theta <- theta
  return.list$Gi.expl <- Gi.expl
  return.list$LR2 <- LR2
  return.list$MRS <- MRS
  return.list$index <- index

  return(return.list)
}

