#' Fits a Lorenz regression
#'
#' \code{Lorenz.Reg} fits a Lorenz regression of a response with respect to several covariates.
#'
#' @param formula An object of class "\code{\link{formula}}" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data An optional data frame, list or environment (or object coercible by \code{\link{as.data.frame}} to a data frame) containing the variables in the model. If not found in \code{data}, the variables are taken from \code{environment(formula)}, typically the environment from which \code{Lorenz.Reg} is called.
#' @param weights An optional vector of sample weights to be used in the fitting process. Should be \code{NULL} or a numeric vector.
#' @param na.action A function which indicates what should happen when the data contain \code{NA}s. The default is set by the \code{na.action} setting of \code{\link{options}}, and is \code{\link{na.fail}} if that is unset. The 'factory-fresh' default is \code{\link{na.omit}}. Another possible value is \code{NULL}, no action. Value \code{\link{na.exclude}} can be useful.
#' @param standardize A logical determining whether the variables should be standardized before the estimation process. Default value is \code{TRUE}.
#' @param penalty A character string specifying the type of penalty on the size of the estimated coefficients of the single-index model.
#' The default value is \code{"none"}, in which case a non-penalized Lorenz regression is fitted using \code{\link{Lorenz.GA}}.
#' Other possible values are \code{"LASSO"} and \code{"SCAD"}, in which case a penalized Lorenz regression is fitted using \code{\link{PLR.wrap}}.
#' @param h.grid Only used if \code{penalty="SCAD"} or \code{penalty="LASSO"}. A vector (grid) of values for the bandwidth of the kernel, determining the smoothness of the approximation of the indicator function. By default, it is set to (0.1,0.2,1,2,5)*n^(-1/5.5), where n is the sample size.
#' @param SCAD.nfwd.grid Only used if \code{penalty="SCAD"}. A vector (grid) of values for the \code{SCAD.nfwd} argument used in the \code{\link{PLR.wrap}} function. Default value is \code{NULL}. If a vector is supplied, the argument \code{h.grid} is modified to only use the first value of the vector.
#' @param eps Only used if \code{penalty="SCAD"} or \code{penalty="LASSO"}. A scalar indicating the step size in the FABS or SCADFABS algorithm. Default value is 0.005.
#' @param ... Additional parameters corresponding to arguments passed in \code{\link{Lorenz.GA}} or \code{\link{PLR.wrap}}, depending on the argument chosen in \code{penalty}.
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
#'    \item{\code{MRS}}{The matrix of estimated marginal rates of substitution. More precisely, if we want the MRS of \eqn{X_1} (numerator) with respect to \eqn{X_2} (denominator),
#'    we should look for row corresponding to \eqn{X_1} and column corresponding to \eqn{X_2}.}
#'    \item{\code{index}}{A vector gathering the estimated index.}
#' }
#' For the Penalized Lorenz Regression, the tuning parameter (i.e. the value on the grid \code{h.grid} or \code{SCAD.nfwd.grid}) and the penalization parameter (lambda) are chosen optimally via the BIC method.
#' The object of class \code{"PLR"} is a list containing the same components as previously, and in addition :
#' \describe{
#'    \item{\code{path}}{A list where the different elements correspond to the values of the tuning parameter. Each element is a matrix where the first line displays the vector of lambda values. The second and third lines display the evolution of the Lorenz-\eqn{R^2} and explained Gini coefficient along that vector. The next lines display the evolution of the BIC score. The remaining lines display the evolution of the estimated coefficients of the single-index model.}
#'    \item{\code{which.lambda}}{the index of the optimal lambda obtained by the BIC method}
#'    \item{\code{which.tuning}}{the index of the optimal tuning parameter obtained by the BIC method.}
#' }
#' In both cases, the list also provides technical information, such as the specified \code{formula}, \code{weights} and \code{call}, as well as the design matrix \code{x} and the response vector \code{y}.
#'
#' @details In the penalized case, the model is fitted for a grid of values of two parameters : the penalty parameter (lambda) and a tuning parameter.
#' By default, the tuning parameter is the bandwidth of the kernel (argument \code{h} in \code{\link{PLR.wrap}}). Its grid is defined by the user via the argument \code{h.grid}.
#' Alternatively, and only if the SCAD penalty is used, the tuning parameter can also be the \eqn{n_{fwd}} parameter described in Jacquemain et al. (argument \code{SCAD.nfwd} in \code{\link{PLR.wrap}}). Its grid is defined by the user via the argument \code{SCAD.nfwd.grid}.
#' To avoid too large computation time, the user must choose between providing a grid for \code{h} or for \code{SCAD.nfwd}. If \code{SCAD.nfwd.grid} is supplied, only the first value of \code{h.grid} is used for the bandwidth.
#' The optimal values of the parameters are selected using the BIC criterion. Other selections methods include bootstrap (obtained with \code{\link{Lorenz.boot}}) and cross-validation (obtained with \code{\link{PLR.CV}}).
#'
#' @seealso \code{\link{Lorenz.GA}}, \code{\link{Lorenz.SCADFABS}}, \code{\link{Lorenz.FABS}}, \code{\link{PLR.wrap}}, \code{\link{Lorenz.boot}}
#'
#' @section References:
#' Heuchenne, C. and A. Jacquemain (2022). Inference for monotone single-index conditional means: A Lorenz regression approach. \emph{Computational Statistics & Data Analysis 167(C)}.
#'
#' Jacquemain, A., C. Heuchenne, and E. Pircalabelu (2024). A penalised bootstrap estimation procedure for the explained Gini coefficient. \emph{Electronic Journal of Statistics 18(1) 247-300}.
#'
#' @examples
#' data(Data.Incomes)
#' set.seed(123)
#' data <- Data.Incomes[sample(1:200,50),]
#' # 1. Non-penalized regression
#' NPLR <- Lorenz.Reg(Income ~ ., data = Data.Incomes, penalty = "none", popSize = 30)
#' # 2. Penalized regression
#' PLR <- Lorenz.Reg(Income ~ ., data = Data.Incomes, penalty = "SCAD",
#'                   eps = 0.02, h.grid=c(0.2,1,2)*nrow(Data.Incomes)^(-1/5.5))
#' # Print method
#' print(NPLR)
#' print(PLR)
#' # Summary method
#' summary(NPLR)
#' summary(PLR)
#' # Coef method
#' coef(NPLR)
#' coef(PLR)
#' # Predict method
#' ## One can predict either the index or the response
#' predict(NPLR,type="response")
#' predict(PLR,type="response")
#' # Plot method
#' plot(NPLR)
#' plot(PLR)
#' ## Traceplot of the penalized coefficients
#' plot(PLR,type="traceplot")
#'
#' @export

Lorenz.Reg <- function(formula,
                       data,
                       weights,
                       na.action,
                       standardize = TRUE,
                       penalty=c("none","SCAD","LASSO"),
                       h.grid=NULL,
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
    x <- as.matrix(rep(1,length(y)))
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
  if(ncol(x)==1){
    penalty <- "none"
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
    if(is.null(h.grid)) h.grid <- c(0.1,0.2,1,2,5)*length(y)^(-1/5.5)
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

  if(!is.null(theta)){
    # Matrix of MRS
    if(penalty == "none"){
      MRS <- outer(theta,theta,"/")
    }else{
      theta_nz <- theta[theta!=0]
      MRS <- outer(theta_nz,theta_nz,"/")
    }
    # Estimated index
    index <- as.vector(theta%*%t(x))
  }else{
    MRS <- index <- NULL
  }

  # Return fitted objects
  return.list$theta <- theta
  return.list$Gi.expl <- Gi.expl
  return.list$LR2 <- LR2
  return.list$MRS <- MRS
  return.list$index <- index

  return(return.list)
}

