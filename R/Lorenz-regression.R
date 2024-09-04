#' Fits a Lorenz regression
#'
#' \code{Lorenz.Reg} fits the Lorenz regression of a response with respect to several covariates.
#'
#' @param formula An object of class "\code{\link{formula}}" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data An optional data frame, list or environment (or object coercible by \code{\link{as.data.frame}} to a data frame) containing the variables in the model. If not found in \code{data}, the variables are taken from \code{environment(formula)}, typically the environment from which \code{Lorenz.Reg} is called.
#' @param weights An optional vector of sample weights to be used in the fitting process. Should be \code{NULL} or a numeric vector.
#' @param na.action A function which indicates what should happen when the data contain \code{NA}s. The default is set by the \code{na.action} setting of \code{\link{options}}, and is \code{\link{na.fail}} if that is unset. The 'factory-fresh' default is \code{\link{na.omit}}. Another possible value is \code{NULL}, no action. Value \code{\link{na.exclude}} can be useful.
#' @param penalty A character string specifying the type of penalty on the size of the estimated coefficients of the single-index model.
#' The default value is \code{"none"}, in which case a non-penalized Lorenz regression is fitted using \code{\link{Lorenz.GA}}.
#' Other possible values are \code{"LASSO"} and \code{"SCAD"}, in which case a penalized Lorenz regression is fitted using \code{\link{Lorenz.FABS}} or \code{\link{Lorenz.SCADFABS}} respectively.
#' @param grid.arg A character string specifying the tuning parameter for which a grid is to be constructed, see Details.
#' @param grid.value A numeric vector specifying the grid values, see Details.
#' @param lambda.list Technical argument used inside the function \code{\link{Lorenz.boot}}.
#' @param ... Additional parameters corresponding to arguments passed in \code{\link{Lorenz.GA}}, \code{\link{Lorenz.FABS}} or \code{\link{Lorenz.SCADFABS}}, depending on the argument chosen in \code{penalty}.
#'
#' @return An object of class \code{"LR"} for the non-penalized Lorenz regression or of class \code{"PLR"} for a penalized Lorenz regression.
#'
#' For both classes, several methods are available. The function \code{summary} is used to summarize the model fits.
#' Information on the coefficient of the single-index model is obtained via \code{coef}.
#' Measures of explained inequality (Gini coefficient and Lorenz-\eqn{R^2}) are extraced with \code{ineqExplained}.
#' The method \code{predict} is used to predict either the response or the index of the model.
#' A visual representation of explained inequality through Lorenz curves is provided with the method \code{plot}.
#'
#' The object of class \code{"LR"} is a list containing the following components:
#' \describe{
#'    \item{\code{theta}}{The estimated vector of parameters.}
#'    \item{\code{Gi.expl}}{The estimated explained Gini coefficient.}
#'    \item{\code{LR2}}{The Lorenz-\eqn{R^2} of the regression.}
#'    \item{\code{index}}{A vector gathering the estimated index.}
#' }
#' The object of class \code{"PLR"} is a list containing the following components:
#' \describe{
#'    \item{\code{path}}{A list where the different elements correspond to the values of the grid parameter. Each element is a matrix where the first line displays the vector of lambda values. The second and third lines display the evolution of the Lorenz-\eqn{R^2} and explained Gini coefficient along that vector. The next lines display the evolution of the BIC score. The remaining lines display the evolution of the estimated coefficients of the single-index model.}
#'    \item{\code{lambda.idx}}{the index of the optimal lambda obtained by the BIC method}
#'    \item{\code{grid.idx}}{the index of the optimal grid parameter obtained by the BIC method.}
#' }
#' In both cases, the list also provides technical information, such as the specified \code{formula}, \code{weights} and \code{call}, as well as the design matrix \code{x} and the response vector \code{y}.
#'
#' @details In the penalized case, the model is fitted for a grid of values of two parameters : the penalty parameter (lambda) and one tuning parameter specified by the arguments \code{grid.arg} and \code{grid.value}.
#' The possibles values for \code{grid.arg} are tuning parameters of the functions \code{\link{Lorenz.FABS}} and \code{\link{Lorenz.SCADFABS}} : \code{''h''} (the default), \code{''SCAD.nfwd''},\code{''eps''}, \code{''kernel''}, \code{''a''} and \code{''gamma''}.
#' The values for the grid are specified with \code{grid.value}. The default is \code{NULL}, in which case no grid is constructed
#'
#' @seealso \code{\link{Lorenz.GA}}, \code{\link{Lorenz.SCADFABS}}, \code{\link{Lorenz.FABS}}, code{\link{Lorenz.boot}}
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
#'                   eps = 0.02, grid.arg = "h",
#'                   grid.value=c(0.5,1,2)*nrow(Data.Incomes)^(-1/5.5))
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
#' @importFrom stats model.response model.weights is.empty.model model.matrix .getXlevels setNames
#'
#'
#' @export

Lorenz.Reg <- function(formula,
                       data,
                       weights,
                       na.action,
                       penalty=c("none","SCAD","LASSO"),
                       grid.arg=c("h","SCAD.nfwd","eps","kernel","a","gamma"),
                       grid.value=NULL,
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

  penalty <- match.arg(penalty)
  grid.arg <- match.arg(grid.arg)

  # 0 > Response and Design ----

  y <- model.response(mf, "numeric")
  w <- as.vector(model.weights(mf))
  if (!is.null(w) && !is.numeric(w))
    stop("'weights' must be a numeric vector")

  if (length(all.vars(mt))==1) {
    x <- as.matrix(rep(1,length(y)))
  } else {
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
    LR <- Lorenz.GA(y, x, weights=w, ...)
  }else{
    if(is.null(grid.value)){
      lth.path <- 1
    }else{
      lth.path <- length(grid.value)
    }
    fun <- switch(penalty,
                  "LASSO" = Lorenz.FABS,
                  "SCAD" = Lorenz.SCADFABS)
    arg.list <- lapply(1:lth.path,function(z)list(y = y, x = x, weights = w))
    for (i in 1:lth.path){
      if(!is.null(lambda.list)) arg.list[[i]]$lambda <- lambda.list[[i]]
      if(!is.null(grid.value)) arg.list[[i]][grid.arg] <- grid.value[i]
    }
    dots <- list(...)
    call.list <- lapply(1:lth.path,function(i)c(arg.list[[i]],dots))
    LR <- lapply(1:lth.path,function(i)do.call(fun,call.list[[i]]))
  }

  # 2. Output of the (P)LR ----

  if(penalty == "none"){
    theta <- LR$theta
    names(theta) <- colnames(x)
    Gi.expl <- LR$Gi.expl
    LR2 <- LR$LR2
    class(return.list) <- "LR"
  }else{
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
    # Optimum grid params for BIC
    # grid refers either to h or to SCAD.nfwd
    grid.idx <- which.max(sapply(1:lth.path,function(i)max(val.BIC[[i]])))
    lambda.idx <- best.BIC[[grid.idx]]
    names(grid.idx) <- names(lambda.idx) <- "BIC"
    return.list$grid.idx <- grid.idx
    return.list$lambda.idx <- lambda.idx
    return.list$grid.value <- grid.value
    class(return.list) <- "PLR"
  }

  # Return fitted objects
  if(penalty == "none"){
    return.list$theta <- theta
    return.list$Gi.expl <- Gi.expl
    return.list$LR2 <- LR2
  }

  return(return.list)
}

