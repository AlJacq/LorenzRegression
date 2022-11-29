#' Undertakes a Lorenz regression
#'
#' \code{Lorenz.Reg} performs the Lorenz regression of a response with respect to several covariates.
#'
#' @param formula A formula object of the form \emph{response} ~ \emph{other_variables}.
#' @param data A data frame containing the variables displayed in the formula.
#' @param standardize Should the variables be standardized before the estimation process? Default value is TRUE.
#' @param weights vector of sample weights. By default, each observation is given the same weight.
#' @param parallel Whether parallel computing should be used to distribute the computations on different CPUs. Either a logical value determining whether parallel computing is used (TRUE) or not (FALSE, the default value). Or a numerical value determining the number of cores to use.
#' @param penalty should the regression include a penalty on the coefficients size.
#' If "none" is chosen, a non-penalized Lorenz regression is computed using function \code{\link{Lorenz.GA}}.
#' If "SCAD" is chosen, a penalized Lorenz regression with SCAD penalty is computed using function \code{\link{Lorenz.SCADFABS}}.
#' IF "LASSO" is chosen, a penalized Lorenz regression with LASSO penalty is computed using function \code{\link{Lorenz.FABS}}.
#' @param h.grid Only used if penalty="SCAD" or penalty="LASSO". Grid of values for the bandwidth of the kernel, determining the smoothness of the approximation of the indicator function. Default value is (0.1,0.2,1,2,5)*n^(-1/5.5), where n is sample size.
#' @param eps Only used if penalty="SCAD" or penalty="LASSO". Step size in the FABS or SCADFABS algorithm. Default value is 0.005.
#' @param sel.choice Only used if penalty="SCAD" or penalty="LASSO". Determines what method is used to determine the optimal regularization parameter. Possibles values are any subvector of c("BIC","CV","Boot"). Default is "BIC". Notice that "Boot" is necessarily added if Boot.inference is set to TRUE.
#' @param nfolds Only used if sel.choice contains "CV". Number of folds in the cross-validation.
#' @param seed.CV Only used if sel.choice contains "CV". Should a specific seed be used in the definition of the folds. Default value is NULL in which case no seed is imposed.
#' @param foldID vector taking value from 1 to nfolds specifying the fold index of each observation. Default value is NULL in which case the folds are defined internally.
#' @param Boot.inference should bootstrap inference be produced ? Default is FALSE. It is automatically turned to TRUE if sel.choice contains "Boot".
#' @param B Only used if Boot.inference is TRUE. Number of bootstrap resamples. Default is 400.
#' @param bootID Only used if Boot.inference is TRUE. matrix where each row provides the ID of the observations selected in each bootstrap resample. Default is NULL, in which case these are defined internally.
#' @param seed.boot Only used if Boot.inference is TRUE. Should a specific seed be used in the definition of the folds. Default value is NULL in which case no seed is imposed.
#' @param LR Estimation on the original sample. Output of a call to \code{\link{Lorenz.GA}} or \code{\link{PLR.wrap}}.
#' @param LR.boot Estimation on the bootstrap resamples. In the non-penalized case, it is the output of a call to \code{\link{Lorenz.boot}}. In the penalized case, it is a list of size length(h.grid), where each element is the output of a call to \code{\link{Lorenz.boot}} and uses a different value of the bandwidth.
#' @param ... Additional parameters corresponding to arguments passed in \code{\link{Lorenz.GA}}, \code{\link{Lorenz.SCADFABS}} or \code{\link{Lorenz.FABS}} depending on the argument chosen in penalty.
#'
#' @return For the Non-penalized Lorenz Regression, a list with the following elements :
#' \describe{
#'    \item{\code{theta}}{the estimated vector of parameters.}
#'    \item{\code{pval.theta}}{Only returned if Boot.inference is TRUE. the pvalues associated to each element of the parameter vector.}
#'    \item{\code{summary}}{a vector including the estimated explained Gini coefficient and the Lorenz-\eqn{R^2}.}
#'    \item{\code{Gi.expl}}{the estimated explained Gini coefficient}
#'    \item{\code{LR2}}{the Lorenz-\eqn{R^2} of the regression.}
#'    \item{\code{MRS}}{the matrix of estimated marginal rates of substitution. More precisely, if we want the MRS of X1 (numerator) with respect to X2 (denominator),
#'    we should look for row corresponding to X1 and column corresponding to X2.}
#'    \item{\code{Fit}}{A data frame containing the response (first column) and the estimated index (second column).}
#'    \item{\code{Gi.star}}{Only returned if Boot.inference is TRUE. A vector gathering the bootstrap estimators of the explained Gini coefficient.}
#'    \item{\code{LR2.star}}{Only returned if Boot.inference is TRUE. A vector gathering the bootstrap estimators of the Lorenz-\eqn{R^2}.}
#'    \item{\code{theta.star}}{Only returned if Boot.inference is TRUE. A matrix gathering the bootstrap estimators of theta (rows refer to bootstrap iterations and columns refer to the different coefficients)}.
#' }
#' For the Penalized Lorenz Regression, a list with the following elements.
#' \describe{
#'    \item{\code{path}}{a list where the different elements correspond to the values of h.grid. Each element is a matrix where the first line displays the path of regularization parameters. The second and third lines display the evolution of the Lorenz-\eqn{R^2} and explained Gini coefficient along that path. The next lines display the evolution of the scores of the methods chosen in sel.choice. The remaining lines display the evolution of the estimated parameter vector.}
#'    \item{\code{theta}}{a matrix where the different lines correspond to the methods chosen in sel.choice. Each line provides the estimated vector of parameters at the optimal value of the regularization parameter.}
#'    \item{\code{summary}}{a matrix where the different lines correspond to the methods chosen in sel.choice. Each line provides the estimated explained Gini coefficient, the Lorenz-\eqn{R^2}, the optimal lambda, the optimal bandwidth, the number of selected variables and the scores at the optimal value of the regularization parameter.}
#'    \item{\code{Gi.expl}}{a vector providing the estimated explained Gini coefficient at the optimal value of the regularization parameter for each method in sel.choice.}
#'    \item{\code{LR2}}{a vector providing the Lorenz-\eqn{R^2} at the optimal value of the regularization parameter for each method in sel.choice.}
#'    \item{\code{MRS}}{a list where the different elements correspond to a method in sel.choice. Each element is a matrix of estimated marginal rates of substitution for non-zero coefficients at the optimal value of the regularization parameter.}
#'    \item{\code{Fit}}{A data frame containing the response (first column). The remaining columns give the estimated index at the optimal value of the regularization parameter, for each method chosen in sel.choice.}
#'    \item{\code{which.h}}{a vector providing the index of the optimal bandwidth for each method in sel.choice.}
#'    \item{\code{which.lambda}}{a vector providing the index of the optimal lambda for each method in sel.choice.}
#'    \item{\code{Gi.star}}{Only returned if Boot.inference is TRUE. A list (each element a different value of the bandwidth h) of lists (each element a different value of the penalty parameter) of vectors (each element a bootstrap iteration) gathering the bootstrap estimators of the explained Gini coefficient.}
#'    \item{\code{LR2.star}}{Only returned if Boot.inference is TRUE. Similarly for the Lorenz-\eqn{R^2}}
#'    \item{\code{theta.star}}{Only returned if Boot.inference is TRUE. A list (each element a different value of the bandwidth h) of lists (each element a different value of the penalty parameter) of matrices (rows are bootstrap iterations and columns refer to the coefficients) gathering the bootstrap estimators of theta.}
#' }
#' In both cases, the list also technical information, namely the formula, data, weights and call.
#'
#' @seealso \code{\link{Lorenz.GA}}, \code{\link{Lorenz.SCADFABS}}, \code{\link{Lorenz.FABS}}, \code{\link{PLR.wrap}}, \code{\link{Lorenz.boot}}
#'
#' @section References:
#' Heuchenne, C. and A. Jacquemain (2022). Inference for monotone single-index conditional means: A Lorenz regression approach. \emph{Computational Statistics & Data Analysis 167(C)}.
#' Jacquemain, A., C. Heuchenne, and E. Pircalabelu (2022). A penalised bootstrap estimation procedure for the explained Gini coefficient.
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
#' @export

Lorenz.Reg <- function(formula,
                       data,
                       standardize=T,
                       weights=NULL,
                       parallel=F,
                       penalty=c("none","SCAD","LASSO"),
                       h.grid=c(0.1,0.2,1,2,5)*nrow(data)^(-1/5.5),
                       eps=0.005,
                       sel.choice=c("BIC","CV","Boot")[1],
                       nfolds=10,
                       seed.CV=NULL,
                       foldID=NULL,
                       Boot.inference=FALSE,
                       B=400,
                       bootID=NULL,
                       seed.boot=NULL,
                       LR=NULL,
                       LR.boot=NULL,
                       ...){

  # Get the call of the function
  call <- deparse(match.call())

  # Check on penalty
  penalty <- match.arg(penalty)

  # Check on sel.choice
  if( !all(sel.choice%in%c("BIC","CV","Boot")) ) stop("sel.choice should be a subvector of c(\"BIC\",\"CV\",\"Boot\")")

  # Check on bootstrap
  # The idea is that if bootstrap is required, we provide all the results that are available anyway
  if( Boot.inference ) sel.choice <- union(sel.choice,"Boot")
  if ("Boot"%in%sel.choice) Boot.inference <- T

  return.list <- list()

  # 0. Obtain YX_mat ----

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

  # 1. (Penalized) Lorenz Regression ----

  n.h <- length(h.grid)

  if(is.null(LR)){
    if(penalty == "none"){
      LR <- Lorenz.GA(YX_mat, standardize=standardize, weights=weights, parallel=parallel, ...)
    }else{
      LR <- lapply(1:n.h,function(i)PLR.wrap(YX_mat, standardize=standardize, weights=weights, penalty=penalty, h = h.grid[i], eps=eps, ...))
    }
  }

  # 2. Output of the PLR ----

  if(penalty == "none"){

    # Estimation of theta
    theta <- LR$theta # Vector of estimated coefficients
    names(theta) <- colnames(YX_mat[,-1])
    # Summary
    summary <- c()
    summary["Explained Gini"] <- LR$Gi.expl
    summary["Lorenz-R2"] <- LR$LR2
    # Matrix of MRS
    MRS <- outer(theta,theta,"/")
    # Estimated index
    Fit <- data.frame(Response = YX_mat[,1], Index = as.vector(theta%*%t(YX_mat[,-1])))
    # Bootstrap
    if(is.null(LR.boot) & Boot.inference){
      LR.boot <- Lorenz.boot(formula, data, standardize = standardize, weights = weights, LR.est = LR, penalty = penalty, B = B, bootID = bootID, seed.boot = seed.boot, parallel = parallel, ...)
      Sigma.hat.star <- n*stats::var(LR.boot$theta.star)
      pval.theta <- sapply(1:p,function(k)2*stats::pnorm(sqrt(n)*abs(theta[k])/sqrt(Sigma.hat.star[k,k]),lower.tail=FALSE))
      names(pval.theta) <- names(theta)
    }
    # Return
    return.list$theta <- theta
    return.list$summary <- summary
    return.list$Gi.expl <- LR$Gi.expl
    return.list$LR2 <- LR$LR2
    return.list$MRS <- MRS
    return.list$Fit <- Fit
    if(Boot.inference){
      return.list$Gi.star <- LR.boot$Gi.star
      return.list$LR2.star <- LR.boot$LR2.star
      return.list$theta.star <- LR.boot$theta.star
      return.list$pval.theta <- pval.theta
    }

    class(return.list) <- "LR"

  }else{

    # Number of variables selected
    n_selected <- lapply(1:n.h,function(i)apply(LR[[i]]$theta,2,function(x)sum(abs(x) > 10^(-10))))
    # Path
    Path <- lapply(1:n.h,function(i)rbind(LR[[i]]$lambda, LR[[i]]$LR2, LR[[i]]$Gi.expl, n_selected[[i]]))
    for(i in 1:n.h) rownames(Path[[i]]) <- c("lambda","Lorenz-R2","Explained Gini", "Number of nonzeroes")
    # BIC and/or CV
    if ("BIC" %in% sel.choice){
      Path_BIC <- lapply(1:n.h,function(i)PLR.BIC(YX_mat, LR[[i]]$theta, weights = weights))
      best.BIC <- lapply(1:n.h,function(i)Path_BIC[[i]]$best)
      val.BIC <- lapply(1:n.h,function(i)Path_BIC[[i]]$val)
      for (i in 1:n.h){
        Path[[i]] <- rbind(Path[[i]], val.BIC[[i]])
        rownames(Path[[i]])[nrow(Path[[i]])] <- "BIC score"
      }

    }
    if ("CV" %in% sel.choice){
      Path_CV <- lapply(1:n.h,function(i)PLR.CV(formula, data, penalty = penalty, h = h.grid[i], PLR.est = LR[[i]], standardize = standardize, weights = weights, eps = eps, nfolds = nfolds, parallel = parallel, seed.CV = seed.CV, foldID = foldID, ...))
      best.CV <- lapply(1:n.h,function(i)Path_CV[[i]]$best)
      val.CV <- lapply(1:n.h,function(i)Path_CV[[i]]$val)
      for (i in 1:n.h){
        Path[[i]] <- rbind(Path[[i]], val.CV[[i]])
        rownames(Path[[i]])[nrow(Path[[i]])] <- "CV score"
      }
    }
    if ("Boot" %in% sel.choice){
      if (is.null(LR.boot)){
        Path_Boot <- lapply(1:n.h,function(i)Lorenz.boot(formula, data, standardize = standardize, weights = weights, LR.est = LR[[i]], penalty = penalty, h = h.grid[i], eps = eps, alpha = alpha, B = B, bootID = bootID, seed.boot = seed.boot, parallel = parallel, ...))
      }else{
        Path_Boot <- LR.boot
      }
      best.Boot <- lapply(1:n.h,function(i)Path_Boot[[i]]$OOB.best)
      val.Boot <- lapply(1:n.h,function(i)Path_Boot[[i]]$OOB.total)
      # In the bootstrap resamples, the algorithm may stop before the end of the lambda path.
      # We decide to give a score of 0 for lambda values where we have at least 5% of missing values.
      # This is already done for the OOB score in function Lorenz.boot. We do it here for BIC
      if ("BIC" %in% sel.choice){
        val.BIC <- lapply(1:n.h,function(i)ifelse( abs(Path_Boot[[i]]$OOB.total ) > 0, 1, -Inf*sign(val.BIC[[i]]))*val.BIC[[i]])
        best.BIC <- lapply(1:n.h,function(i)which.max(val.BIC[[i]]))
        for (i in 1:n.h){
          Path[[i]]["BIC score",] <- val.BIC[[i]]
        }
      }
      for (i in 1:n.h){
        Path[[i]] <- rbind(Path[[i]], val.Boot[[i]])
        rownames(Path[[i]])[nrow(Path[[i]])] <- "Boot score"
      }
    }

    for (i in 1:n.h){
      lth <- nrow(Path[[i]])
      Path[[i]] <- rbind(Path[[i]], LR[[i]]$theta)
      rownames(Path[[i]])[(lth+1):nrow(Path[[i]])] <- colnames(YX_mat[,-1])
    }

    # Output without bootstrap
    theta <- matrix(nrow = length(sel.choice), ncol = ncol(YX_mat)-1)
    colnames(theta) <- colnames(YX_mat[,-1])
    summary <- matrix(nrow = length(sel.choice), ncol = 5 + length(sel.choice))
    colnames(summary) <- c("Explained Gini", "Lorenz-R2", "lambda", "h", "Number of variables", paste0(sel.choice," score"))
    Gi.expl <- rep(NA,length(sel.choice))
    LR2 <- rep(NA,length(sel.choice))
    MRS <- lapply(1:length(sel.choice),function(x)matrix(nrow = ncol(YX_mat)-1, ncol = ncol(YX_mat)-1))
    Fit <- data.frame(Response = YX_mat[,1])
    names(MRS) <- rownames(theta) <- rownames(summary) <- names(Gi.expl) <- names(LR2) <- sel.choice

    which.h <- c()
    which.lambda <- c()

    if( "BIC" %in% sel.choice ){

      i.BIC <- which(sel.choice == "BIC")
      which.h.BIC <- which.max(sapply(1:n.h,function(i)max(val.BIC[[i]])))
      which.h["BIC"] <- which.h.BIC
      which.lambda["BIC"] <- best.BIC[[which.h.BIC]]
      # theta
      theta[i.BIC,] <- LR[[which.h.BIC]]$theta[,best.BIC[[which.h.BIC]]]
      # summary
      summary[i.BIC,1] <- LR[[which.h.BIC]]$Gi.expl[best.BIC[[which.h.BIC]]]
      summary[i.BIC,2] <- LR[[which.h.BIC]]$LR2[best.BIC[[which.h.BIC]]]
      summary[i.BIC,3] <- LR[[which.h.BIC]]$lambda[best.BIC[[which.h.BIC]]]
      summary[i.BIC,4] <- h.grid[which.h.BIC]
      summary[i.BIC,5] <- n_selected[[which.h.BIC]][best.BIC[[which.h.BIC]]]
      summary[i.BIC,"BIC score"] <- val.BIC[[which.h.BIC]][best.BIC[[which.h.BIC]]]
      if ("Boot" %in% sel.choice) summary[i.BIC,"Boot score"] <- val.Boot[[which.h.BIC]][best.BIC[[which.h.BIC]]]
      if ("CV" %in% sel.choice) summary[i.BIC,"CV score"] <- val.CV[[which.h.BIC]][best.BIC[[which.h.BIC]]]
      # Gi.expl and LR2
      Gi.expl[i.BIC] <- summary[i.BIC,1]
      LR2[i.BIC] <- summary[i.BIC,2]
      # Matrix of MRS
      theta.MRS.BIC <- theta[i.BIC,][theta[i.BIC,]!=0]
      MRS$BIC <- outer(theta.MRS.BIC,theta.MRS.BIC,"/")
      # Fit
      Fit$Index.BIC <- as.vector(theta[i.BIC,]%*%t(YX_mat[,-1]))

    }

    if( "Boot" %in% sel.choice ){

      i.Boot <- which(sel.choice == "Boot")
      which.h.Boot <- which.max(sapply(1:n.h,function(i)max(val.Boot[[i]])))
      which.h["Boot"] <- which.h.Boot
      which.lambda["Boot"] <- best.Boot[[which.h.Boot]]
      # theta
      theta[i.Boot,] <- LR[[which.h.Boot]]$theta[,best.Boot[[which.h.Boot]]]
      # summary
      summary[i.Boot,1] <- LR[[which.h.Boot]]$Gi.expl[best.Boot[[which.h.Boot]]]
      summary[i.Boot,2] <- LR[[which.h.Boot]]$LR2[best.Boot[[which.h.Boot]]]
      summary[i.Boot,3] <- LR[[which.h.Boot]]$lambda[best.Boot[[which.h.Boot]]]
      summary[i.Boot,4] <- h.grid[which.h.Boot]
      summary[i.Boot,5] <- n_selected[[which.h.Boot]][best.Boot[[which.h.Boot]]]
      summary[i.Boot,"Boot score"] <- val.Boot[[which.h.Boot]][best.Boot[[which.h.Boot]]]
      if ("BIC" %in% sel.choice) summary[i.Boot,"BIC score"] <- val.BIC[[which.h.Boot]][best.Boot[[which.h.Boot]]]
      if ("CV" %in% sel.choice) summary[i.Boot,"CV score"] <- val.CV[[which.h.Boot]][best.Boot[[which.h.Boot]]]
      # Gi.expl and LR2
      Gi.expl[i.Boot] <- summary[i.Boot,1]
      LR2[i.Boot] <- summary[i.Boot,2]
      # Matrix of MRS
      theta.MRS.Boot <- theta[i.Boot,][theta[i.Boot,]!=0]
      MRS$Boot <- outer(theta.MRS.Boot,theta.MRS.Boot,"/")
      # Fit
      Fit$Index.Boot <- as.vector(theta[i.Boot,]%*%t(YX_mat[,-1]))

    }

    if( "CV" %in% sel.choice ){

      i.CV <- which(sel.choice == "CV")
      which.h.CV <- which.max(sapply(1:n.h,function(i)max(val.CV[[i]])))
      which.h["CV"] <- which.h.CV
      which.lambda["CV"] <- best.CV[[which.h.CV]]
      # theta
      theta[i.CV,] <- LR[[which.h.CV]]$theta[,best.CV[[which.h.CV]]]
      # summary
      summary[i.CV,1] <- LR[[which.h.CV]]$Gi.expl[best.CV[[which.h.CV]]]
      summary[i.CV,2] <- LR[[which.h.CV]]$LR2[best.CV[[which.h.CV]]]
      summary[i.CV,3] <- LR[[which.h.CV]]$lambda[best.CV[[which.h.CV]]]
      summary[i.CV,4] <- h.grid[which.h.CV]
      summary[i.CV,5] <- n_selected[[which.h.CV]][best.CV[[which.h.CV]]]
      summary[i.CV,"CV score"] <- val.CV[[which.h.CV]][best.CV[[which.h.CV]]]
      if ("BIC" %in% sel.choice) summary[i.CV,"BIC score"] <- val.BIC[[which.h.CV]][best.CV[[which.h.CV]]]
      if ("Boot" %in% sel.choice) summary[i.CV,"Boot score"] <- val.Boot[[which.h.CV]][best.CV[[which.h.CV]]]
      # Gi.expl and LR2
      Gi.expl[i.CV] <- summary[i.CV,1]
      LR2[i.CV] <- summary[i.CV,2]
      # Matrix of MRS
      theta.MRS.CV <- theta[i.CV,][theta[i.CV,]!=0]
      MRS$CV <- outer(theta.MRS.CV,theta.MRS.CV,"/")
      # Fit
      Fit$Index.CV <- as.vector(theta[i.CV,]%*%t(YX_mat[,-1]))

    }

    # Return
    return.list$path <- Path
    return.list$theta <- theta
    return.list$summary <- summary
    return.list$Gi.expl <- Gi.expl
    return.list$LR2 <- LR2
    return.list$MRS <- MRS
    return.list$Fit <- Fit
    return.list$which.h <- which.h
    return.list$which.lambda <- which.lambda

    # Output of the bootstrap

    if (Boot.inference){

      return.list$Gi.star <- lapply(1:n.h, function(i)Path_Boot[[i]]$Gi.star)
      return.list$LR2.star <- lapply(1:n.h, function(i)Path_Boot[[i]]$LR2.star)
      return.list$theta.star <- lapply(1:n.h, function(i)Path_Boot[[i]]$theta.star)

    }

    class(return.list) <- "PLR"

  }

  return.list$data <- data
  return.list$formula <- formula
  return.list$weight <- weights
  return.list$call <- call

  return(return.list)
}

