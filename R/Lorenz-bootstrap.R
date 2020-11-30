#' Bootstrap tests or confidence intervals for Lorenz regression
#'
#' \code{Lorenz.boot} produces inference on the vector of parameters of the Lorenz regression.
#' More specifically, this function can be used to obtain individual confidence intervals for each element of the parameter vector.
#' It can also be used to perform a test of joint significance of several covariates.
#' We refer the user to Heuchenne and Jacquemain (2020) for a rigorous exposition of the procedures.
#'
#' @param formula A formula object of the form \emph{response} ~ \emph{other_variables}.
#' @param data A data frame containing the variables displayed in the formula.
#' @param Lorenz.est Output of a Lorenz.Reg on the original dataset. Default value is NULL, in which case the vector of parameters is estimated internally.
#' @param standardize should the variables be standardized before the estimation process? Default value is TRUE
#' @param testorCI whether we should do a bootstrap test of joint significance or a CI or both. Possibles values are "test", "CI" or "both".
#' @param which.CI vector gathering the methods to construct the CI, among Basic Bootstrap ("Basic"), Percentile Bootstrap ("Perc") and Parametric Bootstrap
#' ("Param"). Default value is c("Basic","Perc","Param"), meaning that the three methods are used.
#' @param which.vars.test names of the variables to consider for the test of joint significance.
#' @param alpha level of the CI and/or tests.
#' @param B the number of bootstrap iterations.
#' @param parallel whether parallel computing should be used. Default value is FALSE.
#'
#' @return A list containing at most the following elements
#' \describe{
#'     \item{\code{Basic.theta}}{matrix with p rows and 2 columns (lower and upper bounds) providing the basic bootstrap CI for theta.}
#'     \item{\code{Param.theta}}{idem for the parametric bootstrap.}
#'     \item{\code{Perc.theta}}{idem for the percentile bootstrap}
#'     \item{\code{Basic.Gi}}{vector providing the lower and upper bounds of the basic bootstrap CI for the explained Gini}
#'     \item{\code{Param.Gi}}{idem for the parametric bootstrap.}
#'     \item{\code{Perc.Gi}}{idem for the percentile bootstrap.}
#'     \item{\code{Param.pvalues}}{p-values indicating the significance of each coefficient, computed using parametric bootstrap.}
#'     \item{\code{pvalue}}{p-value of the required test of joint significance}
#'     \item{\code{U}}{observed value of the test statistic}
#' }
#'
#' @seealso \code{\link{Lorenz.Reg}}
#'
#' @section References:
#' Heuchenne, C. and A. Jacquemain. 2020. “Inference for monotone single-index conditional means: a Lorenz regression approach”
#'
#' @examples
#' # # Example A) takes approx 7.8 minutes to complete.
#' # # Example B) takes approx 3.9 minutes to complete.
#' # # In practice, we advise to set B at least 200
#' # # Example A)
#' # a <- Sys.time()
#' # data(Data.Incomes)
#' # Lorenz.boot(Income ~ ., data = Data.Incomes, B = 50, which.vars.test = c("Age","Seniority"))
#' # b <- Sys.time() - a
#' # b
#' # # Computation time can be sped up using parallel computing
#' # c <- Sys.time()
#' # Lorenz.boot(Income ~ ., data = Data.Incomes, B = 50,
#' #             which.vars.test = c("Age","Seniority"), parallel = T)
#' # d <- Sys.time() - c
#' # d
#'
#' @import parallel
#' @import foreach
#' @import doParallel
#'
#' @export

Lorenz.boot<-function(formula,
                      data,
                      Lorenz.est=NULL,
                      standardize=T,
                      testorCI="both",
                      which.CI=c("Basic","Perc","Param"),
                      which.vars.test,
                      alpha=0.05,
                      B,
                      parallel=F
){

  b <- NULL
  # PRE-BOOT ----
  CI.T <- testorCI!="test"
  test.T <- testorCI!="CI"
  Basic.T <- (testorCI!="test")&("Basic"%in%which.CI)
  Perc.T <- (testorCI!="test")&("Perc"%in%which.CI)
  Param.T <- (testorCI!="test")&("Param"%in%which.CI)

  # PRE-BOOT ---- GETTING YX_mat ----

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

  # PRE-BOOT > STANDARDIZE X ----

  if (standardize){
    X <- YX_mat[,-1]
    X.temp <- stats::setNames(data.frame(matrix(ncol = p, nrow = n)), colnames(X))
    which.disc <- which(apply(X,2,function(x)length(unique(x)))<=2)
    n.disc <- length(which.disc)
    X.scale <- rep(NA,p)

    if(n.disc>0){
      X.temp[,which.disc] <- X[,which.disc]
      X.scale[which.disc] <- 1
    }
    if( (n.disc<p) & (n.disc>0) ){
      if (n.disc<(p-1)){
        X.center <- colMeans(X[,-which.disc])
      }else{
        X.center <- mean(X[,-which.disc])
      }
      X[,-which.disc] <- X[,-which.disc] - rep(X.center, rep.int(n,p-n.disc))
      if (n.disc<(p-1)){
        X.scale[-which.disc] <- sqrt(colSums(X[,-which.disc]^2)/(n-1))
      }else{
        X.scale[-which.disc] <- sqrt(sum(X[,-which.disc]^2)/(n-1))
      }
      X[,-which.disc] <- X[,-which.disc] / rep(X.scale[-which.disc], rep.int(n,p-n.disc))
    }
    if(n.disc==0){
      X.center <- colMeans(X)
      X <- X - rep(X.center, rep.int(n,p))
      X.scale <- sqrt(colSums(X^2)/(n-1))
      X <- X / rep(X.scale, rep.int(n,p))
    }
    YX_mat[,-1] <- X
  }

  # Initial Lorenz regression
  if(is.null(Lorenz.est)){
    Lorenz.est <- LorenzRegression::Lorenz.Reg(formula, data, standardize)
  }
  theta.hat <- Lorenz.est$theta
  LR2 <- Lorenz.est$LorenzR2
  Gi.hat <- Lorenz.est$expl.Gini

  # BOOT > PRE-ITER ----

  # BOOT > PRE-ITER > TEST ----

  if(test.T){

    # We need to create the restricted dataset
    col.0 <- unlist(sapply(which.vars.test,function(t)which(grepl(t,colnames(YX_mat)))))
    YX_mat.0 <- YX_mat[,-col.0]

    #We need to estimate theta, LR2 and H under H0. Then we can compute the residuals under H0
    Reg.0 <- LorenzRegression::Lorenz.GA.cpp(YX_mat.0)
    theta.0 <- Reg.0$sol
    LR2.0 <- Reg.0$LR2
    Data.H.Est.0 <- data.frame(Y=YX_mat.0[,1],Index=t(theta.0%*%t(YX_mat.0[,-1])))
    H.hat.0 <- LorenzRegression::Rearrangement.estimation(Data.H.Est.0$Y,Data.H.Est.0$Index)$H
    Res.0 <- YX_mat.0[,1]-H.hat.0
    Res.0 <- Res.0-mean(Res.0)

  }

  # BOOT > INNER ----

  Boot.inner <- function(b){
    Return.list <- list()
    # BOOT > INNER > CI ----
    if(CI.T){
      YX_mat.b <- YX_mat[sample(1:n,replace=T),]
      Lorenz.est.star <- LorenzRegression::Lorenz.GA.cpp(YX_mat.b)
      theta.hat.star <- Lorenz.est.star$sol
      X.b <- as.matrix(YX_mat.b[,-1])
      Y.b <- YX_mat.b[,1]
      Index.b <- X.b%*%theta.hat.star
      Gi.hat.star <- Gini.coef(Y.b,Index.b)
      Return.list$theta.hat.star <- theta.hat.star
      Return.list$Gi.hat.star <- Gi.hat.star
    }

    # BOOT > INNER > TEST ----
    if(test.T){
      Y.b<-H.hat.0+sample(Res.0,replace=T)
      #Compute the bootstrapped LR2 unrestricted
      YX_mat.b<-cbind(Y.b,YX_mat[,-1])
      LR2.b<-LorenzRegression::Lorenz.GA.cpp(YX_mat.b)$LR2
      #Compute the bootstrapped LR2 under H0
      YX_mat.b0<-cbind(Y.b,YX_mat[,-c(1,col.0)])
      LR2.b0<-LorenzRegression::Lorenz.GA.cpp(YX_mat.b0)$LR2
      #We return the bootstrapped test statistic
      U.b<-LR2.b/LR2.b0
      Return.list$U.b <- U.b
    }

    return(Return.list)

  }

  # BOOT > ITERATIONS ----

  if(parallel){
    numCores <- detectCores()
    registerDoParallel(numCores)
    Boot.b <- foreach(b=1:B) %dopar% {
      Boot.inner(b)
    }
    stopImplicitCluster()
  }else{
    Boot.b <- foreach(b=1:B) %do% {
      Boot.inner(b)
    }
  }
  # BOOT > ITERATIONS > CI ----
  if(CI.T){
    theta.hat.star <- t(sapply(1:B,function(b)Boot.b[[b]]$theta.hat.star))
    if (standardize){ # We need the matrix of estimated theta on the original scale !!
      theta.hat.star.orig <- theta.hat.star/rep(X.scale,rep.int(B,p))
      theta.hat.star.orig <- theta.hat.star.orig/apply(theta.hat.star.orig,1,function(x)sum(abs(x)))
      theta.hat.star <- theta.hat.star.orig
    }
    Gi.hat.star <- sapply(1:B,function(b)Boot.b[[b]]$Gi.hat.star) # For the explained Gini, the scale doesn't matter (the ranks of the index are the same)
  }
  # BOOT > ITERATIONS > TEST ----
  if(test.T) U.b.star <- sapply(1:B,function(b)Boot.b[[b]]$U.b)

  # BOOT > POST-ITER ----

  # BOOT > POST-ITER > CI ----
  if(Basic.T|Perc.T){
    l1<-floor(alpha/2*B+1)
    l2<-ceiling((1-alpha/2)*B)
    if(Basic.T){
      Basic.theta <- matrix(nrow=p,ncol=2)
      rownames(Basic.theta) <- colnames(YX_mat[,-1])
      colnames(Basic.theta) <- c("Lower.limit","Upper.limit")
    }
    if(Perc.T){
      Perc.theta <- matrix(nrow=p,ncol=2)
      rownames(Perc.theta) <- colnames(YX_mat[,-1])
      colnames(Perc.theta) <- c("Lower.limit","Upper.limit")
    }

    # BOOT > POST-ITER > CI > THETA ----
    for(k in 1:p){
      theta.hat.star.k <- sort(theta.hat.star[,k])
      q1.k<-theta.hat.star.k[l1]
      q2.k<-theta.hat.star.k[l2]
      # BOOT > POST-ITER > CI > BASIC ----
      if(Basic.T) Basic.theta[k,] <- 2*theta.hat[k] - c(q2.k,q1.k)
      # BOOT > POST-ITER > CI > PERC ----
      if(Perc.T) Perc.theta[k,] <- c(q1.k,q2.k)
    }
    # BOOT > POST-ITER > CI > GINI
    G1.k <- sort(Gi.hat.star)[l1]
    G2.k <- sort(Gi.hat.star)[l2]
    # BOOT > POST-ITER > CI > GINI > BASIC ----
    if(Basic.T){
      Basic.Gi <- 2*Gi.hat - c(G2.k,G1.k)
      names(Basic.Gi) <- c("Lower.limit","Upper.limit")
    }
    # BOOT > POST-ITER > CI > GINI > PERC ----
    if(Perc.T){
      Perc.Gi <- c(G1.k,G2.k)
      names(Perc.Gi) <- c("Lower.limit","Upper.limit")
    }
  }
  # BOOT > POST-ITER > CI > PARAM ----
  if(Param.T){
    # BOOT > POST-ITER > CI > PARAM > THETA ----
    Sigma.hat.star <- n*stats::var(theta.hat.star)
    # if (standardize) Sigma.hat.star.orig <- n*var(theta.hat.star.orig)
    Param.theta <- matrix(nrow=p,ncol=2)
    rownames(Param.theta) <- colnames(YX_mat[,-1])
    colnames(Param.theta) <- c("Lower.limit","Upper.limit")
    Param.pval <- c()
    for(k in 1:p){
      # Computation of the individual p-values
      Param.pval[k] <- 2*stats::pnorm(sqrt(n)*abs(theta.hat[k])/sqrt(Sigma.hat.star[k,k]),lower.tail=FALSE)
      # Computation of the CI per se
      Param.theta[k,] <- theta.hat[k]+c(-1,1)*sqrt(Sigma.hat.star[k,k])/sqrt(n)*stats::qnorm(1-alpha/2)
    }
    names(Param.pval) <- colnames(YX_mat[,-1])
    # BOOT > POST-ITER > CI > PARAM > GINI ----
    Sigma.Gi.star <- n*stats::var(Gi.hat.star)
    Param.Gi <- Gi.hat+c(-1,1)*sqrt(Sigma.Gi.star)/sqrt(n)*stats::qnorm(1-alpha/2)
    names(Param.Gi) <- c("Lower.limit","Upper.limit")
  }
  # BOOT > POST-ITER > TEST ----
  if(test.T){
    U <- LR2/LR2.0
    p.hat <- (sum(U.b.star>=U)+1)/(B+1)
  }

  # RETURN LIST ----

  Return.list <- list()

  if(Basic.T){
    Return.list$Basic.theta <- Basic.theta
    Return.list$Basic.Gi <- Basic.Gi
  }

  if(Perc.T){
    Return.list$Perc.theta <- Perc.theta
    Return.list$Perc.Gi <- Perc.Gi
  }

  if(Param.T){
    Return.list$Param.theta <- Param.theta
    Return.list$Param.Gi <- Param.Gi
    Return.list$Param.pvalues <- Param.pval
  }

  if(test.T){
    Return.list$pvalue <- p.hat
    Return.list$U <- U
  }

  return(Return.list)

}
