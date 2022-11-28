#' Plots for the Penalized Lorenz Regression
#'
#' \code{plot.PLR} provides plots for an object of class \code{PLR}.
#'
#' @param x Output of a call to \code{\link{Lorenz.Reg}}, where \code{penalty!="none"}.
#' @param ... Additional arguments.
#'
#' @return Three types of plots
#' The first is the Lorenz curve of the response and concentration curves of the response with respect to the estimated index (obtained with each selection method).
#' In each of the remaining graphs, the horizontal axis is -log(lambda), lambda being the value of the regularization parameter.
#' The second type of plot is a traceplot, where the vertical axis gives the size of the coefficient attached to each covariate.
#' The third type of plot shows the evolution of the score(s) for each of the selection method chosen in the \code{PLR} object.
#' For comparability reasons, the scores are normalized such that  the larger the better and the optimum is attained in 1.
#' Since the whole path depends on the chosen bandwidth for the kernel, and the optimal bandwidth may depend on the selection method, the plots are produced for each selection method used in the PLR object
#'
#' @details For cross-validation and bootstrap, the scores may not cover the whole path.
#' We run the FABS or SCAD-FABS algorithm on the bootstrap resamples (bootstrap) or folds (cross-validaiton) imposing the lambda grid obtained in the original sample.
#' However, the algorithm may stop before the end of the grid. In fact, it stops as soon as the next iteration would increase the unpenalized loss.
#' We dismiss lambda values for which this situation occurs more than 5\% of the time.
#'
#' @seealso \code{\link{Lorenz.Reg}}
#'
#' @examples
#' data(Data.Incomes)
#' PLR <- Lorenz.Reg(Income ~ ., data = Data.Incomes, penalty = "SCAD", lambda.choice = c("BIC","CV"), eps = 0.005, seed.CV = 123, nfolds = 5)
#' plot(PLR)
#'
#' @import ggplot2
#'
#' @method plot PLR
#' @export

plot.PLR <- function(x, ...){

  PLR <- x
  p0 <- Lorenz.graphs(Response ~ ., PLR$Fit, weights = PLR$weights)
  p0 <- p0 + ggtitle("Observed and explained inequality")

  exclude.trace <- c("lambda","Lorenz-R2","Explained Gini","Number of nonzeroes","BIC score","CV score","Boot score")

  # 1. BIC ----

  if ("BIC" %in% names(PLR$which.h)){

    Path <- PLR$path[[PLR$which.h["BIC"]]]

    # Traceplots

    Path.cov <- Path[!(rownames(Path) %in% exclude.trace),]
    names.var <- rownames(Path.cov)
    n.iter <- ncol(Path.cov)
    Plot.data <- data.frame(
      "Variable" = rep(names.var,n.iter),
      "theta" = as.vector(Path.cov),
      "minloglambda" = rep(-log(Path["lambda",]),each=length(names.var))
    )

    p1 <- ggplot2::ggplot(Plot.data) +
      aes(x = minloglambda, y = theta, colour = Variable) +
      geom_line(size = 1L) +
      scale_color_hue() +
      labs(x = expression(paste("-log(", symbol(lambda), ")",sep="")),
           y = expression(symbol(theta)[k]),
           title = "Traceplot - bandwidth selected by BIC") +
      theme_minimal()

    # Evolution of the score

    Path.score <- Path[grep("score",rownames(Path)),,drop=FALSE]
    Path.score <- apply(Path.score,1,function(x)x/max(x))
    if ("BIC score" %in% colnames(Path.score)) Path.score[,"BIC score"] <- 1/Path.score[,"BIC score"]
    Path.score[Path.score == -Inf] <- 0
    Path.score[Path.score == 0] <- NA

    Score.data <- data.frame(
      "Score" = rep(colnames(Path.score),n.iter),
      "value" = as.vector(t(Path.score)),
      "minloglambda" = rep(-log(Path["lambda",]),each=length(colnames(Path.score)))
    )

    p2 <- ggplot2::ggplot(Score.data) +
      aes(x = minloglambda, y = value, colour = Score) +
      geom_line(size = 1L) +
      scale_color_hue() +
      labs(x = expression(paste("-log(", symbol(lambda), ")",sep="")),
           y = "value",
           title = "Evolution of the scores - bandwidth selected by BIC") +
      theme_minimal()

  }

  # 2. Boot ----

  if ("Boot" %in% names(PLR$which.h)){

    Path <- PLR$path[[PLR$which.h["Boot"]]]

    # Traceplots

    Path.cov <- Path[!(rownames(Path) %in% exclude.trace),]
    names.var <- rownames(Path.cov)
    n.iter <- ncol(Path.cov)
    Plot.data <- data.frame(
      "Variable" = rep(names.var,n.iter),
      "theta" = as.vector(Path.cov),
      "minloglambda" = rep(-log(Path["lambda",]),each=length(names.var))
    )

    p3 <- ggplot2::ggplot(Plot.data) +
      aes(x = minloglambda, y = theta, colour = Variable) +
      geom_line(size = 1L) +
      scale_color_hue() +
      labs(x = expression(paste("-log(", symbol(lambda), ")",sep="")),
           y = expression(symbol(theta)[k]),
           title = "Traceplot - bandwidth selected by bootstrap") +
      theme_minimal()

    # Evolution of the score

    Path.score <- Path[grep("score",rownames(Path)),,drop=FALSE]
    Path.score <- apply(Path.score,1,function(x)x/max(x))
    if ("BIC score" %in% colnames(Path.score)) Path.score[,"BIC score"] <- 1/Path.score[,"BIC score"]
    Path.score[Path.score == -Inf] <- 0
    Path.score[Path.score == 0] <- NA

    Score.data <- data.frame(
      "Score" = rep(colnames(Path.score),n.iter),
      "value" = as.vector(t(Path.score)),
      "minloglambda" = rep(-log(Path["lambda",]),each=length(colnames(Path.score)))
    )

    p4 <- ggplot2::ggplot(Score.data) +
      aes(x = minloglambda, y = value, colour = Score) +
      geom_line(size = 1L) +
      scale_color_hue() +
      labs(x = expression(paste("-log(", symbol(lambda), ")",sep="")),
           y = "value",
           title = "Evolution of the scores - bandwidth selected by bootstrap") +
      theme_minimal()

  }

  # 3. CV ----

  if ("CV" %in% names(PLR$which.h)){

    Path <- PLR$path[[PLR$which.h["CV"]]]

    # Traceplots

    Path.cov <- Path[!(rownames(Path) %in% exclude.trace),]
    names.var <- rownames(Path.cov)
    n.iter <- ncol(Path.cov)
    Plot.data <- data.frame(
      "Variable" = rep(names.var,n.iter),
      "theta" = as.vector(Path.cov),
      "minloglambda" = rep(-log(Path["lambda",]),each=length(names.var))
    )

    p5 <- ggplot2::ggplot(Plot.data) +
      aes(x = minloglambda, y = theta, colour = Variable) +
      geom_line(size = 1L) +
      scale_color_hue() +
      labs(x = expression(paste("-log(", symbol(lambda), ")",sep="")),
           y = expression(symbol(theta)[k]),
           title = "Traceplot - bandwidth selected by cross-validation") +
      theme_minimal()

    # Evolution of the score

    Path.score <- Path[grep("score",rownames(Path)),,drop=FALSE]
    Path.score <- apply(Path.score,1,function(x)x/max(x))
    if ("BIC score" %in% colnames(Path.score)) Path.score[,"BIC score"] <- 1/Path.score[,"BIC score"]
    Path.score[Path.score == -Inf] <- 0
    Path.score[Path.score == 0] <- NA

    Score.data <- data.frame(
      "Score" = rep(colnames(Path.score),n.iter),
      "value" = as.vector(t(Path.score)),
      "minloglambda" = rep(-log(Path["lambda",]),each=length(colnames(Path.score)))
    )

    p6 <- ggplot2::ggplot(Score.data) +
      aes(x = minloglambda, y = value, colour = Score) +
      geom_line(size = 1L) +
      scale_color_hue() +
      labs(x = expression(paste("-log(", symbol(lambda), ")",sep="")),
           y = "value",
           title = "Evolution of the scores - bandwidth selected by cross-validation") +
      theme_minimal()

  }

  # 4. Output ----

  p <- list()
  p[[1]] <- p0
  if ("BIC" %in% names(PLR$which.h)) p <- append(p,list(p1,p2))
  if ("Boot" %in% names(PLR$which.h)) p <- append(p,list(p3,p4))
  if ("CV" %in% names(PLR$which.h)) p <- append(p,list(p5,p6))

  p

}
