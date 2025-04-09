#' Plots for the penalized Lorenz regression
#'
#' \code{autoplot} generates summary plots for an object of class \code{"PLR"} and returns them as \code{ggplot} objects.
#' The \code{plot} method is a wrapper around \code{autoplot} that directly displays the plot,
#' providing a more familiar interface for users accustomed to base R plotting.
#'
#' @aliases plot.PLR autoplot.PLR_boot plot.PLR_boot autoplot.PLR_cv plot.PLR_cv
#' @param x An object of class \code{"PLR"}. The object might also have S3 classes \code{"PLR_boot"} and/or \code{"PLR_cv"} (both inherit from class \code{"PLR"})
#' @param object An object of class \code{"PLR"}. The object might also have S3 classes \code{"PLR_boot"} and/or \code{"PLR_cv"} (both inherit from class \code{"PLR"})
#' @param type A character string indicating the type of plot. Possible values are \code{"explained"}, \code{"traceplot"} and \code{"diagnostic"}.
#' \itemize{
#' \item If \code{"explained"} is selected, the graph displays the Lorenz curve of the response and concentration curve of the response with respect to the estimated index. The grid and penalty parameters used to estimate the index are chosen via the \code{pars.idx} argument.
#' If \code{object} inherits from \code{"PLR_boot"} and \code{LC_store} was set to \code{TRUE} in \code{\link{Lorenz.boot}}, pointwise confidence intervals for the concentration curve are added. Their confidence level is set via the argument \code{band.level}.
#' \item If \code{"traceplot"} is selected, the graph displays a traceplot, where the horizontal axis is -log(lambda), lambda being the value of the penalty parameter. The vertical axis gives the value of the estimated coefficient attached to each covariate.
#' \item If \code{"diagnostic"} is selected, the graph displays a faceted plot, where each facet corresponds to a different value of the grid parameter. Each plot shows the evolution of the scores of each available selection method. For comparability reasons, the scores are normalized such that the larger the better and the optimum is attained in 1.
#' \item If \code{"residuals"} is selected, the graph displays a scatterplot of residuals with respect to the estimated index.
#' The grid and penalty parameters used for estimation are chosen via the \code{pars.idx} argument.
#' Obtaining residuals entail to estimate the link function of the single-index. This is performed via the function \code{\link{Rearrangement.estimation}}, as explained in \code{\link{predict.LR}}.
#' }
#' @param traceplot.which This argument indicates the value of the grid parameter for which the traceplot should be produced (see arguments \code{grid.value} and \code{grid.arg} in function \code{\link{Lorenz.Reg}}).
#' It can be an integer indicating the index in the grid determined via \code{grid.value}.
#' Alternatively, it can be a character string indicating the selection method. In this case the index corresponds to the optimal value according to that selection method.
#' @param pars.idx What grid and penalty parameters should be used for parameter selection. Either a character string specifying the selection method, where the possible values are:
#' \itemize{
#'    \item \code{"BIC"} (default) - Always available.
#'    \item \code{"Boot"} - Available if \code{object} inherits from \code{"PLR_boot"}.
#'    \item \code{"CV"} - Available if \code{object} inherits from \code{"PLR_cv"}.
#' }
#' Or a numeric vector of length 2, where the first element is the index of the grid parameter and the second is the index of the penalty parameter.
#' @param score.df A data.frame providing the scores to be displayed if \code{type} is set to \code{"diagnostic"}. For internal use only.
#' @param band.level Confidence level for the bootstrap confidence intervals.
#' @param ... Additional arguments passed to function \code{\link{Lorenz.graphs}}
#'
#' @return \code{autoplot} returns a \code{ggplot} object representing the desired graph. \code{plot} directly displays this plot.
#'
#' @details The available selection methods depend on the classes of the object: BIC is always available, bootstrap is available if \code{object} inherits from \code{"PLR_boot"}, cross-validation is available if \code{object} inherits from \code{"PLR_cv"}
#'
#' @seealso \code{\link{Lorenz.Reg}}
#'
#' @examples
#' ## For examples see example(Lorenz.Reg), example(Lorenz.boot) and example(PLR.CV)
#'
#' @importFrom ggplot2 ggplot aes geom_line ggtitle scale_color_hue labs theme_minimal facet_wrap labeller theme autoplot geom_point xlab ylab geom_hline
#' @importFrom stats as.formula na.omit predict
#'
#' @method autoplot PLR
#' @export

autoplot.PLR <- function(object, type = c("explained","traceplot","diagnostic","residuals"), traceplot.which = "BIC", pars.idx = "BIC", score.df = NULL, band.level = 0.95, ...){

  type <- match.arg(type)

  if((is.numeric(traceplot.which) & length(traceplot.which)==1)){
    lth <- length(object$path)
    if(!(traceplot.which %in% 1:lth)) stop("The index in traceplot.which is out of bounds.")
  }else if(traceplot.which == "BIC"){
    traceplot.which <- object$grid.idx["BIC"]
  }else if(traceplot.which == "Boot"){
    stop("object is not of class 'PLR_boot'. Therefore traceplot.which cannot be set to 'Boot'.")
  }else if(traceplot.which == "CV"){
    stop("object is not of class 'PLR_cv'. Therefore traceplot.which cannot be set to 'CV'.")
  }else{
    stop("traceplot.which does not have the correct format")
  }

  if((is.numeric(pars.idx) & length(pars.idx)==2)){
    lth1 <- length(object$path)
    if(pars.idx[1] > lth1) stop("Index of grid parameter is out of bond.")
    lth2 <- ncol(object$path[[pars.idx[1]]])
    if(pars.idx[2] > lth2) stop("Index of lambda parameter is out of bond.")
  }else if(pars.idx == "BIC"){
    pars.idx <- c(object$grid.idx["BIC"],object$lambda.idx["BIC"])
  }else if(pars.idx == "Boot"){
    stop("object is not of class 'PLR_boot'. Therefore pars.idx cannot be set to 'Boot'.")
  }else if(pars.idx == "CV"){
    stop("object is not of class 'PLR_cv'. Therefore pars.idx cannot be set to 'CV'.")
  }else{
    stop("pars.idx does not have the correct format")
  }

  # 1. type = "explained" ----

  if(type == "explained"){

    formula <- as.formula(paste(as.character(object$call$formula[[2]]), "~ ."))
    data <- data.frame(object$y,
                       predict.PLR(object, pars.idx = pars.idx))
    names(data) <- c(all.vars(formula)[1],"index")

    g <- Lorenz.graphs(formula, data, weights = object$weights, ...)
    g <- g + ggtitle("Observed and explained inequality")

  }

  # 2. type = "traceplot" ----

  if (type == "traceplot"){

    lambda <- object$path[[traceplot.which]]["lambda",]
    n.iter <- length(lambda)
    path.theta <- sapply(seq_len(n.iter),function(v)coef.PLR(object, renormalize = FALSE, pars.idx = c(traceplot.which,v)))

    df.long <- data.frame(
      "Variable" = rep(rownames(path.theta),n.iter),
      "theta" = as.vector(path.theta),
      "minloglambda" = rep(-log(lambda),each=nrow(path.theta))
    )

    g <- ggplot(df.long) +
      aes(x = minloglambda, y = theta, colour = Variable) +
      geom_line(linewidth = 1L) +
      labs(x = expression(paste("-log(", symbol(lambda), ")",sep="")),
           y = expression(symbol(theta)[k]),
           title = "Traceplot")

  }

  # 3. type = "diagnostic" ----

  if (type == "diagnostic"){

    if(is.null(score.df)){
      score.df <- do.call(rbind, lapply(1:length(object$path), function(i) {
        data.frame(
          grid = i,
          lambda = -log(object$path[[i]]["lambda",]),
          BIC = object$path[[i]]["BIC score",]
        )
      }))
    }else{
      score.df$BIC <- unlist(sapply(1:length(object$path), function(i) object$path[[i]]["BIC score",]))
    }
    score.df$BIC <- max(score.df$BIC)/score.df$BIC
    scores.only <- score.df[,-(1:2),drop=FALSE]

    df.long <- data.frame(
      grid = rep(score.df$grid, ncol(scores.only)), # Repeat 'grid' column values for each method
      lambda = rep(score.df$lambda, ncol(scores.only)),
      method = rep(names(scores.only),each=nrow(scores.only)), # Create method column
      score = unlist(scores.only, use.names = FALSE) # Combine scores
    )

    if(!is.null(object$grid.value)){
      custom_labels <- paste0(object$call$grid.arg," = ",round(object$grid.value,4))
      names(custom_labels) <- 1:length(object$grid.value)
    }else{
      custom_labels <- ""
      names(custom_labels) <- 1
    }

    lambda <- score <- method <- NULL

    g <- ggplot(df.long, aes(x = lambda, y = score, color = method)) +
      geom_line() +
      facet_wrap(~ grid, scales = "free_x", labeller = labeller(grid = custom_labels)) +
      labs(x = expression(paste("-log(", symbol(lambda), ")",sep="")), y = "Score", color = "Selection method") +
      theme(legend.position = "bottom")

  }

  # 4. Residuals plot ----

  if (type == "residuals"){

    data <- data.frame(index = fitted.PLR(object, pars.idx = pars.idx),
                       resid = residuals.PLR(object, pars.idx = pars.idx))

    g <- ggplot(data) +
      aes(x = index, y = resid) +
      geom_point() +
      geom_hline(yintercept = 0) +
      ggtitle("Residuals vs Estimated index") +
      xlab("Estimated index") + ylab("Residuals")

  }

  # 5. Output ----

  return(g)

}

#' @method autoplot PLR_boot
#' @export

autoplot.PLR_boot <- function(object, type = c("explained","traceplot","diagnostic","residuals"), traceplot.which = "BIC", pars.idx = "BIC", score.df = NULL, band.level = 0.95, ...){

  type <- match.arg(type)

  # 1. type = "explained" ----

  if (type == "explained"){

    if(all(pars.idx == "Boot")) pars.idx <- c(object$grid.idx["Boot"],object$lambda.idx["Boot"])
    g <- NextMethod("autoplot")

    if(object$store_LC){

      if(all(pars.idx == "BIC")) pars.idx <- c(object$grid.idx["BIC"],object$lambda.idx["BIC"])
      path.sizes <- sapply(object$path,ncol)
      path.size <- sum(path.sizes)
      lth.path <- length(path.sizes)
      idx <- lapply(1:lth.path,function(i)(cumsum(path.sizes)-path.sizes+1)[i]:cumsum(path.sizes)[i])
      i <- idx[[pars.idx[1]]][pars.idx[2]]
      LC_path_start <- 3*path.size # bootstrap stores first Gi.expl, LR2 and OOB score
      LC_lth <- 100 # We consider a grid of 100 ordinates
      LC_i_start <- LC_path_start + LC_lth * (i-1)
      LC_ordinates <- object$boot_out$t[,LC_i_start + 1:LC_lth]

      g <- Lorenz.bands(g, LC_ordinates, level = band.level, ...)

    }

  }

  # 2. type = "traceplot" ----

  if (type == "traceplot"){

    if(traceplot.which == "Boot") traceplot.which <- object$grid.idx["Boot"]
    g <- NextMethod("autoplot")

  }

  # 3. type = "diagnostic" ----

  if (type == "diagnostic"){

    if(is.null(score.df)){
      score.df <- do.call(rbind, lapply(1:length(object$path), function(i) {
        data.frame(
          grid = i,
          lambda = -log(object$path[[i]]["lambda",]),
          OOB = object$path[[i]]["OOB score",]
        )
      }))
    }else{
      score.df$OOB <- unlist(sapply(1:length(object$path), function(i) object$path[[i]]["OOB score",]))
    }
    score.df$OOB <- score.df$OOB/max(score.df$OOB)

    g <- NextMethod("autoplot", score.df = score.df)

  }

  # 4. type = "residuals" ----

  if (type == "residuals"){

    if(all(pars.idx == "Boot")) pars.idx <- c(object$grid.idx["Boot"],object$lambda.idx["Boot"])
    g <- NextMethod("autoplot")

  }

  return(g)

}

#' @method autoplot PLR_cv
#' @export

autoplot.PLR_cv <- function(object, type = c("explained","traceplot","diagnostic","residuals"), traceplot.which = "BIC", pars.idx = "BIC", score.df = NULL, band.level = 0.95, ...){

  type <- match.arg(type)

  # 1. type = "explained" ----

  if (type %in% c("explained","residuals")){

    if(all(pars.idx == "CV")) pars.idx <- c(object$grid.idx["CV"],object$lambda.idx["CV"])
    g <- NextMethod("autoplot")

  }

  # 2. type = "traceplot" ----

  if (type == "traceplot"){

    if(traceplot.which == "CV") traceplot.which <- object$grid.idx["CV"]
    g <- NextMethod("autoplot")

  }

  # 3. type = "diagnostic" ----

  if (type == "diagnostic"){

    if(is.null(score.df)){
      score.df <- do.call(rbind, lapply(1:length(object$path), function(i) {
        data.frame(
          grid = i,
          lambda = -log(object$path[[i]]["lambda",]),
          CV = object$path[[i]]["CV score",]
        )
      }))
    }else{
      score.df$CV <- unlist(sapply(1:length(object$path), function(i) object$path[[i]]["CV score",]))
    }
    score.df$CV <- score.df$CV/max(score.df$CV)

    g <- NextMethod("autoplot", score.df = score.df)

  }

  return(g)

}

#' @method plot PLR
#' @rdname autoplot.PLR
#' @export
plot.PLR <- function(x, ...) {
  print(autoplot(x, ...))
}

#' @method plot PLR_boot
#' @export
plot.PLR_boot <- function(x, ...) {
  print(autoplot(x, ...))
}

#' @method plot PLR_cv
#' @export
plot.PLR_cv <- function(x, ...) {
  print(autoplot(x, ...))
}


