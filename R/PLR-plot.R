#' Plots for the Penalized Lorenz Regression
#'
#' \code{autoplot} generates summary plots for an object of class \code{"PLR"} and returns them as \code{ggplot} objects.
#' The \code{plot} method is a wrapper around \code{autoplot} that directly displays the plot,
#' providing a more familiar interface for users accustomed to base R plotting.
#'
#' @aliases plot.PLR autoplot.PLR_boot plot.PLR_boot autoplot.PLR_cv plot.PLR_cv
#' @param x An object of class \code{"PLR"}. The object might also have S3 classes \code{"PLR_boot"} and/or \code{"PLR_cv"} (both inherit from class \code{"PLR"})
#' @param type A character string indicating the type of plot. Possible values are \code{"explained"}, \code{"traceplot"} and \code{"diagnostic"}.
#' \itemize{
#' \item If \code{"explained"} is selected, the graph displays the Lorenz curve of the response and concentration curve(s) of the response with respect to the estimated index. More specifically, there is one concentration curve per selection method available.
#' \item If \code{"traceplot"} is selected, the graph displays a traceplot, where the horizontal axis is -log(lambda), lambda being the value of the penalty parameter. The vertical axis gives the value of the estimated coefficient attached to each covariate.
#' \item If \code{"diagnostic"} is selected, the graph displays a faceted plot, where each facet corresponds to a different value of the grid parameter. Each plot shows the evolution of the scores of each available selection method. For comparability reasons, the scores are normalized such that the larger the better and the optimum is attained in 1.
#' }
#' @param traceplot.which This argument indicates the value of the grid parameter for which the traceplot should be produced (see arguments \code{grid.value} and \code{grid.arg} in function \code{\link{Lorenz.Reg}}).
#' It can be an integer indicating the index in the grid determined via \code{grid.value}.
#' Alternatively, it can be a character string indicating the selection method. In this case the index corresponds to the optimal value according to that selection method.
#' @param score.df A data.frame providing the scores to be displayed if \code{type} is set to \code{"diagnostic"}. For internal use only.
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
#' @importFrom ggplot2 ggplot aes geom_line ggtitle scale_color_hue labs theme_minimal facet_wrap labeller theme autoplot
#' @importFrom stats as.formula na.omit
#'
#' @method autoplot PLR
#' @export

autoplot.PLR <- function(x, type = c("explained","traceplot","diagnostic"), traceplot.which = "BIC", score.df = NULL, ...){

  type <- match.arg(type)

  if((is.numeric(traceplot.which) & length(traceplot.which)==1)){
    lth <- length(x$path)
    if(!(traceplot.which %in% 1:lth)) stop("The index in traceplot.which is out of bounds.")
  }else if(traceplot.which == "BIC"){
    traceplot.which <- x$grid.idx["BIC"]
  }else if(traceplot.which == "Boot"){
    stop("x is not of class 'PLR_boot'. Therefore traceplot.which cannot be set to 'Boot'.")
  }else if(traceplot.which == "CV"){
    stop("x is not of class 'PLR_cv'. Therefore traceplot.which cannot be set to 'CV'.")
  }else{
    stop("traceplot.which does not have the correct format")
  }

  # 1. type = "explained" ----

  if(type == "explained"){

    formula <- as.formula(paste(as.character(x$call$formula[[2]]), "~ ."))
    data <- data.frame(x$y,predict.PLR(x))
    names(data) <- c(all.vars(formula)[1],"index (BIC)")

    g <- Lorenz.graphs(formula, data, weights = x$weights, ...)
    g <- g + ggtitle("Observed and explained inequality")

  }

  # 2. type = "traceplot" ----

  if (type == "traceplot"){

    lambda <- x$path[[traceplot.which]]["lambda",]
    n.iter <- length(lambda)
    path.theta <- sapply(seq_len(n.iter),function(v)coef.PLR(x, renormalize = FALSE, pars.idx = c(traceplot.which,v)))

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
      score.df <- do.call(rbind, lapply(1:length(x$path), function(i) {
        data.frame(
          grid = i,
          lambda = -log(x$path[[i]]["lambda",]),
          BIC = x$path[[i]]["BIC score",]
        )
      }))
    }else{
      score.df$BIC <- unlist(sapply(1:length(x$path), function(i) x$path[[i]]["BIC score",]))
    }
    score.df$BIC <- max(score.df$BIC)/score.df$BIC
    scores.only <- score.df[,-(1:2),drop=FALSE]

    df.long <- data.frame(
      grid = rep(score.df$grid, ncol(scores.only)), # Repeat 'grid' column values for each method
      lambda = rep(score.df$lambda, ncol(scores.only)),
      method = rep(names(scores.only),each=nrow(scores.only)), # Create method column
      score = unlist(scores.only, use.names = FALSE) # Combine scores
    )

    if(!is.null(x$grid.value)){
      custom_labels <- paste0(x$call$grid.arg," = ",round(x$grid.value,4))
      names(custom_labels) <- 1:length(x$grid.value)
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

  # 4. Output ----

  return(g)

}

#' @method autoplot PLR_boot
#' @export

autoplot.PLR_boot <- function(x, type = c("explained","traceplot","diagnostic"), traceplot.which = "BIC", score.df = NULL, ...){

  type <- match.arg(type)

  # 1. type = "explained" ----

  if (type == "explained"){

    g <- NextMethod("autoplot")
    y <- x$y
    x <- predict(x, pars.idx = "Boot")
    g <- Lorenz.graphs_add(g, y, x, curve_label = "index (Boot)",...)

  }

  # 2. type = "traceplot" ----

  if (type == "traceplot"){

    if(traceplot.which == "Boot") traceplot.which <- x$grid.idx["Boot"]
    g <- NextMethod("autoplot")

  }

  # 3. type = "diagnostic" ----

  if (type == "diagnostic"){

    if(is.null(score.df)){
      score.df <- do.call(rbind, lapply(1:length(x$path), function(i) {
        data.frame(
          grid = i,
          lambda = -log(x$path[[i]]["lambda",]),
          OOB = x$path[[i]]["OOB score",]
        )
      }))
    }else{
      score.df$OOB <- unlist(sapply(1:length(x$path), function(i) x$path[[i]]["OOB score",]))
    }
    score.df$OOB <- score.df$OOB/max(score.df$OOB)

    g <- NextMethod("autoplot", score.df = score.df)

  }

  return(g)

}

#' @method autoplot PLR_cv
#' @export

autoplot.PLR_cv <- function(x, type = c("explained","traceplot","diagnostic"), traceplot.which = "BIC", score.df = NULL, ...){

  type <- match.arg(type)

  # 1. type = "explained" ----

  if (type == "explained"){

    g <- NextMethod("autoplot")
    y <- x$y
    x <- predict(x, pars.idx = "CV")
    g <- Lorenz.graphs_add(g, y, x, curve_label = "index (CV)",...)

  }

  # 2. type = "traceplot" ----

  if (type == "traceplot"){

    if(traceplot.which == "CV") traceplot.which <- x$grid.idx["CV"]
    g <- NextMethod("autoplot")

  }

  # 3. type = "diagnostic" ----

  if (type == "diagnostic"){

    if(is.null(score.df)){
      score.df <- do.call(rbind, lapply(1:length(x$path), function(i) {
        data.frame(
          grid = i,
          lambda = -log(x$path[[i]]["lambda",]),
          CV = x$path[[i]]["CV score",]
        )
      }))
    }else{
      score.df$CV <- unlist(sapply(1:length(x$path), function(i) x$path[[i]]["CV score",]))
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


