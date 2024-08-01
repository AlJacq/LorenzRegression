#' Plots for the Penalized Lorenz Regression
#'
#' \code{plot.PLR} provides plots for an object of class \code{"PLR"}.
#'
#' @param x An object of class \code{"PLR"}.
#' @param type A character string indicating the type of plot. Possible values are \code{"explained"}, \code{"traceplot"} and \code{"diagnostic"}.
#' If \code{"explained"} is selected, the graph displays the Lorenz curve of the response and concentration curve(s) of the response with respect to the estimated index. More specifically, there is one concentration curve per selection method available.
#' If \code{"traceplot"} is selected, the graph displays a traceplot, where the horizontal axis is -log(lambda), lambda being the value of the penalty parameter. The vertical axis gives the value of the estimated coefficient attached to each covariate.
#' If \code{"diagnostic"} is selected, the graph displays a faceted plot, where each facet corresponds to a different value of the grid parameter. Each plot shows the evolution of the scores of each available selection method. For comparability reasons, the scores are normalized such that the larger the better and the optimum is attained in 1.
#' @param traceplot.which This argument indicates the value of the grid parameter for which the traceplot should be produced (see arguments \code{grid.value} and \code{grid.arg} in function \code{\link{Lorenz.Reg}}).
#' It can be an integer indicating the index in the grid determined via \code{grid.value}.
#' Alternatively, it can be a character string indicating the selection method. In this case the index corresponds to the optimal value according to that selection method.
#' @param ... Additional arguments passed to function \code{\link{Lorenz.graphs}}
#'
#' @return The graph, as an object of class \code{"ggplot"}.
#'
#' @seealso \code{\link{Lorenz.Reg}}
#'
#' @examples
#' ## For examples see example(Lorenz.Reg), example(Lorenz.boot) and example(PLR.CV)
#'
#' @importFrom ggplot2 ggplot aes geom_line ggtitle scale_color_hue labs theme_minimal facet_wrap labeller theme
#' @importFrom stats as.formula na.omit
#'
#' @method plot PLR
#' @export

plot.PLR <- function(x, type = c("explained","traceplot","diagnostic"), traceplot.which = "BIC", ...){

  if (!inherits(x, "PLR")) stop("x must be of class 'PLR'")

  type <- match.arg(type)

  if(is.character(traceplot.which)){
    traceplot.which <- match.arg(tolower(traceplot.which),c("bic","boot","cv"))
  }else{
    if(!(traceplot.which %in% 1:length(x$path))) stop("If traceplot.which is set to an integer, it must correspond to the index of a value on the grid of the grid parameter.")
  }

  # 1. type = explained ----

  if(type == "explained"){

    formula <- as.formula(paste(as.character(x$call$formula[[2]]), "~ ."))

    if( inherits(x, c("PLR_boot","PLR_cv"))){
      data <- data.frame(x$y, t(x$index))
      names(data) <- c(all.vars(formula)[1],paste("index",rownames(x$index),sep="."))
    }else{
      data <- data.frame(x$y,x$index)
      names(data) <- c(all.vars(formula)[1],"index")
    }

    g <- Lorenz.graphs(formula, data, weights = x$weights, ...)
    g <- g + ggtitle("Observed and explained inequality")

  }

  # 2. type = "traceplot" ----

  if (type == "traceplot"){

    if(traceplot.which == "bic"){
      if(!inherits(x,c("PLR_boot","PLR_cv"))){
        traceplot.which <- x$grid.idx
      }else{
        traceplot.which <- x$grid.idx["BIC"]
      }
    }
    if(traceplot.which == "boot"){
      if(!inherits(x,"PLR_boot")) stop("The object must be of class 'PLR_boot'")
      traceplot.which <- x$grid.idx["Boot"]
    }
    if(traceplot.which == "cv"){
      if(!inherits(x,"PLR_cv")) stop("The object must be of class 'PLR_cv'")
      traceplot.which <- x$grid.idx["CV"]
    }

    path <- x$path[[traceplot.which]]

    if(inherits(x,c("PLR_cv","PLR_boot"))){
      var_names <- colnames(x$theta)
    }else{
      var_names <- names(x$theta)
    }

    path.theta <- path[var_names,]
    n.iter <- ncol(path.theta)

    df.long <- data.frame(
      "Variable" = rep(var_names,n.iter),
      "theta" = as.vector(path.theta),
      "minloglambda" = rep(-log(path["lambda",]),each=length(var_names))
    )

    g <- ggplot(df.long) +
      aes(x = minloglambda, y = theta, colour = Variable) +
      geom_line(linewidth = 1L) +
      scale_color_hue() +
      labs(x = expression(paste("-log(", symbol(lambda), ")",sep="")),
           y = expression(symbol(theta)[k]),
           title = "Traceplot") +
      theme_minimal()

  }

  # 3. type = "diagnostic" ----

  if (type == "diagnostic"){

    df.wide <- do.call(rbind, lapply(1:length(x$path), function(i) {
      data.frame(
        grid = i,
        lambda = -log(x$path[[i]]["lambda",]),
        score.BIC = x$path[[i]]["BIC score",],
        score.OOB = if (inherits(x, "PLR_boot")) x$path[[i]]["OOB score",] else NA,
        score.CV = if (inherits(x, "PLR_cv")) x$path[[i]]["CV score",] else NA
      )
    }))

    df.wide$score.BIC <- max(df.wide$score.BIC)/df.wide$score.BIC
    if(inherits(x, "PLR_boot")) df.wide$score.OOB <- df.wide$score.OOB/max(df.wide$score.OOB)
    if (inherits(x, "PLR_cv")) df.wide$score.CV <- df.wide$score.CV/max(df.wide$score.CV)

    df.long <- data.frame(
      grid = rep(df.wide$grid, 3), # Repeat 'grid' column values for each method
      lambda = rep(df.wide$lambda, 3),
      method = c(rep("BIC", nrow(df.wide)), rep("OOB", nrow(df.wide)),rep("CV", nrow(df.wide))), # Create method column
      score = c(df.wide$score.BIC, df.wide$score.OOB, df.wide$score.CV) # Combine scores
    )

    df.long <- na.omit(df.long)

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

  g

}
