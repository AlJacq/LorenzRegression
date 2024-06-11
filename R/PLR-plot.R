#' Plots for the Penalized Lorenz Regression
#'
#' \code{plot.PLR} provides plots for an object of class \code{"PLR"}.
#'
#' @param x An object of class \code{"PLR"}.
#' @param type A character string indicating the type of plot. Possible values are \code{"explained"}, \code{"traceplot"} and \code{"diagnostic"}.
#' If \code{"explained"} is selected, the graph displays the Lorenz curve of the response and concentration curve(s) of the response with respect to the estimated index. More specifically, there is one concentration curve per selection method available.
#' If \code{"traceplot"} is selected, the graph displays a traceplot, where the horizontal axis is -log(lambda), lambda being the value of the penalty parameter. The vertical axis gives the size of the coefficient attached to each covariate.
#' If \code{"diagnostic"} is selected, the graph displays a faceted plot, where each facet corresponds to a different value of the tuning parameter. Each plot shows the evolution of the scores of each available selection method. For comparability reasons, the scores are normalized such that the larger the better and the optimum is attained in 1.
#' @param ... Additional arguments.
#'
#' @return The graph, as an object of class \code{"ggplot"}.
#'
#' @seealso \code{\link{Lorenz.Reg}}
#'
#' @examples
#' data(Data.Incomes)
#' PLR <- Lorenz.Reg(Income ~ ., data = Data.Incomes, penalty = "SCAD",
#'                   sel.choice = c("BIC","CV"), h.grid = nrow(Data.Incomes)^(-1/5.5),
#'                   eps = 0.01, seed.CV = 123, nfolds = 5)
#' plot(PLR)
#'
#' @import ggplot2
#'
#' @method plot PLR
#' @export

plot.PLR <- function(x, type = c("explained","traceplot","diagnostic"), traceplot.which = "BIC", ...){

  if (!inherits(x, "PLR")) stop("x must be of class 'PLR'")

  type <- match.arg(type)

  if(is.character(traceplot.which)){
    traceplot.which <- match.arg(tolower(traceplot.which),c("bic","boot","cv"))
  }else{
    if(!(traceplot.which %in% 1:length(x$path))) stop("If traceplot.which is set to an integer, it must correspond to the index of a value on the grid of the tuning parameter.")
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

    g <- Lorenz.graphs(formula, data, weights = x$weights)
    g <- g + ggtitle("Observed and explained inequality")

  }

  # 2. type = "traceplot" ----

  if (type == "traceplot"){

    if(traceplot.which == "bic"){
      if(!inherits(x,c("PLR_boot","PLR_cv"))){
        traceplot.which <- x$which.tuning
      }else{
        traceplot.which <- x$which.tuning["BIC"]
      }
    }
    if(traceplot.which == "boot"){
      if(!inherits(x,"PLR_boot")) stop("The object must be of class 'PLR_boot'")
      traceplot.which <- x$which.tuning["Boot"]
    }
    if(traceplot.which == "cv"){
      if(!inherits(x,"PLR_cv")) stop("The object must be of class 'PLR_cv'")
      traceplot.which <- x$which.tuning["CV"]
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

    g <- ggplot2::ggplot(df.long) +
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
        tuning = i,
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
      tuning = rep(df.wide$tuning, 3), # Repeat 'tuning' column values for each method
      lambda = rep(df.wide$lambda, 3),
      method = c(rep("BIC", nrow(df.wide)), rep("OOB", nrow(df.wide)),rep("CV", nrow(df.wide))), # Create method column
      score = c(df.wide$score.BIC, df.wide$score.OOB, df.wide$score.CV) # Combine scores
    )

    df.long <- na.omit(df.long)

    g <- ggplot(df.long, aes(x = lambda, y = score, color = method)) +
      geom_line() +
      facet_wrap(~ tuning, scales = "free_x") +
      labs(x = "Lambda", y = "Score", color = "Selection method") +
      theme(legend.position = "bottom")

  }

  # 4. Output ----

  g

}
