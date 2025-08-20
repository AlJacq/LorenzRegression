#' Graphs of concentration curves
#'
#' \code{Lorenz.graphs} traces the Lorenz curve of a response and the concentration curve of the response and each of a series of covariates.
#'
#' @param formula A formula object of the form \emph{response} ~ \emph{other_variables}. The form \emph{response} ~ \emph{1} is used to display only the Lorenz curve of the response.
#' @param data A dataframe containing the variables of interest
#' @param difference A logical determining whether the vertical axis should be expressed in terms of deviation from perfect equality. Default is \code{FALSE}.
#' @param palette A vector of colors. If \code{NULL} (default), the base R
#' palette is used. When provided, the first color is reserved for the line of equality
#' (typically "black"), and the remaining colors are used for the Lorenz and concentration curves.
#' @param ... Further arguments (see Section 'Arguments' in \code{\link{Lorenz.curve}}).
#'
#' @return A plot comprising
#' \itemize{
#'    \item The Lorenz curve of \emph{response}
#'    \item The concentration curves of \emph{response} with respect to each element of \emph{other_variables}
#' }
#'
#' @seealso \code{\link{Lorenz.curve}}, \code{\link{Gini.coef}}
#'
#' @examples
#' data(Data.Incomes)
#' Lorenz.graphs(Income ~ Age + Work.Hours, data = Data.Incomes)
#' # Expressing now the vertical axis as the deviation from perfect equality
#' Lorenz.graphs(Income ~ Age + Work.Hours, data = Data.Incomes, difference = TRUE)
#'
#' @importFrom stats model.response model.weights model.matrix
#' @importFrom ggplot2 ggplot aes scale_color_manual geom_hline stat_function labs
#'
#' @export

Lorenz.graphs <- function(formula, data, difference = FALSE, palette = NULL, ...){

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

  y <- model.response(mf, "numeric")
  w <- as.vector(model.weights(mf))
  if (!is.null(w) && !is.numeric(w))
    stop("'weights' must be a numeric vector")
  x <- model.matrix(mt, mf)[,-1,drop=FALSE]

  p <- NULL

  if(is.null(palette)) palette <- 1:(ncol(mf)+1)

  graph <- ggplot(data.frame(p=c(0,1)),aes(p)) +
    scale_color_manual(values = palette[-1],
                       breaks = colnames(mf))

  # 1 > Lorenz curve ----

  if(difference){
    graph <- graph + geom_hline(yintercept = 0,color=palette[1]) +
      stat_function(fun=function(p)Lorenz.curve(y, ...)(p)-p, geom="line",aes(color=colnames(mf)[1])) +
      labs(x = "Cumulative share of the population",y = "Deviation from perfect equality", color= "Ranking:")
  }else{
    graph <- graph + stat_function(fun=function(p)p, geom="line", color=palette[1]) +
    stat_function(fun=function(p)Lorenz.curve(y, ...)(p), geom="line",aes(color=colnames(mf)[1])) +
    labs(x = "Cumulative share of the population",y = paste0("Cumulative share of ",colnames(mf)[1]), color= "Ranking:")
  }

  # 2 > Concentration curves ----

  if(ncol(x) > 0){

    for (i in 1:ncol(x)){
      if(difference){
        graph <- local({
          j <- i
          graph + stat_function(fun=function(p)Lorenz.curve(y,x[,j], ...)(p)-p, geom="line", aes(color=colnames(mf)[j+1]))
        })
      }else{
        graph <- local({
          j <- i
          graph + stat_function(fun=function(p)Lorenz.curve(y,x[,j], ...)(p), geom="line", aes(color=colnames(mf)[j+1]))
        })
      }

    }

  }

  graph
}

#' @importFrom ggplot2 geom_ribbon aes
#' @importFrom stats quantile
#' @keywords internal

Lorenz.bands <- function(g, LC_ordinates, level, difference = FALSE, palette = NULL) {

  # Determine the upper and lower bounds
  lci <- apply(LC_ordinates, 2, quantile, probs = (1-level)/2)
  uci <- apply(LC_ordinates, 2, quantile, probs = 1-(1-level)/2)
  p <- seq(from = 0, to = 1, length.out = 100)
  if(difference){
    lci <- lci - p
    uci <- uci - p
  }

  # Add the bands
  df_band <- data.frame(p = p, lci = lci, uci = uci)
  if(is.null(palette)) palette <- 1:3
  g <- g + geom_ribbon(data = df_band, aes(x = p, ymin = lci, ymax = uci), fill = palette[3], alpha = 0.3)
  g
}
