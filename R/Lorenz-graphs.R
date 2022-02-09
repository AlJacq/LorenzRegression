#' Graphs of concentration curves
#'
#' \code{Lorenz.graphs} traces the Lorenz curve of a response and the concentration curve of the response and each of a series of covariates.
#'
#' @param formula A formula object of the form \emph{response} ~ \emph{other_variables}.
#' @param data A dataframe containing the variables of interest
#' @param ... other arguments (see Section 'Arguments' in \code{\link{Lorenz.curve}}).
#'
#' @return A plot comprising
#' \describe{
#'    \item{}{The Lorenz curve of \emph{response}}
#'    \item{}{The concentration curves of \emph{response} with respect to each element of \emph{other_variables}}
#' }
#'
#' @seealso \code{\link{Lorenz.curve}}, \code{\link{Gini.coef}}
#'
#' @examples
#' data(Data.Incomes)
#' Lorenz.graphs(Income ~ Age + Work.Hours, data = Data.Incomes)
#'
#' @import ggplot2
#'
#' @export

Lorenz.graphs <- function(formula, data, ...){

  p <- NULL
  Data.temp.X <- as.data.frame(stats::model.matrix(formula,data=data)[,-1])
  Data.temp.Y <- stats::model.frame(formula,data=data)[,1]
  Data.temp <- cbind(Data.temp.Y,Data.temp.X)
  colnames(Data.temp)[1] <- colnames(stats::model.frame(formula,data=data))[1]

  if(length(Data.temp.X[1,])==1) colnames(Data.temp)[2] <- colnames(stats::model.matrix(formula,data=data))[2]

  graph <- ggplot2::ggplot(data.frame(p=c(0,1)),aes(p)) +
    scale_color_manual(values = 2:(length(Data.temp[1,])+1),
                       breaks = colnames(Data.temp)) +
    stat_function(fun=function(p)p, geom="line", color=1) +
    stat_function(fun=function(p)Lorenz.curve(Data.temp.Y, ...)(p), geom="line",aes(color=colnames(Data.temp)[1])) +
    labs(x = "Cumulative share of the population",y = paste0("Cumulative share of ",colnames(Data.temp)[1]), color= "Ranking:")

  for (i in 1:length(Data.temp.X[1,])){
    graph <- local({
    j <- i
    graph + stat_function(fun=function(p)Lorenz.curve(Data.temp.Y,Data.temp.X[,j], ...)(p), geom="line", aes(color=colnames(Data.temp)[j+1]))
    })
  }
  graph
}
