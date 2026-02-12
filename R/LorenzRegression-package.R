#' LorenzRegression : A package to estimate and interpret Lorenz regressions
#'
#' The \code{LorenzRegression} package proposes a toolbox to estimate, produce inference on and interpret Lorenz regressions.
#' As argued in Heuchenne and Jacquemain (2020) and Jacquemain, A., C. Heuchenne, and E. Pircalabelu (2024), these regressions are used to determine the explanatory power of a set of covariates on the inequality of a response variable.
#' In a nutshell, each variable is given a weight in order to maximize the concentration index of the response with respect to a weighted sum of the covariates.
#' The obtained concentration index is called the explained Gini coefficient. If a single-index model with increasing link function is assumed, the explained Gini boils down to the Gini coefficient of the fitted part of the model.
#' This package rests on two main functions: \code{\link{Lorenz.Reg}} for the estimation process and \code{\link{Lorenz.boot}} for more complete inference (tests and confidence intervals).
#'
#' We direct the user to Heuchenne and Jacquemain (2020) and Jacquemain, A., C. Heuchenne, and E. Pircalabelu (2024) for a rigorous exposition of the methodology.
#'
#' @section References:
#' Heuchenne, C. and A. Jacquemain (2022). Inference for monotone single-index conditional means: A Lorenz regression approach. \emph{Computational Statistics & Data Analysis 167(C)}.
#' Jacquemain, A., C. Heuchenne, and E. Pircalabelu (2024). A penalised bootstrap estimation procedure for the explained Gini coefficient. \emph{Electronic Journal of Statistics 18(1) 247-300}.
#'
#' @useDynLib LorenzRegression
#' @importFrom Rcpp evalCpp
#' @keywords internal
"_PACKAGE"
