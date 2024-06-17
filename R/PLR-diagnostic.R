#' Diagnostic for the Penalized Lorenz Regression
#'
#' \code{diagnostic.PLR} provides diagnostic information for an object of class \code{"PLR"}
#' It restricts the path of the PLR to pairs of parameters (tuning, lambda) that satisfy a threshold criterion.
#'
#' @param object An object of class \code{"PLR"}.
#' @param tol A numeric threshold value used to restrict the PLR path. More specifically, we restrict to pairs (tuning,lambda) whose score exceeds \code{tol}. Default value is 0.95.
#' @param method A character string specifying the method used to evaluate the scores.
#'        Options are \code{"union"}, \code{"intersect"}, \code{"BIC"}, \code{"Boot"}, and \code{"CV"}.
#'        \describe{
#'          \item{"BIC"}{The score is the BIC-score.}
#'          \item{"Boot"}{The score is the OOB-score.}
#'          \item{"CV"}{The score is the CV-score.}
#'          \item{"union"}{The threshold requirement must be met for at least one of the selection methods available.}
#'          \item{"intersect"}{The threshold requirement must be met for all selection methods available.}
#'        }
#' @return A list with two elements:
#' \describe{
#'   \item{\code{path}}{The restricted model path, containing only the values of the pair (tuning, lambda) that satisfy the threshold criterion.}
#'   \item{\code{best}}{The best model. It is obtained by considering the pair (tuning, lambda) in the restricted path that maximizes the minimum score across all selection methods available.}
#' }
#'
#' @seealso \code{\link{Lorenz.Reg}}
#'
#' @examples
#' data(Data.Incomes)
#'
#' @export

diagnostic.PLR <- function(object, tol = 0.99, method = c("union","intersect","BIC","Boot","CV")){

  if (!inherits(object, "PLR")) stop("The object must be of class 'PLR'")

  method <- match.arg(method)

  df.wide <- do.call(rbind, lapply(1:length(object$path), function(i) {
    data.frame(
      which.tuning = i,
      which.lambda = 1:ncol(object$path[[i]]),
      minloglambda = -log(object$path[[i]]["lambda",]),
      nnzeroes = object$path[[i]]["Number of nonzeroes",],
      score.BIC = object$path[[i]]["BIC score",],
      score.OOB = if (inherits(object, "PLR_boot")) object$path[[i]]["OOB score",] else NA,
      score.CV = if (inherits(object, "PLR_cv")) object$path[[i]]["CV score",] else NA
    )
  }))

  df.wide$score.BIC <- max(df.wide$score.BIC)/df.wide$score.BIC
  if(inherits(object, "PLR_boot")) df.wide$score.OOB <- df.wide$score.OOB/max(df.wide$score.OOB)
  if (inherits(object, "PLR_cv")) df.wide$score.CV <- df.wide$score.CV/max(df.wide$score.CV)

  df.wide <- df.wide[, colSums(is.na(df.wide)) < nrow(df.wide)]

  if(method == "BIC"){
    chosen <- "score.BIC"
  }else if(method == "Boot"){
    if(!inherits(object,"PLR_boot")) stop("The object must be of class 'PLR_boot'")
    chosen <- "score.OOB"
  }else if(method == "CV"){
    if(!inherits(object,"PLR_cv")) stop("The object must be of class 'PLR_cv'")
    chosen <- "score.CV"
  }else{
    chosen <- grep("score",names(df.wide),value=T)
  }

  exceeds_tol <- sapply(chosen,function(x)df.wide[,x] > tol)
  if(method == "union"){
    to_keep <- apply(exceeds_tol,1,any)
  }else if(method == "intersect"){
    to_keep <- apply(exceeds_tol,1,all)
  }else{
    to_keep <- as.vector(exceeds_tol)
  }

  path.keep <- df.wide[to_keep,]

  if(nrow(path.keep)==0){
    warning("No value of (tuning,lambda) meets the required scores. Consider switching 'method' to another value than 'intersect' or lower the value of 'tol'.")
  }else{
    path.keep2 <- path.keep[path.keep$nnzeroes == min(path.keep$nnzeroes),,drop=FALSE]
    best.keep <- path.keep2[which.max(apply(path.keep2[,grep("score",names(path.keep2),value=T)],1,min)),]
    return(list("path"=path.keep,"best"=best.keep))
  }

}
