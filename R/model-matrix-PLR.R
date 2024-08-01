#' Design matrix in the Penalized Lorenz Regression
#'
#' \code{model_matrix_PLR} is a utilitary function that provides the design matrix for the Penalized Lorenz Regression
#'
#' @param mt Model terms
#' @param mf Model frame
#'
#' @details This function ensures that the design matrix is constructed according to the requirements of the PLR.
#' In PLR, one must exclude the intercept and use one-hot encoding for all variables, except when binary
#'
#' @importFrom parsnip contr_one_hot
#'
#' @return The design matrix
#' @keywords internal

model_matrix_PLR <- function(mt,mf){

  # 1) One-hot encoding
  cat_vars <- all.vars(mt)[sapply(all.vars(mt), function(x) is.factor(mf[[x]]))]
  custom_contrasts <- lapply(cat_vars, function(var) {
    parsnip::contr_one_hot(levels(mf[[var]]))
  })
  x <- model.matrix(mt,mf,contrasts = setNames(custom_contrasts, cat_vars))
  # 2) Exclude intercept and "true" binary variables
  binary_fac <- which(sapply(mf,nlevels)==2)
  to_del <- c("(Intercept)",paste0(names(binary_fac),sapply(mf[binary_fac],function(x)levels(x)[1])))
  x <- x[,!colnames(x)%in%to_del,drop=FALSE]
  return(x)

}
