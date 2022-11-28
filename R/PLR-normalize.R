#' Re-normalizes the estimated coefficients of a penalized Lorenz regression
#'
#' \code{PLR.normalize} transforms the estimated coefficients of a penalized Lorenz regression to match the model where the first category of each categorical variable is omitted.
#'
#' @param PLR Output of a call to \code{\link{Lorenz.Reg}}, where \code{penalty!="none"}.
#'
#' @return A matrix of re-normalized coefficients.
#'
#' @seealso \code{\link{Lorenz.Reg}}
#'
#' @examples
#' data(Data.Incomes)
#' PLR <- Lorenz.Reg(Income ~ ., data = Data.Incomes, penalty = "SCAD", sel.choice = c("BIC","CV"), eps = 0.005, seed.CV = 123, nfolds = 5)
#' PLR.normalize(PLR)
#'
#' @export

PLR.normalize <- function(PLR){

  data <- PLR$data
  formula <- PLR$formula
  theta <- PLR$theta

  data <- stats::model.frame(formula,data)[,-1]
  which.factor <- which(sapply(1:ncol(data),function(i)class(data[,i]))=="factor")
  which.other <- which(sapply(1:ncol(data),function(i)class(data[,i]))!="factor")
  all.cats <- lapply(which.factor,function(i)paste0(colnames(data)[i],levels(data[,i])))
  ref.cats <- sapply(which.factor,function(i)paste0(colnames(data)[i],levels(data[,i])[1]))
  other.vars <- lapply(which.other,function(i)colnames(data)[i])
  all.vars <- c(all.cats,other.vars)

  to.del <- sort(unique(unlist(lapply(ref.cats,function(x)grep(x,colnames(theta))))))

  if(is.null(to.del)){

    theta.new <- theta

  }else{

    names(to.del) <- colnames(theta)[to.del]

    if(nrow(theta)==1){
      theta.new <- t(as.matrix(theta[,-to.del]))
      rownames(theta.new) <- rownames(theta)
    }else{
      theta.new <- as.matrix(theta[,-to.del])
    }

    for (col in to.del){
      names(col) <- colnames(theta)[col]
      # Main effect
      if(length(grep(":",names(col)))==0){

        # We subtract the coeff to the remaining categories of this variable
        which.fac <- which(sapply(all.cats,function(x)names(col)%in%x))
        col.delta <- colnames(theta.new)%in%all.cats[[which.fac]]
        theta.new[,col.delta] <- theta.new[,col.delta] - theta[,col]

        # Interaction
      }else{

        which.ref <- sapply(ref.cats,function(x)grepl(x,names(col)))
        # Interaction between two categorical at their reference level
        if (sum(which.ref)==2){

          which.fac <- which(sapply(all.cats,function(x)length(intersect(x,strsplit(names(col),split=":")[[1]]))>0))
          which.fac.names <- unlist(all.cats[which.fac])
          # Add the coeff at whole remaining interactions implicating these two variables
          col.delta1 <- which(rowSums(sapply(which.fac.names,function(x)grepl(x,colnames(theta.new))))==2)
          theta.new[,col.delta1] <- theta.new[,col.delta1] + theta[,col]
          # Subtract it to the remaining main effects
          col.delta2 <- which(sapply(colnames(theta.new),function(x)any(grepl(x,which.fac.names))))
          theta.new[,col.delta2] <- theta.new[,col.delta2] - theta[,col]

          # Interaction between a categorical at ref level and a continuous
          # Interaction between a categorical at ref level and a categorical not at ref level
          # Both cases are similar
        }else{

          which.var.ref <- names(which(which.ref))
          which.var.other <- setdiff(strsplit(names(col),split=":")[[1]],which.var.ref)
          which.var.ref <- which(sapply(all.cats,function(x)length(intersect(x,which.var.ref))>0))
          which.var.ref.names <- unlist(all.cats[which.var.ref])
          which.vars.names <- c(which.var.other,which.var.ref.names)
          # Subtract the coefficient to remaining interactions implicating these two variables
          col.delta1 <- which(rowSums(sapply(which.vars.names,function(x)grepl(x,colnames(theta.new))))==2)
          theta.new[,col.delta1] <- theta.new[,col.delta1] - theta[,col]
          # Add the coefficient to the main effects of the other variable
          col.delta2 <- which(sapply(colnames(theta.new),function(x)x==which.var.other))
          theta.new[,col.delta2] <- theta.new[,col.delta2] + theta[,col]

        }

      }

    }

    theta.new <- t(apply(theta.new,1,function(x)x/sqrt(sum(x^2))))

  }

  return(theta.new)

}
