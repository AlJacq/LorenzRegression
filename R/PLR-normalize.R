#' Re-normalizes the estimated coefficients of a penalized Lorenz regression
#'
#' \code{PLR.normalize} transforms the estimated coefficients of a penalized Lorenz regression to match the model where the first category of each categorical variable is omitted.
#'
#' @param object An object of S3 class \code{"PLR"}.
#'
#' @return A matrix of re-normalized coefficients.
#'
#' @seealso \code{\link{Lorenz.Reg}}

PLR.normalize <- function(object){

  theta <- object$theta
  if(is.vector(theta)) theta <- t(as.matrix(theta))
  classes_cov <- attr(object$terms,"dataClasses")[-1]
  which.factor <- which(classes_cov=="factor")
  which.other <- which(classes_cov!="factor")
  all.cats <- lapply(names(object$xlevels),function(x)paste0(x,object$xlevels[[x]]))
  ref.cats <- sapply(all.cats,function(x)x[1])
  other.vars <- as.list(names(which.other))
  all.vars <- c(all.cats,other.vars)

  to.del <- sort(unique(unlist(lapply(ref.cats,function(x)grep(x,colnames(theta))))))

  if(length(to.del)==0){

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

  if(nrow(theta.new)==1) theta.new <- setNames(as.vector(theta.new),colnames(theta.new))

  return(theta.new)

}
