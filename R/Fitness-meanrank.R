Fitness_meanrank <- function(u, Y, X, pi){
  theta1 <- c(u,1-sum(abs(u)))
  theta2 <- c(u,-(1-sum(abs(u))))
  index1 <- X%*%theta1
  index2 <- X%*%theta2
  index1_k <- sort(unique(index1))
  pi1_k <- sapply(1:length(index1_k),function(k)sum(pi[index1==index1_k[k]]))
  F1_k <- cumsum(pi1_k) - 0.5*pi1_k
  F1_i <- sapply(1:length(index1),function(i)sum(F1_k[index1_k==index1[i]])) # Ensures that sum(F_i*pi) = 0.5
  index2_k <- sort(unique(index2))
  pi2_k <- sapply(1:length(index2_k),function(k)sum(pi[index2==index2_k[k]]))
  F2_k <- cumsum(pi2_k) - 0.5*pi2_k
  F2_i <- sapply(1:length(index2),function(i)sum(F2_k[index2_k==index2[i]])) # Ensures that sum(F_i*pi) = 0.5
  Obj <- c()
  Obj[1] <- ((pi*Y)%*%F1_i)
  Obj[2] <- ((pi*Y)%*%F2_i)
  Fit <- max(Obj)
  pen <- Fit*abs(sum(abs(theta1))-1)
  return(Fit - pen)
}

