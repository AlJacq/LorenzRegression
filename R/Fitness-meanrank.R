Fitness_meanrank <- function(u, Y, X){
  theta1 <- c(u,1-sum(abs(u)))
  theta2 <- c(u,-(1-sum(abs(u))))
  index1 <- X%*%theta1
  index2 <- X%*%theta2
  Obj <- c()
  Obj[1] <- Y%*%rank(index1, ties.method = "average")
  Obj[2] <- Y%*%rank(index2, ties.method = "average")
  Fit <- max(Obj)
  pen <- Fit*abs(sum(abs(theta1))-1)
  Fit - pen
}
