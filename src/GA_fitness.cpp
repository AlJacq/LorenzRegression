#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// Function to sort x in terms of first y and then if ties occur in terms of z
uvec arma_sort(arma::vec y, arma::vec z) {
  vec y1 = y + ::fabs(min(y));
  vec z1 = z + ::fabs(min(z));
  vec a = nonzeros(diff(y1(arma::sort_index(y1))));
  double bound;
  if (a.n_rows==1){
    bound = a(0) - exp(-10);
  }else{
    bound = min(a) - exp(-10);
  }
  vec b = y1 + z1/max(z1)*bound;
  uvec c = sort_index(b);
  return c;
}

//' @title Computes the fitness used in the GA
//' @description Computes the fitness of a candidate in the genetic algorithm displayed in function Lorenz.GA.cpp
//' @param x vector of size (p-1) giving the proposed candidate, where p is the number of covariates
//' @param Y vector of size n gathering the response, where n is the sample size
//' @param X matrix of dimension (n*p) gathering the covariates
//' @param Z vector of size n gathering iid repetitions of a U[0,1]
//' @param pi vector of size n gathering the observation weights (notice that sum(pi)=1)
//' @return Fitness of candidate x
//' @keywords internal
// [[Rcpp::export(.Fitness_cpp)]]
double Fitness_cpp(arma::vec x, arma::vec Y, arma::mat X, arma::vec Z, arma::vec pi) {
  // 0. Let us define some basic objects
  int nx = x.n_rows;
  // 1. We must acknowledge the fact that the last coefficient is constrained
  vec theta1(nx+1);//But it may be positive
  vec theta2(nx+1);//... or negative
  for (int i=0;i<(nx+1);i++){
    if (i<nx){
      theta1(i) = x(i);
      theta2(i) = x(i);
    }else{
      theta1(i) = 1-accu(abs(x));
      theta2(i) = -(1-accu(abs(x)));
    }
  }
  vec index1 = X*theta1;
  vec index2 = X*theta2;
  // 2. We need to sort the Y's first in terms of index (a) and then in terms of a unif(0,1) (b)
  uvec sort1 = arma_sort(index1,Z);
  uvec sort2 = arma_sort(index2,Z);
  vec Y_sort1 = Y(sort1);
  vec Y_sort2 = Y(sort2);
  // 3. We proceed similarly for the weights
  vec pi_sort1 = pi(sort1);
  vec pi_sort2 = pi(sort2);
  vec Y_pi_sort1 = Y_sort1%pi_sort1;
  vec Y_pi_sort2 = Y_sort2%pi_sort2;
  // 4. We can compute the ranks
  vec rank1 = cumsum(pi_sort1) - pi_sort1/2;
  vec rank2 = cumsum(pi_sort2) - pi_sort2/2;
  vec Obj(2);
  Obj(0) = as_scalar(Y_pi_sort1.t()*rank1);
  Obj(1) = as_scalar(Y_pi_sort2.t()*rank2);
  double Fit = max(Obj);
  double pen = ::fabs(Fit)*::fabs(accu(abs(theta1))-1);
  return Fit-pen;
}
