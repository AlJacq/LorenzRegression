#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// Function to sort x in terms of first y and then if ties occur in terms of z
vec arma_sort(arma::vec x, arma::vec y, arma::vec z) {
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
  return x(sort_index(b));
}

//' @title Computes the fitness used in the GA
//' @description Computes the fitness of a candidate in the genetic algorithm displayed in function Lorenz.GA.cpp
//' @param x vector of size (p-1) giving the proposed candidate, where p is the number of covariates
//' @param Y vector of size n gathering the response, where n is the sample size
//' @param X matrix of dimension (n*p) gathering the covariates
//' @param Z vector of size n gathering iid repetitions of a U[0,1]
//' @return Fitness of candidate x
// [[Rcpp::export]]
double Fitness_cpp(arma::vec x, arma::vec Y, arma::mat X, arma::vec Z) {
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
  vec Y_sort1 = arma_sort(Y,index1,Z);
  vec Y_sort2 = arma_sort(Y,index2,Z);
  vec seq = linspace(1,Y.n_rows,Y.n_rows);
  vec Obj(2);
  Obj(0) = as_scalar(Y_sort1.t()*seq);
  Obj(1) = as_scalar(Y_sort2.t()*seq);
  double Fit = max(Obj);
  double pen = ::fabs(Fit)*::fabs(accu(abs(theta1))-1);
  return Fit-pen;
}
