#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

arma::vec frac_rank_cpp(arma::vec x, arma::vec pi);  // Declaration

// [[Rcpp::export(.Fitness_meanrank)]]

double Fitness_meanrank(arma::vec x, arma::vec Y, arma::mat X, arma::vec pi, double tolerance) {
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
  index1 = round(index1 / tolerance) * tolerance;
  index2 = round(index2 / tolerance) * tolerance;
  // 2. Computation of fractional rank
  vec Fi1 = frac_rank_cpp(index1, pi);
  vec Fi2 = frac_rank_cpp(index2, pi);
  // 3. Computation objective function
  vec Y_pi = Y%pi;
  vec Obj(2);
  Obj(0) = as_scalar(Y_pi.t()*Fi1);
  Obj(1) = as_scalar(Y_pi.t()*Fi2);
  double Fit = max(Obj);
  double pen = ::fabs(Fit)*::fabs(accu(abs(theta1))-1);
  return Fit-pen;
}
