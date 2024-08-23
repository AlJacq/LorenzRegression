#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// [[Rcpp::export(.Fitness_meanrank)]]

double Fitness_meanrank(arma::vec x, arma::vec Y, arma::mat X, arma::vec pi, double tolerance) {
  // 0. Let us define some basic objects
  int nx = x.n_rows;
  int n = Y.n_rows;
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
  vec index1u = index1(arma::find_unique(index1));
  vec index2u = index2(arma::find_unique(index2));
  vec index1k = index1u(arma::sort_index(index1u));
  vec index2k = index2u(arma::sort_index(index2u));
  int n1k = index1k.n_rows;
  int n2k = index2k.n_rows;
  vec pi1k(n1k);
  vec pi2k(n2k);
  vec F1k(n1k);
  vec F2k(n2k);
  for (int k=0; k<n1k;k++){
    for (int j=0; j<n; j++){
     if (index1(j) == index1k(k)) pi1k(k) = pi1k(k) + pi(j);
    }
    for (int l=k; l<n1k; l++){
      F1k(l) = F1k(l) + pi1k(k);
      if (l == k){
        F1k(l) = F1k(l) - pi1k(l)/2;
      }
    }
  }
  for (int k=0; k<n2k;k++){
    for (int j=0; j<n; j++){
      if (index2(j) == index2k(k)) pi2k(k) = pi2k(k) + pi(j);
    }
    for (int l=k; l<n2k; l++){
      F2k(l) = F2k(l) + pi2k(k);
      if (l == k){
        F2k(l) = F2k(l) - pi2k(l)/2;
      }
    }
  }
  vec Fi1(n);
  vec Fi2(n);
  for (int k=0; k<n1k;k++){
    for (int j=0; j<n; j++){
      if (index1(j) == index1k(k)) Fi1(j) = Fi1(j) + F1k(k);
    }
  }
  for (int k=0; k<n2k;k++){
    for (int j=0; j<n; j++){
      if (index2(j) == index2k(k)) Fi2(j) = Fi2(j) + F2k(k);
    }
  }
  // 3. Computation objective function
  vec Y_pi = Y%pi;
  vec Obj(2);
  Obj(0) = as_scalar(Y_pi.t()*Fi1);
  Obj(1) = as_scalar(Y_pi.t()*Fi2);
  double Fit = max(Obj);
  double pen = ::fabs(Fit)*::fabs(accu(abs(theta1))-1);
  return Fit-pen;
}
