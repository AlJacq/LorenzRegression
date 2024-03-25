#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// [[Rcpp::export(.SCAD_derivative_cpp)]]
arma::vec SCAD_derivative_cpp(arma::vec x, double lambda, double a)
{
  int k;
  int p=x.n_rows;
  vec der(p);

  for (k=0; k<p; k++){

    if (x[k] <= lambda){
      der[k] = lambda;
    }
    else if (x[k] > a*lambda){
      der[k] = 0;
    }
    else {
      der[k] = (a*lambda-x[k])/(a-1);
    }

  }

  return der;
}
