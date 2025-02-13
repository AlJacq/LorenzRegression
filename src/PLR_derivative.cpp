#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// [[Rcpp::export(.PLR_derivative_cpp)]]
arma::vec PLR_derivative_cpp(arma::vec y, arma::mat X, arma::vec pi, arma::vec theta, double h, double gamma, int kernel)
{
  int i, j, k;
  double  a0, u=0;
  int n=y.n_rows;
  int p=theta.n_rows;
  vec v(p);
  vec der(p);

  for (k=0; k<p; k++)
    der[k] = 0;

  vec index = X*theta;

  for (i=1; i<n; i++)
  {

    double y_i = y(i);
    double pi_i = pi(i);
    double index_i = index(i);
    arma::rowvec X_i = X.row(i);

    for (j=0; j<i; j++)
    {
      // Avoid computation if increment to der(k) is 0 anyway
      double a = y_i - y(j);
      if(a == 0) continue;

      u =  (index_i - index(j))/h;
      // Avoid computation if increment to der(k) is 0 anyway
      if (u < -1 || u > 1) continue;

      for (k=0; k<p; k++)
        v(k) =  (X_i(k) - X(j,k))/h;

      if (u >= -1 && u <=1 && kernel == 1) a0 = 9.0/8.0 - 15.0/8.0*pow(u,2.0);
      if (u >= -1 && u <=1 && kernel == 2) a0 = 45.0/32.0 - 75.0/16.0*pow(u,2.0) + 105.0/32.0*pow(u,4.0);

      for (k=0; k<p; k++){
        der(k) = der(k) + 1.0 * pi_i*pi(j)*a * a0 * (v(k));
      }

    }
  }

  for (k=0; k<p; k++)
    der[k] = der[k] - 2*gamma*theta[k];

  return der;
}
