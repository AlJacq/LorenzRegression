#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// [[Rcpp::export(.PLR_derivative_cpp)]]
arma::vec PLR_derivative_cpp(arma::vec y, arma::mat X, arma::vec pi, arma::vec theta, double h, double gamma)
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
    for (j=0; j<i; j++)
    {

      if (j != i)
      {
        u =  (index(i) - index(j))/h;

        for (k=0; k<p; k++)
          v(k) =  (X(i,k) - X(j,k))/h;

        if (u < -1 || u > 1) a0=0;
        else a0 = 9.0/8.0 - 15.0/8.0*pow(u,2.0);

        for (k=0; k<p; k++){
          der(k) = der(k) + 1.0 * pi(i)*pi(j)*(y(i)-y(j)) * a0 * (v(k));
        }

      }
    }
  }

  for (k=0; k<p; k++)
    der[k] = der[k] - 2*gamma*theta[k];

  return der;
}
