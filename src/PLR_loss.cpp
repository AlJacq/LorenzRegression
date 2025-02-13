#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// [[Rcpp::export(.PLR_loss_cpp)]]
double PLR_loss_cpp(arma::mat X, arma::vec y, arma::vec pi, arma::vec theta, double h, double gamma, int kernel)
{
  int i, j;
  int k;
  int n=y.n_rows;
  int p=theta.n_rows;
  double sum=0, u=0;
  double pen=0;

  vec index = X*theta;

  for (i=1; i<n; i++)
  {
    double y_i = y(i);
    double pi_i = pi(i);
    double index_i = index(i);

    for (j=0; j<i; j++)
    {
      double a = y_i - y(j);
      // Avoid computation if increment of sum is 0 anyway
      if (a == 0) continue;
      u = (index_i-index(j))/h;
      if (u > -1 && u < 1 && kernel == 1) sum =  sum + 1.0* pi_i*pi(j)*a * (9.0/8.0*u - 5.0/8.0*pow(u,3.0) + 0.5);
      if (u > -1 && u < 1 && kernel == 2) sum =  sum + 1.0* pi_i*pi(j)*a * (45.0/32.0*u - 25.0/16.0*pow(u,3.0) + 21.0/32.0*pow(u,5.0) + 0.5);
      if (u >= 1) sum =  sum + 1.0* pi_i*pi(j)*a;
    }
  }

  for (k=0; k<p; k++)
    pen = pen + pow(theta[k],2.0);

  sum = -sum + gamma*pen;

  return sum;
}
