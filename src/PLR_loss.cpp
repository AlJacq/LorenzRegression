#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// [[Rcpp::export]]
double PLR_loss_cpp(arma::mat X, arma::vec y, arma::vec pi, arma::vec theta, double sigma, double gamma)
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
    for (j=0; j<i; j++)
    {
      u = (index(i)-index(j))/sigma;
      sum =  sum + 1.0* pi(i)*pi(j)*(y(i)-y(j)) / (1 + exp(-1.0*u));
    }
  }

  for (k=0; k<p; k++)
    pen = pen + pow(theta[k],2.0);

  sum = -sum/(n-1) + gamma*pen;

  return sum;
}
