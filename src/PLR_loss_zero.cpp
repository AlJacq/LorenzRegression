#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// [[Rcpp::export(.PLR_loss_cpp_zero)]]
double PLR_loss_cpp_zero(arma::mat X, arma::vec y, arma::vec pi, double h, double gamma, int kernel)
{
  int i, j;
  int k;
  int n=y.n_rows;
  double sum=0, u=0;

  for (i=1; i<n; i++)
  {
    for (j=0; j<i; j++)
    {

      // Loop-skipping 1: if y_i = y_j, contrib = 0
      if(std::abs(y(i)-y(j)) < 1e-12) continue;

      // Remark : there is no "u" here since it is 0 everywhere

      // Computation of loss
      sum = sum + pi(i)*pi(j)*(y(i)-y(j))*0.5;

    }
  }

  sum = -sum;

  return sum;
}
