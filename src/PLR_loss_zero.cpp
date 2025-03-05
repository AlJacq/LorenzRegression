#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// [[Rcpp::export(.PLR_loss_cpp_zero)]]
double PLR_loss_cpp_zero(arma::mat X, arma::vec y, arma::vec ycum, arma::vec pi, double h, double gamma, int kernel)
{
  int i, j;
  int n=y.n_rows;
  double sum=0;

  // Convert y and pi for speed
  std::vector<double> y_std(y.begin(), y.end());
  std::vector<double> pi_std(pi.begin(), pi.end());

  for (i=1; i<n; i++)
  {
    // Loop-skipping 1: if y_i = y_j, contrib = 0
    int j_end = i - ycum[i];
    if (j_end < 0) j_end = 0;
    for (j = 0; j <= j_end; j++)
    {

      // Remark : there is no "u" here since it is 0 everywhere

      // Computation of loss
      sum = sum + pi_std[i]*pi_std[j]*(y_std[i]-y_std[j])*0.5;

    }
  }

  sum = -sum;

  return sum;
}
