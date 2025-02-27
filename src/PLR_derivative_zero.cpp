#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// [[Rcpp::export(.PLR_derivative_cpp_zero)]]
arma::vec PLR_derivative_cpp_zero(arma::vec y, arma::vec ycum, arma::mat X, arma::vec pi, arma::vec theta, double h, double gamma, int kernel)
{
  int i, j, k;
  double  kerh;
  int n=y.n_rows;
  int p=theta.n_rows;
  std::vector<double> der(p);

  // Initializing der(.)
  for (k=0; k<p; k++)
    der[k] = 0;

  // Convert y, X and pi for speed
  std::vector<double> y_std(y.begin(), y.end());
  std::vector<double> pi_std(pi.begin(), pi.end());
  std::vector<std::vector<double>> X_std(n, std::vector<double>(p));
  for (i = 0; i < n; i++) {
    for (k = 0; k < p; k++) {
      X_std[i][k] = X(i, k);
    }
  }

  // Computation of k(u)
  if (kernel == 1) kerh = 9.0/8.0 / h;
  if (kernel == 2) kerh = 45.0/32.0 / h;

  for (i=1; i<n; i++)
  {
    // Loop-skipping 1: if y_i = y_j, contrib = 0
    int j_end = i - ycum[i];
    if (j_end < 0) j_end = 0;
    for (j = 0; j <= j_end; j++)
    {
      // Remark : there is no "u" here since it is 0 everywhere

      // Computation of der(k)
      double contrib = pi_std[i] * pi_std[j] * (y_std[i] - y_std[j]) * kerh;
      const std::vector<double>& Xi = X_std[i];
      const std::vector<double>& Xj = X_std[j];
      for (k = 0; k < p; k++) {
        der[k] += contrib * (Xi[k] - Xj[k]);
      }

    }
  }

  for (k=0; k<p; k++)
    der[k] = der[k] - 2*gamma*theta[k];

  return der;
}

