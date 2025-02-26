#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// [[Rcpp::export(.PLR_derivative_cpp_zero)]]
arma::vec PLR_derivative_cpp_zero(arma::vec y, arma::vec ycum, arma::mat X, arma::vec pi, arma::vec theta, double h, double gamma, int kernel)
{
  int i, j, k;
  double  ker, u=0, xd=0;
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

  for (i=1; i<n; i++)
  {
    // Loop-skipping 1: if y_i = y_j, contrib = 0
    int j_end = i - ycum[i];
    if (j_end < 0) j_end = 0;
    for (j = 0; j <= j_end; j++)
    {
      // Remark : there is no "u" here since it is 0 everywhere

      // Computation of k(u)
      if (kernel == 1) ker = 9.0/8.0;
      if (kernel == 2) ker = 45.0/32.0;

      // Computation of der(k)
      double contrib = pi_std[i] * pi_std[j] * (y_std[i] - y_std[j]) * ker;
      for (k = 0; k < p; k++) {
        xd = (X_std[i][k] - X_std[j][k]) / h;
        // Loop-skipping 4: if x_ik = x_jk, contrib = 0
        if (std::abs(xd) > 1e-12) {
          der[k] += contrib * xd;
        }
      }

    }
  }

  for (k=0; k<p; k++)
    der[k] = der[k] - 2*gamma*theta[k];

  return der;
}

