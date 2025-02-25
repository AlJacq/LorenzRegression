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
  vec der(p);

  // Initializing der(.)
  for (k=0; k<p; k++)
    der[k] = 0;

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
      double contrib = pi(i) * pi(j) * (y(i) - y(j)) * ker;
      for (k = 0; k < p; k++) {
        xd = (X(i, k) - X(j, k)) / h;
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

