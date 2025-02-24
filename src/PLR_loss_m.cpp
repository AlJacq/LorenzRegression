#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// [[Rcpp::export(.PLR_loss_cpp_m)]]
double PLR_loss_cpp_m(double lossz, arma::mat X, arma::vec y, arma::vec pi, arma::vec theta, double h, double gamma, int kernel)
{
  int i, j;
  int k;
  int n=y.n_rows;
  int p=theta.n_rows;
  double sum=-lossz, u=0, kerd;
  double pen=0;
  vec index = X*theta;

  for (i=1; i<n; i++)
  {
    for (j=0; j<i; j++)
    {

      // Loop-skipping 1: if y_i = y_j, contrib = 0
      if(std::abs(y(i)-y(j)) < 1e-12) continue;

      // Computation of u_{ij}
      u =  (index(i) - index(j))/h;
      // Loop-skipping 2: if u_{ij}<-1, contrib = 0
      // Loop-skipping 3: if u_{ij}=0, the contrib is the same as loss_0 and therefore we can skip
      if (u < -1 || std::abs(u) < 1e-12) continue;

      // Computation of difference k(u)-k(0)
      if (u < 1 && kernel == 1) kerd = 9.0/8.0*u - 5.0/8.0*pow(u,3.0);
      if (u < 1 && kernel == 2) kerd = 45.0/32.0*u - 25.0/16.0*pow(u,3.0) + 21.0/32.0*pow(u,5.0);
      if (u >= 1) kerd = 0.5;

      // Computation of loss
      sum = sum + pi(i)*pi(j)*(y(i)-y(j))*kerd;

    }
  }

  for (k=0; k<p; k++)
    pen = pen + pow(theta[k],2.0);

  sum = -sum + gamma*pen;

  return sum;
}
