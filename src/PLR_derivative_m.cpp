#include <RcppArmadillo.h>
#include <chrono>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace std::chrono;

// Global timing variables
static double time_all = 0.0;
static int call_count = 0;

// [[Rcpp::export(.PLR_derivative_cpp_m)]]
arma::vec PLR_derivative_cpp_m(arma::vec derz,arma::vec y, arma::vec ycum, arma::mat X, arma::vec pi, arma::vec theta, double h, double gamma, int kernel)
{

  auto start_time = high_resolution_clock::now();

  int i, j, k;
  double  kerd, u=0, xd=0;
  int n=y.n_rows;
  int p=theta.n_rows;
  vec der = derz;
  vec index = X*theta;

  for (i=1; i<n; i++)
  {
    // Loop-skipping 1: if y_i = y_j, contrib = 0
    int j_end = i - ycum[i];
    if (j_end < 0) j_end = 0;
    for (j = 0; j <= j_end; j++)
    {

      // Computation of u_{ij}
      u =  (index(i) - index(j))/h;
      // Loop-skipping 2: if out of kernel bonds, contrib = 0
      // Loop-skipping 3: if u_{ij}=0, the contrib is the same as der_0 and therefore we can skip
      if (std::abs(u) < 1e-12) continue;

      // Computation of difference k(u)-k(0)
      if (kernel == 1){
        if(u < -1 || u > 1){
          kerd = - 9.0/8.0;
        } else {
          kerd = - 15.0/8.0*pow(u,2.0);
        }
      } else if (kernel == 2){
        if(u < -1 || u > 1){
          kerd = - 45.0/32.0;
        } else {
          kerd = - 75.0/16.0*pow(u,2.0) + 105.0/32.0*pow(u,4.0);
        }
      }

      // Computation of der(k)
      double contrib = pi(i) * pi(j) * (y(i) - y(j)) * kerd;
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

  auto end_time = high_resolution_clock::now();
  time_all += duration<double>(end_time - start_time).count();
  call_count++;
  if (call_count % 10 == 0) {
    Rcpp::Rcout << "Timing summary after " << call_count << " calls:" << std::endl;
    Rcpp::Rcout << "Complete time: " << time_all << " seconds" << std::endl;
  }

  return der;
}
