#include <RcppArmadillo.h>
#include <chrono>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// [[Rcpp::export(.PLR_derivative_cpp_m)]]
arma::vec PLR_derivative_cpp_m(arma::vec derz,arma::vec y, arma::vec ycum, int y_skipped, arma::mat X, arma::vec pi, arma::vec theta, double h, double gamma, int kernel)
{

  int i, j, k;
  double  kerd, u=0;
  int n=y.n_rows;
  int p=theta.n_rows;
  std::vector<double> der(derz.begin(), derz.end());
  arma::vec index = X * theta;

  // Determine potential for index skipping
  arma::uvec o_idx = arma::sort_index(index);
  arma::vec index_idx = index.elem(o_idx);
  arma::vec y_idx = y.elem(o_idx);
  arma::mat X_idx = X.rows(o_idx);
  arma::vec pi_idx = pi.elem(o_idx);
  // Cumulative counts of sorted unique values
  std::unordered_map<double, int> count_map;
  arma::vec icum(index_idx.n_elem);
  for (size_t j = 0; j < index_idx.n_elem; j++) {
    count_map[index_idx[j]]++; // Increase count
    icum[j] = count_map[index_idx[j]]; // Store cumulative count
  }
  // Compute index_skipped (sum of k*(k-1)/2 for each unique value count)
  int index_skipped = 0;
  for (const auto& pair : count_map) {
    int k = pair.second;
    index_skipped += (k * (k - 1)) / 2;
  }

  // The loop is way faster with std::vectors
  // For y-skipping loop
  std::vector<double> y_std(y.begin(), y.end());
  std::vector<double> index_std(index.begin(), index.end());
  std::vector<double> pi_std(pi.begin(), pi.end());
  std::vector<std::vector<double>> X_std(n, std::vector<double>(p));
  for (int i = 0; i < n; i++) {
    for (int k = 0; k < p; k++) {
      X_std[i][k] = X(i, k);
    }
  }
  // For index-skipping loop
  std::vector<double> y_idx_std(y_idx.begin(), y_idx.end());
  std::vector<double> index_idx_std(index_idx.begin(), index_idx.end());
  std::vector<double> pi_idx_std(pi_idx.begin(), pi_idx.end());
  std::vector<std::vector<double>> X_idx_std(n, std::vector<double>(p));
  for (int i = 0; i < n; i++) {
    for (int k = 0; k < p; k++) {
      X_idx_std[i][k] = X_idx(i, k);
    }
  }

  // Divide index by h outside the loop
  for (double& val : index_std) {
    val /= h;
  }
  for (double& val : index_idx_std) {
    val /= h;
  }

  // Store kernel constants (and externalize the h that would go in contrib)
  double k1_a = -9.0 / 8.0 / h;
  double k1_b = -15.0 / 8.0 / h;
  double k2_a = -45.0 / 32.0 / h;
  double k2_b = -75.0 / 16.0 / h;
  double k2_c = 105.0 / 32.0 / h;

  if(index_skipped > y_skipped){

    for (i=1; i<n; i++)
    {
      // Loop-skipping 1: if index_i = index_j, contrib = 0
      int j_end = i - icum[i];
      if (j_end < 0) j_end = 0;
      for (j = 0; j <= j_end; j++)
      {

        // Loop-skipping 2: if y_i = y_j, contrib = 0
        if(std::abs(y_idx_std[i]-y_idx_std[j]) < 1e-12) continue;

        // Computation of u_{ij}
        u =  (index_idx_std[i] - index_idx_std[j]);

        // Computation of difference k(u)-k(0)
        if (kernel == 1){
          if(u < -1 || u > 1){
            kerd = k1_a;
          } else {
            kerd = k1_b * u * u;
          }
        } else if (kernel == 2){
          if(u < -1 || u > 1){
            kerd = k2_a;
          } else {
            double u2 = u * u;
            kerd = k2_b * u2 + k2_c * u2 * u2;
          }
        }

        // Computation of der(k)
        double contrib = pi_idx_std[i] * pi_idx_std[j] * (y_idx_std[i] - y_idx_std[j]) * kerd;
        const std::vector<double>& Xi = X_idx_std[i];
        const std::vector<double>& Xj = X_idx_std[j];
        for (k = 0; k < p; k++) {
          der[k] += contrib * (Xi[k] - Xj[k]);
        }

      }
    }

  }else{

    for (i=1; i<n; i++)
    {
      // Loop-skipping 1: if y_i = y_j, contrib = 0
      int j_end = i - ycum[i];
      if (j_end < 0) j_end = 0;
      for (j = 0; j <= j_end; j++)
      {

        // Computation of u_{ij}
        u =  (index_std[i] - index_std[j]);
        // Loop-skipping 2: if index_i = index_j, contrib = 0
        if(std::abs(u) < 1e-12) continue;

        // Computation of difference k(u)-k(0)
        if (kernel == 1){
          if(u < -1 || u > 1){
            kerd = k1_a;
          } else {
            kerd = k1_b * u * u;
          }
        } else if (kernel == 2){
          if(u < -1 || u > 1){
            kerd = k2_a;
          } else {
            double u2 = u * u;
            kerd = k2_b * u2 + k2_c * u2 * u2;
          }
        }

        // Computation of der(k)
        double contrib = pi_std[i] * pi_std[j] * (y_std[i] - y_std[j]) * kerd;
        const std::vector<double>& Xi = X_std[i];
        const std::vector<double>& Xj = X_std[j];
        for (k = 0; k < p; k++) {
          der[k] += contrib * (Xi[k] - Xj[k]);
        }

      }
    }

  }

  for (k=0; k<p; k++)
    der[k] = der[k] - 2*gamma*theta[k];

  return der;
}
