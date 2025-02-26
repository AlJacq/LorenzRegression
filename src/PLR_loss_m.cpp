#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// [[Rcpp::export(.PLR_loss_cpp_m)]]
double PLR_loss_cpp_m(double lossz, arma::mat X, arma::vec y, arma::vec ycum, int y_skipped, arma::vec pi, arma::vec theta, double h, double gamma, int kernel)
{
  int i, j;
  int k;
  int n=y.n_rows;
  int p=theta.n_rows;
  double sum=-lossz, u=0, kerd;
  double pen=0;
  vec index = X*theta;

  // Determine potential for index skipping
  arma::uvec o_idx = arma::sort_index(index);
  arma::vec index_idx = index.elem(o_idx);
  arma::vec y_idx = y.elem(o_idx);
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
  // For index-skipping loop
  std::vector<double> y_idx_std(y_idx.begin(), y_idx.end());
  std::vector<double> index_idx_std(index_idx.begin(), index_idx.end());
  std::vector<double> pi_idx_std(pi_idx.begin(), pi_idx.end());

  // Divide index by h outside the loop
  for (double& val : index_std) {
    val /= h;
  }
  for (double& val : index_idx_std) {
    val /= h;
  }

  constexpr double k1_a = 9.0/8.0;
  constexpr double k1_b = 5.0/8.0;
  constexpr double k2_a = 45.0/32.0;
  constexpr double k2_b = 25.0/16.0;
  constexpr double k2_c = 21.0/32.0;

  if (index_skipped > y_skipped){

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
        if(u < -1){
          kerd = -0.5;
        }else if (u>= 1){
          kerd = 0.5;
        }else{
          if(kernel == 1){
            kerd = k1_a * u - k1_b * u * u * u;
          }else if(kernel == 2){
            double u3 = u * u * u;
            kerd = k2_a*u - k2_b*u3 + k2_c * u3 * u * u;
          }
        }

        // Computation of loss
        sum = sum + pi_idx_std[i]*pi_idx_std[j]*(y_idx_std[i]-y_idx_std[j])*kerd;

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
        if(u < -1){
          kerd = -0.5;
        }else if (u>= 1){
          kerd = 0.5;
        }else{
          if(kernel == 1){
            kerd = k1_a * u - k1_b * u * u * u;
          }else if(kernel == 2){
            double u3 = u * u * u;
            kerd = k2_a*u - k2_b*u3 + k2_c * u3 * u * u;
          }
        }

        // Computation of loss
        sum = sum + pi_std[i]*pi_std[j]*(y_std[i]-y_std[j])*kerd;

      }
    }

  }



  for (k=0; k<p; k++)
    pen = pen + pow(theta[k],2.0);

  sum = -sum + gamma*pen;

  return sum;
}
