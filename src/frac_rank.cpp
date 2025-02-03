#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

//' @title Computes fractional ranks
//' @description Computes the vector of fractional ranks related to a given vector
//' @param x vector of size n gathering the values for which the fractional rank should be computed
//' @param pi vector of size n gathering the observation weights (notice that sum(pi)=1)
//' @return Fractional rank related to vector x
//' @keywords internal
// [[Rcpp::export(.frac_rank_cpp)]]
arma::vec frac_rank_cpp(arma::vec x, arma::vec pi) {

  // 1. Sort x and reorder pi accordingly
  uvec sorted_indices = sort_index(x);  // Get sorted indices
  arma::vec x_sorted = x(sorted_indices); // Sort x
  arma::vec pi_sorted = pi(sorted_indices); // Reorder pi accordingly

  // 2. Compute cumulative probabilities
  arma::vec cum_pi = cumsum(pi_sorted);  // Cumulative sum of pi
  arma::vec ranks(x.n_elem, fill::zeros);  // Initialize ranks

  // 3. Handle ties
  int n = x.n_elem;
  int i = 0;

  while (i < n) {
    int j = i;
    // Group tied values
    while (j < n - 1 && x_sorted(j) == x_sorted(j + 1)) {
      j++;
    }

    // Assign the mean rank to all tied values
    double mean_rank = cum_pi(j) - sum(pi_sorted.subvec(i, j)) / 2;
    ranks.subvec(i, j).fill(mean_rank);

    // Move to the next group
    i = j + 1;
  }

  // 4. Reorder ranks to match the original order
  arma::vec final_ranks(x.n_elem);
  final_ranks(sorted_indices) = ranks;

  return final_ranks;
}
