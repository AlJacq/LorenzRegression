#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

#include <chrono>
static double total_time_test = 0.0, total_matmul_test = 0.0, total_sort_test = 0.0, total_penalty_test = 0.0;
static int call_count_test = 0;

//' @title Computes the fitness used in the GA
//' @description Computes the fitness of a candidate in the genetic algorithm displayed in function Lorenz.GA.cpp
//' @param x vector of size (p-1) giving the proposed candidate, where p is the number of covariates
//' @param Y vector of size n gathering the response, where n is the sample size
//' @param X matrix of dimension (n*p) gathering the covariates
//' @param Z vector of size n gathering iid repetitions of a U[0,1]
//' @param pi vector of size n gathering the observation weights (notice that sum(pi)=1)
//' @param tolerance A small positive number used to determine the threshold for considering two floating-point numbers as equal. This is primarily used to
//' address issues with floating-point precision when comparing values that should theoretically be identical but may differ slightly due to numerical inaccuracies.
//' @return Fitness of candidate x
//' @keywords internal
// [[Rcpp::export(.Fitness_cpp_test)]]
double Fitness_cpp_test(arma::vec x, arma::vec Y, arma::mat X, arma::vec Z, arma::vec pi, double tolerance) {
  auto start_total_test = std::chrono::high_resolution_clock::now();
  // 0. Let us define some basic objects
  int nx = x.n_rows;
  auto start_matmul_test = std::chrono::high_resolution_clock::now();
  // 1. We must acknowledge the fact that the last coefficient is constrained
  vec theta1(nx+1);//But it may be positive
  vec theta2(nx+1);//... or negative
  for (int i=0;i<(nx+1);i++){
    if (i<nx){
      theta1(i) = x(i);
      theta2(i) = x(i);
    }else{
      theta1(i) = 1-accu(abs(x));
      theta2(i) = -(1-accu(abs(x)));
    }
  }
  vec index1 = X*theta1;
  vec index2 = X*theta2;
  auto end_matmul_test = std::chrono::high_resolution_clock::now();
  total_matmul_test += std::chrono::duration<double>(end_matmul_test - start_matmul_test).count();
  auto start_sort_test = std::chrono::high_resolution_clock::now();
  index1 = round(index1 / tolerance) + Z/2;
  index2 = round(index2 / tolerance) + Z/2;
  // 2. We need to sort the Y's first in terms of index (a) and then in terms of a unif(0,1) (b)
  uvec sort1 = sort_index(index1);
  uvec sort2 = sort_index(index2);
  auto end_sort_test = std::chrono::high_resolution_clock::now();
  total_sort_test += std::chrono::duration<double>(end_sort_test - start_sort_test).count();
  vec Y_sort1 = Y(sort1);
  vec Y_sort2 = Y(sort2);
  // 3. We proceed similarly for the weights
  vec pi_sort1 = pi(sort1);
  vec pi_sort2 = pi(sort2);
  vec Y_pi_sort1 = Y_sort1%pi_sort1;
  vec Y_pi_sort2 = Y_sort2%pi_sort2;
  // 4. We can compute the ranks
  vec rank1 = cumsum(pi_sort1) - pi_sort1/2;
  vec rank2 = cumsum(pi_sort2) - pi_sort2/2;
  vec Obj(2);
  Obj(0) = as_scalar(Y_pi_sort1.t()*rank1);
  Obj(1) = as_scalar(Y_pi_sort2.t()*rank2);
  double Fit = max(Obj);
  double pen = ::fabs(Fit*(accu(abs(theta1))-1));

  auto end_total_test = std::chrono::high_resolution_clock::now();
  total_time_test += std::chrono::duration<double>(end_total_test - start_total_test).count();
  call_count_test++;

  if (call_count_test % 5000 == 0) {
    Rcpp::Rcout << "Fitness_cpp avg time: " << (total_time_test / call_count_test)
                << " sec | MatMul: " << (total_matmul_test / call_count_test)
                << " | Sort: " << (total_sort_test / call_count_test)
                << " over " << call_count_test << " calls.\n";
  }

  return Fit-pen;
}
