#include <Rcpp.h>
using namespace Rcpp;

// Function to calculate w_ij
// [[Rcpp::export]]
double w_ij(double x, int i, int j) {
  return pow(i*x+1.0,1.0/j);
}
// Modify this function to change the behavior of the model it must be vectorized on x


// Function to calculate w_i
// [[Rcpp::export]]
double w_i(double x, int i, int ntype) {
  double sum_wij = 0.0;
  for (int j = 1; j <= ntype; j++) {
    sum_wij += w_ij(x, i, j);
  }
  return sum_wij;
}

// Function to calculate I_i
// [[Rcpp::export]]
double I_i(double x, double theta, int i, int ntype) {
  if (x > 0) {
    double sum_log_w = 0.0;
    double sum_log_theta_w = 0.0;
    for (int k = 0; k < x; k++) {
      double w = w_i(k, i, ntype);
      sum_log_w += log(w);
      sum_log_theta_w += log(theta + w);
    }
    double log_out = -log(theta + w_i(x, i, ntype)) + sum_log_w - sum_log_theta_w;
    return exp(log_out);
  } else if (x == 0) {
    return exp(-log(theta + w_i(0, i, ntype)));
  }
  return 0.0; // Default return value, although x < 0 is not expected
}

// [[Rcpp::export]]
double S_i(double x, double theta, int i, int ntype){
  double sum_log_w = 0.0;
  double sum_log_theta_w = 0.0;
  for (int k = 0; k <= x; k++) {
    double w = w_i(k, i, ntype);
    sum_log_w += log(w);
    sum_log_theta_w += log(theta + w);
  }
  double log_out = sum_log_w - sum_log_theta_w;
  return exp(log_out);
  
}


// Function to compute the sum from x = 0 to x_max
// [[Rcpp::export]]
double compute_sum(double theta, int i, int j, int ntype, int x_max) {
  double total_sum = 0.0;
  for (int x = 0; x <= x_max; x++) {
    double w = w_ij(x, i, j);
    double I = I_i(x, theta, i, ntype);
    total_sum += w * I;
  }
  return total_sum;
}