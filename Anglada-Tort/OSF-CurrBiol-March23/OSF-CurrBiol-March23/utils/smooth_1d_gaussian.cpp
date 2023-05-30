#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector smooth_1d_gaussian(
  const NumericVector data_x,
  const NumericVector data_val,
  const NumericVector probe_x,
  const double sigma_x
) {
  const int data_N = data_x.size();
  const int probe_N = probe_x.size();
  NumericVector probe_val = NumericVector(probe_N);
  
  for (int i = 0; i < probe_N; i ++) {
    double weighted_sum = 0.0;
    double sum_of_weights = 0.0;
    for (int j = 0; j < data_N; j ++) {
      const double dist = fabs(data_x[j] - probe_x[i])/sigma_x;
      const double weight = R::dnorm(dist, 0.0, 1.0, false);
      weighted_sum += weight * data_val[j];
      sum_of_weights += weight;
    }
    probe_val[i] = weighted_sum / sum_of_weights;
  }
  
  return probe_val;
}
