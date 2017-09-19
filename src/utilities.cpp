#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix cumulativeRows(NumericMatrix M) {
  for(int j = 1; j < M.ncol(); ++j) {
    M(_, j) = M(_, j-1) + M(_, j);
  }
  return M;
}

// [[Rcpp::export]]
NumericVector decideEvents(NumericVector events, LogicalMatrix probs_random) {
  for(int e = probs_random.ncol(); e >= 1; --e) {
    events[probs_random(_, e-1)] = e-1;
  }
  return events;
}
