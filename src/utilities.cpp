// Copyright 2017 Venelin Mitov
//
// This file is part of patherit.
// 
// patherit is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// patherit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

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
