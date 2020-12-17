#include <Rcpp.h>
using namespace Rcpp;
//' @title Mean function
//' @description calculate mean of vcector x
//' @param x vector
//' @return mean(double)\code{n}
//' @export
// [[Rcpp::export]]
double meanC(NumericVector x) {
  int n = x.size();
  double total = 0;
  for(int i = 0; i < n; ++i) {
    total += x[i];
  }
  return total / n;
}
