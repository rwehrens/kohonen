#include <Rcpp.h>
typedef double (*DistanceFunctionPtr)(double *, double *, int, int);

double brayCurtisDissim(double *data, double *codes, int n, int nNA) {
  if (nNA > 0) return NA_REAL;
  
  double num = 0.0, denom = 0.0;
  for (int i = 0; i < n; i++) {
    num += std::abs(data[i] - codes[i]);
    denom += data[i] + codes[i];
  }
  
  return num/denom;
}

// [[Rcpp::export]]
Rcpp::XPtr<DistanceFunctionPtr> BrayCurtis() {
  return Rcpp::XPtr<DistanceFunctionPtr>(new DistanceFunctionPtr(&brayCurtisDissim));
}
