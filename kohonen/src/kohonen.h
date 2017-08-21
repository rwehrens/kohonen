#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::List RcppSupersom(
  Rcpp::NumericMatrix data,
  Rcpp::NumericMatrix codes,
  Rcpp::IntegerVector numVars,
  Rcpp::NumericVector weights,
  Rcpp::ExpressionVector distanceFunctions,
  Rcpp::IntegerMatrix numNAs,
  Rcpp::NumericMatrix neighbourhoodDistances,
  int neighbourhoodFct,
  Rcpp::NumericVector alphas,
  Rcpp::NumericVector radii,
  int numEpochs
  );

// [[Rcpp::export]]
Rcpp::List RcppBatchSupersom(
  Rcpp::NumericMatrix data,
  Rcpp::NumericMatrix codes,
  Rcpp::IntegerVector numVars,
  Rcpp::NumericVector weights,
  Rcpp::ExpressionVector distanceFunctions,
  Rcpp::IntegerMatrix numNAs,
  Rcpp::NumericMatrix neighbourhoodDistances,
  int neighbourhoodFct,
  Rcpp::NumericVector radii,
  int numEpochs
  );
  
// [[Rcpp::export]]
Rcpp::List RcppParallelBatchSupersom(
  Rcpp::NumericMatrix data,
  Rcpp::NumericMatrix codes,
  Rcpp::IntegerVector numVars,
  Rcpp::NumericVector weights,
  Rcpp::ExpressionVector distanceFunctions,
  Rcpp::IntegerMatrix numNAs,
  Rcpp::NumericMatrix neighbourhoodDistances,
  int neighbourhoodFct,
  Rcpp::NumericVector radii,
  int numEpochs,
  int numCores
  );  

// [[Rcpp::export]]
Rcpp::List RcppMap(
    Rcpp::NumericMatrix data,   /* objects to be mapped */
    Rcpp::IntegerVector numVars,
    Rcpp::IntegerMatrix numNAs,
    Rcpp::NumericMatrix codes,
    Rcpp::NumericVector weights,
    Rcpp::ExpressionVector distanceFunctions);
    
