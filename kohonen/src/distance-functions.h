#ifndef DISTANCE_H
#define DISTANCE_H

#include <R.h>
#include <Rcpp.h>
#include <cfloat>

#include "kohonen_types.h"

#define EPS 1e-8 /* relative test of equality of distances */

Rcpp::XPtr<DistanceFunctionPtr> CreateStdDistancePointer(int type, bool considerNaNs);
Rcpp::XPtr<DistanceFunctionPtr> CreateNaNDistanceFunctionXPtr(int type);
Rcpp::XPtr<DistanceFunctionPtr> CreateNonNaNDistanceFunctionXPtr(int type);

Rcpp::ExpressionVector CreateStdDistancePointers(const Rcpp::IntegerVector &types);

std::vector<DistanceFunctionPtr> GetDistanceFunctions(const Rcpp::ExpressionVector &distanceFunctionXPtrs);

inline DistanceFunctionPtr AsDistanceFunctionPtr(const Rcpp::XPtr<DistanceFunctionPtr> &xptr) {
  return (*xptr);
}

void FindBestMatchingUnit(
  double *object,
  double *codes,
  int *layerOffsets,
  int *numNAs,
  int numCodes,
  int numLayers,
  int *numVars,
  int totalVars,
  const std::vector<DistanceFunctionPtr> &distanceFunctions,
  double *weights,
  int &index,
  double &distance);
  
// [[Rcpp::export]]
Rcpp::NumericVector ObjectDistances(Rcpp::NumericMatrix data,
				    Rcpp::IntegerVector numVars,
				    Rcpp::IntegerMatrix numNAs,
				    Rcpp::ExpressionVector distanceFunctions,
				    Rcpp::NumericVector weights);
// [[Rcpp::export]]
Rcpp::NumericVector LayerDistances(Rcpp::NumericMatrix data,
				   Rcpp::NumericMatrix codes,
				   Rcpp::IntegerVector uclassif,
				   Rcpp::IntegerVector numVars,
				   Rcpp::IntegerMatrix numNAs,
				   Rcpp::ExpressionVector distanceFunctions,
				   Rcpp::NumericVector weights);
            
double EuclideanDistance(double *dataVector, double *codeVector, int n, int nNA);
double SumOfSquaresDistance(double *dataVector, double *codeVector, int n, int nNA);
double TanimotoDistance(double *dataVector, double *codeVector, int n, int nNA);
double ManhattanDistance(double *dataVector, double *codeVector, int n, int nNA);

double EuclideanDistanceNaN(double *dataVector, double *codeVector, int n, int nNA);
double SumOfSquaresDistanceNaN(double *dataVector, double *codeVector, int n, int nNA);
double TanimotoDistanceNaN(double *dataVector, double *codeVector, int n, int nNA);
double ManhattanDistanceNaN(double *dataVector, double *codeVector, int n, int nNA);

#endif
