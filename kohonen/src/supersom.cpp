/* supersom.c: supervised SOMs for data fusion, started Sept 15, 2006
   Revision, July 2016

   Rcpp version of supersoms

   Addition of handling NAs, following code in R_euclidean
   (distance.c). We ignore individual NAs and rescale the distance
   according to the number of non-NA variables. The R-function
   supersom should make sure that every object at least has one non-NA
   value for every matrix. There should be no NAs in the codebook
   vectors, nor in the unit.distances that are calculated here: if an
   object has more NAs than allowed (indicated as a fraction of the
   total number), it is removed from the data set in the pre-C R
   code. It can be mapped later using the whatmaps argument in the
   map.kohonen function. 
   NOTE: this means that our distance functions never should return NAs...

   Author: Ron Wehrens and Johannes Kruisselbrink
*/

/* Oct. 18 2007: copied the check for equality of distances from BDR's
   class library */

#include <Rcpp.h>
#include <Rmath.h>

#include "kohonen.h"
#include "distance-functions.h"
#include "neighbourhood-functions.h"

#define RANDIN  GetRNGstate()
#define RANDOUT PutRNGstate()
#define UNIF unif_rand()

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
  )
{
  int
    numObjects = data.ncol(),     /* number of objects */
    numLayers = numVars.size(),   /* number of layers */
    numCodes = codes.ncol(),      /* number of units in the map */
    totalVars = data.nrow(),      /* total number of variables sum(numVars) */
    cd,                           /* counter over units */
    i,                            /* randomly drawn object */
    j,                            /* counter over variables */
    k,                            /* counter over iterations */
    l,                            /* counter over layers */
    m,                            /* counter over epochs */
    nearest, totalIters, curIter = 0;
  double distance, tmp, threshold, alpha;

  Rcpp::IntegerVector offsets(numLayers);
  Rcpp::NumericMatrix changes(numLayers, numEpochs);

  double
    *pCodes = REAL(codes),
    *pWeights = REAL(weights),
    *pChanges = REAL(changes),
    *pData = REAL(data),
    *pNeighbourhoodDistances = REAL(neighbourhoodDistances),
    *pObject;
  int
    *pOffsets = INTEGER(offsets),
    *pNumVars = INTEGER(numVars),
    *pNumNAs = INTEGER(numNAs);
  
  /* Get the distance function pointers. */
  std::vector<DistanceFunctionPtr> distanceFunctionPtrs =
    GetDistanceFunctions(distanceFunctions);
  
  /* Create the neighborhood influence function pointer. */
  NeighbourhoodFunctionPtr neighbourhoodFunctionPtr =
    CreateNeighbourhoodFunction((NeighbourhoodFunctionType)neighbourhoodFct);

  /* Compute the layer data offsets and the total object length. */
  totalVars = 0;
  for (l = 0; l < numLayers; l++) {
    offsets[l] = totalVars;
    totalVars += numVars[l];
  }

  totalIters = numEpochs * numObjects;

  RANDIN;
  
  /* Outer loop: number of iterations */
  for (m = 0; m < numEpochs; m++) {

    /* Inner loop: loop over (bootstrapped) objects */
    for (k = 0; k < numObjects; k++) {

      /* Select random object */
      i = (int)(numObjects * UNIF);

      /* Find best matching unit index and distance */
      pObject = &pData[i * totalVars];

      /* Find best matching unit index and distance */
      FindBestMatchingUnit(
        pObject,
        pCodes,
        pOffsets,
        &pNumNAs[i * numLayers],
        numCodes,
        numLayers,
        pNumVars,
        totalVars,
        distanceFunctionPtrs,
        pWeights,
        nearest,
        distance);

      if (nearest < 0) {
        ::Rf_error("No nearest neighbour found...");
      }

      /* Linear decays for radius and learning parameter */
      tmp = (double)(curIter) / (double)(totalIters);
      threshold = radii[0] - (radii[0] - radii[1]) * tmp;
      if (threshold < 1.0) {
        threshold = 0.5;
      }
      alpha = alphas[0] - (alphas[0] - alphas[1]) * tmp;

      /* Update changes */
      for (l = 0; l < numLayers; l++) {
        distance = 0.0;
        for (j = pOffsets[l]; j < pOffsets[l] + pNumVars[l]; j++) {
          if (!std::isnan(pObject[j])) {
            tmp = pObject[j] - pCodes[nearest * totalVars + j];
            distance += tmp * tmp;
          }
        }
        if (pNumNAs[i * numLayers + l] > 0) {
          distance = distance * pNumVars[l] / (pNumVars[l] - pNumNAs[i * numLayers + l]);
        }
        pChanges[m * numLayers + l] += distance;
      }

      /* Update all maps */
      for (cd = 0; cd < numCodes; cd++) {
        tmp = neighbourhoodFunctionPtr(pNeighbourhoodDistances[numCodes * nearest + cd], threshold);
        if (tmp > 0) {
          for (j = 0; j < totalVars; j++) {
            if (!std::isnan(pObject[j])) {
              pCodes[cd * totalVars + j] += tmp * alpha * (pObject[j] - pCodes[cd * totalVars + j]);
            }
          }
        }
      }

      ++curIter;
    }

    /* Mean of the nearest layer distances of this iteration */
    for (l = 0; l < numLayers; l++) {
      pChanges[m * numLayers + l] =
        sqrt(pChanges[m * numLayers + l] / pNumVars[l]) / numObjects;
    }
  }

  RANDOUT;

  return Rcpp::List::create(
    Rcpp::Named("codes") = codes,
    Rcpp::Named("changes") = changes);
}
