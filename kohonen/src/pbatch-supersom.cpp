/* batch-supersom.c
   Batch supervised SOMs for data fusion.
   Authors: Ron Wehrens and Johannes Kruisselbrink
*/

#include <Rcpp.h>
#include <Rmath.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "kohonen.h"
#include "distance-functions.h"
#include "neighbourhood-functions.h"

#define RANDIN  GetRNGstate()
#define RANDOUT PutRNGstate()
#define UNIF unif_rand()

#ifdef _OPENMP
// [[Rcpp::plugins(openmp)]]
#endif
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
  )
{
  #ifdef _OPENMP

  int numObjects = data.ncol(),   /* number of objects */
    numLayers = numVars.size(),   /* number of layers */
    numCodes = codes.ncol(),      /* number of units in the map */
    totalVars = data.nrow(),      /* total number of variables sum(numVars) */
    cd,                           /* counter over units */
    i,                            /* randomly drawn object */
    j,                            /* counter over variables */
    l,                            /* counter over layers */
    m,                            /* counter over epochs */
    nearest;
  double dist, tmp, radius;

  Rcpp::IntegerVector offsets(numLayers);
  Rcpp::NumericMatrix changes(numLayers, numEpochs);
  Rcpp::NumericMatrix codeSums(totalVars, numCodes);
  Rcpp::NumericVector codeWeights(numCodes);

  double
    *pCodes = REAL(codes),
    *pWeights = REAL(weights),
    *pCodeSums = REAL(codeSums),
    *pCodeWeights = REAL(codeWeights),
    *pChanges = REAL(changes),
    *pData = REAL(data),
    *pNeighbourhoodDistances = REAL(neighbourhoodDistances),
    *pObject;
  int
    *pOffsets = INTEGER(offsets),
    *pNumVars = INTEGER(numVars),
    *pNumNAs = INTEGER(numNAs);

  std::vector<DistanceFunctionPtr> distanceFunctionPtrs =
    GetDistanceFunctions(distanceFunctions);

  /* Create the neighborhood influence function pointer. */
  NeighbourhoodFunctionPtr neighbourhoodFunctionPtr =
    CreateNeighbourhoodFunction((NeighbourhoodFunctionType)neighbourhoodFct);;

  /* Compute the layer data offsets and the total object length. */
  totalVars = 0;
  for (l = 0; l < numLayers; l++) {
    offsets[l] = totalVars;
    totalVars += numVars[l];
  }

  /* Set number of cores for parallel code blocks */
  int availablecores = omp_get_num_procs();
  if (numCores <= 0 || numCores > availablecores) {
     numCores = availablecores;
  }
  omp_set_num_threads(numCores);

  RANDIN;

  #pragma omp parallel private(cd, l, j, nearest, dist, tmp)
  {
    #pragma omp single
    {
      m = 0;
    }

    std::vector<double> threadCodeWeights(numCodes);
    std::vector<double> threadCodeSums(numCodes * totalVars);
    std::vector<double> threadChanges(numLayers);

    /* Outer loop: number of epochs */
    while (m < numEpochs) {

      #pragma omp single
      {
        /* Linear decays for radius */
        radius =
          radii[0] - (radii[0] - radii[1]) * ((double)m / (double)numEpochs);
        if (radius < EPS) {
          radius = EPS;
        }

        /* Initialize new codes vectors and weights */
        std::fill(codeWeights.begin(), codeWeights.end(), 0.0);
        std::fill(codeSums.begin(), codeSums.end(), 0.0);
      }

      /* Initialize new codes vectors and weights */
      std::fill(threadCodeWeights.begin(), threadCodeWeights.end(), 0.0);
      std::fill(threadCodeSums.begin(), threadCodeSums.end(), 0.0);
      std::fill(threadChanges.begin(), threadChanges.end(), 0.0);

      /* Loop: number of objects */
      #pragma omp for
      for (i = 0; i < numObjects; i++) {

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
          dist);

        if (nearest < 0) {
          ::Rf_error("No nearest neighbour found...");
        }

        /* Update changes */
        for (l = 0; l < numLayers; l++) {
          dist = 0.0;
          for (j = 0; j < numVars[l]; j++) {
            if (!std::isnan(data[i * totalVars + offsets[l] + j])) {
              tmp = data[i * totalVars + offsets[l] + j] -
                codes[nearest * totalVars + offsets[l] + j];
              dist += tmp * tmp;
            }
          }
          if (numNAs[i * numLayers + l] > 0) {
            dist = dist * numVars[l] / (numVars[l] - numNAs[i * numLayers + l]);
          }
          threadChanges[l] += dist;
        }

        /* Accumulate sums and weights */
        for (cd = 0; cd < numCodes; cd++) {
          tmp = neighbourhoodFunctionPtr(
            pNeighbourhoodDistances[numCodes * nearest + cd], radius);
          if (tmp > 0) {
            threadCodeWeights[cd] += tmp;
            for (j = 0; j < totalVars; j++) {
              if (!std::isnan(data[i * totalVars + j])) {
                threadCodeSums[cd * totalVars + j] +=
                  tmp * data[i * totalVars + j];
              }
            }
          }
        }
      }

      #pragma omp critical
      {
        for (cd = 0; cd < numCodes; cd++) {
          if (threadCodeWeights[cd] > 0) {
            codeWeights[cd] += threadCodeWeights[cd];
            for (j = 0; j < totalVars; j++) {
              pCodeSums[cd * totalVars + j] +=
                threadCodeSums[cd * totalVars + j];
            }
          }
        }
        for (l = 0; l < numLayers; l++) {
          pChanges[m * numLayers + l] += threadChanges[l];
        }
      }

      #pragma omp barrier
      #pragma omp single
      {
        /* Update all maps */
        for (cd = 0; cd < numCodes; cd++) {
          if (codeWeights[cd] > 0) {
            for (j = 0; j < totalVars; j++) {
              pCodes[cd * totalVars + j] =
                pCodeSums[cd * totalVars + j] / pCodeWeights[cd];
            }
          }
        }

        /* Mean of the nearest layer distances of this iteration */
        for (l = 0; l < numLayers; l++) {
          pChanges[m * numLayers + l] =
            sqrt(pChanges[m * numLayers + l] / numVars[l]) / numObjects;
        }

        m++;
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("codes") = codes,
    Rcpp::Named("changes") = changes);

  RANDOUT;
  
  #else
    ::Rf_warning("OpenMP not available: cannot run in parallel mode");
    return RcppBatchSupersom(data, codes, numVars, weights, distanceFunctions,
			     numNAs, neighbourhoodDistances, neighbourhoodFct,
			     radii, numEpochs);  
  #endif
}
