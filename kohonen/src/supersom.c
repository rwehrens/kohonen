/* supersom.c: supervised SOMs for data fusion, started Sept 15, 2006

   This version: April 13, 2007. 

   Addition of handling NAs, following
   code in R_euclidean (distance.c). We ignore individual NAs and
   rescale the distance according to the number of non-NA
   variables. The R-function supersom should make sure that every object at
   least has one non-NA value for every matrix. There should be no NAs
   in the codebook vectors, nor in the unit.distances that are
   calculated here: if an object has more NAs than allowed (indicated
   as a fraction of the total number), it is removed from the data set
   in the pre-C R code. It can be mapped later using the whatmaps
   argument in the map.kohonen function.

   Addition April 2007: R code counts the number of NAs before calling
   this function, and provides it as an argument.

   Author: Ron Wehrens
*/

/* Oct. 18 2007: copied the check for equality of distances from BDR's
   class library */

#include <R_ext/Arith.h>

#include <R.h>
#include <Rmath.h>

#define RANDIN  GetRNGstate()
#define RANDOUT PutRNGstate()
#define UNIF unif_rand()

#define EPS 1e-4                /* relative test of equality of distances */


void supersom(double *data, double *codes,
	      double *nhbrdist,
	      double *alphas, double *radii,
	      double *weights, double *changes,
	      double *unitdistances, double *maxdists, /* working arrays */
	      Sint *pn, Sint *pnmat, Sint *nvar, Sint *nNA,
	      Sint *pncodes, Sint *prlen)
{
  int n = *pn,         /* number of objects */
    nmat = *pnmat,     /* number of data matrices */
    ncodes = *pncodes, /* number of units in the map */
    rlen = *prlen,     /* number of epochs */
    cd,                /* counter over units */
    i,                 /* randomly drawn object */
    j,                 /* counter over variables */ 
    k,                 /* counter over iterations */
    l,                 /* counter over nmat */
    m,                 /* counter over rlen */
    nearest, niter, dataoffset, nind;
  double dist, tmp, threshold, alpha;
  
  RANDIN;

  niter = rlen * n;
  for (k = 0; k < niter; k++) {
    /* select random object */
    i = (int)(n * UNIF);
    /*    fprintf(stderr, "\nObject %d:", i);
	  i = k % n;
	  if (i > 3) return; */
    
    /* calculate distances in all spaces, and keep track of the
       largest values */
    dataoffset = 0;
    /*    fprintf(stderr, "\nMaximal distances: "); */
     for (l = 0; l < nmat; l++) {
      maxdists[l] = 0.0;
      if (l > 0) dataoffset += nvar[l-1];
      
      for (cd = 0; cd < ncodes; cd++) {
	unitdistances[cd + l*ncodes] = 0.0;

	for (j = 0; j < nvar[l]; j++) {
	  if (!ISNAN(data[i + j*n + n*dataoffset])) {
	    tmp = data[i + j*n + n*dataoffset] - 
	      codes[cd + j*ncodes + ncodes*dataoffset];
	    unitdistances[cd + l*ncodes] += tmp * tmp;
	  }
	}
	
	if (nNA[i + l*n] > 0) 
	  unitdistances[cd + l*ncodes] *= nvar[l] / (nvar[l] - nNA[i + l*n]);
	
	unitdistances[cd + l*ncodes] = sqrt(unitdistances[cd + l*ncodes]);

	if (unitdistances[cd + l*ncodes] > maxdists[l]) 
	  maxdists[l] = unitdistances[cd + l*ncodes];
      }

      /*      fprintf(stderr, "%.5lf ", maxdists[l]); */

    }
    
    
    /* scale all distances so that largest value equals 1, and 
       sum to take total distance. Find smallest distance.  */
    for (l = 0; l < nmat; l++) {
      for (cd = 0; cd < ncodes; cd++) {
	/* sqrt is necessary because of the interpretation of the
	   weights: it should be a general distance, not a squared
	   distance */
	unitdistances[cd + l*ncodes] /=  maxdists[l];
      }
    }
    
    /* find overall smallest distance */
    nind = 0; dist = DOUBLE_XMAX; nearest = -1;
    for (cd = 0; cd < ncodes; cd++) {
      tmp = 0.0;
      for (l = 0; l < nmat; l++) 
	tmp += weights[l] * unitdistances[cd + l*ncodes];      
      
      if (tmp <= dist * (1 + EPS)) {
	if (tmp < dist * (1 - EPS)) {
	  nind = 0;
	  nearest = cd;
	} else {
	  if(++nind * UNIF < 1.0) nearest = cd;
	}
	dist = tmp;
      }
    }
    
    if (nearest < 0)
      error("No nearest neighbour found...");
    
    /* linear decays for radius and learning parameter */
    threshold = radii[0] - (radii[0] - radii[1]) * (double) k / (double) niter;
    if (threshold < 1.0) threshold = 0.5;    
    alpha = alphas[0] - (alphas[0] - alphas[1]) * (double)k/(double)niter;
    
    m = (int)(k/n);
        
    /* update all maps */
    dataoffset = 0;
    for (l = 0; l < nmat; l++) {
      if (l > 0) dataoffset += nvar[l-1];

      for (cd = 0; cd < ncodes; cd++) {
	if(nhbrdist[cd + ncodes*nearest] > threshold) continue;

	if (cd == nearest) dist = 0.0;
	
	for(j = 0; j < nvar[l]; j++) {
	  if (!ISNAN(data[i + j*n + n*dataoffset])) {
	    tmp = data[i + j*n + n*dataoffset] - 
	      codes[cd + j*ncodes + ncodes*dataoffset];
	    codes[cd + j*ncodes + ncodes*dataoffset] += tmp * alpha;

	    if (cd == nearest) dist += tmp * tmp;
	  } 
	}
	
	if ((cd == nearest) & (nNA[i + l*n] > 0))
	  dist = dist * nvar[l] / (nvar[l] - nNA[i + l*n]);
      }
      
      changes[m + l*rlen] += dist;
    }
  }

    /* mean difference per variable per object */
  for (m = 0; m < rlen; m++) 
    for (l = 0; l< nmat; l++)
      changes[m + l*rlen] = 
	sqrt(changes[m + l*rlen] / nvar[l]) / n;
  
  RANDOUT;
}
