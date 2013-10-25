#include <R.h>

#define RANDIN  GetRNGstate()
#define RANDOUT PutRNGstate()
#define UNIF unif_rand()

#define EPS 1e-4                /* relative test of equality of distances */

/* SOM function for the kohonen package, modelled after BDR's function
   in the class package. */
/* Oct. 18 2007: copied the check for equality of distances from BDR's
   class library */

void SOM_online(double *data, double *codes, double *nhbrdist,
		double *alphas, double *radii, double *changes,
		Sint *pn, Sint *pp, Sint *pncodes, Sint *prlen)
{
  int n = *pn, p = *pp, ncodes = *pncodes, rlen = *prlen;
  int cd, i, j, k, l, nearest, niter, nind;
  double dm, dist, tmp, alpha, threshold;
  
  RANDIN;

  niter = rlen * n;

  for (k = 0; k < niter; k++) {
    /* pick a random data point */
    i = (int)(n * UNIF);

    /* find the nearest code 'near' */
    nind = 0; dm = DOUBLE_XMAX; nearest = -1;
    for (cd = 0; cd < ncodes; cd++) {
      dist = 0.0;
      for (j = 0; j < p; j++) {
	tmp = data[i + j*n] - codes[cd + j*ncodes];
	dist += tmp * tmp;
      }

      if (dist <= dm * (1 + EPS)) {
	if (dist < dm * (1 - EPS)) {
	  nind = 0;
	  nearest = cd;
	} else {
	  if(++nind * UNIF < 1.0) nearest = cd;
	}
	dm = dist;
      }
    }

    if (nearest < 0)
      error("No nearest neighbour found...");
    
    /* update all codes within threshold of 'nearest'. Linear decrease
       for both radius and learning parameter. */ 
    threshold = radii[0] - (radii[0] - radii[1]) * (double)k/(double)niter;
    if (threshold < 1.0) threshold = 0.5;
    alpha = alphas[0] - (alphas[0] - alphas[1]) * (double)k/(double)niter;

    l = (int)(k/n);

    for (cd = 0; cd < ncodes; cd++) {
      if(nhbrdist[cd + ncodes*nearest] > threshold) continue;

      for(j = 0; j < p; j++) {
	tmp = data[i + j*n] - codes[cd + j*ncodes];
	codes[cd + j*ncodes] += tmp * alpha;

	if (cd == nearest) changes[l] += tmp * tmp;
      }
    }
  }

  for (l = 0; l < rlen; l++) {
    /* mean difference per variable per object */
    changes[l] = sqrt(changes[l]/p)/n;  
  }

  RANDOUT;
}
