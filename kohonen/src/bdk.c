/* bdk.c: supervised SOMs, started Oct 10, 2005
   This version: Oct 10, 2005
   Author: Ron Wehrens

   BDK differs from XYF in that two separate updates are performed for
   X and Y: X is updated according to the distance in the Y-map, and
   vice versa.
*/

/* Oct. 18 2007: copied the check for equality of distances from BDR's
   class library */


#include <R.h>

#define RANDIN  GetRNGstate()
#define RANDOUT PutRNGstate()
#define UNIF unif_rand()

#define EPS 1e-4                /* relative test of equality of distances */



void BDK_Eucl(double *data, double *Ys, 
	      double *codes, double *codeYs,
	      double *nhbrdist,
	      double *alphas, double *radii, 
	      double *xweight,
	      double *changes,
	      double *xdists, double *ydists, /* working arrays */
	      Sint *pn, Sint *ppx, Sint *ppy, 
	      Sint *pncodes, Sint *prlen)
{
  int n = *pn, py = *ppy, px = *ppx, ncodes = *pncodes, rlen = *prlen;
  int cd, i, j, k, l, xnearest, ynearest, niter, nind;
  double dm, xdist, ydist, dist, tmp, maxx, maxy, distwght, 
    alpha, threshold;

  RANDIN;

  niter = rlen * n;
  for (k = 0; k < niter; k++) {
    /* i is a counter over objects in data, cd is a counter over units
       in the map, and j is a counter over variables */
    i = (int)(n * UNIF);
    dm = DOUBLE_XMAX;

    maxx = maxy = 0.0;
    /* calculate distances in x and y spaces, and keep track of the
       largest values */
    for (cd = 0; cd < ncodes; cd++) {
      xdist = ydist = 0.0;
      for (j = 0; j < px; j++) {
	tmp = data[i + j*n] - codes[cd + j*ncodes];
	xdist += tmp * tmp;
      }
      xdists[cd] = sqrt(xdist);
      if (xdists[cd] > maxx) maxx = xdists[cd];

      for (j = 0; j < py; j++) {
	tmp = Ys[i + j*n] - codeYs[cd + j*ncodes]; 
	ydist += tmp * tmp;
      }
      ydists[cd] = sqrt(ydist);
      if (ydists[cd] > maxy) maxy = ydists[cd];
    }
    
    /* distwght goes linearly to .5 during training */
    distwght = *xweight - (*xweight - 0.5) * k / niter; 
    
    /* scale x and y distances so that largest value in both cases is 1;
       find nearest in x for updating Y */
    nind = 0; dist = DOUBLE_XMAX; xnearest = -1;
    for (cd = 0; cd < ncodes; cd++) {
      xdists[cd] /= maxx;
      ydists[cd] /= maxy;
      tmp = distwght * xdists[cd] + (1.0 - distwght) * ydists[cd];

      if (tmp <= dist * (1 + EPS)) {
	if (tmp < dist * (1 - EPS)) {
	  nind = 0;
	  xnearest = cd;
	} else {
	  if(++nind * UNIF < 1.0) xnearest = cd;
	}
	dist = tmp;
      }
    }

    /* find nearest in y for updating X */
    nind = 0; dist = DOUBLE_XMAX; ynearest = -1;
    for (cd = 0; cd < ncodes; cd++) {
      tmp = (1.0 - distwght) * xdists[cd] + distwght * ydists[cd];

      if (tmp <= dist * (1 + EPS)) {
	if (tmp < dist * (1 - EPS)) {
	  nind = 0;
	  ynearest = cd;
	} else {
	  if(++nind * UNIF < 1.0) ynearest = cd;
	}
	dist = tmp;
      }
    }

    if ((xnearest < 0) | (ynearest < 0))
      error("No nearest neighbour found...");

    /* linear decays for radius and learning parameter */
    threshold = radii[0] - (radii[0] - radii[1]) * (double) k / (double) niter;
    if (threshold < 1.0) threshold = 0.5;
    alpha = alphas[0] - (alphas[0] - alphas[1]) * (double)k/(double)niter;

    l = (int)(k/n);

    /* update all X codes within threshold of 'ynearest' */
    for (cd = 0; cd < ncodes; cd++) {
      if(nhbrdist[cd + ncodes*ynearest] > threshold) continue;

      for(j = 0; j < px; j++) {
	tmp = data[i + j*n] - codes[cd + j*ncodes];
	codes[cd + j*ncodes] += tmp * alpha;

	if (cd == ynearest) changes[l] += tmp*tmp;
      }
    }

    /* update all Y codes within threshold of 'xnearest' */
    for (cd = 0; cd < ncodes; cd++) {
      if(nhbrdist[cd + ncodes*xnearest] > threshold) continue;
      
      for(j = 0; j < py; j++) {
	tmp = Ys[i + j*n] - codeYs[cd + j*ncodes];	
	codeYs[cd + j*ncodes] += tmp * alpha;
	
	if (cd == xnearest) changes[l+rlen] += tmp*tmp;
      }
    }
  }

  for (l = 0; l < rlen; l++) {
    changes[l] /= n*px*px;  /* mean difference per variable per object */
    changes[l+rlen] /= n*py*py;
  }

  RANDOUT;
}


void BDK_Tani(double *data, double *Ys, 
	      double *codes, double *codeYs,
	      double *nhbrdist,
	      double *alphas, double *radii, 
	      double *xweight,
	      double *changes,
	      double *xdists, double *ydists, /* working arrays */
	      Sint *pn, Sint *ppx, Sint *ppy, 
	      Sint *pncodes, Sint *prlen)
{
  int n = *pn, py = *ppy, px = *ppx, ncodes = *pncodes, rlen = *prlen;
  int cd, i, j, k, l, xnearest, ynearest, niter, nind;
  double dm, xdist, ydist, dist, tmp, maxx, maxy, distwght, 
    alpha, threshold;

  RANDIN;

  niter = rlen * n;
  for (k = 0; k < niter; k++) {
    /* i is a counter over objects in data, cd is a counter over units
       in the map, and j is a counter over variables */
    i = (int)(n * UNIF);
    dm = DOUBLE_XMAX;

    maxx = maxy = 0.0;
    /* calculate distances in x and y spaces, and keep track of the
       largest values */
    for (cd = 0; cd < ncodes; cd++) {
      xdist = ydist = 0.0;
      for (j = 0; j < px; j++) {
	tmp = data[i + j*n] - codes[cd + j*ncodes];
	xdist += tmp * tmp;
      }
      xdists[cd] = sqrt(xdist);
      if (xdists[cd] > maxx) maxx = xdists[cd];

      /* Tanimoto distance */
      tmp = 0;
      for (j = 0; j < py; j++) {
	if ((Ys[i + j*n] > .5 && codeYs[cd + j*ncodes] < .5) ||
	    (Ys[i + j*n] <= .5 && codeYs[cd + j*ncodes] >= .5))
	  tmp += 1.0;
      }
      ydists[cd] = tmp/(double)py;
    }
    
    /* distwght goes linearly to .5 during training */
    distwght = *xweight - (*xweight - 0.5) * k / niter; 
    
    /* scale x distance so that largest value is 1; automatically true
       for Tanimoto distances; find nearest in x for updating Y */
    nind = 0; dist = DOUBLE_XMAX; xnearest = -1;
    for (cd = 0; cd < ncodes; cd++) {
      xdists[cd] /= maxx;
      tmp = distwght * xdists[cd] + (1.0 - distwght) * ydists[cd];

      if (tmp < dist * (1 + EPS)) {
	if (tmp < dist * (1 - EPS)) {
	  nind = 0;
	  xnearest = cd;
	} else {
	  if(++nind * UNIF < 1.0) xnearest = cd;
	}
	dist = tmp;
      }
    }
    
    /* find nearest in y for updating X */
    nind = 0; dist = DOUBLE_XMAX; ynearest = -1;
    for (cd = 0; cd < ncodes; cd++) {
      tmp = (1.0 - distwght) * xdists[cd] + distwght * ydists[cd];
      
      if (tmp <= dist * (1 + EPS)) {
	if (tmp < dist * (1 - EPS)) {
	  nind = 0;
	  ynearest = cd;
	} else {
	  if(++nind * UNIF < 1.0) ynearest = cd;
	}
	dist = tmp;
      }
    }
    
    if ((xnearest < 0) | (ynearest < 0))
      error("No nearest neighbour found...");

    /* linear decays for radius and learning parameter */
    threshold = radii[0] - (radii[0] - radii[1]) * (double) k / (double) niter;
    if (threshold < 1.0) threshold = 0.5;
    alpha = alphas[0] - (alphas[0] - alphas[1]) * (double)k/(double)niter;

    l = (int)(k/n);

    /* update all X codes within threshold of 'ynearest' */
    for (cd = 0; cd < ncodes; cd++) {
      if(nhbrdist[cd + ncodes*ynearest] > threshold) continue;

      for(j = 0; j < px; j++) {
	tmp = data[i + j*n] - codes[cd + j*ncodes];
	codes[cd + j*ncodes] += tmp * alpha;

	if (cd == ynearest) changes[l] += tmp*tmp;
      }
    }

    /* update all Y codes within threshold of 'xnearest' */
    for (cd = 0; cd < ncodes; cd++) {
      if(nhbrdist[cd + ncodes*xnearest] > threshold) continue;
      
      for(j = 0; j < py; j++) {
	tmp = Ys[i + j*n] - codeYs[cd + j*ncodes];	
	codeYs[cd + j*ncodes] += tmp * alpha;
	
	if (cd == xnearest) changes[l+rlen] += tmp*tmp;
      }
    }
  }
  
  for (l = 0; l < rlen; l++) {
    changes[l] /= n*px*px;  /* mean difference per variable per object */
    changes[l+rlen] /= n*py*py;
  }

  RANDOUT;
}

