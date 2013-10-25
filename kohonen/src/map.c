/* map.c: calculate distances of objects to unis
   This version: Sept 26, 2006
   Author: Ron Wehrens
*/

#include <R.h>

void mapKohonen(double *data, double *codes,
		Sint *pncodes, Sint *pnd, Sint *pp, double *dists)
{
  int ncodes = *pncodes, nd = *pnd, p = *pp;
  int i, j, cd, distcounter, NAcounter;
  double tmp;
  
  /* i is a counter over objects in data, cd  is a counter over SOM
     units, and j is a counter over variables. */
  distcounter = -1;
  for (i = 0; i < nd; i++) {
    for (cd = 0; cd < ncodes; cd++) {

      dists[++distcounter] = 0.0;
      NAcounter = 0;

      for (j = 0; j < p; j++) {
	if (!ISNAN(data[i + j*nd])) {
	  tmp = data[i + j*nd] - codes[cd + j*ncodes];
	  dists[distcounter] += tmp * tmp;
	} else {
	  NAcounter++;
	}
      }

      if (NAcounter == p) {
	dists[distcounter] = NA_REAL;
      } else { /* correct distance for lower number of non-NA variables */
	if (NAcounter > 0) dists[distcounter] *= p/(p - NAcounter);
      }
    }
  }
}
