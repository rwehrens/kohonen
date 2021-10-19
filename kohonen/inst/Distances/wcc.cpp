/* C++ version of the WCCd dissimilarity, originally written in C for
   the wccsom package */
#include <Rcpp.h>

typedef double (*DistanceFunctionPtr)(double *, double *, int, int);
double wcc_corr(double *, double *, int);
double wcc_autocorr(double *, int);

double WCCd20(double *p1, double *p2, int np, int nNA) {
  int i;
  double crosscov, wght, acor1, acor2, trwdth = 20.0;

  acor1 = wcc_autocorr(p1, np);
  acor2 = wcc_autocorr(p2, np);
  crosscov = wcc_corr(p1, p2, np);
  
  for (i=1; i<trwdth; i++) {
    wght = 1.0 - (double) i / trwdth;
    crosscov += (wcc_corr(p1, p2+i, np-i)*wght);
    crosscov += (wcc_corr(p1+i, p2, np-i)*wght);
  }

  crosscov = 1.0 - crosscov / (acor1 * acor2);
    
  return crosscov;
}

double wcc_corr(double *p1, double *p2, int npoints)
{
  int i;
  double anum;
  
  anum = 0.0;
  for (i=0; i<npoints; i++)
    anum += p1[i]*p2[i];
  
  return anum;
}

double wcc_autocorr(double *p1, int np) {
  int i;
  double autocov, wght, trwdth = 20.0;

  autocov = wcc_corr(p1, p1, np);
  for (i=1; i<trwdth; i++) {
    wght = 1.0 - (double) i / trwdth;
    autocov += 2.0*(wcc_corr(p1, p1+i, np-i)*wght);
  }

  return(std::sqrt(autocov));
}

/* If you want to use wccd inside R, you can call the function with this
   applyDF wrapper:
   wccvalue <- applyDF(WCCd() pattern1, pattern2, numNA) */

// [[Rcpp::export]]
double applyDF(SEXP xpsexp, Rcpp::NumericVector data,
	       Rcpp::NumericVector codes, int nNA) {
  int n = data.size();

  Rcpp::XPtr<DistanceFunctionPtr> xpfun(xpsexp);
  DistanceFunctionPtr fun = *xpfun;
  return fun(REAL(data), REAL(codes), n, nNA);
}

// [[Rcpp::export]]
Rcpp::XPtr<DistanceFunctionPtr> WCCd() {
  return (Rcpp::XPtr<DistanceFunctionPtr>(new DistanceFunctionPtr(&WCCd20)));
}

