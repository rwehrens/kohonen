#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP kohonen_CreateStdDistancePointer(SEXP, SEXP);
extern SEXP kohonen_CreateStdDistancePointers(SEXP, SEXP);
extern SEXP kohonen_ObjectDistances(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP kohonen_RcppBatchSupersom(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP kohonen_RcppMap(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP kohonen_RcppParallelBatchSupersom(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP kohonen_RcppSupersom(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"kohonen_CreateStdDistancePointer",  (DL_FUNC) &kohonen_CreateStdDistancePointer,   2},
    {"kohonen_CreateStdDistancePointers", (DL_FUNC) &kohonen_CreateStdDistancePointers,  2},
    {"kohonen_ObjectDistances",           (DL_FUNC) &kohonen_ObjectDistances,            5},
    {"kohonen_RcppBatchSupersom",         (DL_FUNC) &kohonen_RcppBatchSupersom,         10},
    {"kohonen_RcppMap",                   (DL_FUNC) &kohonen_RcppMap,                    6},
    {"kohonen_RcppParallelBatchSupersom", (DL_FUNC) &kohonen_RcppParallelBatchSupersom, 11},
    {"kohonen_RcppSupersom",              (DL_FUNC) &kohonen_RcppSupersom,              11},
    {NULL, NULL, 0}
};

void R_init_kohonen(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
