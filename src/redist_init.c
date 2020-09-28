#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/*
 The following symbols/expressions for .NAME have been omitted
 
 _redist_cppGeneratePartitions
 _redist_countpartitions
 _redist_calcPWDh
 _redist_segregationcalc
 _redist_rsg
 _redist_sample_partition
 _redist_swMH
 _redist_genAlConn
 _redist_findBoundary
 
 Most likely possible values need to be added below.
 */

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP _redist_rsg(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_redist_rsg", (DL_FUNC) &_redist_rsg, 6},
  {NULL, NULL, 0}
};

void R_init_redist(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

