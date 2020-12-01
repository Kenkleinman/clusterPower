#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran calls */
extern void F77_NAME(computeparameter)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(der_likelihood_time)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(legendre_handle)(void *, void *, void *, void *, void *);
extern void F77_NAME(linearpower_notime_subroutine)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(syminverse)(void *, void *, void *);
extern void F77_NAME(vectorsquare)(void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
  {"computeparameter",              (DL_FUNC) &F77_NAME(computeparameter),               8},
  {"der_likelihood_time",           (DL_FUNC) &F77_NAME(der_likelihood_time),           18},
  {"legendre_handle",               (DL_FUNC) &F77_NAME(legendre_handle),                5},
  {"linearpower_notime_subroutine", (DL_FUNC) &F77_NAME(linearpower_notime_subroutine), 12},
  {"syminverse",                    (DL_FUNC) &F77_NAME(syminverse),                     3},
  {"vectorsquare",                  (DL_FUNC) &F77_NAME(vectorsquare),                   3},
  {NULL, NULL, 0}
};

void R_init_clusterPower(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}