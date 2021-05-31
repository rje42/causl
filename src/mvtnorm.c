#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <Rmath.h>
#include <R_ext/Rdynload.h>

void F77_NAME(mvtdst)(int *n, int *nu, double *lower, double *upper,
              int *infin, double *corr, double *delta,
              int *maxpts, double *abseps, double *releps,
              double *error, double *value, int *inform);

extern void C_mvtdst(int *n, int *nu, double *lower, double *upper,
              int *infin, double *corr, double *delta,
              int *maxpts, double *abseps, double *releps,
              double *error, double *value, int *inform, int *rnd)
{

  if (rnd[0]) GetRNGstate();

  /* call FORTRAN subroutine */
  F77_CALL(mvtdst)(n, nu, lower, upper,
           infin, corr, delta,
           maxpts, abseps, releps,
           error, value, inform);

  if (rnd[0]) PutRNGstate();

}

