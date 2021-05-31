#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

  void F77_NAME(mvtdst)(int *n, int *nu, double *lower, double *upper,
                int *infin, double *corr, double *delta,
                int *maxpts, double *abseps, double *releps,
                double *error, double *value, int *inform);

  void C_mvtdst(int *n, int *nu, double *lower, double *upper,
                int *infin, double *corr, double *delta,
                int *maxpts, double *abseps, double *releps,
                double *error, double *value, int *inform, int *rnd);

#ifdef __cplusplus
}
#endif
