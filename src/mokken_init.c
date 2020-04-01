#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern int GeneticAlgorithm(int *population, int *itercount, int *popsize, int *nitem,
                     int *npers, int *iter, double *pxover, double *pmutation, double *critval,
                     double *alpha, double *VAR, double *MAXVAR, double *SijMatrix, double *fitness);

static const R_CMethodDef CEntries[] = {
  {"GeneticAlgorithm", (DL_FUNC) &GeneticAlgorithm, 14},
  {NULL, NULL, 0}
};

void R_init_mokken(DllInfo* info) {
	R_registerRoutines(info, CEntries, NULL, NULL, NULL);
	R_useDynamicSymbols(info, TRUE);
	R_forceSymbols(info, TRUE);
}
