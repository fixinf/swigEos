#include "setconst.h"

extern double E_rho_free(double* n, int dimN,  set_const * C, double * inplace = NULL, int dim_inplace = 0);
extern void step_free(double n, double * init, int dimInit, double * f_init, int dimF_init, double * out, int dimOut, int iter, set_const* C);
