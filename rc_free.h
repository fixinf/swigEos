#include "setconst.h"

extern double E_rho_free(double* n, int dimN,  double nc, set_const * C, double * inplace = NULL, int dim_inplace = 0);
extern void step_free(double n, double * init, int dimInit, double * f_init, int dimF_init, double * out, int dimOut, int iter, set_const* C);
extern double wrap_eqf_rc_free(double f, double * n, int dimN, double rc2, set_const * C);
void wrap_fun_free(double n_tot,
                   double * init, int dimInit, 
                   set_const * C, 
                   double * out, int dimOut);