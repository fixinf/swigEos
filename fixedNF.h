#include "setconst.h"

void stepE_rho_nf(double n, double f, double * init, int dimInit, 
                double * out, int dimOut, int iter, double mu_init, set_const* C);

extern void fun_n_eq_rho_nf(double * p, double * hx, int m, int n, void * adata);

void stepE_nf(double n, double f, double * init, int dimInit, 
                double * out, int dimOut, int iter, double mu_init, set_const* C);

extern void fun_n_eq_nf(double * p, double * hx, int m, int n, void * adata);