/*
 * aux.h
 *
 *  Created on: 06 июля 2014 г.
 *      Author: const
 */

#ifndef AUX_H_
#define AUX_H_
#include "setconst.h"
#include "eos.h"

struct func_f_eq_params{
    double * n;
    int dimN;
    double df;
    set_const * C;
};

struct func_f_eq_params_rho{
    double * n;
    int dimN;
    double df;
    set_const * C;
    double mu_c;
};

extern void f_eq(double * n, int dimN, double * init, int dimInit, double * res, int dimRes, set_const * C);
//extern double f_eq(double * n, int dimN, set_const * C, double init = 1e-1);
extern void func_f_eq(double * p, double * hx, int m, int _n, void * adata);


extern void f_eq_rho(double * n, int dimN, double * init, int dimInit, double * res, int dimRes, double mu_c, set_const * C);
extern void func_f_eq_rho(double * p, double * hx, int m, int _n, void * adata);

extern double K(double n, set_const *C);
extern double EBind(double * n, int dimN, set_const *C);
extern double J(double n, set_const * C);
extern double J(double n, set_const * C, double f);


#endif /* AUX_H_ */
