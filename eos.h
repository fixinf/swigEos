/*
 * eos.h
 *
 *  Created on: 09 �������� 2014 ��.
 *      Author: const
 */



#ifndef EOS_H_
#define EOS_H_

#include "setconst.h"
#include <vector>
#include "gsl/gsl_vector_double.h"
#include <cmath>
#include "constants.h"


extern void stepE(double n, double * init, int dimInit, double * f_init, int dimF_init, double * out, int dimOut, int iter, set_const *);

extern void stepE_rho(double n, double * init, int dimInit, double * f_init, int dimF_init,
		double * out, int dimOut, int iter, double mu_init, set_const *);


void potentials(double * n, int dimN, double * out, int dimOut, set_const * C);

//extern double stepF(var, set_const *);
extern double _E(double * n, int dimN, set_const *, double * inplace = NULL, int dim_inplace = 0);
extern double E(double * n, int dimN, set_const *, double * inplace = NULL, int dim_inplace = 0);
extern double E_rho(double * n, int dimN, double mu_c, set_const * C, double * inplace = NULL, int dim_inplace = 0);
extern double sum(std::vector<double> x);
extern double kineticInt(double n, double m, double f);
extern double p_f(double n, double gamma);



namespace calc{
	extern double mu(double * n, int dimN, int i, set_const * C);
	extern double mu_deriv(double *n, int dimN, int i, double mu_c, set_const *C);
	extern double mu_rho(double * n,  int dimN, int i, double mu_c, set_const * C);
	extern void fun_n_eq_rho_anal(double * p, double * hx, int m, int n, void * adata);
}

extern float sumTest(double * in, int n);
extern float sumTest2(double * in, int n, double * in2, int n2);
extern void solveF(double n, double E, double P, double * init, int dimInit, double * out, int dimOut, set_const * C);
#endif /* EOS_H_ */
