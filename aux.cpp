/*
 * aux.cpp
 *
 *  Created on: 06 июля 2014 г.
 *      Author: const
 */

#include "aux.h"
#include "levmar.h"

struct func_f_eq_params{
	double * n;
	int dimN;
	double df;
	set_const * C;
};

void func_f_eq(double * p, double * hx, int m, int _n, void * adata){
	bool debug = 1;
	func_f_eq_params * params = (func_f_eq_params *) adata;
	bool sprime = (params->C->sprime and (m > 1));
	double * n = new double[params->dimN + m];
	if (debug){
		printf("sprime = %i \n", sprime);
	}
	for (int i = 0; i < m; i++){
		n[i] = p[i];
	}
	for (int i = m; i <= params->dimN; i++){
		n[i] = params->n[i-m];
	}
	if (debug) {
		printf("f_eq: n = ");
		for (int i = 0; i < params->dimN + m; i++){
			printf("%f ", n[i]);
		}
		printf("\n");
	}
	double dE;
	for (int i = 0; i < m; i++){
		n[i] += params->df;
		dE = _E(n, params->dimN + m, params->C);
		n[i] -= params->df;
		dE -= _E(n, params->dimN + m, params->C);
		dE /= params->df;
		hx[i] = dE;
		if (debug){
			printf("dE[%i] = %f  ", i, dE);
		}
	}
	if (debug){
		printf("\n");
	}
	delete [] n;
}

double f_eq(double * n, int dimN, set_const * C, double init){
	double opts[5];
	printf("I'm in \n");
	func_f_eq_params p = {n, dimN, 1e-4, C};
	int m = 1;
	double * x = new double[m];
	double * lb = new double[m];
	double * ub = new double[m];
	double * fun = new double[m];
	double info[LM_INFO_SZ];
	//double x[3] = {v.n[0], v.n[1], v.f};
	x[0] = init;
	lb[0] = 0.0;
	ub[0] = 1.0;
	opts[0]= LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-25; opts[3]=1E-20;
	opts[4]= -1e-5;
	int iter = 300;
	dlevmar_bc_dif(func_f_eq, x, NULL, m, m, lb, ub, NULL, iter, opts, info, NULL, NULL, &p);
	double res = x[0];
	delete[] x;
	delete[] fun;
	delete[] lb;
	delete[] ub;
	return res;
}

void f_eq(double * n, int dimN, double * init, int dimInit, double * out, int dimOut, set_const * C){
	double opts[5];
	func_f_eq_params p = {n, dimN, 1e-4, C};
	int m = 1 + C->sprime;
	double * x = new double[m];
	double * lb = new double[m];
	double * ub = new double[m];
	double * fun = new double[m];
	double info[LM_INFO_SZ];
	//double x[3] = {v.n[0], v.n[1], v.f};
	x[0] = init[0];
	if (C->sprime){
		x[1] = init[1];
	}
	lb[0] = 0.0;
	ub[0] = 1.0;
	opts[0]= LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-25; opts[3]=1E-20;
		opts[4]= -1e-5;
	int iter = 300;
	dlevmar_bc_dif(func_f_eq, x, NULL, m, m, lb, ub, NULL, iter, opts, info, NULL, NULL, &p);
	double res = x[0];
	out[0] = x[0];
	if (C->sprime){
		out[1] = x[1];
	}
	delete[] x;
	delete[] fun;
	delete[] lb;
	delete[] ub;
}

double EBind(double * n, int dimN, set_const *C){
	double sum = 0;
	for (int i = 1; i < dimN; ++i) {
		sum += n[i];
	}
	return 135.0*( E(n, dimN, C)/sum - C->M[0]);
}

double K(double n, set_const *C){
	double dn = 1e-3;
	double _n[2] = {(n+dn)/2, (n+dn)/2};
	double f = f_eq(_n, 2, C);
	double n_E[3] = {f, (n+dn)/2, (n + dn)/2};
	double d2E = EBind(n_E, 3, C);

	_n[0] -= dn/2;
	_n[1] -= dn/2;
	n_E[0] = f_eq(_n, 2, C);
	n_E[1] -= dn/2;
	n_E[2] -= dn/2;

	d2E -= 2*EBind(n_E, 3, C);

	_n[0] -= dn/2;
	_n[1] -= dn/2;
	n_E[0] = f_eq(_n, 2, C);
	n_E[1] -= dn/2;
	n_E[2] -= dn/2;

	d2E += EBind(n_E, 3, C);
	d2E /= dn*dn;
	return 9*n*n*d2E;
}



double J(double n, set_const * C){
	double dn = 1e-3;
	double _n[2] = {(n-dn)/2, (n+dn)/2};
	double f = f_eq(_n, 2, C);
	double n_E[3] = {f, (n - dn)/2, (n + dn)/2};
	double d2E = E(n_E, 3, C);

	_n[0] += dn/2;
	_n[1] -= dn/2;
	n_E[0] = f_eq(_n, 2, C);
	n_E[1] += dn/2;
	n_E[2] -= dn/2;

	d2E -= 2*E(n_E, 3, C);

	_n[0] += dn/2;
	_n[1] -= dn/2;
	n_E[0] = f_eq(_n, 2, C);
	n_E[1] += dn/2;
	n_E[2] -= dn/2;

	d2E += E(n_E, 3, C);
	d2E /= (dn/2)*(dn/2);
	return 135.0*n*d2E/8;
}


