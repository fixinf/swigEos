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

double K_f(double n, double m){
	return 0.;
}

void func_f_eq(double * p, double * hx, int m, int _n, void * adata){
	bool debug = 0;
	func_f_eq_params * params = (func_f_eq_params *) adata;
	bool sprime = (params->C->sprime and (m > 1));
	double * n = new double[params->dimN + m];
	if (debug){
		printf("sprime = %i \n", sprime);
	}
	for (int i = 0; i < m; i++){
		n[i] = p[i];
	}
	for (int i = m; i < m + params->dimN; i++){
		n[i] = params->n[i-m];
	}
	if (debug) {
		printf("f_eq: n = ");
		for (int i = 0; i < params->dimN + m; i++){
			printf("%f ", n[i]);
		}
		printf("\n");
	}
	double df = params->df;
	bool anal = 0;
	if (!anal){
	double dE;
	for (int i = 0; i < m; i++){
//		n[i] += params->df;
//
//		dE = _E(n, params->dimN + m, params->C);
//		n[i] -= params->df;
//		dE -= _E(n, params->dimN + m, params->C);
//		dE /= params->df;
//		hx[i] = dE;
//		n[i] += params->df;

//Increased precision:
//		n[i] += 2*params->df;
//		double dE = -0.25*_E(n, params->dimN + m, params->C);
//		n[i] -= params->df;
//		dE += 2*_E(n, params->dimN + m, params->C);
//		n[i] -= 2*params->df;
//		dE += -2*_E(n, params->dimN + m, params->C);
//		n[i] -= params->df;
//		dE += 0.25*_E(n, params->dimN + m, params->C);
//		n[i] += 2*params->df;
//		dE /= 3*params->df;
//		hx[i] = dE;

//WE NEED MORE!
		n[i] += 3*df;
		double dE = _E(n, params->dimN + m, params->C);
		n[i] -= df;
		dE += -9 * _E(n, params->dimN + m, params->C);
		n[i] -= df;
		dE += 45 * _E(n, params->dimN + m, params->C);
		n[i] -= 2*df;
		dE += -45 * _E(n, params->dimN + m, params->C);
		n[i] -= df;
		dE += 9 * _E(n, params->dimN + m, params->C);
		n[i] -= df;
		dE += -1 * _E(n, params->dimN + m, params->C);
		n[i] += 3*df;
		dE /= 60 * df;
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
	else{

	}
}

double f_eq(double * n, int dimN, set_const * C, double init){
	double opts[5];
	printf("I'm in \n");
	func_f_eq_params p = {n, dimN, 1e-3, C};
	int m = 1;
	double * x = new double[m];
	double * lb = new double[m];
	double * ub = new double[m];
	double * fun = new double[m];
	double info[LM_INFO_SZ];
	//double x[3] = {v.n[0], v.n[1], v.f};
	x[0] = init;
	lb[0] = 0.0;
	ub[0] = C->fmax;
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

//
//	double opts[5];
//	func_f_eq_params p = {n, dimN, 1e-4, C};
//	int m = 1;
//	double * x = new double[m];
//	double * lb = new double[m];
//	double * ub = new double[m];
//	double * fun = new double[m];
//	double info[LM_INFO_SZ];
//	//double x[3] = {v.n[0], v.n[1], v.f};
//	x[0] = init;
//	lb[0] = 0.0;
//	ub[0] = 1.0;
//	opts[0]= LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-25; opts[3]=1E-20;
//		opts[4]= -1e-5;
//	int iter = 300;
//	dlevmar_bc_dif(func_f_eq, x, NULL, m, m, lb, ub, NULL, iter, opts, info, NULL, NULL, &p);
//	double res = x[0];
//	delete[] x;
//	delete[] fun;
//	delete[] lb;
//	delete[] ub;
//	return res;
}

void f_eq(double * n, int dimN, double * init, int dimInit, double * res, int dimRes, set_const * C){
	double opts[5];
	func_f_eq_params p = {n, dimN, 1e-5, C};
	int m = 1 + C->sprime;
	double * x = new double[m];
	double * lb = new double[m];
	double * ub = new double[m];
	double * fun = new double[m];
	double info[LM_INFO_SZ];
	//double x[3] = {v.n[0], v.n[1], v.f};

	for (int i = 0; i < m; i++){
		x[i] = init[i];
		lb[i] = 0.0;
		ub[i] = 1.0;
	}

	//debug only
	double sum = 0;
	bool debug = 0;
	for (int i = 0; i < dimN; i++){
		sum += n[i];
	}

	if (sum > 2.93){
		debug = 0;
	}

	ub[0] = C->fmax;
	ub[1] = 1.;
	opts[0]= LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-25; opts[3]=1E-24;
		opts[4]= -1e-5;
	int iter = 300;
	dlevmar_bc_dif(func_f_eq, x, NULL, m, m, lb, ub, NULL, iter, opts, info, NULL, NULL, &p);

	if (debug){
		printf("info: ");
		for (int i = 0; i < LM_INFO_SZ; i++){
			printf("[%i] = %f ", i, info[i]);
		}
		printf("\n");
	}

	for (int i = 0; i < m; i++){
		res[i] = x[i];
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
return 135.0*( _E(n, dimN, C)/sum - C->M[0]);
}

double K(double n, set_const *C){
	double dn = 1e-2;
	double _n[2] = {(n+2*dn)/2, (n+2*dn)/2};
	double init[1] = {C->f0};
	double out[1];
	int old_sprime = C->sprime;
	C->sprime = 0;
	f_eq(_n, 2, init, 1, out, 1, C);
	double f = out[0];
	double n_E[3] = {f, (n+2*dn)/2, (n + 2*dn)/2};
	double d2E = -EBind(n_E, 3, C);

	_n[0] -= dn/2;
	_n[1] -= dn/2;
	f_eq(_n, 2, init, 1, out, 1, C);
	n_E[0] = out[0];
	n_E[1] -= dn/2;
	n_E[2] -= dn/2;

	d2E += 16.*EBind(n_E, 3, C);

	_n[0] -= dn/2;
	_n[1] -= dn/2;
	f_eq(_n, 2, init, 1, out, 1, C);
	n_E[0] = out[0];
	n_E[1] -= dn/2;
	n_E[2] -= dn/2;

	d2E += -30.*EBind(n_E, 3, C);

	_n[0] -= dn/2;
	_n[1] -= dn/2;
	f_eq(_n, 2, init, 1, out, 1, C);
	n_E[0] = out[0];
	n_E[1] -= dn/2;
	n_E[2] -= dn/2;

	d2E += 16.*EBind(n_E, 3, C);

	_n[0] -= dn/2;
	_n[1] -= dn/2;
	f_eq(_n, 2, init, 1, out, 1, C);
	n_E[0] = out[0];
	n_E[1] -= dn/2;
	n_E[2] -= dn/2;

	d2E += -1.*EBind(n_E, 3, C);
	d2E /= 12*dn*dn;
	C->sprime = old_sprime;
	return 9*n*n*d2E;
}



double J(double n, set_const * C){
	double dn = 1e-3;
	int old_sprime = C->sprime;
	C->sprime = 0;
	double _n[2] = {(n-2*dn)/2, (n+2*dn)/2};
	double out[1];
	double init[1] = {C->f0};
	f_eq(_n, 2, init, 1, out, 1, C);
	double n_E[3] = {out[0], (n - 2*dn)/2, (n + 2*dn)/2};
	double d2E = -1.*_E(n_E, 3, C);

	_n[0] += dn/2;
	_n[1] -= dn/2;

	f_eq(_n, 2, init, 1, out, 1, C);
	n_E[0] = out[0];
	n_E[1] += dn/2;
	n_E[2] -= dn/2;

	d2E += 16*_E(n_E, 3, C);

	_n[0] += dn/2;
	_n[1] -= dn/2;
	f_eq(_n, 2, init, 1, out, 1, C);
	n_E[0] = out[0];
	n_E[1] += dn/2;
	n_E[2] -= dn/2;

	d2E += -30.*_E(n_E, 3, C);

	_n[0] += dn/2;
	_n[1] -= dn/2;
	f_eq(_n, 2, init, 1, out, 1, C);
	n_E[0] = out[0];
	n_E[1] += dn/2;
	n_E[2] -= dn/2;

	d2E += 16.*_E(n_E, 3, C);

	_n[0] += dn/2;
	_n[1] -= dn/2;
	f_eq(_n, 2, init, 1, out, 1, C);
	n_E[0] = out[0];
	n_E[1] += dn/2;
	n_E[2] -= dn/2;

	d2E += -1.*_E(n_E, 3, C);

	d2E /= 12*(dn/2)*(dn/2);
	C->sprime = old_sprime;
	return 135.0*n*d2E/8;
}

double J(double n, set_const * C, double f){
	double dn = 1e-3;
	int old_sprime = C->sprime;
	C->sprime = 0;
	double _n[2] = {(n-dn)/2, (n+dn)/2};
	double out[1];
	double init[1] = {f};
	f_eq(_n, 2, init, 1, out, 1, C);
	double n_E[3] = {out[0], (n - dn)/2, (n + dn)/2};
	double d2E = _E(n_E, 3, C);

	_n[0] += dn/2;
	_n[1] -= dn/2;

	f_eq(_n, 2, init, 1, out, 1, C);
	n_E[0] = out[0];
	n_E[1] += dn/2;
	n_E[2] -= dn/2;

	d2E -= 2*_E(n_E, 3, C);

	_n[0] += dn/2;
	_n[1] -= dn/2;
	f_eq(_n, 2, init, 1, out, 1, C);
	n_E[0] = out[0];
	n_E[1] += dn/2;
	n_E[2] -= dn/2;

	d2E += _E(n_E, 3, C);
	d2E /= (dn/2)*(dn/2);
	C->sprime = old_sprime;
	return 135.0*n*d2E/8;
}

