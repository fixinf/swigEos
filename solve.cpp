/*
 * solve.cpp
 *
 *  Created on: 06 июля 2014 г.
 *      Author: const
 */

#include "solve.h"
#include "levmar.h"
#include "eos.h"
#include "aux.h"

struct solve_params{
	double n0;
	double E0;
	double K0;
	double J0;
	double f0;
	set_const * C;
};

void func_solve(double * x, double * hx, int m, int _n, void * adata){
	bool debug = 0;
	solve_params * p = (solve_params *) adata;
	if (debug){
			for (int i = 0; i < m; i++){
			printf("x[%i] = %.10f ",i, x[i]);
		}
		printf("\n");
	}

	p->C->set(x, 5);
	p->C->f0 = p->f0;
	double n[3];
	double init[1] = {p->C->f0};
	double out[1];
	double n_f[2] = {p->n0/2, p->n0/2};
	f_eq(n_f, 2, init, 1, out, 1, p->C);
	n[0] = out[0];
	if (debug){
		printf("f = %f \n", n[0]);
	}
	n[1] = n_f[0];
	n[2] = n_f[1];
	if (debug){
		printf("n = %f %f %f\n", n[0], n[1], n[2]);
	}
	hx[0] = EBind(n, 3, p->C) - p->E0; //E(n0) = e0
	if (debug){
		printf("EBind = %f \n", EBind(n, 3, p->C));
	}
	hx[1] = n[0] - p->f0;//m_eff(n0) = M0
	double dn = 1e-5;
	n_f[0] += dn/2;
	n_f[1] += dn/2;
	f_eq(n_f, 2, init, 1, out, 1, p->C);
	n[0] = out[0];
	n[1] = n_f[0];
	n[2] = n_f[1];
//	printf("n1 = %f, n2 = %f \n", n[1], n[2]);

	double dEBind = EBind(n, 3, p->C);
//	printf("f0 = %f dEbind = %f \n",n[0], dEBind);
	n_f[0] -= dn;
	n_f[1] -= dn;
	f_eq(n_f, 2, init, 1, out, 1, p->C);
	n[0] = out[0];
	n[1] = n_f[0];
	n[2] = n_f[1];
//	printf("n1 = %f, n2 = %f \n", n[1], n[2]);
	dEBind -= EBind(n, 3, p->C);
//	printf("f0 = %f dEbind = %f \n",n[0], dEBind);
	dEBind /= 2*dn;
//	printf("f0 = %f dEbind = %f \n",n[0], dEBind);

	hx[2] = dEBind; //dEBind/dn (n0) = 0;
	hx[3] = p->K0 - K(p->n0, p->C); //K(n0) - K0 = 0
	hx[4] = p->J0 - J(p->n0, p->C);
}

int solve(double n0, double E0, double f0, double K0, double J0, set_const* C, int iter, int mu_scale) {
	double opts[5];
	solve_params p = {n0, E0, K0, J0, f0, C};
	int m = 5;
	double * x = new double[m];
	double lb[5] = {0.0, 0.0, 0.0, -1., -1.};
	double ub[5] = {1000. ,1000. , 1000., 1. ,1.};
	double scale[5] = {1e-2, 1e-2, 1e-2, 1e2, 1e2};
	double * fun = new double[m];
	double info[LM_INFO_SZ];
	//double x[3] = {v.n[0], v.n[1], v.f};
	x[0] = C->Cs;
	x[1] = C->Co;
	x[2] = C->Cr;
	x[3] = C->b;
	x[4] = C->c;

	func_solve(x, fun, m, m, &p);
	for (int i = 0; i < m; i++){
		printf("f%i = %f  ", i, fun[i]);
	}
	printf("\n");

	opts[0]= mu_scale*LM_INIT_MU; opts[1]=1E-15; opts[2]=0.0; opts[3]=1E-14;
		opts[4]= -1e-4;

	dlevmar_bc_dif(func_solve, x, NULL, m, m, lb, ub,scale,iter, opts, info, NULL, NULL, &p);
//	dlevmar_dif(func_solve, x, NULL, m, m, iter, opts, info, NULL, NULL, &p);

	printf("info: ");
	for (int i = 0; i < LM_INFO_SZ; i++){
		printf(" %i : %f ", i, info[i]);
	}
	printf("\n");

	printf("Cs = %.10f, Co = %.10f, Cr = %.10f \n b = %.10f, c = %.10f \n", x[0], x[1], x[2], x[3], x[4]);

	func_solve(x, fun, m, m, &p);
	for (int i = 0; i < m; i++){
		printf("f%i = %f  ", i, fun[i]);
	}
	printf("\n");
	delete[] x;
	delete[] fun;
//	delete[] lb;
//	delete[] ub;
//	delete[] scale;
	return info[6];
}





