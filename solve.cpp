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
	p->C->set(x, 5);
	p->C->f0 = p->f0;
	double n[3];
	double n_f[2] = {p->n0/2, p->n0/2};
	n[0] = f_eq(n_f, 2, p->C, p->f0);
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
	double dn = 1e-3;
	n_f[0] += dn/2;
	n_f[1] += dn/2;
	n[0] = f_eq(n_f, 2, p->C);
	n[1] = n_f[0];
	n[2] = n_f[1];
//	printf("n1 = %f, n2 = %f \n", n[1], n[2]);
	double dEBind = EBind(n, 3, p->C);
//	printf("f0 = %f dEbind = %f \n",n[0], dEBind);
	n_f[0] -= dn/2;
	n_f[1] -= dn/2;
	n[0] = f_eq(n_f, 2, p->C);
	n[1] = n_f[0];
	n[2] = n_f[1];
//	printf("n1 = %f, n2 = %f \n", n[1], n[2]);
	dEBind -= EBind(n, 3, p->C);
//	printf("f0 = %f dEbind = %f \n",n[0], dEBind);
	dEBind /= dn;
//	printf("f0 = %f dEbind = %f \n",n[0], dEBind);

	hx[2] = dEBind; //dEBind/dn (n0) = 0;
	hx[3] = p->K0 - K(p->n0, p->C); //K(n0) - K0 = 0
	hx[4] = p->J0 - J(p->n0, p->C);
}

void solve(double n0, double E0, double f0, double K0, double J0, set_const* C) {
	double opts[5];
	solve_params p = {n0, E0, K0, J0, f0, C};
	int m = 5;
	double * x = new double[m];
	double * lb = new double[m];
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

	opts[0]= LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-30; opts[3]=1E-20;
		opts[4]= -1e-5;

	int iter = 1000;
	dlevmar_dif(func_solve, x, NULL, m, m, iter, opts, info, NULL, NULL, &p);

	printf("info: ");
	for (int i = 0; i < LM_INFO_SZ; i++){
		printf(" %i : %f ", i, info[i]);
	}
	printf("\n");

	for (int i = 1; i < m; i++){
		printf(",n%i = %f ", i, x[i]);
	}
	printf("\n");

	func_solve(x, fun, m, m, &p);
	for (int i = 0; i < m; i++){
		printf("f%i = %f  ", i, fun[i]);
	}
	printf("\n");
	delete[] x;
	delete[] fun;
	delete[] lb;
}





