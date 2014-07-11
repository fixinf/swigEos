/*
 * eos.cpp
 *
 *  Created on: 09 июня 2014 г.
 *      Author: const
 */


#include "eos.h"
#include "aux.h"
//#include <gsl/gsl_vector_double.h>
#include <levmar.h>
#include <algorithm>
//#include <cmath>
#include <cstdio>
#include <iterator>
//#include <vector>

//#include "constants.h"
//#include "setconst.h"

namespace calc{
	struct fun_n_eq_params{
		set_const * C;
		double n;
		double f_init;
	};


	double p_f(double n) {
		return pow(3.0 * M_PI * M_PI * D * n, 1.0 / 3.0);
	}

	double mu(double * n, int dimN, int i, set_const * C){
		double dn = 1e-4;
		n[i] += dn;
		double dE = _E(n, dimN, C);
		n[i] -= 2*dn;
		dE -= _E(n, dimN, C);
		n[i] += dn;
//		printf("mu: i = %i, res = %f \n", i, dE/dn);
		return dE/(2*dn);
	}

	void fun_n_eq(double * p, double * hx, int m, int n, void * adata){

		fun_n_eq_params * par = (fun_n_eq_params *) adata;
		set_const * C = par->C;

		double n_sum = 0.0;
		double n_n = par->n;
		double * n_in = new double [m+2]; //input set for _E and mu
		double * n_f = new double [m+1]; //input set for f_eq; actually is {n_n,n_p,...,n_X0}
		for (int i = 0; i < m; i++){
			n_n -= p[i];
			n_in[i+2] = p[i];
			n_f[i+1] = p[i];
		}

		n_f[0] = n_n;

		double f = f_eq(p, m, C, par->f_init);
		par->f_init = f;

		n_in[0] = f;
		n_in[1] = n_n;

		double sum=0, sum_ch=0, sum_o = 0.0, sum_rho = 0.0;
		for (int i = 0; i < m; i++){
			sum += n_f[i];
			sum_ch += n_f[i]*C->Q[i];
			sum_o += n_f[i]*C->X_o[i];
			sum_rho += n_f[i]*C->X_r[i]*C->T[i];
		}

		double mu_n = mu(n_in, m, 1, C);
		double mu_p = mu(n_in, m, 2, C);
		double mu_e = mu_n - mu_p;
		double n_e = 0, n_mu = 0;
		if (mu_e > m_e){
			n_e += pow(mu_e*mu_e - m_e*m_e,1.5)/(3*M_PI*M_PI);
		}
		if (mu_e > m_mu){
			n_mu += pow(mu_e*mu_e - m_mu*m_mu,1.5)/(3*M_PI*M_PI);
		}

		hx[0] = sum_ch - n_e - n_mu;

		for (int i = 1; i < m; i++){
			hx[i] = calc::p_f(p[i]);

			double m_eff = C->M[i+1]*C->phi_n(C->X_s[i+1] * (C->M[0]/C->M[i+1]) * f);

			double res = pow(
					mu_n - C->Q[i+1]*mu_e - C->Co/pow(C->M[0],2) * C->X_o[i+1] * sum_o
					- C->Cr/pow(C->M[0],2) * C->X_r[i+1]*C->T[i+1] * sum_rho,
					2.0);

			res -= m_eff*m_eff;

			if (res > 0){
				hx[i] -= sqrt(res);
			}
		}
		delete[] n_in;
		delete[] n_f;
	}

}//namespace calc

/**\brief Energy density functional itself.
 * Provides the energy density functional evaluated at some point v = {n, f}
 */
double kineticInt(double n, double m, double f){
	double p_f = calc::p_f(n);
	double result2 = p_f*sqrt(m*m + p_f*p_f)*(m*m + 2*p_f*p_f);
	if (m > 0.0){
		result2 -= pow(m,4)*asinh(p_f/m);
	}
	result2 = result2/(8*M_PI*M_PI);
	return result2;
}

double _E(double * n, int dimN, set_const * C){
	double f = n[0];
	double res = pow(m_n, 4.0)*f*f*C->eta_s(f)/(2*C->Cs);
	double sum = 0;
	double sum_t3 = 0;
	for (int i = 1; i < dimN; ++i){
		res += kineticInt(n[i], (C->M)[i-1] * C->phi_n(C->X_s[i-1] * (C->M[0]/C->M[i-1]) * f), f);
//		printf("i = %i, n[i] = %f, pf(n[i]) = %f \n", i, v.n[i], calc::p_f(v.n[i]));
//		printf("K = %f \n", kineticInt(v.n[i], (C->M)[i] * C->phi_n(v.f), v.f));
//		printf("M_PI = %f \n", M_PI);
//		printf("asinh(1) = %f\n", asinh(1.0));
		sum += n[i] * C->X_o[i-1];
		sum_t3 += n[i]*(C->T)[i-1] * C->X_r[i-1];
	}
	res += C->Co * sum*sum/(2.0*m_n*m_n*C->eta_o(f));
//	printf("sum_t3 = %f  \n", sum_t3);
	res += C->Cr * pow(sum_t3/m_n, 2.0) / (2 * C->eta_r(f));
	res += C->U(f);
	return res;
}

double stepF(var v, set_const *C){
	return 0.0;
}
//Stepper function for E
void stepE(double n, double * init, int initN, double * out, int dim_Out, int iter, set_const* C) {
	double opts[5];
	bool debug = 1;
	calc::fun_n_eq_params p = {C, n, 0.1};
	int m = initN;
	double * x = new double[m];
	double * lb = new double[m];
	double * fun = new double[m];
	double info[LM_INFO_SZ];
	//double x[3] = {v.n[0], v.n[1], v.f};
	x[0] = init[0];
	lb[0] = 0.0;
	for (int i = 1; i < m; i++){
		x[i] = init[i];
		lb[i] = 0.0;
//		if (i > 2) lb[i] = -100500.0;
	}

	opts[0]= LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-25; opts[3]=1E-12;
		opts[4]= -1e-5;
	dlevmar_bc_dif(calc::fun_n_eq, x, NULL, m, m, lb, NULL, NULL, iter, opts, info, NULL, NULL, &p);
//	dlevmar_dif(calc::fun_n_eq, x, NULL, m, m, 2000, opts, NULL, NULL, NULL, &p);

	if (debug) {
		printf("info: ");
		for (int i = 0; i < LM_INFO_SZ; i++){
			printf(" %i : %f ", i, info[i]);
		}
		printf("\n");

		printf("n = %f, f_eq = %f", n, x[0]);
		for (int i = 1; i < m; i++){
			printf(",n%i = %f ", i, x[i]);
		}
		printf("\n");

		calc::fun_n_eq(x, fun, m, m, &p);
		for (int i = 0; i < m; i++){
			printf("f%i = %f  ", i, fun[i]);
		}
		printf("\n");
	}
	for (int i = 0; i < m; i++){
		out[i] = x[i];
	}
	delete[] x;
	delete[] fun;
	delete[] lb;
}

//TEST FUNCTION
double sum(std::vector<double> x){
	double res = 0;
	for (int i = 0; i < x.size(); i++){
		res += x[i];
	}
	return res;
}


float sumTest(double * in, int n){
	float res = 0;
	for (int i = 0; i < n; i++){
		res += in[i];
	}
	return res;
}

double E(double* n, int dimN, set_const* C) {
	double f = n[0];
	double res = pow(m_n, 4.0)*f*f*C->eta_s(f)/(2*C->Cs);
	double sum = 0;
	double sum_t3 = 0;
	for (int i = 1; i < dimN; ++i){
		res += kineticInt(n[i], (C->M)[i-1] * C->phi_n(C->X_s[i-1] * (C->M[0]/C->M[i-1]) * f), f);
//		printf("i = %i, n[i] = %f, pf(n[i]) = %f \n", i, v.n[i], calc::p_f(v.n[i]));
//		printf("K = %f \n", kineticInt(v.n[i], (C->M)[i] * C->phi_n(v.f), v.f));
//		printf("M_PI = %f \n", M_PI);
//		printf("asinh(1) = %f\n", asinh(1.0));
		sum += n[i] * C->X_o[i-1];
		sum_t3 += n[i]*(C->T)[i-1] * C->X_r[i-1];
	}
	res += C->Co * sum*sum/(2.0*m_n*m_n*C->eta_o(f));
//	printf("sum_t3 = %f  \n", sum_t3);
	res += C->Cr * pow(sum_t3/m_n, 2.0) / (2 * C->eta_r(f));
	res += C->U(f);
	double mu_n = calc::mu(n, dimN, 1, C);
	double mu_p = calc::mu(n, dimN, 2, C);
	double mu_e = mu_n - mu_p;
	double n_e = 0, n_mu = 0;
	if (mu_e > m_e){
		n_e += pow(mu_e*mu_e - m_e*m_e,1.5)/(3*M_PI*M_PI);
	}
	if (mu_e > m_mu){
		n_mu += pow(mu_e*mu_e - m_mu*m_mu,1.5)/(3*M_PI*M_PI);
	}
	res += kineticInt(n_e, m_e, f);
	res += kineticInt(n_mu, m_mu, f);
	return res;
}

float sumTest2(double * in, int n, double * in2, int n2){
	float res = 0;
	for (int i = 0; i < n; i++){
		res += in[i];
	}
	for (int i = 0; i < n2; i++){
			res += in2[i];
	}
	return res;
}


