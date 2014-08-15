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
		bool debug = false;
		double dn = 1e-3;
		n[i] += dn;
		double dE = _E(n, dimN, C);
		n[i] -= 2.0*dn;
		dE -= _E(n, dimN, C);
		n[i] += dn;
//		printf("mu: i = %i, res = %f \n", i, dE/dn);
		if (debug) {
			printf("mu: n[0] = %f, n[1] = %f, n[2] = %f", n[0], n[1], n[2]);
			for (int j = 3; j < dimN; j++){
				printf("n[%i] = %e",j, n[j]);
			}
			printf(" res=%f", dE/(2*dn));
			printf("\n");
		}
		return dE/(2.0*dn);
	}

	void fun_n_eq(double * p, double * hx, int m, int n, void * adata){
		bool debug = false;
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
		if (debug) {
			printf("n_f = ");
			for (int i = 0; i < m+1; i++){
				printf("%e ", n_f[i]);
			}
			printf("\n");
		}
		double f = f_eq(n_f, m+1, C, par->f_init);//m -> m+1 fixed 14.07.2014
		par->f_init = f;
		if (debug) {
			printf("f = %f \n", f);
		}

		n_in[0] = f;
		n_in[1] = n_n;

		double sum=0, sum_ch=0, sum_o = 0.0, sum_rho = 0.0, sum_p = 0;
		for (int i = 0; i < m + 1; i++){
			sum += n_f[i];
			sum_ch += n_f[i]*C->Q[i];
			sum_o += n_f[i]*C->X_o[i];
			sum_rho += n_f[i]*C->X_r[i]*C->T[i];
			if (C->phi_meson){
				sum_p += n_f[i]*C->X_p[i];
			}
//			printf("sum %f sum_ch %f sum_o %f sum_rho %f \n", sum, sum_ch, sum_o, sum_rho);
		}

		double mu_n = mu(n_in, m+2, 1, C);
		double mu_p = mu(n_in, m+2, 2, C);
		double mu_e = mu_n - mu_p;
//		printf("n = %f %f %f \n", n_in[0], n_in[1], n_in[2]);
//		printf("%f %f %f\n" ,mu_n, mu_p, sum_ch);
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
					mu_n - C->Q[i+1]*mu_e - C->Co/pow(C->M[0],2) * C->X_o[i+1] * sum_o / C->eta_o(f)
					- C->Cr/pow(C->M[0],2) * C->X_r[i+1]*C->T[i+1] * sum_rho / C->eta_r(f)
					- C->Co/pow(C->M[0], 2) * C->X_p[i+1] * sum_p * pow(m_o / m_p,2.0) / C->eta_p(f),
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
	bool debug = 0;
	if (debug){
		printf("n = ");
		for (int i = 0; i < dimN; i++){
			printf("%f ", n[i] );
		}
		printf("\n");
	}
	double f = n[0];
	double res = pow(m_n, 4.0)*f*f*C->eta_s(f)/(2*C->Cs);
	double sum = 0;
	double sum_t3 = 0;
	double sum_p = 0;
	if (debug){
		printf("res_f : %f \n", res);
	}
	for (int i = 1; i < dimN; ++i){
		res += kineticInt(n[i], (C->M)[i-1] * C->phi_n(C->X_s[i-1] * (C->M[0]/C->M[i-1]) * f), f);
//		printf("i = %i, n[i] = %f, pf(n[i]) = %f \n", i, v.n[i], calc::p_f(v.n[i]));
//		printf("K = %f \n", kineticInt(v.n[i], (C->M)[i] * C->phi_n(v.f), v.f));
//		printf("M_PI = %f \n", M_PI);
//		printf("asinh(1) = %f\n", asinh(1.0));
		sum += n[i] * C->X_o[i-1];
		sum_t3 += n[i]*(C->T)[i-1] * C->X_r[i-1];
		if (C->phi_meson){
			sum_p += n[i]*C->X_p[i-1];
		}
	}
	//omega
	res += C->Co * sum*sum/(2.0*m_n*m_n*C->eta_o(f));
	if (debug){
		printf("res_om : %f \n", res);
	}
	//phi
	res += pow(m_o/m_p ,2.0)*C->Co * sum_p*sum_p/(2.0*m_n*m_n*C->eta_p(f));
	if (debug){
		printf("res_phi : %f \n", res);
	}
	//rho
	if (debug){
		printf("res_rho : %f \n", res);
	}
	res += C->Cr * pow(sum_t3/m_n, 2.0) / (2 * C->eta_r(f));

	res += C->U(f);
	return res;
}

double E(double* n, int dimN, set_const* C) {
	double f = n[0];
	double res = _E(n, dimN, C);

//	double f = n[0];
//	double res = pow(m_n, 4.0)*f*f*C->eta_s(f)/(2*C->Cs);
//	double sum = 0;
//	double sum_t3 = 0;
//	for (int i = 1; i < dimN; ++i){
//		res += kineticInt(n[i], (C->M)[i-1] * C->phi_n(C->X_s[i-1] * (C->M[0]/C->M[i-1]) * f), f);
////		printf("i = %i, n[i] = %f, pf(n[i]) = %f \n", i, v.n[i], calc::p_f(v.n[i]));
////		printf("K = %f \n", kineticInt(v.n[i], (C->M)[i] * C->phi_n(v.f), v.f));
////		printf("M_PI = %f \n", M_PI);
////		printf("asinh(1) = %f\n", asinh(1.0));
//		sum += n[i] * C->X_o[i-1];
//		sum_t3 += n[i]*(C->T)[i-1] * C->X_r[i-1];
//	}
//	res += C->Co * sum*sum/(2.0*m_n*m_n*C->eta_o(f));
////	printf("sum_t3 = %f  \n", sum_t3);
//	res += C->Cr * pow(sum_t3/m_n, 2.0) / (2 * C->eta_r(f));
//	res += C->U(f);

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

double stepF(var v, set_const *C){
	return 0.0;
}
//Stepper function for E
void stepE(double n, double * init, int initN, double f_init, double * out, int dim_Out, int iter, set_const* C) {
	double opts[5];
	bool debug = 0;
	calc::fun_n_eq_params p = {C, n, f_init};
	int m = initN;
	double * x = new double[m];
	double * lb = new double[m];
	double * fun = new double[m];
	double info[LM_INFO_SZ];
	//double x[3] = {v.n[0], v.n[1], v.f};
	x[0] = init[0];
	if (C->sprime){
		x[1] = init[1];
	}
	lb[0] = 0.0;
	for (int i = 1+C->sprime; i < m; i++){
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

		printf("n = %f, n_p = %e", n, x[0]);
		printf(",n_L = %e ", x[1]);
		printf(",n_S- = %e ", x[2]);
		printf(",n_S0 = %e ", x[3]);
		printf(",n_S+ = %e ", x[4]);
		printf(",n_X- = %e ", x[5]);
		printf(",n_X0 = %e ", x[6]);
		printf("\n");

		calc::fun_n_eq(x, fun, m, m, &p);
		for (int i = 0; i < m; i++){
			printf("f%i = %e  ", i, fun[i]);
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
