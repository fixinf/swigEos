/*
 * eos.cpp
 *
 *  Created on: 09 июня 2014 г.
 *      Author: const
 */


#include "eos.h"

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
	};

	double mu(double * n, int dimN, int i, set_const * C){
		double dn = 1e-4;
		n[i] += dn;
		double dE = E(n, dimN, C);
		n[i] -= dn;
		dE -= E(n, dimN, C);
//		printf("mu: i = %i, res = %f \n", i, dE/dn);
		return dE/dn;
	}

	void fun_n_eq(double * p, double * hx, int m, int n, void * adata){
		double f_eq = p[0];
//		printf("i'm in\n");
		fun_n_eq_params * par = (fun_n_eq_params *) adata;
		set_const * C = par->C;
//		printf("i'm in\n");
//		vec _n;
		double sum=0, sum_ch=0;
//		printf("n = ");
//		printf("p = ");
		for (int i = 1; i < m; i++){
//			printf("%f ", p[i]);
			//_n.push_back(p[i]);
			sum += p[i];
			sum_ch += p[i]*C->Q[i-1];
		}
//		printf("\n");
//		printf("i'm in\n");
		double df = 1e-4;
//		var v;
//		v.n = _n;
//		v.f = f_eq + df;
		double dE = E(p, m, C);
		p[0] -= df;
		dE -= E(p, m ,C);
		dE /= df;
		p[0] += df; //returns to the given value
		hx[0] = dE;
		double mu_e = mu(p, m, 1, C) - mu(p, m, 2, C);
		double n_e = 0, n_mu = 0;
		if (mu_e > m_e){
			n_e += pow(mu_e*mu_e - m_e*m_e,1.5)/(3*M_PI*M_PI);
		}
		if (mu_e > m_mu){
			n_mu += pow(mu_e*mu_e - m_mu*m_mu,1.5)/(3*M_PI*M_PI);
		}
		hx[1] = par->n - sum;
		hx[2] = sum_ch - n_e - n_mu;
		//var() for testing the existence of hyperons
		double * p_temp = new double[m];
		p_temp[0] = f_eq;
		for (int i = 1; i < m; i++){
			p_temp[i] = p[i];
		}
		for (int i = 3; i < m; i++){
//			if (v.n[i] < 1e-6){
//				hx[i] = 0.0;
//			}
//			else{
//			printf("i = %i, %i\n",i,v_temp.n.size());

//			printf("mu = %f, m = %f \n", mu(v, 0, C), C->M[i-1]);
			p_temp[i] = 0;
			if (mu(p_temp, m, 1, C) - C->Q[i-1]*mu_e >= mu(p_temp, m, i, C)){//C->M[i-1]*C->phi_n(C->X_s[i] * (C->M[0]/C->M[i]) * v.f)){
				hx[i] = mu(p, m, 1, C) - C->Q[i-1]*mu_e - mu(p, m, i, C);
//				printf("! eq = %f \n control = %f \n", hx[i], mu(p_temp, m, 1, C) -
//						C->Q[i-1]*mu_e - mu(p_temp, m, i, C));
			}
			else{
				hx[i] = 0;
				sum -= p[i];
				hx[1] = par->n - sum;
				sum_ch -= p[i]*C->Q[i-1];
				hx[2] = sum_ch - n_e - n_mu;
				p[i] = 0;
			}
//			}
			p_temp[i] = p[i];
		}
//		printf("f = ");
		for (int i = 0; i < m; i++){
//			printf("%f ", hx[i]);
		}
//		printf("\n");
	}

	double p_f(double n) {
		return pow(3.0 * M_PI * M_PI * D * n, 1.0 / 3.0);
	}
}

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

double E(double * n, int dimN, set_const * C){
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
	calc::fun_n_eq_params p = {C, n};
	int m = initN;
	double * x = new double[m];
	double * lb = new double[m];
	double * fun = new double[m];
	double info[LM_INFO_SZ];
	//double x[3] = {v.n[0], v.n[1], v.f};
	x[0] = init[0];
	for (int i = 1; i < m; i++){
		x[i] = init[i];
		lb[i] = 0.0;
//		if (i > 2) lb[i] = -100500.0;
	}

	opts[0]= LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-12;
		opts[4]= -1e-3;
	dlevmar_bc_dif(calc::fun_n_eq, x, NULL, m, m, lb, NULL, NULL, iter, opts, info, NULL, NULL, &p);
//	dlevmar_dif(calc::fun_n_eq, x, NULL, m, m, 2000, opts, NULL, NULL, NULL, &p);

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

	for (int i = 0; i < m; i++){
		out[i] = x[i];
	}
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



