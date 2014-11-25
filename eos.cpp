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
double p_f(double n) {
	return pow(3.0 * M_PI * M_PI * D * n, 1.0 / 3.0);
}

double kineticInt(double n, double m, double f){
	double pf = p_f(n);
	double result2 = pf*sqrt(m*m + pf*pf)*(m*m + 2*pf*pf);
	if (m > 0.0){
		result2 -= pow(m,4)*asinh(pf/m);
	}
	result2 = result2/(8*M_PI*M_PI);
	return result2;
}


void potentials(double * n, int dimN, double * result, int dimResult, set_const * C){
	double sum_o = 0;
	double sum_r = 0;
	double sum_phi = 0;
	int sc = 1 + C->sprime;
	for (int i = sc; i < dimN; i++){
		sum_o += n[i] * C->X_o[i-sc];
		sum_r += n[i]*(C->T)[i-sc] * C->X_r[i-sc];
		if (C->phi_meson){
			sum_phi += n[i]*C->X_p[i-sc];
		}
	}
	double f = n[0];
	double fp = 0;
	if (C->sprime){
		fp = n[1];
	}
	result[0] = -f*C->M[0]; ///(sqrt(C->Cs/C->eta_s(f)));
	result[1] = -fp*C->M[0]/(sqrt(C->Csp));
	result[2] = (C->Co/(C->eta_o(f)))*sum_o/pow(C->M[0],2);
	result[3] = (C->Cr/(C->eta_r(f)))*sum_r/pow(C->M[0],2);
	result[4] = (C->Co/C->eta_p(f)) * pow(m_o/m_p, 2) * sum_phi / pow(C->M[0], 2);
}

namespace calc{
	struct fun_n_eq_params{
		set_const * C;
		double n;
		double * f_init;
		int dimF_init;
	};



	double mu(double * n, int dimN, int i, set_const * C){
//		bool debug = false;
//		double dn = 1e-3;
//
//		n[i] += dn;
//		double dE = _E(n, dimN, C);
//		n[i] -= 2.0*dn;
//		dE -= _E(n, dimN, C);
//		n[i] += dn;
//
////		Increased precision: f'(x) = (1/3h)(2(f(1) - f(-1)) - 0.25 (f(2) - f(-2)) )
////
////		n[i] += 2*dn;
////		double dE = -0.25*_E(n, dimN, C);
////		n[i] -= dn;
////		dE += 2*_E(n, dimN, C);
////		n[i] -= 2*dn;
////		dE += -2*_E(n, dimN, C);
////		n[i] -= dn;
////		dE += 0.25*_E(n, dimN, C);
////		n[i] += 2*dn;
//
//		if (debug) {
//			printf("mu: n[0] = %f, n[1] = %f, n[2] = %f", n[0], n[1], n[2]);
//			for (int j = 3; j < dimN; j++){
//				printf("n[%i] = %e",j, n[j]);
//			}
//			printf(" res=%f", dE/(2*dn));
//			printf("\n");
//		}
//		return dE/(2.0*dn);
//		return dE/(3.0*dn);

		int sp = 1 + C->sprime;
		i = i - sp;
		double out[5];
		potentials(n, dimN, out, 5, C);
		double f = n[0];
		double fp = 0.0;
		if (C->sprime){
			fp = n[1];
		}
		double xs = 0.;
		if (C->sigma_kind == 0){
			xs = C->X_s[i];
		}
		else{
			xs = C->Xs(i, f);
		}
		double m_eff_arg = xs*(C->M[0]/C->M[i])*f + C->X_sp[i]*(C->M[0]/C->M[i])*fp;
//		printf("f = %f, m_eff_arg = %f \n", f, m_eff_arg);
		double m_eff = C->M[i] * C->phi_n(i,m_eff_arg);
		double res = sqrt(pow(p_f(n[i+sp]), 2.0) + pow(m_eff, 2.0));
//		double res = 0.0;
		res += C->X_o[i]*out[2];
		res += C->X_r[i]*C->T[i]*out[3];
		res += C->X_p[i]*out[4];
//		printf("mu[%i] res = %f \n", i, res);
		return res;
	}


	void fun_n_eq(double * p, double * hx, int m, int n, void * adata){
		bool debug = 0;
		fun_n_eq_params * par = (fun_n_eq_params *) adata;
		set_const * C = par->C;
		int sc = 1 + C->sprime;
		double n_sum = 0.0;
		double n_n = par->n;
		double * n_in = new double [m+sc+1]; //input set for _E and mu
		double * n_f = new double [m+sc]; //input set for f_eq; actually is {n_n,n_p,...,n_X0}
		for (int i = 0; i < m; i++){
			n_n -= p[i];
			n_in[i + 1 + sc] = p[i]; //scalar + neutron(1) offset
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
		double * out = new double[sc];
		f_eq(n_f, m+1, par->f_init, sc, out, sc, C);//m -> m+1 fixed


		for (int i = 0; i < sc; i++){
			n_in[i] = out[i];
		}
		n_in[sc] = n_n;
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

		double mu_n = mu(n_in, m + sc + 1, sc + 0, C);
		double mu_p = mu(n_in, m + sc + 1, sc + 1, C);
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

		double fp = 0;
		double f = out[0];
		if (C->sprime){
			fp = out[1];
		}

		for (int i = 1; i < m; i++){
			hx[i] = p_f(p[i]);
			double xs = 0.;
			if (C->sigma_kind == 0){
				xs = C->X_s[i+1];
			}
			else{
				xs = C->Xs(i+1, f);
			}
			double m_eff = C->M[i+1]*C->phi_n(i+1, xs * (C->M[0]/C->M[i+1]) * f  + C->X_sp[i+1] * (C->M[0]/C->M[i+1]) * fp);

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
	int sc = 1 + C->sprime;
	double fp = 0;
	if (C->sprime){
		fp = n[1];
	}
	double res = pow(C->M[0], 4.0)*f*f*C->eta_s(f)/(2*C->Cs);
	if (debug){
		printf("res_f : %f\n", res);
	}

	res += C->U(f);

	res += pow(C->M[0], 4.0)*fp*fp/(2*C->Csp);

	double sum = 0;
	double sum_t3 = 0;
	double sum_p = 0;
	double meff_arg = 0;
	if (debug){
		printf("res_Uf : %f \n", res);
	}
	for (int i = sc; i < dimN; ++i){
		double xs = 0.;
		if (C->sigma_kind == 0){
			xs = C->X_s[i-sc];
		}
		else{
			xs = C->Xs(i-sc, f);
		}
//		printf("xs = %f \n", xs);
		meff_arg = xs * (C->M[0]/C->M[i-sc]) * f + C->X_sp[i-sc] * (C->M[0]/C->M[i-sc])*fp;
		res += kineticInt(n[i], (C->M)[i-sc] * C->phi_n(i-sc,meff_arg), f);
//		printf("i = %i, n[i] = %f, pf(n[i]) = %f \n", i, v.n[i], calc::p_f(v.n[i]));
//		printf("K = %f \n", kineticInt(v.n[i], (C->M)[i] * C->phi_n(v.f), v.f));
//		printf("M_PI = %f \n", M_PI);
//		printf("asinh(1) = %f\n", asinh(1.0));
		sum += n[i] * C->X_o[i-sc];
		sum_t3 += n[i]*(C->T)[i-sc] * C->X_r[i-sc];
		if (C->phi_meson){
			sum_p += n[i]*C->X_p[i-sc];
		}
	}
	//omega
	res += C->Co * sum*sum/(2.0*C->M[0]*C->M[0]*C->eta_o(f));
	if (debug){
		printf("res_om : %f \n", res);
	}
	//phi
	res += pow(m_o/m_p ,2.0)*C->Co * sum_p*sum_p/(2.0*C->M[0]*C->M[0]*C->eta_p(f));
	if (debug){
		printf("res_phi : %f \n", res);
	}
	//rho

	res += C->Cr * pow(sum_t3/C->M[0], 2.0) / (2 * C->eta_r(f));
	if (debug){
		printf("res_rho : %f \n", res);
	}

	return res;
}

double E(double* n, int dimN, set_const* C) {
	double f = n[0];
	double res = _E(n, dimN, C);

//	double f = n[0];
//	double res = pow(C->M[0], 4.0)*f*f*C->eta_s(f)/(2*C->Cs);
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
//	res += C->Co * sum*sum/(2.0*C->M[0]*C->M[0]*C->eta_o(f));
////	printf("sum_t3 = %f  \n", sum_t3);
//	res += C->Cr * pow(sum_t3/C->M[0], 2.0) / (2 * C->eta_r(f));
//	res += C->U(f);
	int sp = 1 + C->sprime;
	double mu_n = calc::mu(n, dimN, sp, C);
	double mu_p = calc::mu(n, dimN, sp+1, C);
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
void stepE(double n, double * init, int initN, double * f_init, int dimF_init, double * out, int dim_Out, int iter, set_const* C) {
	double opts[5];
	bool debug = 0;
	calc::fun_n_eq_params p = {C, n, f_init};
	int m = initN;
	double * x = new double[m];
	double * lb = new double[m];
	double * fun = new double[m];
	double info[LM_INFO_SZ];
	//double x[3] = {v.n[0], v.n[1], v.f};

	for (int i = 0; i < m; i++){
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

struct fun_eq2_params{
	double E;
	double P;
	double n;
	double Co;
	double Cr;
	double mn;
};

void fun_n_eq(double * p, double * hx, int m, int n, void * adata){
	fun_eq2_params * params = (fun_eq2_params *) adata;
	double f = p[0];
	double np = p[1];
	double nn = params->n - np;
	double mn = params->mn;
	printf("f = %f, np = %f \n", f, np);
	hx[0] = (params->E + params->P -
			 params->Co * pow(params->n/mn, 2) -
			 params->Cr * pow((np - nn)/mn, 2)/4);
	hx[0] -= nn*pow(pow(p_f(nn),2) + pow(mn*(1-f),2), 0.5) +
			 np*pow(pow(p_f(np),2) + pow(mn*(1-f),2), 0.5);
	double ne = 0.;
	double nmu = 0.;
	double mu_e = pow(pow(p_f(nn),2) + pow(mn*(1-f),2), 0.5) -
				  pow(pow(p_f(np),2) + pow(mn*(1-f),2), 0.5) +
				  params->Cr * (nn - np)/(2*mn*mn);
	printf("mu_e = %f \n", mu_e);
	if (mu_e > m_e){
		ne += pow(mu_e*mu_e - m_e*m_e, 1.5)/(3*M_PI*M_PI);
	}

	if (mu_e > m_mu){
		nmu += pow(mu_e*mu_e - m_mu*m_mu, 1.5)/(3*M_PI*M_PI);
	}
	hx[0] -= ne*mu_e + nmu*mu_e;
	hx[1] = ne + nmu - np;
	printf("out[0] = %f, out[1] = %f \n", hx[0], hx[1]);
}

void solveF(double n, double E, double P, double * init, int initN, double * out, int dim_Out, set_const * C){
	fun_eq2_params par = {E, P, n, C->Co, C->Cr, C->M[0]};
	int m = 2;
	double opts[5];
	double * x = new double[m];
	double * lb = new double[m];
	double * fun = new double[m];
	double info[LM_INFO_SZ];
	int iter = 1000;
	//double x[3] = {v.n[0], v.n[1], v.f};

	for (int i = 0; i < m; i++){
		x[i] = init[i];
		lb[i] = 0.0;
//		if (i > 2) lb[i] = -100500.0;
	}

	opts[0]= LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-25; opts[3]=1E-20;
		opts[4]= -1e-5;

	printf("Invoking solver \n");
	dlevmar_bc_dif(fun_n_eq, x, NULL, m, m, lb, NULL, NULL, iter, opts, info, NULL, NULL, &par);

	printf("info: ");
	for (int i = 0; i < LM_INFO_SZ; i++){
		printf(" %i : %f ", i, info[i]);
	}
	printf("\n");

	printf("n = %f, f = %e", n, x[0]);
	printf(",n_p = %e ", x[1]);
	printf("\n");

	fun_n_eq(x, fun, m, m, &par);
	for (int i = 0; i < m; i++){
		printf("f%i = %e  ", i, fun[i]);
	}
	printf("\n");
	out[0] = x[0];
	out[1] = x[1];
}

