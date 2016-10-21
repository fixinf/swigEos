/*
 * eos.cpp
 *
 *  Created on: 09 �������� 2014 ��.
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
double p_f(double n, double gamma) {
	if (n < 0)
		return 0.;
	return pow(6. * M_PI * M_PI * D * n / gamma, 1.0 / 3.0);
}

double kineticInt(double n, double m, double gamma){
	double pf = p_f(n, gamma);
	double result2 = pf*sqrt(m*m + pf*pf)*(m*m + 2*pf*pf);
	if (m > 0.0){
		result2 -= pow(m,4)*asinh(pf/m);
	}
	result2 = gamma * result2/(16*M_PI*M_PI);
	return result2;
}


void potentials(double * n, int dimN, double * result, int dimResult, set_const * C){
	double sum_o = 0;
	double sum_r = 0;
	double sum_phi = 0;
	int sc = 1 + C->sprime;
	for (int i = sc; i < dimN; i++){
		sum_o += n[i] * C->X_o[i-sc];
		sum_r += n[i]* C->T[i-sc] * C->X_r[i-sc];
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
	(*C).X_s == C->X_s;
}

namespace calc{
	struct fun_n_eq_params{
		set_const * C;
		double n;
		double * f_init;
		int dimF_init;
		double misc;
	};

	struct fun_n_eq_f_params{
			set_const * C;
			double f;
			double misc;
		};

    double mu_deriv(double *n, int dimN, int i, double mu_c, set_const *C){
	    bool debug = false;
		double dn = 1e-4;

//		n[i] += dn;
//		double dE = _E(n, dimN, C);
//		n[i] -= 2.0*dn;
//		dE -= _E(n, dimN, C);
//		n[i] += dn;

		//Increased precision: f'(x) = (1/3h)(2(f(1) - f(-1)) - 0.25 (f(2) - f(-2)) )

		n[i] += 2*dn;
		double dE = -0.25*E_rho(n, dimN, mu_c, C);
		n[i] -= dn;
		dE += 2*E_rho(n, dimN, mu_c, C);
		n[i] -= 2*dn;
		dE += -2*E_rho(n, dimN, mu_c, C);
		n[i] -= dn;
		dE += 0.25*E_rho(n, dimN, mu_c, C);
		n[i] += 2*dn;

		if (debug) {
			printf("mu: n[0] = %f, n[1] = %f, n[2] = %f", n[0], n[1], n[2]);
			for (int j = 3; j < dimN; j++){
				printf("n[%i] = %e",j, n[j]);
			}
			printf(" res=%f", dE/(2*dn));
			printf("\n");
		}
//		return dE/(2.0*dn);
		return dE/(3.0*dn);

	}

    double mu_rho(double * n,  int dimN, int i, double mu_c, set_const * C){
		int sp = 1 + C->sprime;
		i = i - sp;
		double out[5];
		potentials(n, dimN, out, 5, C);

		double sum_o = 0.;
		double sum_r = 0.;
		double sum_phi = 0.;

		for (int j = sp; j < dimN; j++){
//			printf("n[%i]=%.6f, Xo[%i]= %.6f", j, n[j], j, C->X_o[j-sp]);
			sum_o += n[j] * C->X_o[j-sp];
			sum_r += n[j]*(C->T)[j-sp] * C->X_r[j-sp];
			if (C->phi_meson){
				sum_phi += n[j]*C->X_p[j-sp];
			}
		}
//		printf("\n sum_o=%.6f, sum_r=%.6f, sum_p=%.6f \n" ,sum_o, sum_r,
//				sum_phi);
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
		double res = sqrt(pow(p_f(n[i + sp], 2*C->S[i] + 1), 2.0) + pow(m_eff, 2.0));
//		double res = 0.0;
		res += C->X_o[i]*out[2];
		res += C->X_r[i]*C->T[i]*out[3];
		res += C->X_p[i]*out[4];
//		printf("mu[%i] res = %f \n", i, res);

		///////////////Rho-condensate contribution//////////////
    	double n_rho = 2 * C->m_rho * pow(C->M[0]*C->phi_n(0, f),2.)* sqrt(C->eta_r(f)) / (C->Cr * C->chi_prime
      (f)) *
    		(1 - mu_c/(C->m_rho * C->phi_n(0, f)));

    	if (fabs(sum_r) > n_rho/2){
    		res -= C->Cr / (pow(C->M[0],2.) * C->eta_r(f)) * (fabs(sum_r) - n_rho/2) * C->X_r[i] * C->T[i] *
    				((sum_r > 0) - (sum_r < 0)); //sign(sum_r) = - sign (n_n - n_p)
    	}

		return res;

    }

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
		double res = sqrt(pow(p_f(n[i + sp], 2*C->S[i] + 1), 2.0) + pow(m_eff, 2.0));
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
			hx[i] = p_f(p[i], 2*C->S[i+1]+1);
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



	void fun_n_eq_f(double * p, double * hx, int m, int n, void * adata){
			bool debug = 0;
			fun_n_eq_f_params * par = (fun_n_eq_f_params *) adata;
			set_const * C = par->C;
			int sc = 1 + C->sprime;
			double * n_in = new double [m+1]; //input set for _E and mu
			n_in[0] = par->f;
			for (int i = 0; i < m; i++){
				n_in[i + 1] = p[i];
			}
//
			double f = par->f;
			double out[1];
		    func_f_eq_params params = {p, m, 1e-5, C};
		    double f_in[1] = {par->f};
		    func_f_eq(f_in, out, 1, 1, &params);


			double sum=0, sum_ch=0, sum_o = 0.0, sum_rho = 0.0, sum_p = 0;
			for (int i = 0; i < m; i++){
				sum += p[i];
				sum_ch += p[i]*C->Q[i];
				sum_o += p[i]*C->X_o[i];
				sum_rho += p[i]*C->X_r[i]*C->T[i];
				if (C->phi_meson){
					sum_p += p[i]*C->X_p[i];
				}
				if (debug)
				printf("sum %f sum_ch %f sum_o %f sum_rho %f \n", sum, sum_ch, sum_o, sum_rho);
			}

			double mu_n = mu(n_in, m + sc + 1, sc + 0, C);
			double mu_p = mu(n_in, m + sc + 1, sc + 1, C);
			double mu_e = mu_n - mu_p;
			if (debug){
			printf("n = %f %f %f \n", n_in[0], n_in[1], n_in[2]);
			printf("%f %f %f\n" ,mu_n, mu_p, sum_ch);
			}
			double n_e = 0, n_mu = 0;
			if (mu_e > m_e){
				n_e += pow(mu_e*mu_e - m_e*m_e,1.5)/(3*M_PI*M_PI);
			}
			if (mu_e > m_mu){
				n_mu += pow(mu_e*mu_e - m_mu*m_mu,1.5)/(3*M_PI*M_PI);
			}

			hx[0] = sum_ch - n_e - n_mu;

			hx[1] = out[0];

			for (int i = 2; i < m; i++){
				hx[i] = p_f(p[i], 2*C->S[i]+1);
				double xs = 0.;
				if (C->sigma_kind == 0){
					xs = C->X_s[i];
				}
				else{
					xs = C->Xs(i, f);
				}
				double m_eff = C->M[i]*C->phi_n(i, xs * (C->M[0]/C->M[i]) * f);

				double res = pow(
						mu_n - C->Q[i]*mu_e - C->Co/pow(C->M[0],2) * C->X_o[i] * sum_o / C->eta_o(f)
						- C->Cr/pow(C->M[0],2) * C->X_r[i]*C->T[i+1] * sum_rho / C->eta_r(f)
						- C->Co/pow(C->M[0], 2) * C->X_p[i] * sum_p * pow(m_o / m_p,2.0) / C->eta_p(f),
						2.0);

				res -= m_eff*m_eff;

				if (res > 0){
					hx[i] -= sqrt(res);
				}
			}
			if (debug){
				for (int i = 0; i < m; i++) printf("hx[%i] = %.6f ",i, hx[i]);
				printf("\n");
			}
			delete[] n_in;
		}

	void fun_n_eq_f_rho(double * p, double * hx, int m, int n, void * adata){
		bool debug = 0;
		fun_n_eq_f_params * par = (fun_n_eq_f_params *) adata;
		set_const * C = par->C;
		int sc = 1 + C->sprime;
		double * n_in = new double [m+1-1]; //input set for _E and mu
		n_in[0] = par->f;
		for (int i = 0; i < m-1; i++){
			n_in[i + 1] = p[i];
		}

		double mu_c = p[m-1];
//
		double f = par->f;
		double out[1];
		func_f_eq_params_rho params = {p, m-1, 1e-5, C, mu_c};
		double f_in[1] = {par->f};
		func_f_eq_rho(f_in, out, 1, 1, &params);


		double sum=0, sum_ch=0, sum_o = 0.0, sum_rho = 0.0, sum_p = 0;
		for (int i = 0; i < m-1; i++){
			sum += p[i];
			sum_ch += p[i]*C->Q[i];
			sum_o += p[i]*C->X_o[i];
			sum_rho += p[i]*C->X_r[i]*C->T[i];
			if (C->phi_meson){
				sum_p += p[i]*C->X_p[i];
			}
		}

		if (debug)
			printf("sum %f sum_ch %f sum_o %f sum_rho %f \n", sum, sum_ch, sum_o, sum_rho);

		double mu_n = mu_rho(n_in, m + sc, sc + 0, mu_c, C);
		double mu_p = mu_rho(n_in, m + sc, sc + 1, mu_c, C);
		double mu_e = mu_n - mu_p;
		if (debug){
			for (int i = 0; i < m + sc; i++){
				printf("n_in[%i] = %.6f ", i, n_in[i]);
			}
		printf("\n");
//		printf("n = %f %f %f , mu_c = %f \n", n_in[0], n_in[1], n_in[2], mu_c);
		printf("%f %f %f\n" ,mu_n, mu_p, sum_ch);
		}
		double n_e = 0, n_mu = 0;
		if (mu_e > m_e){
			n_e += pow(mu_e*mu_e - m_e*m_e,1.5)/(3*M_PI*M_PI);
		}
		if (mu_e > m_mu){
			n_mu += pow(mu_e*mu_e - m_mu*m_mu,1.5)/(3*M_PI*M_PI);
		}

		hx[0] = sum_ch - n_e - n_mu;

		hx[1] = out[0];


    	double n_rho = 2 * C->m_rho * pow(C->M[0]*C->phi_n(0, f),2.)* sqrt(C->eta_r(f)) / (C->Cr * C->chi_prime
      (f)) *
    		(1 - mu_c/(C->m_rho * C->phi_n(0, f)));
    	par->misc = 0.;

    	if (debug)
    	printf("n_rho = %.6f, sum_rho = %.6f \n", n_rho, sum_rho);

		double  n_c = 0.;

		if(fabs(sum_rho) > n_rho/2){
			double r_c2 = (fabs(sum_rho) - n_rho/2) / (2 * C->m_rho * C->chi_prime(f) * sqrt(C->eta_r(f)));
			n_c = 2 * C->m_rho * C->phi_n(0, f) * r_c2;

			hx[0] -= n_c;

			par->misc = n_c;

//			printf("sum_ch = %.6f, n_e = %.6f, n_mu = %.6f, n_c = %.6f \n", sum_ch, n_e, n_mu, n_c);
		}

		for (int i = 2; i < m-1; i++){
			hx[i] = p_f(p[i], 2*C->S[i]+1);
			double xs = 0.;
			if (C->sigma_kind == 0){
				xs = C->X_s[i];
			}
			else{
				xs = C->Xs(i, f);
			}
			double m_eff = C->M[i]*C->phi_n(i+1, xs * (C->M[0]/C->M[i]) * f);

			double res = mu_n - C->Q[i]*mu_e - C->Co/pow(C->M[0],2) * C->X_o[i] * sum_o / C->eta_o(f)
					- C->Cr/pow(C->M[0],2) * C->X_r[i]*C->T[i] * sum_rho / C->eta_r(f)
					- C->Co/pow(C->M[0], 2) * C->X_p[i] * sum_p * pow(m_o / m_p,2.0) / C->eta_p(f);

	    	if (fabs(sum_rho) > n_rho/2){
	    		res += C->Cr / (pow(C->M[0],2.) * C->eta_r(f)) * (fabs(sum_rho) - n_rho/2) * C->X_r[i] * C->T[i] *
	    				((sum_rho > 0) - (sum_rho < 0));
	    	}
	    	res = pow(res, 2.0);

			res -= m_eff*m_eff;

			if (res > 0){
				hx[i] -= sqrt(res);
			}
		}


		hx[m-1] = mu_c - mu_e;

		if (debug){
			for (int i = 0; i < m; i++) printf("hx[%i] = %.6f ",i, hx[i]);
			printf("\n");
		}


		delete[] n_in;
	}
  
  void fun_n_eq_dsym(double * p, double * hx, int m, int k, void * adata){
      double n_d = p[0];
      fun_n_eq_params * par = (fun_n_eq_params *) adata;
      set_const * C = par->C;
      double n = par->n;
      double n_n = (n - n_d)/2; //Nucleon concentration

      double n_f[6] = {n_n, n_n, n_d/4, n_d/4, n_d/4, n_d/4};
      double  * out = new double[1];
      f_eq(n_f, 6, par->f_init, 1, out, 1, C);
      double f = out[0];

      double n_in[7];
      n_in[0] = f;
      for (int i = 0; i < 6; i++){
        n_in[i+1] = n_f[i];
      }
      int i = 2;
      double sum_o = C->X_o[0] * 2 * n_n + C->X_o[2] * n_d;
      double mu_n = mu(n_in, 7, 1, C);
      hx[0] = p_f(n_d/4, 2*C->S[2]+1);
			double xs = 0.;
			if (C->sigma_kind == 0){
				xs = C->X_s[3];
			}
			else{
				xs = C->Xs(3, f);
			}
			double m_eff = C->M[i+1]*C->phi_n(i+1, xs * (C->M[0]/C->M[i+1]) * f);

			double res = mu_n - C->Co/pow(C->M[0],2) * C->X_o[i+1] * sum_o / C->eta_o(f);
					

	    res = pow(res, 2.0);

			res -= m_eff*m_eff;

			if (res > 0){
				hx[0] -= sqrt(res);
			}
	}

  void fun_n_eq_dsym_f(double * p, double * hx, int m, int k, void * adata){
        double n_d = p[0];
        fun_n_eq_params * par = (fun_n_eq_params *) adata;
        set_const * C = par->C;
        double n = par->n;
        double n_n = (n - n_d)/2; //Nucleon concentration

        double n_f[6] = {n_n, n_n, n_d/4, n_d/4, n_d/4, n_d/4};
        double  * out = new double[1];
//        f_eq(n_f, 6, par->f_init, 1, out, 1, C);
        double f = par->f_init[0];

        double n_in[7];
        n_in[0] = f;
        for (int i = 0; i < 6; i++){
          n_in[i+1] = n_f[i];
        }
        int i = 2;
        double sum_o = C->X_o[0] * 2 * n_n + C->X_o[2] * n_d;
        double mu_n = mu(n_in, 7, 1, C);
        hx[0] = p_f(n_d/4, 2*C->S[2]+1);
  			double xs = 0.;
  			if (C->sigma_kind == 0){
  				xs = C->X_s[3];
  			}
  			else{
  				xs = C->Xs(3, f);
  			}
  			double m_eff = C->M[i+1]*C->phi_n(i+1, xs * (C->M[0]/C->M[i+1]) * f);

  			double res = mu_n - C->Co/pow(C->M[0],2) * C->X_o[i+1] * sum_o / C->eta_o(f);


  	    res = pow(res, 2.0);

  			res -= m_eff*m_eff;

  			if (res > 0){
  				hx[0] -= sqrt(res);
  			}
  	}

	void fun_n_eq_rho(double * p, double * hx, int m, int n, void * adata){
		bool debug = 0;
		fun_n_eq_params * par = (fun_n_eq_params *) adata;
		set_const * C = par->C;
		int sc = 1 + C->sprime;
		double n_sum = 0.0;
		double n_n = par->n;
		double * n_in = new double [m+sc+1]; //input set for _E and mu
		double * n_f = new double [m+sc]; //input set for f_eq; actually is {n_n,n_p,...,n_X0}
		for (int i = 0; i < m-1; i++){
			n_n -= p[i];
			n_in[i + 1 + sc] = p[i]; //scalar + neutron(1) offset
			n_f[i+1] = p[i];
		}

		double mu_c = p[m-1];

		n_f[0] = n_n;
		if (debug) {
			printf("n_f = ");
			for (int i = 0; i < m+1; i++){
				printf("%e ", n_f[i]);
			}
			printf("\n");
		}
		double * out = new double[sc];
		f_eq_rho(n_f, m, par->f_init, sc, out, sc, mu_c, C);//m -> m+1 fixed


		for (int i = 0; i < sc; i++){
			n_in[i] = out[i];
			printf("f[%i] = %.3f ", i, out[i]);
		}
		printf("\n");

		n_in[sc] = n_n;
		double sum=0, sum_ch=0, sum_o = 0.0, sum_rho = 0.0, sum_p = 0;
		for (int i = 0; i < m; i++){
			sum += n_f[i];
			sum_ch += n_f[i]*C->Q[i];
			sum_o += n_f[i]*C->X_o[i];
			sum_rho += n_f[i]*C->X_r[i]*C->T[i];
			if (C->phi_meson){
				sum_p += n_f[i]*C->X_p[i];
			}
//			printf("sum %f sum_ch %f sum_o %f sum_rho %f \n", sum, sum_ch, sum_o, sum_rho);
		}

		double mu_n = mu_deriv(n_in, m + sc + 1, sc + 0, mu_c, C);
		double mu_p = mu_deriv(n_in, m + sc + 1, sc + 1, mu_c, C);
		double mu_e = mu_n - mu_p;

		printf("mu_n = %.3f, mu_p = %.3f, mu_e = %.3f", mu_n, mu_p, mu_e);

//		printf("n = %f %f %f \n", n_in[0], n_in[1], n_in[2]);
//		printf("%f %f %f\n" ,mu_n, mu_p, sum_ch);
		double n_e = 0, n_mu = 0;
		if (mu_e > m_e){
			n_e += pow(mu_e*mu_e - m_e*m_e,1.5)/(3*M_PI*M_PI);
		}
		if (mu_e > m_mu){
			n_mu += pow(mu_e*mu_e - m_mu*m_mu,1.5)/(3*M_PI*M_PI);
		}

		double * eparts = new double[9];
		E_rho(n_in, m+sc+1, mu_c, C, eparts, 9);
		double n_rho = eparts[8];


		hx[0] = sum_ch - n_e - n_mu - n_rho;

		double fp = 0;
		double f = out[0];
		if (C->sprime){
			fp = out[1];
		}

		for (int i = 1; i < m-1; i++){
			double mu_i = mu_deriv(n_in, m + sc + 1, sc + i, mu_c, C);
			printf("mu_%i = %.3f", i, mu_i);
			hx[i] = mu_i - (mu_n - C->Q[i] * mu_e);
		}

		hx[m] = mu_c - mu_e;

		for (int i = 0; i < m; i++){
			printf("hx[%i] = %.6e ", i, hx[i]);
		}
		printf("\n");

		delete[] n_in;
		delete[] n_f;
	}


	void fun_n_eq_rho_anal(double * p, double * hx, int m, int n, void * adata){
		bool debug = 0;
		fun_n_eq_params * par = (fun_n_eq_params *) adata;
		set_const * C = par->C;
		int sc = 1 + C->sprime;
    int num_part = m - 1; //actual num. of particles. p = [n_1, n_2, ..., n_num, mu_c]
		double n_sum = 0.0;
		double n_n = par->n;
		double * n_in = new double [m+sc]; //input set for _E and mu
		double * n_f = new double [m+sc]; //input set for f_eq; actually is {n_n,n_p,...,n_X0}
		for (int i = 0; i < num_part; i++){
			n_n -= p[i];
			n_in[i + 1 + sc] = p[i]; //scalar + neutron(1) offset
			n_f[i+1] = p[i];
		}

		if (debug){
			for (int i = 0; i < m; i++){
				printf("p[%i]=%.6f,", i, p[i]);
			}
			printf("\n");
		}

		double mu_c = p[m-1];

		n_f[0] = n_n;
		if (debug) {
			printf("n_f = ");
			for (int i = 0; i < m; i++){
				printf("%e ", n_f[i]);
			}
			printf("\n");
		}
		double * out = new double[sc];
		f_eq_rho(n_f, m, par->f_init, sc, out, sc, mu_c, C);
    //par->f_init[0] = out[0];


		for (int i = 0; i < sc; i++){
			n_in[i] = out[i];
			if (debug)
				printf("f[%i] = %.3f ", i, out[i]);
		}
		if (debug) printf("\n");


		n_in[sc] = n_n;

        if (debug){
            printf("n_in = [");
            for (int i = 0; i < num_part + 1 + sc; i++){
                printf("%.6f ", n_in[i]);
            }
            printf("] \n");
        }
		double sum=0, sum_ch=0, sum_o = 0.0, sum_rho = 0.0, sum_p = 0;
		for (int i = 0; i < num_part + 1; i++){
			sum += n_f[i];
			sum_ch += n_f[i]*C->Q[i];
			sum_o += n_f[i]*C->X_o[i];
			sum_rho += n_f[i]*C->X_r[i]*C->T[i];
			if (C->phi_meson){
				sum_p += n_f[i]*C->X_p[i];
			}
//			printf("sum %f sum_ch %f sum_o %f sum_rho %f \n", sum, sum_ch, sum_o, sum_rho);
		}

		double mu_n = mu_rho(n_in, num_part + 1 + sc, sc + 0, mu_c, C);
		double mu_p = mu_rho(n_in, num_part + 1 + sc, sc + 1, mu_c, C);
		double mu_e = mu_n - mu_p;

		if (debug)
		printf("mu_n = %.3f, mu_p = %.3f, mu_e = %.3f \n", mu_n, mu_p, mu_e);

//		printf("n = %f %f %f \n", n_in[0], n_in[1], n_in[2]);
//		printf("%f %f %f\n" ,mu_n, mu_p, sum_ch);
		double n_e = 0, n_mu = 0;
		if (mu_e > m_e){
			n_e += pow(mu_e*mu_e - m_e*m_e,1.5)/(3*M_PI*M_PI);
		}
		if (mu_e > m_mu){
			n_mu += pow(mu_e*mu_e - m_mu*m_mu,1.5)/(3*M_PI*M_PI);
		}

//		double * eparts = new double[9];
//		E_rho(n_in, m+sc+1, mu_c, C, eparts, 9);
//		double n_rho = eparts[8];



		double fp = 0;
		double f = out[0];
		if (C->sprime){
			fp = out[1];
		}



    	double n_rho = 2 * C->m_rho * pow(C->M[0]*C->phi_n(0, f),2.)* sqrt(C->eta_r(f)) / (C->Cr * C->chi_prime
      (f)) *
    		(1 - mu_c/(C->m_rho * C->phi_n(0, f)));
    	par->misc = 0.;

    	if (debug)
    	printf("n_rho = %.6f, sum_rho = %.6f \n", n_rho, sum_rho);

		hx[0] = sum_ch - n_e - n_mu;

		double  n_c = 0.;

		if(fabs(sum_rho) > n_rho/2){
			double r_c2 = (fabs(sum_rho) - n_rho/2) / (2 * C->m_rho * C->chi_prime(f) * sqrt(C->eta_r(f)));
			n_c = 2 * C->m_rho * C->phi_n(0, f) * r_c2;

			hx[0] -= n_c;

			par->misc = n_c;

//			printf("sum_ch = %.6f, n_e = %.6f, n_mu = %.6f, n_c = %.6f \n", sum_ch, n_e, n_mu, n_c);
		}

		for (int i = 1; i < num_part; i++){
			hx[i] = p_f(p[i], 2*C->S[i+1]+1);
			double xs = 0.;
			if (C->sigma_kind == 0){
				xs = C->X_s[i+1];
			}
			else{
				xs = C->Xs(i+1, f);
			}
			double m_eff = C->M[i+1]*C->phi_n(i+1, xs * (C->M[0]/C->M[i+1]) * f  + C->X_sp[i+1] * (C->M[0]/C->M[i+1]) * fp);

			double res = mu_n - C->Q[i+1]*mu_e - C->Co/pow(C->M[0],2) * C->X_o[i+1] * sum_o / C->eta_o(f)
					- C->Cr/pow(C->M[0],2) * C->X_r[i+1]*C->T[i+1] * sum_rho / C->eta_r(f)
					- C->Co/pow(C->M[0], 2) * C->X_p[i+1] * sum_p * pow(m_o / m_p,2.0) / C->eta_p(f);

	    	if (fabs(sum_rho) > n_rho/2){
	    		res += C->Cr / (pow(C->M[0],2.) * C->eta_r(f)) * (fabs(sum_rho) - n_rho/2) * C->X_r[i+1] * C->T[i+1] *
	    				((sum_rho > 0) - (sum_rho < 0));
	    	}
	    	res = pow(res, 2.0);

			res -= m_eff*m_eff;

			if (res > 0){
				hx[i] -= sqrt(res);
			}
		}

		hx[m-1] = mu_c - mu_e;
		if (debug)
		printf("mu_c: %.6f \n", mu_c);

		if (debug){
		for (int i = 0; i < m; i++){
			printf("hx[%i] = %.6e ", i, hx[i]);
		}
		printf("\n");
		}

		delete[] n_in;
		delete[] n_f;
	}


void fun_n_eq_rho_anal2(double * p, double * hx, int m, int n, void * adata){
		bool debug = 0;
		fun_n_eq_params * par = (fun_n_eq_params *) adata;
		set_const * C = par->C;
		int sc = 1 + C->sprime;
    int num_part = m - 1; //actual num. of particles. p = [n_1, n_2, ..., n_num, mu_c]
		double n_sum = 0.0;
		double n_n = par->n;
		double * n_in = new double [m+sc - 1]; //input set for _E and mu
		double * n_f = new double [m+sc]; //input set for f_eq; actually is {n_n,n_p,...,n_X0}
		for (int i = 0; i < num_part; i++){
			n_n -= p[i];
			n_in[i + 1 + sc] = p[i]; //scalar + neutron(1) offset
			n_f[i+1] = p[i];
		}

		if (debug){
			for (int i = 0; i < m; i++){
				printf("p[%i]=%.6f,", i, p[i]);
			}
			printf("\n");
		}

		double mu_c = p[m-2];
    double f = p[m-1];


		n_f[0] = n_n;
		if (debug) {
			printf("n_f = ");
			for (int i = 0; i < m; i++){
				printf("%e ", n_f[i]);
			}
			printf("\n");
		}
//		double * out = new double[sc];
//		f_eq_rho(n_f, m, par->f_init, sc, out, sc, mu_c, C);


//		for (int i = 0; i < sc; i++){
//			n_in[i] = out[i];
//			if (debug)
//				printf("f[%i] = %.3f ", i, out[i]);
//		}
//		if (debug) printf("\n");


    n_in[0] = f;
		n_in[sc] = n_n;

        if (debug){
            printf("n_in = [");
            for (int i = 0; i < num_part + 1 + sc; i++){
                printf("%.6f ", n_in[i]);
            }
            printf("] \n");
        }
		double sum=0, sum_ch=0, sum_o = 0.0, sum_rho = 0.0, sum_p = 0;
		for (int i = 0; i < num_part + 1; i++){
			sum += n_f[i];
			sum_ch += n_f[i]*C->Q[i];
			sum_o += n_f[i]*C->X_o[i];
			sum_rho += n_f[i]*C->X_r[i]*C->T[i];
			if (C->phi_meson){
				sum_p += n_f[i]*C->X_p[i];
			}
//			printf("sum %f sum_ch %f sum_o %f sum_rho %f \n", sum, sum_ch, sum_o, sum_rho);
		}

		double mu_n = mu_rho(n_in, num_part + sc, sc + 0, mu_c, C);
		double mu_p = mu_rho(n_in, num_part + sc, sc + 1, mu_c, C);
		double mu_e = mu_n - mu_p;

		if (debug)
		printf("mu_n = %.3f, mu_p = %.3f, mu_e = %.3f \n", mu_n, mu_p, mu_e);

//		printf("n = %f %f %f \n", n_in[0], n_in[1], n_in[2]);
//		printf("%f %f %f\n" ,mu_n, mu_p, sum_ch);
		double n_e = 0, n_mu = 0;
		if (mu_e > m_e){
			n_e += pow(mu_e*mu_e - m_e*m_e,1.5)/(3*M_PI*M_PI);
		}
		if (mu_e > m_mu){
			n_mu += pow(mu_e*mu_e - m_mu*m_mu,1.5)/(3*M_PI*M_PI);
		}

//		double * eparts = new double[9];
//		E_rho(n_in, m+sc+1, mu_c, C, eparts, 9);
//		double n_rho = eparts[8];



		double fp = 0;
//		double f = out[0];
//		if (C->sprime){
//			fp = out[1];
//		}



    	double n_rho = 2 * C->m_rho * pow(C->M[0]*C->phi_n(0, f),2.)* sqrt(C->eta_r(f)) / (C->Cr * C->chi_prime
      (f)) *
    		(1 - mu_c/(C->m_rho * C->phi_n(0, f)));
    	par->misc = 0.;

    	if (debug)
    	printf("n_rho = %.6f, sum_rho = %.6f \n", n_rho, sum_rho);

		hx[0] = sum_ch - n_e - n_mu;

		double  n_c = 0.;

		if (fabs(sum_rho) > n_rho/2){
			double r_c2 = (fabs(sum_rho) - n_rho/2) / (2 * C->m_rho * C->chi_prime(f) * sqrt(C->eta_r(f)));
			n_c = 2 * C->m_rho * C->phi_n(0, f) * r_c2;

			hx[0] -= n_c;

			par->misc = n_c;

//			printf("sum_ch = %.6f, n_e = %.6f, n_mu = %.6f, n_c = %.6f \n", sum_ch, n_e, n_mu, n_c);
		}

		for (int i = 1; i < num_part; i++){
			hx[i] = p_f(p[i], 2*C->S[i+1]+1);
			double xs = 0.;
			if (C->sigma_kind == 0){
				xs = C->X_s[i+1];
			}
			else{
				xs = C->Xs(i+1, f);
			}
			double m_eff = C->M[i+1]*C->phi_n(i+1, xs * (C->M[0]/C->M[i+1]) * f  + C->X_sp[i+1] * (C->M[0]/C->M[i+1]) * fp);

			double res = mu_n - C->Q[i+1]*mu_e - C->Co/pow(C->M[0],2) * C->X_o[i+1] * sum_o / C->eta_o(f)
					- C->Cr/pow(C->M[0],2) * C->X_r[i+1]*C->T[i+1] * sum_rho / C->eta_r(f)
					- C->Co/pow(C->M[0], 2) * C->X_p[i+1] * sum_p * pow(m_o / m_p,2.0) / C->eta_p(f);

	    	if (fabs(sum_rho) > n_rho/2){
	    		res += C->Cr / (pow(C->M[0],2.) * C->eta_r(f)) * (fabs(sum_rho) - n_rho/2) * C->X_r[i+1] * C->T[i+1] *
	    				((sum_rho > 0) - (sum_rho < 0));
	    	}
	    	res = pow(res, 2.0);

			res -= m_eff*m_eff;

			if (res > 0){
				hx[i] -= sqrt(res);
			}
		}

		hx[m-2] = mu_c - mu_e;
    hx[m-1] = wrap_func_feq_rho(f, n_in, m + sc - 1, mu_c, C);
    
		if (debug)
		printf("mu_c: %.6f \n", mu_c);

		if (debug){
		for (int i = 0; i < m; i++){
			printf("hx[%i] = %.6e ", i, hx[i]);
		}
		printf("\n");
		}

		delete[] n_in;
		delete[] n_f;
	}


}//namespace calc

/**\brief Energy density functional itself.
 * Provides the energy density functional evaluated at some point v = {n, f}
 *
 * Out consists of:
 * 0: m_n^4 f^2/ ...
 * 1: U(f)
 * 2 -- (dimN-1): Kinetic integrals
 * dimN -- dimN + 3:
 *
 */


double _E(double * n, int dimN, set_const * C, double * out, int dim_Out){
	bool debug = 0;
	bool ret_parts = (out) && (dim_Out == 9);
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
		fp = 0.;
	}



	double res = 0.;
	double part_s = pow(C->M[0], 4.0)*f*f*C->eta_s(f)/(2*C->Cs);

	res += part_s;

	if (debug){
		printf("res_f : %f\n", part_s);
	}
	if (ret_parts){
		out[0] = part_s;
	}

	double part_U = C->U(f);
	res += part_U;
	if (ret_parts){
		out[1] = part_U;
	}

	res += pow(C->M[0], 4.0)*fp*fp/(2*C->Csp);

	double sum = 0;
	double sum_t3 = 0;
	double sum_p = 0;
	double meff_arg = 0;
	if (debug){
		printf("res_Uf : %f \n", res);
	}
	double kin = 0.;
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
		kin = kineticInt(n[i], (C->M)[i-sc] * C->phi_n(i-sc,meff_arg), 2*C->S[i-sc] + 1);
		res += kin;
		if (i == 1 or i == 2){
			if (ret_parts){
				out[i+1] = kin;
			}
		}
//		printf("i = %i, n[i] = %f, pf(n[i]) = %f, meff_arg = %f, spin = %f \n",
//			   i, n[i], p_f(n[i], 2*C->S[i]+1), meff_arg, 2*C->S[i-sc] + 1);
//		printf("K = %f \n", kineticInt(n[i], (C->M)[i-sc] * C->phi_n(i-sc,meff_arg), 2*C->S[i-sc] + 1));
//		printf("M_PI = %f \n", M_PI);
//		printf("asinh(1) = %f\n", asinh(1.0));
		sum += n[i] * C->X_o[i-sc];
		sum_t3 += n[i]*(C->T)[i-sc] * C->X_r[i-sc];
		if (C->phi_meson){
			sum_p += n[i]*C->X_p[i-sc];
		}
	}
	//omega
	double part_om = C->Co * sum*sum/(2.0*C->M[0]*C->M[0]*C->eta_o(f));
	res += part_om;
	if (debug){
		printf("res_om : %f \n", part_om);
	}
	if (ret_parts){
		out[4] = part_om;
	}
	//rho
	double part_rho = C->Cr * pow(sum_t3/C->M[0], 2.0) / (2 * C->eta_r(f));
	res += part_rho;
	if (debug){
		printf("res_rho : %f \n", res);
	}
	if (ret_parts){
		out[5] = part_rho;
	}
	//phi
	double part_phi = pow(m_o/m_p ,2.0)*C->Co * sum_p*sum_p/(2.0*C->M[0]*C->M[0]*C->eta_p(f));
	res += part_phi;
	if (debug){
		printf("res_phi : %f \n", part_phi);
	}
	if (ret_parts){
		out[6] = part_phi;
	}


	return res;
}


double E_rho(double * n, int dimN, double mu_c, set_const * C, double *inplace, int dim_inplace){
	bool debug = 0;
	bool ret_parts = (inplace) && (dim_inplace == 10);
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
		fp = 0.;
	}

	double res = 0.;
	double part_s = pow(C->M[0], 4.0)*f*f*C->eta_s(f)/(2*C->Cs);

	res += part_s;

	if (debug){
		printf("res_f : %f\n", part_s);
	}
	if (ret_parts){
		inplace[0] = part_s;
	}

	double part_U = C->U(f);
	res += part_U;
	if (ret_parts){
		inplace[1] = part_U;
	}

	res += pow(C->M[0], 4.0)*fp*fp/(2*C->Csp);

	double sum = 0;
	double sum_t3 = 0;
	double sum_p = 0;
	double meff_arg = 0;
	if (debug){
		printf("res_Uf : %f \n", res);
	}
	double kin = 0.;
	for (int i = sc; i < dimN; ++i){
		double xs = 0.;
		if (C->sigma_kind == 0){
			xs = C->X_s[i-sc];
		}
		else{
			xs = C->Xs(i-sc, f);
		}
		meff_arg = xs * (C->M[0]/C->M[i-sc]) * f + C->X_sp[i-sc] * (C->M[0]/C->M[i-sc])*fp;
		kin = kineticInt(n[i], (C->M)[i-sc] * C->phi_n(i-sc,meff_arg), 2*C->S[i-sc] + 1);
		res += kin;
		if (i == 1 or i == 2){
			if (ret_parts){
				inplace[i+1] = kin;
			}
		}
		sum += n[i] * C->X_o[i-sc];
		sum_t3 += n[i]*(C->T)[i-sc] * C->X_r[i-sc];
		if (C->phi_meson){
			sum_p += n[i]*C->X_p[i-sc];
		}
	}
	//omega
	double part_om = C->Co * sum*sum/(2.0*C->M[0]*C->M[0]*C->eta_o(f));
	res += part_om;
	if (debug){
		printf("res_om : %f \n", part_om); 
	}
	if (ret_parts){
		inplace[4] = part_om;
	}
	//rho
	double gr = sqrt(C->Cr/C->eta_r((f)))*(C->m_rho* (1-f)) / C->M[0];
	double n_rho = 2 * C->m_rho * pow(C->M[0]*C->phi_n(0, f),2.)* sqrt(C->eta_r(f)) / (C->Cr * C->chi_prime(f)) *
		(1 - mu_c/(C->m_rho * C->phi_n(0, f)));

	double E_r =  C->Cr * sum_t3*sum_t3/(2.0*C->M[0]*C->M[0]*C->eta_r(f));
	if (debug){
		printf("sum_t3 = %f, n_rho = %f, E_r = %f \n", sum_t3, n_rho, E_r);
		printf("Check fabs(sum_t3) > n_rho/2: %i \n", fabs(sum_t3) > n_rho/2);
	}
	if (fabs(sum_t3) > n_rho/2){
	    E_r -= C->Cr / (2*pow(C->M[0],2.) * C->eta_r(f)) * pow(fabs(-sum_t3) - n_rho/2, 2.);
	    if (debug)
	    printf("rho-condensate on! E_c = %.6f \n", C->Cr / (2*pow(C->M[0],2.) * C->eta_r(f)) * pow(fabs(sum_t3) - n_rho/2, 2.));
	};

	//phi
	double part_phi = pow(m_o/m_p ,2.0)*C->Co * sum_p*sum_p/(2.0*C->M[0]*C->M[0]*C->eta_p(f));
	res += part_phi;
	if (debug){
		printf("res_phi : %f \n", part_phi);
	}

	if (ret_parts) {
		inplace[6] = part_phi;
		inplace[8] = n_rho;
	}

	return res + E_r;
}

double E(double* n, int dimN, set_const* C,  double * out, int dim_Out) {
	double f = n[0];
	bool ret_parts = (out) && (dim_Out == 9);
	double res = _E(n, dimN, C, out, dim_Out);


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
	double kin_e = kineticInt(n_e, m_e, 2.);
	res += kin_e;
	double kin_mu = kineticInt(n_mu, m_mu, 2.);
	res += kin_mu;
	if (ret_parts){
		out[7] = kin_e;
		out[8] = kin_mu;
	}
	return res;
}

double stepF(var v, set_const *C){
	return 0.0;
}
//Stepper function for E
void stepE(double n, double * init, int initN, double * f_init, int dimF_init, double * out, int dim_Out, int iter, set_const* C) {
	double opts[5];
	bool debug = 1;
	calc::fun_n_eq_params p = {C, n, f_init, dimF_init, 0.0};
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

void stepE_f(double f, double * init, int initN, double * out, int dim_Out, int iter, set_const* C) {
	double opts[5];
	bool debug = 1;
	calc::fun_n_eq_f_params p = {C, f, 0.0};
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
	printf("Entering levmar: \n");
	dlevmar_bc_dif(calc::fun_n_eq_f, x, NULL, m, m, lb, NULL, NULL, iter, opts, info, NULL, NULL, &p);
//	dlevmar_dif(calc::fun_n_eq, x, NULL, m, m, 2000, opts, NULL, NULL, NULL, &p);

	if (debug) {
		printf("info: ");
		for (int i = 0; i < LM_INFO_SZ; i++){
			printf(" %i : %f ", i, info[i]);
		}
		printf("\n");

		printf("f = %f, n_p = %e", f, x[0]);
		printf(",n_L = %e ", x[1]);
		printf(",n_S- = %e ", x[2]);
		printf(",n_S0 = %e ", x[3]);
		printf(",n_S+ = %e ", x[4]);
		printf(",n_X- = %e ", x[5]);
		printf(",n_X0 = %e ", x[6]);
		printf("\n");

		calc::fun_n_eq_f(x, fun, m, m, &p);
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


void stepE_dsym(double n, double * init, int initN, double * f_init, int dimF_init, double * out, int dim_Out, int iter, set_const* C) {
	double opts[5];
	bool debug = 1;
	calc::fun_n_eq_params p = {C, n, f_init, dimF_init, 0.0};
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
	dlevmar_bc_dif(calc::fun_n_eq_dsym, x, NULL, m, m, lb, NULL, NULL, iter, opts, info, NULL, NULL, &p);
//	dlevmar_dif(calc::fun_n_eq, x, NULL, m, m, 2000, opts, NULL, NULL, NULL, &p);

	if (debug) {
		printf("info: ");
		for (int i = 0; i < LM_INFO_SZ; i++){
			printf(" %i : %f ", i, info[i]);
		}
		printf("\n");

		printf("n = %f, n_d = %e", n, x[0]);
		printf("\n");

		calc::fun_n_eq_dsym(x, fun, m, m, &p);
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
	hx[0] -= nn*pow(pow(p_f(nn, 2.),2) + pow(mn*(1-f),2), 0.5) +
			 np*pow(pow(p_f(np, 2.),2) + pow(mn*(1-f),2), 0.5);
	double ne = 0.;
	double nmu = 0.;
	double mu_e = pow(pow(p_f(nn, 2.),2) + pow(mn*(1-f),2), 0.5) -
				  pow(pow(p_f(np, 2.),2) + pow(mn*(1-f),2), 0.5) +
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


void stepE_rho(double n, double * init, int initN, double * f_init, int dimF_init, double * out, int dim_Out, int iter, double mu_init, set_const* C) {
	double opts[5];
	bool debug = 1;
	calc::fun_n_eq_params p = {C, n, f_init, dimF_init, 0.0};
	int m = initN;
	double * x = new double[m];
	double * lb = new double[m];
	double * ub = new double[m];
	double * fun = new double[m];
	double * scale = new double [m];
  double * D = new double[m]; //linear inequality constraints matrix
  double * b = new double[1]; //linear inequality rhs
  for (int i = 0; i < m; i++){
    D[i] = -1;
  }
  D[m-1] = 0;
  
  b[0] = -n;

	double info[LM_INFO_SZ];
	//double x[3] = {v.n[0], v.n[1], v.f};
  

	for (int i = 0; i < m-1; i++){
		x[i] = init[i];
		lb[i] = 0.0;
    ub[i] = n;
		scale[i] = 0.1;
//		if (i > 2) lb[i] = -100500.0;
	}
	lb[m-1] = 0.;
  ub[m-1] = C->M[0]/2;
	x[m-1] = init[m-1];
	scale[m-1] = 1.;

	opts[0]= LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-25; opts[3]=1E-12;
		opts[4]= -1e-5;
  info[6] = 3;
//  while ((int)info[6] == 3){
//
//

	dlevmar_bc_dif(calc::fun_n_eq_rho_anal, x, NULL, m, m, lb, ub, scale, iter, opts, info, NULL, NULL, &p);


//	dlevmar_bleic_dif(calc::fun_n_eq_rho_anal, x, NULL, m, m, lb, ub, NULL, NULL, 0, D, b, 1,iter, opts, info, NULL, NULL, &p);



//	dlevmar_dif(calc::fun_n_eq, x, NULL, m, m, 2000, opts, NULL, NULL, NULL, &p);
//  if ((int)info[6] == 3){
    //if (p.f_init[0] + 0.1 < 1.){
    //  p.f_init[0] = p.f_init[0] + 0.1;
    //}
    //
    
//    p.f_init[0] += 0.05;
//  }
//  }

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

		calc::fun_n_eq_rho_anal(x, fun, m, m, &p);
		for (int i = 0; i < m; i++){
			printf("f%i = %e  ", i, fun[i]);
		}
		printf("\n");
	}
	for (int i = 0; i < m; i++){
		out[i] = x[i];
	}

	out[m] = p.misc;
	delete[] x;
	delete[] fun;
	delete[] lb;
}

void stepE_rho_f(double f, double * init, int initN, double * out, int dim_Out, int iter, set_const* C) {
	double opts[5];
	bool debug = 1;
	calc::fun_n_eq_f_params p = {C, f, 0.0};
	int m = initN;
	double * x = new double[m];
	double * lb = new double[m];
	double * ub = new double[m];
	double * fun = new double[m];
	double * scale = new double [m];

	double info[LM_INFO_SZ];

	for (int i = 0; i < m-1; i++){
		x[i] = init[i];
		lb[i] = 0.0;
		scale[i] = 1.;
	}
	lb[m-1] = 0.;
	x[m-1] = init[m-1];
	scale[m-1] = 1.;

	opts[0]= LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-25; opts[3]=1E-12;
		opts[4]= -1e-5;
  info[6] = 3;

	dlevmar_bc_dif(calc::fun_n_eq_f_rho, x, NULL, m, m, lb, NULL, scale, iter, opts, info, NULL, NULL, &p);

	if (debug) {
		printf("info: ");
		for (int i = 0; i < LM_INFO_SZ; i++){
			printf(" %i : %f ", i, info[i]);
		}
		printf("\n");

		printf("f = %f, n_n = %e", f, x[0]);
		printf(",n_p = %e ", x[1]);
		printf(",mu = %e ", x[2]);
		printf("\n");

		calc::fun_n_eq_f_rho(x, fun, m, m, &p);
		for (int i = 0; i < m; i++){
			printf("f%i = %e  ", i, fun[i]);
		}
		printf("\n");
	}
	for (int i = 0; i < m; i++){
		out[i] = x[i];
	}
	if (dim_Out > m) out[m] = p.misc;
//	printf("return\n");
	delete[] x;
//	printf("deleted x\n");
	delete[] fun;
//	printf("deleted fun\n");
	delete[] lb;
	delete[] ub;
	delete[] scale;
//	printf("deleted lb\n");
}


void wrap_fun(double * n, int dimN, set_const * C, double * out, int dim_Out){
  double *n_in = new double[dimN-2];
  for (int i = 0; i < dimN-2; i++){
    n_in[i] = n[i+2];
  }
  double ntot = 0;
  for (int i = 0; i < dimN-1; i++){
    ntot += n[i+1];
  }
//  printf("ntot = %.6f \n", ntot);
  double finit[1] = {n[0]};
  calc::fun_n_eq_params p = {C, ntot, finit, 1,  0.};
  calc::fun_n_eq(n_in, out, dimN-2, dimN-2, &p);
  delete[] n_in;
  return;
}

void wrap_fun_rho(double * n, int dimN, set_const * C, double * out, int dim_Out){
  double *n_in = new double[dimN-2];
  for (int i = 0; i < dimN-2; i++){
    n_in[i] = n[i+2];
  //  printf("%.6f ", n_in[i]);
  }
  //printf("\n");
  double ntot = 0;
  for (int i = 1; i < dimN-1; i++){
    ntot += n[i];
  }
  //printf("ntot = %.6f \n", ntot);
  double finit[1] = {n[0]};
  calc::fun_n_eq_params p = {C, ntot, finit, 1,  0.};
  calc::fun_n_eq_rho_anal(n_in, out, dimN-2, dimN-2, &p);
  delete[] n_in;
  return;
}

double wrap_fun_dsym(double n, double nd, double f_init, set_const * C){
  double n_n = (n - nd)/2;
  double n_in[6] = {n_n, n_n, nd/4, nd/4, nd/4, nd/4};
  double * vn = new double[1];
//  for (int i = 0; i < 6; i++){
//	  vn[i] = n_in[i];
//  }
  vn[0] = nd;

  //printf("\n");
  double * finit = new double[1];
  finit[0] = f_init;
  calc::fun_n_eq_params p = {C, n, finit, 1,  0.};
  double hx[1];
  void* vp = &p;

  double * res = new double[1];
  calc::fun_n_eq_dsym(vn, res, 1, 1, vp);
  double out = res[0];
  delete[] res;
  return out;
}

double wrap_fun_dsym_f(double n, double nd, double f_init, set_const * C){
  double n_n = (n - nd)/2;
  double n_in[6] = {n_n, n_n, nd/4, nd/4, nd/4, nd/4};
  double * vn = new double[1];
//  for (int i = 0; i < 6; i++){
//	  vn[i] = n_in[i];
//  }
  vn[0] = nd;

  //printf("\n");
  double * finit = new double[1];
  finit[0] = f_init;
  calc::fun_n_eq_params p = {C, n, finit, 1,  0.};
  double hx[1];
  void* vp = &p;

  double * res = new double[1];
  calc::fun_n_eq_dsym_f(vn, res, 1, 1, vp);
  double out = res[0];
  delete[] res;
  return out;
}


void wrap_fun_np_f(double nn, double np, double f, set_const * C, double * out, int dimOut){
  double n_in[2] = {nn, np};

  calc::fun_n_eq_f_params p = {C, f};
  void* vp = &p;

  double * res = new double[2];
  calc::fun_n_eq_f(n_in, res, 1, 1, vp);
  out[0] = res[0];
  out[1] = res[1];
  delete[] res;
}

	
void stepE_rho2(double n, double * init, int initN, double * f_init, int dimF_init, double * out, int dim_Out, int iter, double mu_init, set_const* C) {
	double opts[5];
	bool debug = 1;
	calc::fun_n_eq_params p = {C, n, f_init, dimF_init, 0.0};
	int m = initN;
	double * x = new double[m];
	double * lb = new double[m];
	double * ub = new double[m];
	double * fun = new double[m];
	double * scale = new double [m];
  double * D = new double[m]; //linear inequality constraints matrix
  double * b = new double[1]; //linear inequality rhs
  for (int i = 0; i < m; i++){
    D[i] = -1;
  }
  D[m-1] = 0;
  
  b[0] = -n;

	double info[LM_INFO_SZ];
	//double x[3] = {v.n[0], v.n[1], v.f};
  

	for (int i = 0; i < m-1; i++){
		x[i] = init[i];
		lb[i] = 0.0;
    ub[i] = n;
		scale[i] = 0.1;
//		if (i > 2) lb[i] = -100500.0;
	}
	lb[m-2] = 0.;
  ub[m-2] = C->M[0]/2;
	x[m-2] = init[m-2];
	scale[m-2] = 1.;
  lb[m-1] = 0.;
  ub[m-1] = 1.;

  x[m-1] = init[m-1];


	opts[0]= LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-25; opts[3]=1E-12;
		opts[4]= -1e-5;
  info[6] = 3;
//  while ((int)info[6] == 3){
//
//

	dlevmar_bc_dif(calc::fun_n_eq_rho_anal2, x, NULL, m, m, lb, ub, scale, iter, opts, info, NULL, NULL, &p);


//	dlevmar_bleic_dif(calc::fun_n_eq_rho_anal, x, NULL, m, m, lb, ub, NULL, NULL, 0, D, b, 1,iter, opts, info, NULL, NULL, &p);



//	dlevmar_dif(calc::fun_n_eq, x, NULL, m, m, 2000, opts, NULL, NULL, NULL, &p);
//  if ((int)info[6] == 3){
    //if (p.f_init[0] + 0.1 < 1.){
    //  p.f_init[0] = p.f_init[0] + 0.1;
    //}
    //
    
//    p.f_init[0] += 0.05;
//  }
//  }

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

		calc::fun_n_eq_rho_anal(x, fun, m, m, &p);
		for (int i = 0; i < m; i++){
			printf("f%i = %e  ", i, fun[i]);
		}
		printf("\n");
	}
	for (int i = 0; i < m; i++){
		out[i] = x[i];
	}

	out[m] = p.misc;
	delete[] x;
	delete[] fun;
	delete[] lb;
}



