#include "fixedNF.h"
#include "eos.h"
#include "aux.h"
#include "levmar.h"

using namespace calc;

struct fun_n_eq_params_nf : fun_n_eq_params {
    double f;
};

void fun_n_eq_rho_nf(double * p, double * hx, int m, int n, void * adata){
    bool debug = 0;
    fun_n_eq_params_nf * par = (fun_n_eq_params_nf *) adata;
    set_const * C = par->C;
    double f = par->f;
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

    double n_rho = 2 * C->m_rho * pow(C->M[0]*C->phi_rho(f),2.)* sqrt(C->eta_r(f)) / (C->Cr * C->chi_prime
    (f)) *
        (1 - mu_c/(C->m_rho * C->phi_rho(f)));
    par->misc = 0.;

    if (debug)
    printf("n_rho = %.6f, sum_rho = %.6f \n", n_rho, sum_rho);

    hx[0] = sum_ch - n_e - n_mu;

    double  n_c = 0.;

    if(fabs(sum_rho) > n_rho/2){
        double r_c2 = (fabs(sum_rho) - n_rho/2) / (2 * C->m_rho * C->chi_prime(f) * sqrt(C->eta_r(f)));
        n_c = 2 * C->m_rho * C->phi_rho(f) * r_c2;

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

void fun_n_eq_nf(double * p, double * hx, int m, int n, void * adata){
    bool debug = 0;
    fun_n_eq_params_nf * par = (fun_n_eq_params_nf *) adata;
    set_const * C = par->C;
    double f = par->f;
    int sc = 1 + C->sprime;
    int num_part = m; //actual num. of particles. p = [n_1, n_2, ..., n_num, mu_c]
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

    n_f[0] = n_n;
    if (debug) {
        printf("n_f = ");
        for (int i = 0; i < m; i++){
            printf("%e ", n_f[i]);
        }
        printf("\n");
    }

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

    double mu_n = mu(n_in, num_part + 1 + sc, sc + 0, C);
    double mu_p = mu(n_in, num_part + 1 + sc, sc + 1, C);
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

    // double n_rho = 2 * C->m_rho * pow(C->M[0]*C->phi_/rho(f),2.)* sqrt(C->eta_r(f)) / (C->Cr * C->chi_prime
    // (f)) *
    //     (1 - mu_c/(C->m_rho * C->phi_rho(f)));
    par->misc = 0.;

    // if (debug)
    // printf("n_rho = %.6f, sum_rho = %.6f \n", n_rho, sum_rho);

    hx[0] = sum_ch - n_e - n_mu;

    // double  n_c = 0.;

//     if(fabs(sum_rho) > n_rho/2){
//         double r_c2 = (fabs(sum_rho) - n_rho/2) / (2 * C->m_rho * C->chi_prime(f) * sqrt(C->eta_r(f)));
//         n_c = 2 * C->m_rho * C->phi_rho(f) * r_c2;

//         hx[0] -= n_c;

//         par->misc = n_c;

// //			printf("sum_ch = %.6f, n_e = %.6f, n_mu = %.6f, n_c = %.6f \n", sum_ch, n_e, n_mu, n_c);
//     }

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

        // if (fabs(sum_rho) > n_rho/2){
        //     res += C->Cr / (pow(C->M[0],2.) * C->eta_r(f)) * (fabs(sum_rho) - n_rho/2) * C->X_r[i+1] * C->T[i+1] *
        //             ((sum_rho > 0) - (sum_rho < 0));
        // }
        res = pow(res, 2.0);

        res -= m_eff*m_eff;

        if (res > 0){
            hx[i] -= sqrt(res);
        }
    }

    // hx[m-1] = mu_c - mu_e;
    if (debug)
    // printf("mu_c: %.6f \n", mu_c);

    if (debug){
    for (int i = 0; i < m; i++){
        printf("hx[%i] = %.6e ", i, hx[i]);
    }
    printf("\n");
    }

    delete[] n_in;
    delete[] n_f;
}

void stepE_rho_nf(double n, double f, double * init, int dimInit, double * out, int dimOut, int iter, double mu_init, set_const* C) {
	double opts[5];
	bool debug = 1;
	// fun_n_eq_params_nf p = {C, n, NULL, 1, 0.0, f};
    fun_n_eq_params_nf p;
    p.C = C;
    p.n = n;
    p.f_init = NULL;
    p.dimF_init = 0;
    p.f = f;
	int m = dimInit;
	double * x = new double[m];
	double * lb = new double[m];
	double * ub = new double[m];
	double * fun = new double[m];
	double * scale = new double [m];

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

	dlevmar_bc_dif(&fun_n_eq_rho_nf, x, NULL, m, m, lb, ub, scale, iter, opts, info, NULL, NULL, &p);

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

		fun_n_eq_rho_nf(x, fun, m, m, &p);
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

void stepE_nf(double n, double f, double * init, int dimInit, double * out, int dimOut, int iter, double mu_init, set_const* C) {
	double opts[5];
	bool debug = 1;
	// fun_n_eq_params_nf p = {C, n, NULL, 1, 0.0, f};
    fun_n_eq_params_nf p;
    p.C = C;
    p.n = n;
    p.f_init = NULL;
    p.dimF_init = 0;
    p.f = f;
	int m = dimInit;
	double * x = new double[m];
	double * lb = new double[m];
	double * ub = new double[m];
	double * fun = new double[m];
	double * scale = new double [m];

	double info[LM_INFO_SZ];
	//double x[3] = {v.n[0], v.n[1], v.f};
  

	for (int i = 0; i < m; i++){
		x[i] = init[i];
		lb[i] = 0.0;
    ub[i] = n;
		scale[i] = 0.1;
//		if (i > 2) lb[i] = -100500.0;
	}
// 	lb[m-1] = 0.;
//   ub[m-1] = C->M[0]/2;
// 	x[m-1] = init[m-1];
	// scale[m-1] = 1.;

	opts[0]= LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-25; opts[3]=1E-12;
		opts[4]= -1e-5;
  info[6] = 3;

	dlevmar_bc_dif(&fun_n_eq_nf, x, NULL, m, m, lb, ub, scale, iter, opts, info, NULL, NULL, &p);

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

		fun_n_eq_rho_nf(x, fun, m, m, &p);
		for (int i = 0; i < m; i++){
			printf("f%i = %e  ", i, fun[i]);
		}
		printf("\n");
	}
	for (int i = 0; i < m; i++){
		out[i] = x[i];
	}

	// out[m] = p.misc;
	delete[] x;
	delete[] fun;
	delete[] lb;
}