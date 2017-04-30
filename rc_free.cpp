#include "rc_free.h"
#include "eos.h"
#include "aux.h"
#include "levmar.h"

using namespace calc;

double E_rho_free(double * n, int dimN, double rc2, set_const * C, double *inplace, int dim_inplace){
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
    double sum_ch = 0.;
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
        sum_ch += n[i] * C->Q[i];
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
	

	double E_r =  C->Cr * sum_t3*sum_t3/(2.0*C->M[0]*C->M[0]*C->eta_r(f));
   
    double n_e = 0;
    double n_mu = 0;
    
    double m_rho_eff = C->m_rho * C->phi_n(0, f);


    E_r += m_rho_eff * m_rho_eff * rc2;


	//phi
	double part_phi = pow(m_o/m_p ,2.0)*C->Co * sum_p*sum_p/(2.0*C->M[0]*C->M[0]*C->eta_p(f));
	res += part_phi;
	if (debug){
		printf("res_phi : %f \n", part_phi);
	}

	if (ret_parts) {
		inplace[6] = part_phi;
		//inplace[8] = n_rho;
	}

	return res + E_r;
}

void func_f_eq_rc_free(double * p, double * hx, int m, int _n, void * adata){
	bool debug = 0;
	func_f_eq_params_rho * params = (func_f_eq_params_rho *) adata;
	bool sprime = (params->C->sprime and (m > 1));
	double * n = new double[params->dimN + m];
    double rc2 = params->mu_c;
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
		n[i] += 3*df;
		double dE = E_rho_free(n, params->dimN + m, rc2, params->C);
		n[i] -= df;
		dE += -9 * E_rho_free(n, params->dimN + m, rc2, params->C);
		n[i] -= df;
		dE += 45 * E_rho_free(n, params->dimN + m, rc2, params->C);
		n[i] -= 2*df;
		dE += -45 * E_rho_free(n, params->dimN + m, rc2, params->C);
		n[i] -= df;
		dE += 9 * E_rho_free(n, params->dimN + m, rc2, params->C);
		n[i] -= df;
		dE += -1 * E_rho_free(n, params->dimN + m, rc2, params->C);
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


double wrap_eqf_rc_free(double f, double * n, int dimN, double rc2, set_const * C){
    func_f_eq_params_rho params = {n, dimN, 1e-6, C, rc2};
    double res[1];
    double f_in[1] = {f};
    func_f_eq_rc_free(f_in, res, 1, 1, &params);
    return res[0];
}

void fun_n_eq_free(double * p, double * hx, int m, int n, void * adata){
    bool debug = 0;
    fun_n_eq_params * par = (fun_n_eq_params *) adata;
    set_const * C = par->C;
    double * n_in = new double [3];
    double * n_f = new double [2];
    double f = p[0];

    n_in[0] = f;
    n_in[1] = par->n - p[1];
    n_in[2] = p[1];

    n_f[0] = n_in[1]; 
    n_f[1] = n_in[2];

    double rc2 = p[2];
    func_f_eq_params_rho f_par = {n_f, 2, 1e-6, C, rc2};
    double res[1];
    double f_in[1] = {f};
    func_f_eq_rc_free(f_in, res, 1, 1, &f_par);

    hx[0] = res[0];

    double mu_n = mu (n_in, m + 1, 1, C);
    double mu_p = mu(n_in, m + 1, 2, C);
    double mu_e = mu_n - mu_p;

    double n_e = 0;
    double n_mu = 0;
    
    //double mu_e = C->m_rho * C->phi_n(0, f);

    if (mu_e * mu_e > m_e * m_e){
        n_e += pow (mu_e * mu_e - m_e * m_e, 1.5) / (3 * pow(M_PI,2));
    }

    if (mu_e * mu_e > m_mu * m_mu){
        n_mu += pow (mu_e * mu_e - m_mu * m_mu, 1.5) / (3 * pow(M_PI,2));
    }

    double Q = n_in[2] - n_e - n_mu - 2 * C->m_rho * C->phi_n(0, f) * rc2;

	// if (Q > 0){
    //     printf("hey \n");
    //     hx[1] = mu_e - C->m_rho * C->phi_n(0, f);
	// }
    // else{
    //     printf("hey2 \n");
    //     hx[1] = Q;
    // }
    hx[1] = Q;

    //if (rc2 > 0){
    hx[2] = (C->m_rho * C->phi_n(0, f) - mu_e) * rc2;
    // printf("m_rho^* = %.6f, mu_e = %.6f, hx[2] = %.6f \n", C->m_rho * C->phi_n(0, f), mu_e, hx[2]);
    //}
}

void wrap_fun_free(double n_tot,
                   double * init, int initN, 
                   set_const * C, 
                   double * out, int dimOut)
{
    double f[1] = {0.};
    calc::fun_n_eq_params p = {C, n_tot, f, 1, 0.};
    fun_n_eq_free(init, out, initN, initN, &p);
}

void step_free(double n, double * init, int initN, double * f_init, int dimF_init, double * out, int dim_Out, int iter, set_const* C) {
	double opts[5];
	bool debug = 1;
	calc::fun_n_eq_params p = {C, n, f_init, dimF_init, 0.0};
	int m = initN;
	double * x = new double[m];
	double * lb = new double[m];
    double * ub = new double[m];
    double * scale = new double[m];
	double * fun = new double[m];
	double info[LM_INFO_SZ];
	//double x[3] = {v.n[0], v.n[1], v.f};

	for (int i = 0; i < m; i++){
		x[i] = init[i];
		lb[i] = 0.0;
        ub[i] = n;
        scale[i] = 1.;
//		if (i > 2) lb[i] = -100500.0;
	}

    scale[0] = 1;
    //Upper limit for the f variable
    ub[0] = 1.;
    ub[2] = 10;


	opts[0]= LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-55; opts[3]=1E-15;
		opts[4]= -1e-5;
	dlevmar_bc_dif(fun_n_eq_free, x, NULL, m, m, lb, NULL, NULL, iter, opts, info, NULL, NULL, &p);
//	dlevmar_dif(calc::fun_n_eq, x, NULL, m, m, 2000, opts, NULL, NULL, NULL, &p);

	if (debug) {
		printf("info: ");
		for (int i = 0; i < LM_INFO_SZ; i++){
			printf(" %i : %f ", i, info[i]);
		}
		printf("\n");

		printf("n = %f, f = %e", n, x[0]);
		printf(",n_p = %e ", x[1]);
		printf(",rc2 = %e ", x[2]);
		printf("\n");

		fun_n_eq_free(x, fun, m, m, &p);
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
