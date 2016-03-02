//
// Created by const on 11/19/15.
//

#include "RhoCond.h"
//#include "eos.h"
#include <gsl/gsl_deriv.h>

namespace RhoCond{
    struct mu_params{
        int i;
        double * n;
        int dimN;
        set_const * C;
        double rho_0, rho_c;
        double mu_c;
    };

    double mu_func(double x, void * params){
        mu_params * p = (mu_params *) params;
        double * n = new double[p->dimN];
        for (int i = 0; i < p->dimN; i++){
            n[i] = p->n[i];
        }
        n[p->i] = x;
        return E_rho(n, p->dimN, p->rho_0, p->rho_c, p->mu_c, p->C);
    }

    double mu_deriv(int i, double * n, int dimN, double rho_0, double rho_c, double mu_c, set_const * C){
        gsl_function f;
        f.function = &mu_func;
        mu_params params{i, n, dimN, C, rho_0, rho_c, mu_c};
        f.params = &params;
        double res, abserr;
        gsl_deriv_central(&f, n[i], 1e-8, &res, &abserr);
        return res;
    }

    double mu_sol(double * n, int dimN, double rho_0, double rho_c, set_const * C){

    }

    double E(double *n, int dimN, double rho_0, double rho_c, double mu_c, set_const *C, double *inplace,
                 int dim_inplace) {
        bool debug = 0;
        bool ret_parts = (inplace) && (dim_inplace == 9);
        if (debug) {
            printf("n = ");
            for (int i = 0; i < dimN; i++) {
                printf("%f ", n[i]);
            }
            printf("\n");
        }
        double f = n[0];
        int sc = 1 + C->sprime;
        double fp = 0;
        if (C->sprime) {
            fp = 0.;
        }

        double res = 0.;
        double part_s = pow(C->M[0], 4.0) * f * f * C->eta_s(f) / (2 * C->Cs);

        res += part_s;

        if (debug) {
            printf("res_f : %f\n", part_s);
        }
        if (ret_parts) {
            inplace[0] = part_s;
        }

        double part_U = C->U(f);
        res += part_U;
        if (ret_parts) {
            inplace[1] = part_U;
        }

        res += pow(C->M[0], 4.0) * fp * fp / (2 * C->Csp);

        double sum = 0;
        double sum_t3 = 0;
        double sum_p = 0;
        double meff_arg = 0;
        if (debug) {
            printf("res_Uf : %f \n", res);
        }
        double kin = 0.;
        for (int i = sc; i < dimN; ++i) {
            double xs = 0.;
            if (C->sigma_kind == 0) {
                xs = C->X_s[i - sc];
            }
            else {
                xs = C->Xs(i - sc, f);
            }
            meff_arg = xs * (C->M[0] / C->M[i - sc]) * f + C->X_sp[i - sc] * (C->M[0] / C->M[i - sc]) * fp;
            kin = kineticInt(n[i], (C->M)[i - sc] * C->phi_n(i - sc, meff_arg), 2 * C->S[i - sc] + 1);
            res += kin;
            if (i == 1 or i == 2) {
                if (ret_parts) {
                    inplace[i + 1] = kin;
                }
            }
            sum += n[i] * C->X_o[i - sc];
            sum_t3 += n[i] * (C->T)[i - sc] * C->X_r[i - sc];
            if (C->phi_meson) {
                sum_p += n[i] * C->X_p[i - sc];
            }
        }
        //omega
        double part_om = C->Co * sum * sum / (2.0 * C->M[0] * C->M[0] * C->eta_o(f));
        res += part_om;
        if (debug) {
            printf("res_om : %f \n", part_om);
        }
        if (ret_parts) {
            inplace[4] = part_om;
        }
        //rho
        double gr = sqrt(C->Cr / C->eta_r((f))) * (C->m_rho * (1 - f)) / C->M[0];
        double E_r = gr * sum_t3 * rho_0 - .5 * pow(rho_0 * C->m_rho * C->phi_n(0, f), 2);
        E_r -= pow(rho_c, 2) * (pow(gr * rho_0, 2) - pow(C->m_rho * C->phi_n(0, f), 2));
        //phi
        double part_phi = pow(m_o / m_p, 2.0) * C->Co * sum_p * sum_p / (2.0 * C->M[0] * C->M[0] * C->eta_p(f));
        res += part_phi;
        if (debug) {
            printf("res_phi : %f \n", part_phi);
        }

        if (ret_parts) {
            inplace[6] = part_phi;
        }

        return res + E_r;
    }
}