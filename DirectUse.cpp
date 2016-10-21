//
// Created by const on 8/4/15.
//
#include "KVOR.h"
#include "TOV.h"
#include "KVDriver.h"
#include "SCDelta.h"
#include "eos.h"
#include "aux.h"
#include <iostream>
#include <sstream>
using namespace std;

KVDriver * parseEos(string eosfile){
    ifstream f(eosfile);
    int num_lines = 0;
    string line;
    while (getline(f, line))
        if (line[0] != '#')
            num_lines++;
    f.clear();
    f.seekg(0, ios::beg);

    double * N = new double[num_lines];
    double * E = new double[num_lines];
    double * P = new double[num_lines];
    int i = 0;
    while (getline(f, line)){
        if (line[0] != '#') {
//            cout << line << endl;
            stringstream line_ss(line);
            line_ss >> N[i] >> E[i] >> P[i];
            i++;
        }
    }
    for (int i = 1; i < num_lines; i++){
        if (E[i] < E[i-1] || N[i] < N[i-1] || (P[i] < P[i-1])){
            cout << i << " " << N[i] << " " << E[i] << " " << P[i] << endl;
            cout << " " << " " << N[i-1] << " " << E[i-1] << " " << P[i-1] << endl;
        }
    }
    return new KVDriver(E, num_lines, P, num_lines, N, num_lines);
}

void testMass(){
    KVOR * m = new KVOR();
    cout << m->b << endl;
    double result[3];
    KVDriver * dr = parseEos("/home/const/workspace2/swigEosWrapper/Cuts/EPmpi4_KVOR06.dat");
    star_crust2(1., result, 3, dr, 1e-10);
}

void testDelta(){
    set_const * m = new KVOR_cut();
    set_const * m2 = new KVORcut_d();
    m->Hyper = 1;
    m->phi_meson = 0;
    m->sigma_kind = 0;
    m->SetHyperConstants(2);

    m2->Hyper = 1;
    m2->phi_meson = 0;
    m2->sigma_kind = 0;
//    m2->SetHyperConstants(2);
    double n= 2.833225;
    double dn = .8;

    int dim = 8;
    double init[dim] = {6.748172e-01,
                      0.000000e+00,
                      0.000000e+00,
                      0.000000e+00,
                      0.000000e+00,
                      5.304112e-01,
                      0.000000e+00,
                        0.
    };
    double finit[1] = {.52};
    double out[dim];
    int iter = 100;
    stepE(n+dn, init, dim, finit, 1, out, dim, 100, m);
    for (int i = 0; i < dim; i++){
        printf("%f  ", out[i]);
    }
    printf("\n");

    stepE(n+dn, init, dim, finit, 1, out, dim, 100, m2);
    for (int i = 0; i < dim; i++){
        printf("%f  ", out[i]);
    }
    printf("\n");

    double res[1];
    int dimF = 7;
    double n_f[8] = {
            n - 6.748172e-01 - 5.304112e-01,
            6.748172e-01,
            0.000000e+00,
            0.000000e+00,
            0.000000e+00,
            0.000000e+00,
            5.304112e-01,
            0.000000e+00,
    };
    f_eq(n_f, dimF, finit, 1, res, 1, m);
    cout << res[0] << endl;
    f_eq(n_f, dimF, finit, 1, res, 1, m2);
    cout << res[0] << endl;

    for (int i = 0; i < 8; i++){
        cout << m->X_s[i] <<  " " << m2->X_s[i] << endl;
    }

    for (int i = 0; i < 8; i++){
        cout << m->X_o[i] <<  " " << m2->X_o[i] << endl;
    }

    for (int i = 0; i < 8; i++){
        cout << m->X_r[i] <<  " " << m2->X_r[i] << endl;
    }

    for (int i = 0; i < 8; i++){
        cout << m->X_p[i] <<  " " << m2->X_p[i] << endl;
    }
}

void testDeltaWalecka(){
    Walecka_d * C = new Walecka_d();
    double n0 = 2*C->n0;
    cout << C->n0 << endl;
    double Eparts[9];
    cout << C->Co << endl;
    C->X_s[8] = 1.4;
    C->X_s[9] = 1.4;
    C->X_s[10] = 1.4;
    C->X_s[11] = 1.4;
    for (int i = 0; i < 12; i++){
        printf("%f ", C->X_s[i]);
    }
    printf("\n");

    for (int i = 0; i < 12; i++){
        printf("%f ", C->X_o[i]);
    }
    printf("\n");

    for (int i = 0; i < 12; i++){
        printf("%f ", C->M[i]);
    }
    printf("\n");

    for (int i = 0; i < 12; i++){
        printf("%f ", C->S[i]);
    }
    printf("\n");

    for (int i = 0; i < 12; i++){
        printf("%f ", C->T[i]);
    }
    printf("\n");
//    C->X_s[12] = 1.4;
//    C->X_s[13] = 1.4;
    double n_in[13] = {.6, n0/4, n0/4, 0., 0., 0., 0., 0., 0.,
                        n0/8, n0/8, n0/8, n0/8};
    for (int i = 0; i < 12; i++){
        printf("%f ", n_in[i]);
    }
    printf("\n");
    double E = _E(n_in, 13, C, Eparts, 9);
    cout << E << endl;
    for (int i = 0; i < 9; i++){
        cout << Eparts[i] << endl;
    }
    cout << endl;
    double K = kineticInt(n0/8, C->M[10]*C->phi_n(10, C->M[0]/C->M[10]*C->X_s[10]*.6), 4);
    cout << n0/8 << endl;
    cout << K << endl;
    cout << C->M[0]/C->M[10]*C->X_s[10]*.6 << endl;
}

void testStepRho(){
	set_const * C = new KVOR();
	double n = C->n0;
	double init[1] = {0.};
	double finit[1] = {0.};
	double res[1];

	stepE_rho(n, init, 2, finit, 1, res, 1, 30, 0., C);
	cout << res[0] << "  " << res[1] << endl;

}

void testStepRhoInv(){
	set_const * C = new KVOR();
	C->SetHyperConstants(2);
	printf("%.6f\n", C->M[0]);
	double f_init = C->f0;
	double init[3] = {C->n0, 0.1, 0.5};
	double res[3];
	for (int j = 0; j < 10; j++){
	stepE_rho_f(C->f0 + 0.01*j, init, 3, res, 3, 30, C);
	for (int i = 0; i < 3; i++){
		printf("res[%i] = %.6f, ", i, res[i]);
	}
	printf("\n");
	}
}

int main(){
    testStepRhoInv();
}
