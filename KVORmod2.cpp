/*
 * KVORmod2.cpp
 *
 *  Created on: 19 сент. 2014 г.
 *      Author: const
 */

#include "KVORmod2.h"

#include <gsl/gsl_odeiv2.h>
#include <cstdio>
KVOR_mod2::KVOR_mod2() : KVOR_mod() {
	rho_kind = 0;
	omega_kind = 0;
	e = 0.;
	d = 0.;
	phi_a = 0.;
	phi_gamma = 0.0;
	phi_z = 0.0;
	beta = 1.0;
	gamma = 2.0;
	alpha = 1.0;
	omega_c = 0.;
	omega_a = 0.;
	omega_f = 1.;
	rho_a = 0.;
	rho_f = 1.;
	phi_f = 1.;
	this->exp_alpha = 0;
	this->Csp = 1.;
	this->rho_power = 2.;
}

KVOR_mod2::~KVOR_mod2() {

}

double KVOR_mod2::eta_o(double f){
//	f = (1.0-f0)*this->func(f - f0)/this->func(1-f0) + f0;
	double res = pow(KVOR::eta_o(f), alpha);
//	double res = pow(KVOR::eta_o(f), alpha) - 1;
//	res *= (tanh(1e5*pow(f - f0,3)) + 1)/2 ;
//	res += 1;
	if (f < this->omega_f){
		return res;
	}
	else{
//		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		switch (omega_kind){
		case 0 :
			return res + omega_a*pow(f - omega_f, 3) + omega_c*pow(f - omega_f, 11);
			break;
		case 1:
			return res / cosh(omega_a * pow(f - omega_f, 2));
			break;
		}
	}
}

double KVOR_mod2::eta_r(double f){
	double res = pow((1 + beta * pow(f, rho_power))/(1 + beta*pow(f0, rho_power)), gamma);
//	printf("%i \n", rho_kind);
	if (f < this->rho_f){
		return res;
	}
	else{
		switch (rho_kind){
		case 0 :
			return res + rho_a*pow(f - rho_f, 3);
			break;
		case 1 :
			return res / cosh(rho_a * pow(f - rho_f,2));
			break;
		}
	}
}

double KVOR_mod2::eta_s(double f){
	double c1 = this->c - 8*Cs*pow(this->b, 2.0) / 9;
	return pow(1.0 - 2*Cs*b*f/3 - Cs*pow(f,2)*c1/2 + d*pow(f,3)/3 + e*pow(f,4)/4, -1);
}

double KVOR_mod2::U(double f){
	return 0.0;
}

double KVOR_mod2::phi_n(int i, double f){
	double res = 1-f;
	if (i > 1){
		return 1 - f;
	}
	if (f < this->phi_f){
		return res;
	}
	else{
//		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
//		return res + phi_a*pow(f - phi_f, 3)/cosh(f - phi_f);
		return res + phi_a*pow(f - phi_f, 3);
	}
}

double KVOR_mod2::eta_p(double f){
	return pow((1 + phi_z*this->f0)/(1 + phi_z*f), this->phi_gamma);
}
