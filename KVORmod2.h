/*
 * KVORmod2.h
 *
 *  Created on: 19 сент. 2014 г.
 *      Author: const
 */

#ifndef KVORMOD2_H_
#define KVORMOD2_H_

#include "KVORmod.h"

class KVOR_mod2: public KVOR_mod {
public:
	KVOR_mod2();
	virtual ~KVOR_mod2();
	double omega_f;
	double omega_a;
	double d;
	double e;
	double rho_f;
	double rho_a;
	double rho_val;
	double beta;
	double alpha;
	double gamma;
	double beta2;
	double gamma2;
	double phi_f;
	double phi_a;
	double phi_gamma;
	double phi_z;
	double rho_sat_val;
	double rho_sat_a;
	double rho_sat_f1;
	double rho_sat_f2;
	double omega_f_low;
	double rho_f_low;
	double omega_a_low;
	double rho_a_low;
	double rho_a1;
	double rho_a2;
	double rho_a0;
	double rho_a3;
	double rho_gamma_low;
	double rho_d;
	double rho_e;
	double rho_width_f;
	double rho_width_power;
	double rho;
	double drho;
	double d2rho;
	double rho_tan_a;
	double rho_tan_b;
	double rho_tan_c;
	double beta1;
	double c1;
	double dom;
	double d2om;
	double rho_b_low;
	double rho_a4;
	double Delta;
	double omega_c;
	int rho_kind;
	double rho_power;
	int omega_kind;
	double omega_b;
	bool rho_ld;
	virtual double eta_o(double);
	virtual double eta_r(double);
	virtual double eta_s(double);
	virtual double eta_p(double);
	virtual double phi_n(int,double);
	virtual double U(double);
	double f(double x);
	double s(double f);
	double s_exp(double f);
	int umode;
};

class KVOR_cut: public KVOR_mod2{
public:
	KVOR_cut(): KVOR_mod2(){

	}

	double eta_r(double f){
		double res = KVOR::eta_o(f)/(KVOR::eta_o(f) +
			4*(this->Co/this->Cr)*(KVOR::eta_o(f)-1.0));

//		if (f > this->rho_f){
////			printf("!");
//			res /= cosh(this->rho_a*pow(f-rho_f, 2.));
//		}

		if (f > this->rho_f){
			res += rho_a * pow(f-rho_f, 3.);
		}

		return res;
	}

	double eta_s(double f){
		return 1.;
	}

	double U(double f){
		return KVOR::U(f);
	}
};

//class Improved : KVOR_mod2{
//public:
//	double fRhoSat;
//	double fASat;
//	Improved() : KVOR_mod2(){
//		fRhoSat = 1.;
//	}
//
//	double eta_r(double f){
//		double res = KVOR_mod2::eta_r(f);
//
//	}
//
//};

class ImprovedLDParam: public KVOR_mod2{
public:
	ImprovedLDParam() : KVOR_mod2(){
		ar = 0.;
		br = 0;
		cr = 0;
		a_sigma = 0;
		power_0 = 0;
		power_1 = 0;
		f_stop = -1.;
	}

	double ar, br, cr;
	double a_sigma;
	double f_stop_denom;
	double power_0, power_1;
	double f_stop;

	double eta_r(double f){
		double res = KVOR_mod2::eta_r(f);
		if (f < f_stop){
			res -= ar*pow(f/f_stop_denom, power_0)*pow(1 - f/f_stop_denom, power_1);
			res += br*pow(f/f_stop_denom, power_0+1)*pow(1 - f/f_stop_denom, power_1);
			res += cr*pow(f/f_stop_denom, power_0)*pow(1 - f/f_stop_denom, power_1+1);
		}
		return res;
	}

	double eta_s(double f){
		double res = KVOR_mod2::eta_s(f);
		if (f < f_stop){
			res += a_sigma * pow(1 - f/f_stop, 4);
		}
		return res;
	}

};


#endif /* KVORMOD2_H_ */
