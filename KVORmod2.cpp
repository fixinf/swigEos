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
	omega_f_low = 1.;
	omega_a_low = 0.;
	rho_f_low = 0.;
	rho_a_low = 0.;
	rho_gamma_low = 3.;
	this->Delta = 0.;
	this->umode = 0;
	this->rho_val = 0.;
	rho_d = 0;
	rho_ld = 0;
}

KVOR_mod2::~KVOR_mod2() {

}

double KVOR_mod2::eta_o(double f){
//	f = (1.0-f0)*this->func(f - f0)/this->func(1-f0) + f0;
	double res = pow(KVOR::eta_o(f), alpha);
//	double res = pow(KVOR::eta_o(f), alpha) - 1;
//	res *= (tanh(1e5*pow(f - f0,3)) + 1)/2 ;
//	res += 1;

	if (omega_kind ==2){
		double s = 0.5*(1 + tanh(omega_b*(f - omega_f)));
		return res +  omega_a * s;
	}

	if (f < this->omega_f){
		if (f < this->omega_f_low)
		return res + omega_a_low * pow(f - omega_f_low, 3.);
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

double KVOR_mod2::f(double x){
	if (x > 0){
		return exp(-1./x);
	}
	else{
		return 0;
	}
}

double KVOR_mod2::eta_r(double f){
	double res = pow((1 + beta * pow(f, rho_power))/(1 + beta*pow(f0, rho_power)), gamma);
	double f1, fs, s;

	if (rho_kind == 6){
		res = rho_a + rho_a0 * f/(1 + rho_a4 * f) + rho_a1*f*f/(1 + rho_a2*f*f + rho_a3*f*f*f) +
				beta*exp(-rho_sat_val*pow(pow(f-rho_f,2)/(1+rho_d*pow(f - rho_width_f,1)/(1 + rho_e * pow(f - rho_width_f,2))), rho_width_power)) +
				beta1*exp(-this->c1 * pow(f - rho_f, 2));

		if (f > this->rho_sat_f1){
			res += rho + drho*(f - rho_sat_f1) + 0.5*d2rho*pow(f - rho_sat_f1,2 )
			+ rho_tan_a* tanh(rho_tan_b*pow(f/rho_sat_f1 - 1, 4));
		}

		if (f > rho_sat_f2){
			res += rho_tan_a*pow(f/rho_sat_f2 - 1, 3) + rho_tan_b*pow(f/rho_sat_f2 - 1, 5)
			+ rho_tan_c*pow(f/rho_sat_f2 - 1, 7);
		}

//		if (f < f0){
//			res += rho_a_low*pow(f/f0-1,3)/pow(f0,3)
//				+ rho_b_low*pow(f/f0-1,5)/pow(f0,5);
//		}
//
//		if (rho_ld){
//			if (f < f0)
//			return 1 + drho*(f-f0) + d2rho*pow(f-f0,2)/2
//					+ rho_a_low*pow(f-f0,3)/pow(f0,3)
//					+ rho_b_low*pow(f-f0,5)/pow(f0,5);
//		}

		return res;
	}

	if (rho_kind == 7){
		s = 0.5*(1 + tanh(rho_sat_a*(f - rho_sat_f1)));
		res = rho_a + gamma*s + beta*exp(-rho_sat_val*pow(f-rho_f,2));
//		if (f < this->rho_f_low){
//			res = res + rho_a_low*pow(f - rho_f_low, 3.);
//		}
		if (f < f0){
			return 1 + drho*(f-f0) + d2rho*pow(f-f0,2)/2
					+ rho_a_low*pow(f-f0,3)/pow(f0,3)
					+ rho_b_low*pow(f-f0,5)/pow(f0,5);
		}
		return res;
	}

	if (f < this->rho_f){
		if (f < this->rho_f_low){
//			double R = 2 * pow(this->M[0], 2) * Delta + Cs - Co/this->eta_o(0.);
//			R = 4*Cr/R;
//			double etar0 = pow((1 + beta * pow(0., rho_power))/(1 + beta*pow(f0, rho_power)), gamma);
//			double mult = (1 - (1 - R) / etar0 * pow((1 - f/rho_f_low), rho_a_low));
//			printf("f = %f, mult = %f \n", f, mult);
//			return res * mult;
//			res = res + rho_a_low*pow(f - rho_f_low, 3.);
		}

		res = res;
	}
	else{
		switch (rho_kind){
		case 0 :
			res = res + rho_a*pow(f - rho_f, 3);
			break;
		case 1 :
			res =  res / cosh(rho_a * pow(f - rho_f,2));
			break;
		case 2:
			res = res / cosh(rho_a * pow(f - rho_f,2));
//			double s = 0.5*(tanh(rho_sat_a * pow(f-rho_sat_f, 1)) + 1.);
			f1 = f - rho_sat_f1;
//			s = this->f(f1*rho_sat_a)/(this->f(f1*rho_sat_a) +
//					this->f((rho_sat_f2-f)*rho_sat_a));
			s = this->s_exp(f);
//			printf("f = %f S = %f \n", f, s);
			return res*(1-s) + rho_sat_val * s;
			fs = f*(1-s) + rho_sat_f2*s;
//			printf("f = %f fs=%f \n", f, fs);
			res = pow((1 + beta * pow(fs, rho_power))/(1 + beta*pow(f0, rho_power)), gamma);
			res /= cosh(rho_a * pow(fs - rho_f,2));
			return res;
		case 3:
			res = res / cosh(rho_a * pow(f - rho_f,2));
////			double s = 0.5*(tanh(rho_sat_a * pow(f-rho_sat_f, 1)) + 1.);
//			f1 = f - rho_sat_f1;
////			double s = this->f(f1*rho_sat_a)/(this->f(f1*rho_sat_a) +
////					this->f((rho_sat_f2-f)*rho_sat_a));
//			s = this->s(f);
////			printf("f = %f S = %f \n", f, s);
////			res = res*(1-s) + rho_sat_val * s;
//			break;
//			fs = f*(1-s) + rho_sat_f2*s;
////			printf("f = %f fs=%f \n", f, fs);
////			res = pow((1 + beta * pow(fs, rho_power))/(1 + beta*pow(f0, rho_power)), gamma);
////			res /= cosh(rho_a * pow(fs - rho_f,2));
////			return res;
			break;
		case 4:
			res = res *(rho_sat_val + (1-rho_sat_val) / cosh(rho_a * pow(f - rho_f,2)));
			break;
		}
	}
	if (rho_kind == 3){
//		s = this->s(f);
		res = res + rho_sat_val * 0.5*(1+tanh(rho_sat_a*pow(f - rho_sat_f1, 3)));
	}
	return res;
}

double KVOR_mod2::s(double f){
	double f1 = f - rho_sat_f1;
	double s = this->f(f1*rho_sat_a)/(this->f(f1*rho_sat_a) +
			this->f((rho_sat_f2-f)*rho_sat_a));
//	return s;
	return 0.5*(1 + tanh(rho_sat_a*(f - rho_sat_f1)));
}

double KVOR_mod2::s_exp(double f){
	double f1 = f - rho_sat_f1;
	double s = this->f(f1*rho_sat_a)/(this->f(f1*rho_sat_a) +
			this->f((rho_sat_f2-f)*rho_sat_a));
	return s;
//	return 0.5*(1 + tanh(rho_sat_a*(f - rho_sat_f1)));
}

double KVOR_mod2::eta_s(double f){
	if (this->umode == 0){
		double c1 = this->c - 8*Cs*pow(this->b, 2.0) / 9;
		return pow(1.0 - 2*Cs*b*f/3 - Cs*pow(f,2)*c1/2 + d*pow(f,3)/3 + e*pow(f,4)/4, -1);
	}
	else{
		return 1.;
	}
}

double KVOR_mod2::U(double f){
	if (this->umode == 0){
		return 0.0;
	}
	else{
		return pow(this->M[0], 4.)*(b*pow(f,3)/3. + c*pow(f, 4)/4);
	}
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
