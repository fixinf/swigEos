/*
 * KVOR.h
 *
 *  Created on: 09 июня 2014 г.
 *      Author: const
 */

#ifndef KVOR_H_
#define KVOR_H_

#include "setconst.h"

class KVOR: public set_const {
public:
	KVOR();
	virtual ~KVOR();
	double z;
	virtual double eta_o(double);
	virtual double eta_r(double);
	virtual double eta_s(double);
	virtual double phi_n(int, double);
	virtual double eta_p(double);
	virtual double U(double);
};

class KVOR_MD: public KVOR{
public:
	KVOR_MD(): KVOR(){

	}
	double eta_o(double f){
		double res = KVOR::eta_o(f);
		if (f < f0){
			return 1 + om_prime*(f-f0) + b_om * pow(f - f0, 2)/pow(f0,2);
		}
		return res;
	}

	double eta_r(double f){
		double res = KVOR::eta_r(f);
		if (f < f0){
			return 1 + rho_prime*(f-f0) + b_rho * pow(f - f0, 2)/pow(f0,2);
		}
		return res;
	}

	double eta_s(double f){
		if (f < f0)	return 1 + b_sigma * pow(f-f0,2)/pow(f0, 2);
		else return 1.;
	}



	double b_sigma;
	double om_prime;
	double rho_prime;
	double b_om;
	double b_rho;
};

#endif /* KVOR_H_ */
