/*
 * KVORmod.h
 *
 *  Created on: 08 июля 2014 г.
 *      Author: const
 */

#ifndef KVORMOD_H_
#define KVORMOD_H_
#include "KVOR.h"


class KVOR_mod: public KVOR {
public:
	KVOR_mod();
	virtual ~KVOR_mod();
	double omega_f;
	double omega_a;
	double d;
	double e;
	double rho_f;
	double rho_a;
	double beta;
	double alpha;
	double gamma;
	double phi_f;
	double phi_a;
	double phi_gamma;
	double phi_z;
	double omega_c;
	virtual double eta_o(double);
	virtual double eta_r(double);
	virtual double eta_s(double);
	virtual double eta_p(double);
	virtual double phi_n(int,double);
	virtual double U(double);
};

#endif /* KVORMOD_H_ */
