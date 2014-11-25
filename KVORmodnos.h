/*
 * KVORmodnos.h
 *
 *  Created on: 29 сент. 2014 г.
 *      Author: const
 */

#ifndef KVORMODNOS_H_
#define KVORMODNOS_H_

class KVORmod_nos {
public:
	KVORmod_nos();
	virtual ~KVORmod_nos();
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

#endif /* KVORMODNOS_H_ */
