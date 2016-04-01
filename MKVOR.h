/*
 * MKVOR.h
 *
 *  Created on: Mar 30, 2016
 *      Author: const
 */

#ifndef MKVOR_H_
#define MKVOR_H_

#include "setconst.h"


class MKVOR: virtual public set_const {
public:
	MKVOR() : set_const(){
		this->z = 0.65;
	};

	double phi_n(int i, double f){
		double res = 1-f;
			if (i > 1){
				return 1 - f;
			}
			else{
				return 1 - f;
			}
	}

	double eta_s(double f){
		double c1 = c - 8./9. * Cs *  b*b;
		return 1./(1 - 2./3. *Cs * b * f - 0.5 * Cs * c1 * f*f +
				d * f*f*f /3.);
	};

	double eta_o(double f){
		double res = pow((1 + z * f0)/(1 + z*f), alpha);
		return res + .5 * a_om * (1 + tanh(b_om * (f - f_om)));
	}

	double eta_r(double f){
		double res = a_r0 + a_r1 * f + a_r2 * f * f/(1 + a_r3 * f * f) +
				beta * exp(-gamma * (pow(f - f_r, 2)*(1 + e_r * pow(f - f0, 2)))/
						(1 + d_r * (f - f0) + e_r * pow(f - f0, 2)));
		return res;
	}

	double eta_p(double f){
		if (phi_kind == 1){
			return pow(1 - f, 2);
		}
		else{
			return 1.;
		}
	}

	double U(double f){
		return 0.;
	}

	virtual ~MKVOR();

	int phi_kind;

	double alpha;
	double a_om;
	double b_om;
	double f_om;
	double z;

	double d;

	double a_r0;
	double a_r1;
	double a_r2;
	double a_r3;
	double beta;
	double gamma;
	double f_r;
	double e_r;
	double d_r;
};

#endif /* MKVOR_H_ */
