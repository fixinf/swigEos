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

#endif /* KVOR_H_ */
