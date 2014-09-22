/*
 * Walecka.h
 *
 *  Created on: 10 июня 2014 г.
 *      Author: const
 */

#ifndef WALECKA_H_
#define WALECKA_H_

#include "setconst.h"

class Walecka: public set_const {
public:
	Walecka();
	virtual ~Walecka();
	virtual double eta_o(double);
	virtual double eta_r(double);
	virtual double eta_s(double);
	virtual double phi_n(int, double);
	virtual double U(double);
	virtual double eta_p(double);
};

#endif /* WALECKA_H_ */
