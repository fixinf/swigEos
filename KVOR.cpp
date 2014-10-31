/*
 * KVOR.cpp
 *
 *  Created on: 09 июня 2014 г.
 *      Author: const
 */

#include "KVOR.h"

#include <cmath>
#include <cstdio>

#include "constants.h"

KVOR::KVOR() : set_const() {

	this->f0 = 0.195;
	this->z = 0.65;
}

KVOR::~KVOR() {

}

double KVOR::eta_o(double f){
//	printf("%f \n", f);
	return (1 + z*f0)/(1 + z*f);
}

double KVOR::eta_s(double f){
	return 1.0;
}

double KVOR::eta_r(double f){
	return this->eta_o(f)/(this->eta_o(f) +
		4*(this->Co/this->Cr)*(this->eta_o(f)-1.0));
}

double KVOR::phi_n(int i, double f){
	return 1.0 - f;
}

double KVOR::eta_p(double f) {
	return 1.0;
}

double KVOR::U(double f){
	return pow(M[0],4)*(b * pow(f, 3)/3 + c * pow(f,4)/4 );
}
