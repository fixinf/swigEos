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
	// TODO Auto-generated constructor stub
	this->f0 = 0.195;
	this->z = 0.65;
}

KVOR::~KVOR() {
	// TODO Auto-generated destructor stub
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
		4*(this->Co/this->Cr)*(this->eta_o(f)-1.0));;
}

double KVOR::phi_n(double f){
	return 1.0 -f;
}

double KVOR::U(double f){
	return pow(m_n,4)*(b * pow(f, 3)/3 + c * pow(f,4)/4 );
}
