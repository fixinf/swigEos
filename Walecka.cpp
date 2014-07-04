/*
 * Walecka.cpp
 *
 *  Created on: 10 июня 2014 г.
 *      Author: const
 */

#include "Walecka.h"
#include <cmath>
#include "constants.h"
Walecka::Walecka() : set_const() {
	// TODO Auto-generated constructor stub
	Co = 54.6041;
	Cs = 164.462;
	Cr = 121.69;
	b = 0.0202832;
	c = 0.0471633;
}

Walecka::~Walecka() {
	// TODO Auto-generated destructor stub
}

double Walecka::eta_o(double f){
	return 1.0;
}

double Walecka::eta_s(double f){
	return 1.0;
}

double Walecka::eta_r(double f){
	return 1.0;
}

double Walecka::phi_n(double f){
	return 1.0 - f;
}

double Walecka::U(double f){
	return pow(m_n,4)*(b * pow(f, 3)/3 + c * pow(f,4)/4 );
}
