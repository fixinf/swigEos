/*
 * KVORmod.cpp
 *
 *  Created on: 08 июля 2014 г.
 *      Author: const
 */

#include "KVORmod.h"

#include <cmath>
#include <cstdio>

KVOR_mod::KVOR_mod() : KVOR() {
	// TODO Auto-generated constructor stub
	this->phi_a = 0.0;
	this->omega_a = 0.0;
	this->rho_a = 0.0;
	this->d = 0.0;
	this->e = 0.0;
	this->beta = 2.8;
	this->gamma = 2.0;
	this->alpha = 1.0;
	this->omega_f = 1.0;
	this->rho_f = 1.0;
	this->phi_f = 1.0;
}

double KVOR_mod::eta_o(double f){
	if (f < this->omega_f){
		return pow(KVOR::eta_o(f), alpha);
	}
	else{
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		return pow(KVOR::eta_o(f), alpha) + omega_a*pow(f - omega_f, 3);
	}
}

double KVOR_mod::eta_r(double f){
	double res = pow(1 + beta * f/(1 + beta*f0), gamma);
	if (f < this->rho_f){
		return res;
	}
	else{
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		return res + rho_a*pow(f - rho_f, 3);
	}
}

double KVOR_mod::eta_s(double f){
	double c1 = this->c - 8*Cs*pow(this->b, 2.0) / 9;
	return pow(1.0 - 2*Cs*b*f/3 - Cs*pow(f,2)*c1/2 + d*pow(f,3)/3 + e*pow(f,4)/4, -1);
}

double KVOR_mod::U(double f){
	return 0.0;
}

double KVOR_mod::phi_n(double f){
	double res = 1-f;
	if (f < this->phi_f){
		return res;
	}
	else{
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		return res + phi_a*pow(f - phi_f, 3);
	}
}

KVOR_mod::~KVOR_mod() {
	// TODO Auto-generated destructor stub
}

