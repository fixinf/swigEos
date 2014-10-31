/*
 * KVDriver.cpp
 *
 *  Created on: 25 июня 2014 г.
 *      Author: const
 */

#include "KVDriver.h"

#include <cmath>
#include <cstdio>
//#include <fstream>


using namespace std;
KVDriver::KVDriver() {

}

KVDriver::KVDriver(double * E, int dimE, double * P, int dimP, double * n, int dimN){
	for (int i = 0; i < dimE; i++){
		this->E[i] = E[i];
		this->P[i] = P[i];
		this->n[i] = n[i];
	}
	this->count = dimN;
}

KVDriver::KVDriver(set_const* C, string fname) {
	this->C = C;
	this->fname = fname;
}

double KVDriver::PofE(double _E){
	double min = 1e45;
	int iMin;
	for (int i = 0; i < count; i++){
		if (abs(E[i] - _E) < min){
			min = abs(E[i] - _E);
			iMin = i;
		}
	}
	return P[iMin];
}

double KVDriver::EofN(double _n) {
	double min = 1e45;
	int iMin;
	for (int i = 0; i < count; i++){
		if (abs(n[i] - _n) < min){
			min = abs(n[i] - _n);
			iMin = i;
		}
	}
	return E[iMin];
}

double KVDriver::PofN(double _n) {
	double min = 1e45;
	int iMin;
	for (int i = 0; i < count; i++){
		if (abs(n[i] - _n) < min){
			min = abs(n[i] - _n);
			iMin = i;
		}
	}
	return P[iMin];
}

double KVDriver::EofP(double _P) {
	double min = 1e45;
	int iMin;
	for (int i = 0; i < count; i++){
		if (abs(P[i] - _P) < min){
			min = abs(P[i] - _P);
			iMin = i;
		}
	}
	return E[iMin];
}

double KVDriver::NofP(double _P){
	double min = 1e45;
	int iMin;
	for (int i = 0; i < count; i++){
		if (abs(P[i] - _P) < min){
			min = abs(P[i] - _P);
			iMin = i;
		}
	}
	return n[iMin];
}

double KVDriver::NofE(double _E){
	double min = 1e45;
	int iMin;
	for (int i = 0; i < count; i++){
		if (abs(E[i] - _E) < min){
			min = abs(E[i] - _E);
			iMin = i;
		}
	}
	return n[iMin];
}

int KVDriver::lookFor(double* src, int dim_src, double what) {
	double max_diff = 10e42;
	bool got = 0;
	int a = 0;
	int b = dim_src - 1;
	while (a != b){
		if (src[a] < what){
			if (src[b] > what){
				b = floor((b - a)/2);
			}
			else{
				b = a + b/2;
				if (a == ceil(b - a/2)){
					return a;
				}
				a = ceil((b - a)/2);
			}
		}
	}
	return a;
}

void KVDriver::set(double * E, int dimE, double * P, int dimP, double * n, int dimN){
	this->E = new double[dimE];
	this->P = new double[dimP];
	this->n = new double[dimN];
	for (int i = 0; i < dimE; i++){
//		printf("%f %f %f \n", E[i], P[i], n[i]);
		this->E[i] = E[i];
		this->P[i] = P[i];
		this->n[i] = n[i];
	}
	this->count = dimN;
}

KVDriver::~KVDriver() {

}

