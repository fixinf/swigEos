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
#include <gsl/gsl_interp.h>
#include <gsl/gsl_errno.h>
#include <iostream>


using namespace std;
KVDriver::KVDriver() {

}

KVDriver::KVDriver(double * E, int dimE, double * P, int dimP, double * n, int dimN){
	if (dimE != dimP or dimE!=dimN or dimP != dimN){
		throw 20;
	}
	this->isSet = 0;
	this->set(E, dimE, P, dimP, n, dimN);
}

double KVDriver::NofP(double P){
	return gsl_spline_eval(iNofP, P, accNofP);
}

double KVDriver::NofE(double E){
	return gsl_spline_eval(iNofE, E, accNofE);
}

double KVDriver::EofP(double P){
	return gsl_spline_eval(iEofP, P, accEofP);
}

void interp_error_handler(const char * reason,
      const char * file,
	  int line,
	  int gsl_errno){
	cout << "GSL error, doing nothing yet!" << endl;
}

double KVDriver::EofN(double N){
	gsl_error_handler_t * old_handler;
	old_handler = gsl_set_error_handler_off();
	try{
		return gsl_spline_eval(iEofN, N, accEofN);
	}
	catch (int e){
		cout << "Error caught: " << e << endl;
		return 0.;
	}
	gsl_set_error_handler(old_handler);
}

double KVDriver::PofN(double N){
	return gsl_spline_eval(iPofN, N, accPofN);
}

double KVDriver::dEofN(double N){
	return gsl_spline_eval_deriv(iEofN, N, accEofN);
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
	cout << "KVDriver::set" << endl;
	if (this->isSet) {
		delete iEofN;
		delete accEofN;
		delete iNofE;
		delete accNofE;
		delete iNofP;
		delete accNofP;
		delete iEofP;
		delete accEofP;
		delete iPofN;
		delete accPofN;
	}
	accEofN = gsl_interp_accel_alloc();
	iEofN = gsl_spline_alloc(gsl_interp_cspline, dimN);
	gsl_spline_init(iEofN, n, E, dimN);
	cout << "E spline init" << endl;
	accNofE = gsl_interp_accel_alloc();
	iNofE = gsl_spline_alloc(gsl_interp_cspline, dimN);
	gsl_spline_init(iNofE, E, n, dimN);


	accNofP = gsl_interp_accel_alloc();
	iNofP = gsl_spline_alloc(gsl_interp_cspline, dimN);
	gsl_spline_init(iNofP, P, n, dimN);


	accEofP = gsl_interp_accel_alloc();
	iEofP = gsl_spline_alloc(gsl_interp_cspline, dimN);
	gsl_spline_init(iEofP, P, E, dimN);


	accPofN = gsl_interp_accel_alloc();
	iPofN = gsl_spline_alloc(gsl_interp_cspline, dimN);
	gsl_spline_init(iPofN, n, P, dimN);
	this->isSet = 1;
	cout << "Set!" << endl;
}

KVDriver::~KVDriver() {

}

