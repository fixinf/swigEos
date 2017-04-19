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
#include <gsl/gsl_integration.h>
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

double KVDriver::PofE(double E){
	return gsl_spline_eval(iPofE, E, accPofE);
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
		delete iPofE;
		delete accPofE;
		delete proper_N;
	}
	const gsl_interp_type * kind = gsl_interp_linear;

	count = dimN;
	this->E = E;
	this->P = P;
	this->n = n;
	accEofP = gsl_interp_accel_alloc();
	iEofP = gsl_spline_alloc(kind, dimN);
	gsl_spline_init(iEofP, P, E, dimN);

	accPofE = gsl_interp_accel_alloc();
	iPofE = gsl_spline_alloc(kind, dimN);
	gsl_spline_init(iPofE, E, P, dimN);

	// cout << "Set!" << endl;

	calculateN();

	accPofN = gsl_interp_accel_alloc();
	iPofN = gsl_spline_alloc(kind, dimN);
	gsl_spline_init(iPofN, proper_N, P, dimN);

	accEofN = gsl_interp_accel_alloc();
	iEofN = gsl_spline_alloc(kind, dimN);
	gsl_spline_init(iEofN, proper_N, E, dimN);

	accNofE = gsl_interp_accel_alloc();
	iNofE = gsl_spline_alloc(kind, dimN);
	gsl_spline_init(iNofE, E, proper_N, dimN);


	accNofP = gsl_interp_accel_alloc();
	iNofP = gsl_spline_alloc(kind, dimN);
	gsl_spline_init(iNofP, P, proper_N, dimN);
	this->isSet = 1;
}

double properN_func(double x, void * params){
	DriverBase * dr = (DriverBase*) params;
	return 1/(x + dr->PofE(x));
}

void KVDriver::calculateN(){
	gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (1000);
  
  double result, error;
  double expected = -4.0;
  double alpha = 1.0;

  gsl_function F;
  F.function = &properN_func;
  F.params = this;

// The following is stupid. but I don't want to intepolate
  int i1 = 0;
  for (int i = 0; i < count; i++){
	  if (this->n[i] > 0.75){
		  i1 = i;
		  break;
	  }
  }

  double E1 = E[i1];
  printf("%i, %.6f \n", i1, E1);
  this->proper_N = new double[count];
  proper_N[0] = 0.;
  for (int i = 1; i < count; i++){
	printf("%i \n", i);
	gsl_integration_qags (&F, E1, this->E[i], 1e-7, 1e-5, 1000,
                        w, &result, &error);
	proper_N[i] = this->n[i1] * exp(result);
  }
}

KVDriver::~KVDriver() {

}
