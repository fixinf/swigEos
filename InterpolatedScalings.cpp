/*
 * InterpolatedScalings.cpp
 *
 *  Created on: 26 нояб. 2014 г.
 *      Author: const
 */

#include "InterpolatedScalings.h"

InterpolatedScalings::InterpolatedScalings() : set_const() {
	this->splineU = 0;
	this->accU = 0;
	this->splineO = 0;
	this->accO = 0;
	this->splineR = 0;
	this->accR = 0;
	this->splineS = 0;
	this->accS = 0;
	this->rho_akima = 0;
}

InterpolatedScalings::~InterpolatedScalings() {

}

double InterpolatedScalings::U(double f) {
	if (debug) printf("U: f = %f ", f);
	double res = gsl_spline_eval(splineU, f, accU);
	if (debug) printf(" res = %f \n", res);
	return res;
}

double InterpolatedScalings::phi_n(int int1, double f) {
	return 1-f;
}

double InterpolatedScalings::eta_s(double f) {
	if (debug) printf("eta_o: f = %f ", f);
	double res = gsl_spline_eval(splineS, f, accS);
	if (debug) printf(" res = %f \n", res);
	return res;
}

double InterpolatedScalings::eta_o(double f) {
	if (debug) printf("eta_o: f = %f ", f);
	double res = gsl_spline_eval(splineO, f, accO);
	if (debug) printf(" res = %f \n", res);
	return res;
}

double InterpolatedScalings::eta_r(double f) {
	if (debug) printf("eta_r: f = %f ", f);
	double res = gsl_spline_eval(splineR, f, accR);
	if (debug) printf(" res = %f \n", res);
	return res;
}

double InterpolatedScalings::eta_p(double double1) {
	return 1.;
}

void InterpolatedScalings::set_U(double* f_in, int dimF_in, double* y_in,
		int dimY_in) {
	if (this->splineU != 0){
//		printf("U spline is already set!");
//		return;
		delete splineU;
	}

	if (dimF_in != dimY_in){
		printf("F and Y shapes must match!");
		return;
	}

	this->accU = gsl_interp_accel_alloc();
	this->splineU = gsl_spline_alloc(gsl_interp_cspline, dimF_in);
	gsl_spline_init(this->splineU, f_in, y_in, dimF_in);
}

void InterpolatedScalings::set_eta_s(double* f_in, int dimF_in, double* y_in,
		int dimY_in) {
	if (this->splineS != 0){
//		printf("U spline is already set!");
//		return;
		delete splineS;
	}

	if (dimF_in != dimY_in){
		printf("F and Y shapes must match!");
		return;
	}

	this->accS = gsl_interp_accel_alloc();
	this->splineS = gsl_spline_alloc(gsl_interp_cspline, dimF_in);
	gsl_spline_init(this->splineS, f_in, y_in, dimF_in);
}

void InterpolatedScalings::set_eta_o(double* f_in, int dimF_in, double* y_in,
		int dimY_in) {
	if (this->splineO != 0){
		printf("Omega spline is already set!");
		return;
	}

	if (dimF_in != dimY_in){
		printf("F and Y shapes must match!");
		return;
	}

	this->accO = gsl_interp_accel_alloc();
	this->splineO = gsl_spline_alloc(gsl_interp_cspline, dimF_in);
	gsl_spline_init(this->splineO, f_in, y_in, dimF_in);
}

void InterpolatedScalings::set_eta_r(double* f_in, int dimF_in, double* y_in,
		int dimY_in) {
	if (this->splineR != 0){
//		printf("Omega spline is already set!");
//		return;
		delete this->splineR;
		delete this->accR;
	}

	if (dimF_in != dimY_in){
		printf("F and Y shapes must match!");
		return;
	}

	this->accR = gsl_interp_accel_alloc();
	if (rho_akima){
		this->splineR = gsl_spline_alloc(gsl_interp_akima, dimF_in);
	}
	else{
		this->splineR = gsl_spline_alloc(gsl_interp_cspline, dimF_in);
	}

	gsl_spline_init(this->splineR, f_in, y_in, dimF_in);
}
