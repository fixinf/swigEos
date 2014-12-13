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
	this->accU = 0;
}

InterpolatedScalings::~InterpolatedScalings() {

}

double InterpolatedScalings::U(double f) {
//	printf("f = %f \n", f);
	return gsl_spline_eval(splineU, f, accU);
}

double InterpolatedScalings::phi_n(int int1, double f) {
	return 1-f;
}

double InterpolatedScalings::eta_s(double double1) {
	return 1.;
}

double InterpolatedScalings::eta_o(double f) {
	printf("[eta_0] f = %f \n", f);
	double res = gsl_spline_eval(splineO, f, accO);
	printf("[eta_o] res = %f \n", res);
	return res;
}

double InterpolatedScalings::eta_r(double double1) {
	return 1.;
}

double InterpolatedScalings::eta_p(double double1) {
	return 1.;
}

void InterpolatedScalings::set_eta_s(double* f_in, int dimF_in, double* y_in,
		int dimY_in) {
	if (this->splineU != 0){
		printf("U spline is already set!");
		return;
	}

	if (dimF_in != dimY_in){
		printf("F and Y shapes must match!");
		return;
	}

	this->accU = gsl_interp_accel_alloc();
	this->splineU = gsl_spline_alloc(gsl_interp_cspline, dimF_in);
	gsl_spline_init(this->splineU, f_in, y_in, dimF_in);
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
