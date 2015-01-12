/*
 * InterpolatedScalings.h
 *
 *  Created on: 26 нояб. 2014 г.
 *      Author: const
 */

#ifndef INTERPOLATEDSCALINGS_H_
#define INTERPOLATEDSCALINGS_H_

#include "setconst.h"
//#include "interpolation.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

class InterpolatedScalings: public set_const {
public:
	InterpolatedScalings();
	virtual ~InterpolatedScalings();

	virtual double U(double);
	virtual double phi_n(int, double);
	virtual double eta_s(double);
	virtual double eta_o(double);
	virtual double eta_r(double);
	virtual double eta_p(double);
	bool debug;
	void set_eta_s(double * f_in, int dimF_in, double * y_in, int dimY_in);
	void set_eta_o(double * f_in, int dimF_in, double * y_in, int dimY_in);
	void set_eta_r(double * f_in, int dimF_in, double * y_in, int dimY_in);
	void set_U(double * f_in, int dimF_in, double * y_in, int dimY_in);
	bool rho_akima;
//	alglib::spline1dinterpolant * splineU;

private:
	gsl_spline * splineU;
	gsl_spline * splineO;
	gsl_spline * splineR;
	gsl_spline * splineS;
	gsl_interp_accel * accU;
	gsl_interp_accel * accO;
	gsl_interp_accel * accR;
	gsl_interp_accel * accS;

};

#endif /* INTERPOLATEDSCALINGS_H_ */
