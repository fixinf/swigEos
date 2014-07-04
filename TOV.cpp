/*
 * TOV.cpp
 *
 *  Created on: 25 июня 2014 г.
 *      Author: const
 */


#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <cmath>
#include <cstdio>
#include <iostream>

#include "DriverBase.h"

class set_const;

int eq_volkov(double r, const double y[], double f[], void *params) {
	bool debug = false;
	double A = 1.4766;
	double B = 1.0;
	double C = 1.47934585e-12;
	double D = 2.9532;
	double E1 = 1.479345850e-12;
	if (debug)
		printf("r :: %f P :: %f M :: %.e E :: %f \n", r, y[0], y[1], y[2]);
	f[0] = -A * y[1] * y[2] / (r * r);
	if (debug)
		printf("f[0]:1 :: %f \n", f[0]);
	f[0] *= (1.0 + B * y[0] / y[2]);
	if (debug)
		printf("f[0]:2 :: %f \n", f[0]);
	f[0] *= (1.0 + C * y[0] * pow(r, 3.0) / y[1]);
	if (debug)
		printf("f[0]:3 :: %f \n", f[0]);
	f[0] /= (1.0 - D * y[1] / r);
	if (debug)
		printf("f[0]:4 :: %f \n", f[0]);
	f[1] = r * r * E1 * y[2];
	if (debug)
		printf("f[1] :: %f \n", f[1]);
	f[2] = 0.0;
	return GSL_SUCCESS;
}

void star(double rho_init, double * result, int dimResult, DriverBase* D) {
	gsl_odeiv2_system sys = { eq_volkov, NULL, 3, NULL};
	double delta = 5e-3;
	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys,
		gsl_odeiv2_step_rkf45, 1e-3, 1e-3, 1e-3);
	int i;
	double r_init = 1e-6;
	double t = r_init, t1 = 15.0 + r_init;
	double P_init = D->PofN(rho_init);
	double E_init = D->EofN(rho_init);
	double y[3] = { D->PofN(rho_init), 1.333 * M_PI * pow(r_init, 3.0) * 1.4793,
		E_init };
	double f[3];
	int status = eq_volkov(r_init, y, f, NULL);
	for (i = 1; i <= 10000; i++) {
		printf("%f %f %f %f \n \r", t, y[0], y[1], y[2]);
		double ti = i * t1 / 1000.0;
		if (y[0] > delta * P_init) {
			status = gsl_odeiv2_driver_apply(d, &t, ti, y);

			if (status != GSL_SUCCESS) {
				printf("error, return value=%d\n", status);
				break;
			}
		} else {
			std::cout << "RES2 " << y[1] << "      " << t << std::endl;
			result[0] = y[1];
			result[1] = t;
			gsl_odeiv2_driver_free(d);
			return;
		}
		y[2] = D->EofP(y[0]);
	}
	return;
}
