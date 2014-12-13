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
	double bmass = 0.;
	double f[3];
	int status = eq_volkov(r_init, y, f, NULL);
	for (i = 1; i <= 10000; i++) {
//		printf("%f %f %f %f \n \r", t, y[0], y[1], y[2]);
		double ti = i * t1 / 1000.0;
		if (y[0] > delta * P_init) {
			status = gsl_odeiv2_driver_apply(d, &t, ti, y);

			if (status != GSL_SUCCESS) {
				printf("error, return value=%d\n", status);
				break;
			}
		} else {
//			std::cout << "RES2 " << y[1] << "      " << t << std::endl;
			result[0] = y[1];
			result[1] = t;
//			if (dimResult > 2){
//				y[2] = bmass;
//			}
			gsl_odeiv2_driver_free(d);
			return;
		}
		y[2] = D->EofP(y[0]);
//		bMass += (1.479345850e-12)*D->NofE(y[2])*pow((1.0 - (2.0*y[1])/t),-0.5)*t*t*(t1/n_points);
	}
	return;
}

//void star_crust(double rho_init, double* result, int dimResult, DriverBase* D) {
//	gsl_odeiv2_system sys = { eq_volkov, NULL, 3, NULL};
//		double delta = 0.;
//		gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys,
//			gsl_odeiv2_step_rkf45, 1e-3, 1e-3, 1e-3);
//		int i;
//		double step = 0.0;
//		double r_init = 1e-6;
//		double t = r_init, t1 = 30.0 + r_init;
//		double P_init = D->PofN(rho_init);
//		double E_init = D->EofN(rho_init);
//		double y[3] = { D->PofN(rho_init), 1.333 * M_PI * pow(r_init, 3.0) * 1.4793,
//			E_init };
//		double f[3];
//		double bMass = 0.;
//		double Dim = pow(135.,4);
//		int status = eq_volkov(r_init, y, f, NULL);
//		int n_points = 10000;
//		double nmin = 1e-11;
//		for (i = 1; i <= n_points; i++) {
////			printf("%f %f %f %f \n \r", t, y[0], y[1], y[2]);
//			double ti = i * t1 / n_points;
//			if ((y[0] > delta * P_init)&&(D->NofP(y[0]) > nmin)) {
//				status = gsl_odeiv2_driver_apply(d, &t, ti, y);
//				printf("r = %f n = %f nmin= %f P = %f \n", ti, D->NofP(y[0]), nmin, y[0]);
//				if (status != GSL_SUCCESS) {
//					printf("error, return value=%d\n", status);
//					break;
//				}
//			} else {
//	//			std::cout << "RES2 " << y[1] << "      " << t << std::endl;
//				result[0] = y[1];
//				result[1] = t;
//				if (dimResult > 2){
//					result[2] = bMass;
//				}
//				gsl_odeiv2_driver_free(d);
//				return;
//			}
//			y[2] = D->EofP(y[0]);
//			bMass += (1.479345850e-12)*D->NofE(y[2])*pow((1.0 - (2.0*y[1])/t),-0.5)*t*t*(t1/n_points);
//		}
//		return;
//}

void star_crust(double rho_init, double* result, int dimResult, DriverBase* D, double nmin) {
	gsl_odeiv2_system sys = { eq_volkov, NULL, 3, NULL};
		double delta = 0.;
		gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys,
			gsl_odeiv2_step_rkf45, 1e-3, 1e-3, 1e-3);
		int i;
		double step = 0.0;
		double r_init = 1e-6;
		double t = r_init, t1 = 30.0 + r_init;
		double P_init = D->PofN(rho_init);
		double E_init = D->EofN(rho_init);
		double y[3] = { D->PofN(rho_init), 1.333 * M_PI * pow(r_init, 3.0) * 1.4793,
			E_init };
		double f[3];
		double bMass = 0.;
		double Dim = pow(135.,4);
		int status = eq_volkov(r_init, y, f, NULL);
		int n_points = 10000;
		for (i = 1; i <= n_points; i++) {
//			printf("%f %f %f %f \n \r", t, y[0], y[1], y[2]);
			double ti = i * t1 / n_points;
			if ((y[0] > delta * P_init)&&(D->NofP(y[0]) > nmin)) {
				status = gsl_odeiv2_driver_apply(d, &t, ti, y);
//				printf("r = %f n = %f nmin= %f P = %f \n", ti, D->NofP(y[0]), nmin, y[0]);
				if (status != GSL_SUCCESS) {
					printf("error, return value=%d\n", status);
					break;
				}
			} else {
	//			std::cout << "RES2 " << y[1] << "      " << t << std::endl;
				result[0] = y[1];
				result[1] = t;
				if (dimResult > 2){
					result[2] = bMass;
				}
				gsl_odeiv2_driver_free(d);
				return;
			}
			y[2] = D->EofP(y[0]);
			bMass += (1.479345850e-12)*D->NofE(y[2])*pow((1.0 - (2.0*y[1]*1.4766)/t),-0.5)*t*t*(t1/n_points);
		}
		return;
}

int eq_volkov2(double r, const double y[], double f[], void *params) {
	bool debug = false;
	double A = 1.4766;
	double B = 1.0;
	double C = 0.0004898007281478712;
	double D = 2.9532;
	double E1 = 0.0004898007281478712;
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

void star2(double rho_init, double * result, int dimResult, DriverBase* D) {
	gsl_odeiv2_system sys = { eq_volkov2, NULL, 3, NULL};
	double delta = 5e-3;
	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys,
		gsl_odeiv2_step_rkf45, 1e-3, 1e-3, 1e-3);
	int i;
	double r_init = 1e-6;
	double t = r_init, t1 = 15.0 + r_init;
	double P_init = D->PofN(rho_init);
	double E_init = D->EofN(rho_init);
	double y[3] = { D->PofN(rho_init), 1.333 * M_PI * pow(r_init, 3.0) * 1.4766,
		E_init };
	double bmass = 0.;
	double f[3];
	int status = eq_volkov(r_init, y, f, NULL);
	for (i = 1; i <= 10000; i++) {
//		printf("%f %f %f %f \n \r", t, y[0], y[1], y[2]);
		double ti = i * t1 / 1000.0;
		if (y[0] > delta * P_init) {
			status = gsl_odeiv2_driver_apply(d, &t, ti, y);

			if (status != GSL_SUCCESS) {
				printf("error, return value=%d\n", status);
				break;
			}
		} else {
//			std::cout << "RES2 " << y[1] << "      " << t << std::endl;
			result[0] = y[1];
			result[1] = t;
//			if (dimResult > 2){
//				y[2] = bmass;
//			}
			gsl_odeiv2_driver_free(d);
			return;
		}
		y[2] = D->EofP(y[0]);
//		bMass += (1.479345850e-12)*D->NofE(y[2])*pow((1.0 - (2.0*y[1])/t),-0.5)*t*t*(t1/n_points);
	}
	return;
}

void star_crust2(double rho_init, double* result, int dimResult, DriverBase* D, double nmin) {
	gsl_odeiv2_system sys = { eq_volkov2, NULL, 3, NULL};
	double delta = 0.;
	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys,
		gsl_odeiv2_step_rkf45, 1e-3, 1e-3, 1e-3);
	int i;
	double step = 0.0;
	double r_init = 1e-6;
	double t = r_init, t1 = 50.0 + r_init;
	double P_init = D->PofN(rho_init);
	double E_init = D->EofN(rho_init);
	double y[3] = { D->PofN(rho_init), 1.333 * M_PI * pow(r_init, 3.0) * 1.4766,
		E_init };
	double f[3];
	double bMass = 0.;
	double Dim = pow(135.,4);
	int status = eq_volkov(r_init, y, f, NULL);
	int n_points = 10000;

	printf("Deleting... ");

	if (D->lastNstar != 0){
		printf("N... ");
		delete D->lastNstar;
	}
	if (D->lastRstar != 0){
		printf("N... ");
		delete D->lastRstar;
	}

	if (D->lastMstar != 0){
		printf("N... ");
		delete D->lastMstar;
	}

	if (D->lastEstar != 0){
		printf("N... ");
		delete D->lastEstar;
	}

	if (D->lastPstar != 0){
		printf("N... ");
		delete D->lastPstar;
	}

	printf(" done.\n");


//	D->nSize = n_points - 1;
	printf("Allocating... ");
	printf("N... ");
	D->lastNstar = new double[n_points-1];
	printf("R... ");
	D->lastRstar = new double[n_points-1];
	printf("M... ");
	D->lastMstar = new double[n_points-1];
	printf("E... ");
	D->lastEstar = new double[n_points-1];
	printf("P... ");
	D->lastPstar = new double[n_points-1];
	printf(" done.\n");
	for (i = 1; i <= n_points; i++) {
//			printf("%f %f %f %f \n \r", t, y[0], y[1], y[2]);
		double ti = i * t1 / n_points;
		if ((y[0] > delta * P_init)&&(D->NofP(y[0]) > nmin)) {
			status = gsl_odeiv2_driver_apply(d, &t, ti, y);
//				printf("r = %f n = %f nmin= %f P = %f \n", ti, D->NofP(y[0]), nmin, y[0]);
			if (status != GSL_SUCCESS) {
				printf("error, return value=%d\n", status);
				break;
			}
		} else {
//			std::cout << "RES2 " << y[1] << "      " << t << std::endl;
			result[0] = y[1];
			result[1] = t;
			if (dimResult > 2){
				result[2] = bMass;
			}
			gsl_odeiv2_driver_free(d);
			D->nSize = i;
			return;
		}
		y[2] = D->EofP(y[0]);
		bMass += (0.0004898007281478712)*D->NofE(y[2])*pow((1.0 - (2.0*y[1]*1.4766)/t),-0.5)*t*t*(t1/n_points);
		D->lastNstar[i-1] = D->NofE(y[2]);
		D->lastRstar[i-1] = ti;
		D->lastMstar[i-1] = y[1];
		D->lastEstar[i-1] = y[2];
		D->lastPstar[i-1] = y[0];
	}
	return;
}


