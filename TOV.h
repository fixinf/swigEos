/*
 * TOV.h
 *
 *  Created on: 25 июня 2014 г.
 *      Author: const
 */

#ifndef TOV_H_
#define TOV_H_

#include <gsl/gsl_odeiv2.h>
#include "setconst.h"
#include "DriverBase.h"

void star(double rho_init, double * result, int dimResult, DriverBase* D);
void star2(double rho_init, double * result, int dimResult, DriverBase* D);

//void star_crust(double rho_init, double * result, int dimResult, DriverBase* D);
void star_crust(double rho_init, double * result, int dimResult, DriverBase* D, double nmin);
void star_crust2(double rho_init, double * result, int dimResult, DriverBase* D, double nmin);

#endif /* TOV_H_ */
