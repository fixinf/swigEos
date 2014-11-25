/*
 * EnthalpyTOV.cpp
 *
 *  Created on: 21 нояб. 2014 г.
 *      Author: const
 */

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <cmath>
#include <cstdio>
#include <iostream>

int equation(double r, const double y[], double f[], void *params) {
	double m_sun = 1.4766;
	double D = 5.7207e-5;
	double E_const = 4.898007e-4;
}
