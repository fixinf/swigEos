/*
 * aux.h
 *
 *  Created on: 06 июля 2014 г.
 *      Author: const
 */

#ifndef AUX_H_
#define AUX_H_
#include "setconst.h"
#include "eos.h"

extern void f_eq(double * n, int dimN, double * init, int dimInit, double * res, int dimRes, set_const * C);
//extern double f_eq(double * n, int dimN, set_const * C, double init = 1e-1);

extern double K(double n, set_const *C);
extern double EBind(double * n, int dimN, set_const *C);
extern double J(double n, set_const * C);
extern double J(double n, set_const * C, double f);


#endif /* AUX_H_ */
