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

extern double f_eq(double * n, int dimN, set_const * C, double init = 1e-1);
extern double K(double n, set_const *C);
extern double EBind(double * n, int dimN, set_const *C);
extern double J(double n, set_const * C);


#endif /* AUX_H_ */
