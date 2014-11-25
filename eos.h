/*
 * eos.h
 *
 *  Created on: 09 июня 2014 г.
 *      Author: const
 */



#ifndef EOS_H_
#define EOS_H_

#include "setconst.h"
#include <vector>
#include "gsl/gsl_vector_double.h"
#include <cmath>
#include "constants.h"


extern void stepE(double n, double * init, int dimInit, double * f_init, int dimF_init, double * out, int dimOut, int iter, set_const *);
void potentials(double * n, int dimN, double * out, int dimOut, set_const * C);
//extern double stepF(var, set_const *);
extern double _E(double * n, int dimN, set_const *);
extern double E(double * n, int dimN, set_const *);
extern double sum(std::vector<double> x);
extern double kineticInt(double n, double m, double f);
extern double p_f(double n);
namespace calc{
extern double mu(double * n, int dimN, int i, set_const * C);
}

extern float sumTest(double * in, int n);
extern float sumTest2(double * in, int n, double * in2, int n2);
extern void solveF(double n, double E, double P, double * init, int dimInit, double * out, int dimOut, set_const * C);
#endif /* EOS_H_ */
