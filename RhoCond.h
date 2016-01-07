//
// Created by const on 11/19/15.
//

#ifndef EOSWRAP_RHOCOND_H
#define EOSWRAP_RHOCOND_H

#include "setconst.h"

extern double E_rho(double * n, int dimN, double rho0, double rho_c, double mu_c,
                    set_const * C, double * inplace = NULL, int dim_inplace = 0);



#endif //EOSWRAP_RHOCOND_H
