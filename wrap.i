%module eosWrap

%{
#define SWIG_FILE_WITH_INIT
#include "KVOR.h"
#include "Walecka.h"
#include "eos.h"
#include "setconst.h"
#include "DriverBase.h"
#include "KVDriver.h"
#include "TOV.h"
#include "aux.h"
#include "solve.h"
#include "KVORmod.h"
%}

%include "std_vector.i"

%include "carrays.i"
%include "numpy.i"

%init %{
import_array();
%}



namespace std{
%template(vec) vector <double>;
}

%array_class(double, dArray)

%apply (double* INPLACE_ARRAY1, int DIM1) {(double* in, int n)};
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* in, int n),(double* in2, int n2)};
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* n, int dimN)};
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* init, int dimInit)};
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* out, int dimOut)};
%include "eos.h" 
%include "setconst.h"
%include "KVOR.h"
%include "KVORmod.h"
%include "Walecka.h"
%include "DriverBase.h"
%apply (double * INPLACE_ARRAY1, int DIM1) {(double * E, int dimE), (double * P, int dimP), (double * n, int dimN)};
%include "KVDriver.h"
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* result, int dimResult)};
%include "TOV.h"
%apply (double * INPLACE_ARRAY1, int DIM1) {(double * n, int dimN)};
%include "aux.h"
%include "solve.h"