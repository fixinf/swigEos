///*
// * setconst.cpp
// *
// *  Created on: 09 июня 2014 г.
// *      Author: const
// */
//
//#include "setconst.h"
//#include "constants.h"
//
//set_const::set_const() {

//	this->M = new vec();
//	M->push_back(m_n);
//	M->push_back(m_n);
//	this->Q = new vec();
//	Q->push_back(0.0);
//	Q->push_back(1.0);
//	this->T3 = new vec();
//	T3->push_back(0.5);
//	T3->push_back(-0.5);
//}
//
//set_const::~set_const() {

//}


/*
 * set_const.cpp
 *
 *  Created on: 16.05.2013
 *      Author: fixinf
 */

#include "setconst.h"

#include <math.h>
#include <cstdio>
//#include <iostream>
//#include <string>

#include "constants.h"

//double set_const::phi_n(double f){
//	return 1.0 - f;
//	//return 1.0/(1.0 + f);
//}

double set_const::diff_phi_n(double f){
	double df = 0.0001;
	return (this->phi_n(0,f+df) - this->phi_n(0,f))/(df);
}

//double set_const::eta_r(double f){
//	return 1.0;
//	//return this->eta_o(f)/(this->eta_o(f) +
//		//4*pow(this->C_o/this->C_r,2.0)*(this->eta_o(f)-1.0));
//}

double set_const::dU(double f){
	return pow(M[0],4.0)*(b*pow(f,2.0) + c*pow(f,3.0));
}


//double set_const::eta_o(double f){
//	return 1;
//	//return (1 + this->z*0.195)/( 1 + this->z*f);
//}
//
//
//double set_const::eta_s(double f){
//	return 1.0;
//}
//
//double set_const::U(double f){
//
//	//std::cout << "Cs = " << C_s << "B = " << this->b << " C = " << this->c << std::endl;
//	return pow(m_n,4.0)*(b * pow(f,3.0)/3.0 + c*pow(f,4.0)/4.0);
//}

std::string set_const::repr(){
	return this->name;
}


set_const::set_const(double C_s, double C_o, double C_r, double b, double c, double z) {
	this->init(C_s,C_o,C_r,b,c,z);
}

set_const::set_const(std::string name, double C_s, double C_o , double C_r, double b, double c, double z){
	this->init(C_s,C_o,C_r,b,c,z);
	this->name = name;
}

void set_const::set_name(std::string name){
	this->name = name;
}

void set_const::init(double C_s, double C_o, double C_r, double b, double c, double z){
	this->Hyper = true;
	this->SetHyperConstants(0);
	this->Co = C_o;
	this->Cr = C_r;
	this->Cs = C_s;
	this->b = b;
	this->c = c;
	this->z = z;
	this->fmax = 1.;
	this->exp_alpha = 1.0;
}

int set_const::SetHyperConstants(int type){
	//Порядок следования барионов:
	// n p L0 S- S0 S+ X- X0

	this->X_o.clear();
	this->X_s.clear();
	this->X_r.clear();
	this->X_p.clear();
	this->X_sp.clear();
	this->T.clear();
	this->Q.clear();
	this->M.clear();
	double xo[8];
	double xr[8];
	double ebind[8];
	double xp[8];
	double xsp[8];
	double sq2 = sqrt(2.0);
	if (true){
		switch (type) {
			case 0://Quark counting
				double _xo[8] = { 1.0, 1.0, 2.0 / 3.0, 2.0/3.0, 2.0/3.0, 2.0/3.0, 1.0 / 3.0, 1.0 / 3.0 };
				double _xr[8] = { 1.0, 1.0, 0, 2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0 };
				double _ebind[8] = {0, 0, -30, 30, 30, 30, -18, -18};
				double _xp[8] = {0,0,0,0,0,0,0,0};
				for (int i =0; i<8; i++){
					xo[i] = _xo[i];
					xr[i] = _xr[i];
					xp[i] = _xp[i];
					ebind[i] = _ebind[i];
				}
				break;
			case 1://SU(3)
				double _2xo[8] = { 1.0, 1.0, 2.0 / 3.0, 2.0/3.0, 2.0/3.0, 2.0/3.0, 1.0 / 3.0, 1.0 / 3.0 };
				double _2xr[8] = { 1.0, 1.0, 0, 1.0, 1.0, 1.0, 1.0, 1.0};
				double _2ebind[8]= {0, 0, -30, 30, 30, 30, -18, -18};
				double _2xp[8] = {0,0,0,0,0,0,0,0};
				for (int i =0; i<8; i++){
					xo[i] = _2xo[i];
					xr[i] = _2xr[i];
					xp[i] = _2xp[i];
					ebind[i] = _2ebind[i];
				}
				break;

			case 2: //SU(6)
				double _3xo[8] = { 1.0, 1.0, 2.0 / 3.0, 2.0/3.0, 2.0/3.0, 2.0/3.0, 1.0 / 3.0, 1.0 / 3.0 };
				double _3xr[8] = { 1.0, 1.0, 0, 2.0, 2.0, 2.0, 1.0, 1.0};
				double _3ebind[8] = {0, 0, -30, 30, 30, 30, -18, -18};
				double _3xp[8] = {0.0, 0.0, -sq2/3,-sq2/3,-sq2/3,-sq2/3,-2*sq2/3,-2*sq2/3};
				for (int i =0; i<8; i++){
					xo[i] = _3xo[i];
					xr[i] = _3xr[i];
					xp[i] = _3xp[i];
					ebind[i] = _3ebind[i];
				}
				break;

			default:
				printf("Something Wrong with hyper constants label!!!\n");
				return -2;
				break;
		}
//		double ebind[8] = {0, 0, -30, 50, 50, 50, -18, -18};
		double m[8] = { 938/135.0, 938/135.0, 1116/135.0, 1195/135.0,
				1195/135.0, 1195/135.0 , 1317/135.0, 1317/135.0};
		//	double xs[8] = { 1, 1, 0, 0, 0, 0, 0, 0 };
		//	double xo[8] = { 1, 1, 0, 0, 0, 0, 0, 0 };
		//	double xr[8] = { 1, 1, 0, 0, 0, 0, 0, 0 };
		double xs[8];
		xs[0] = 1.0;
		xs[1] = 1.0;
		for (int i = 2; i < 8; i++){
//			xs[i] = (80.73*xo[i] - ebind[i]) / 140.70;
			xs[i] = ((Co*n0*135/(m[0]*m[0]))*xo[i] - ebind[i]) / (135.0*m[0]*f0);
		}

		double ebindLambda = -5.0;
		xsp[0] = 0.0;
		xsp[1] = 0.0;
		double xsp_l = ((Co*xo[2]*n0*135/(m[0]*m[0]))*xo[2] - ebindLambda) / (135.0*m[0]*f0);
		xsp[2] = xsp_l;
		xsp[3] = xsp_l;
		xsp[4] = xsp_l;
		xsp[5] = xsp_l;
		xsp[6] = 2*xsp_l;
		xsp[7] = 2*xsp_l;

		double t[8] = { -0.5, 0.5, 0.0, -1.0, 0.0, 1.0, -0.5, 0.5 };
		double q[8] = { 0, 1, 0, -1, 0, 1, -1, 0 };
		for (int i = 0; i < 8; i++){
				this->X_s.push_back(xs[i]);
				this->X_o.push_back(xo[i]);
				this->X_r.push_back(xr[i]);
				this->X_p.push_back(xp[i]);
				this->X_sp.push_back(xsp[i]);
				this->T.push_back(t[i]);
				this->Q.push_back(q[i]);
				this->M.push_back(m[i]);
		}
	}

	return 0;
}

void set_const::set(double* p, int dimP) {
	if(dimP < 5){
		printf("Array elements must be >= 5!");
		return;
	}
	this->Cs = p[0];
	this->Co = p[1];
	this->Cr = p[2];
	this->b = p[3];
	this->c = p[4];
}

void set_const::set_xo(double* x, int dimX) {
	this->X_o.clear();
	for (int i = 0; i < dimX; i++){
		X_o.push_back(x[i]);
	}
}

void set_const::set_xr(double* x, int dimX) {
	this->X_r.clear();
	for (int i = 0; i < dimX; i++){
		X_r.push_back(x[i]);
	}
}

void set_const::set_xp(double* x, int dimX) {
	this->X_p.clear();
	for (int i = 0; i < dimX; i++){
		X_p.push_back(x[i]);
	}
}

void set_const::set_xs(double* x, int dimX) {
	this->X_s.clear();
	X_s.push_back(1.0);
	X_s.push_back(1.0);
	for (int i = 2; i <dimX; i++){
		printf("x[i] = %f \n", x[i]);
		X_s.push_back(((Co*n0*135/(M[0]*M[0]))*X_o[i] - x[i])
				/ (135.0*M[0]*f0));
	}
}

void set_const::set_hs_z(double* x, int dimX) {
	printf("hey! \n");
	this->hs_z.clear();
	for (int i = 0; i < dimX; i++){
		printf("%f \n", x[i]);
		this->hs_z.push_back(x[i]);
	}
}

void set_const::set_hs_alpha(double* x, int dimX) {
	this->hs_alpha.clear();
	for (int i = 0; i < dimX; i++){
		this->hs_alpha.push_back(x[i]);
	}
}


double set_const::func(double x){
	if (x <= 0){
		return 0.0;
	}
	else{
		return exp(-exp_alpha/x);
	}
}

double set_const::Xs(int i, double f) {
	double res = this->X_s[i];
//	printf("i = %i \n", i);
	if (i > 1){
		res *= pow((1 + hs_z[i] * f0)/ (1 + hs_z[i]*f), hs_alpha[i]);
	}
	return res;
}
