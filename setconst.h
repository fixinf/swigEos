///*
// * setconst.h
// *
// *  Created on: 09 июня 2014 г.
// *      Author: const
// */
//
//#ifndef SETCONST_H_
//#define SETCONST_H_
//#include <vector>
//
//typedef std::vector<double> vec;
//
//struct var{
//	std::vector<double> n;
//	double f;
//};
//
//class set_const {
//public:
//	set_const();
//	virtual ~set_const();
//	virtual double eta_s(double) = 0;
//	virtual double eta_o(double) = 0;
//	virtual double eta_r(double) = 0;
//	virtual double phi_n(double) = 0;
//	virtual double U(double) = 0;
//	vec * M;
//	vec * Q;
//	vec * T3;
//	double Cs;
//	double Co;
//	double Cr;
//	double b;
//	double c;
//	double f0;
//};
//
//
//#endif /* SETCONST_H_ */


#include <string>
#include <functional>
#include <vector>
#include <math.h>
/*
 * set_const.h
 *
 *  Created on: 16.05.2013
 *      Author: fixinf
 */

#ifndef SET_CONST_H_
#define SET_CONST_H_

typedef std::vector<double> vec;

struct var{
	std::vector<double> n;
	double f;
};


class set_const {
public:

	//			Cs		Co		Cr		b		c		z
	set_const(){
		init(179.56,87.6,100.64,7.7346e-3,3.4462e-4, 0.65);
		this->n0 = pow(197.33/135.0, 3) * 0.16;
		this->phi_meson = 0;
		this->sprime = 0;
		this->Csp = 1.0;
		this->sigma_kind = 0;
	}
	set_const(double, double, double, double, double, double);
	//			Cs		Co		Cr		b		c		z

	set_const(std::string name, double Cs, double Co , double Cr, double b, double c, double z);
	void init(double, double, double, double, double, double);
	double diff_phi_n(double);
	void set(double * p, int dimP);
	virtual double U(double) = 0;
	double dU(double);
	void set_name(std::string name);
	virtual double phi_n(int, double) = 0;
	virtual double eta_s(double) = 0;
	virtual double eta_o(double) = 0;
	virtual double eta_r(double) = 0;
	virtual double eta_p(double) = 0;
	double Cs;
	double Co;
	double Cr;
	double Csp;
	std::string name;
	double b;
	double c;
	double z;
	double f0;
	double exp_alpha;
	std::vector<double> X_s;
	std::vector<double> X_o;
	std::vector<double> X_p;
	std::vector<double> hs_z;
	std::vector<double> hs_alpha;
	vec X_r;
	vec X_sp;
	std::vector<double> Q;
	std::vector<double> T;
	vec M;
	bool Hyper;
	int SetHyperConstants(int);
	std::string repr();
	bool phi_meson;
	double n0;
	bool sprime;
	void set_xo(double * x, int dimX);
	void set_xr(double * x, int dimX);
	void set_xp(double * x, int dimX);
	void set_xs(double * x, int dimX);
	double func(double x);
	double Xs(int i, double f);
	void set_hs_z(double* x, int dimX);
	void set_hs_alpha(double* x, int dimX);
	int sigma_kind;
	bool muons;
	double fmax;
	//std::function<double(double)> phi_n;
	//std::function<double(double)> eta_s;
	//std::function<double(double)> eta_o;
	//std::function<double(double)> eta_r;
};

#endif /* SET_CONST_H_ */
