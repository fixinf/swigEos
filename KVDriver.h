/*
 * KVDriver.h
 *
 *  Created on: 25 июня 2014 г.
 *      Author: const
 */

#ifndef KVDRIVER_H_
#define KVDRIVER_H_

#include <string>

#include "DriverBase.h"
#include <gsl/gsl_spline.h>

using namespace std;
class KVDriver: public DriverBase {
public:
	KVDriver();
	KVDriver(double * E, int dimE, double * P, int dimP, double * n, int dimN);
	virtual ~KVDriver();


	const std::string& getName() const {
		return name;
	}

	void setName(const std::string& name) {
		this->name = name;
	}

	const string& getFname() const {
		return fname;
	}

	void setFname(const string& fname) {
		this->fname = fname;
	}

	double NofP(double);
	double NofE(double);
	double EofP(double);
	double EofN(double);
	double PofN(double);
	int lookFor(double * src, int dim_src, double what);
	void set(double * E, int dimE, double * P, int dimP, double * n, int dimN);
private:
	gsl_spline * iNofE;
	gsl_spline * iNofP;
	gsl_spline * iEofP;
	gsl_spline * iEofN;
	gsl_spline * iPofN;

	gsl_interp_accel * accNofE;
	gsl_interp_accel * accNofP;
	gsl_interp_accel * accEofP;
	gsl_interp_accel * accEofN;
	gsl_interp_accel * accPofN;

	std::string name;
	string fname;
	bool isSet;
};

#endif /* KVDRIVER_H_ */
