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
#include "setconst.h"

using namespace std;
class KVDriver: public DriverBase {
public:
	KVDriver();
	KVDriver(double * E, int dimE, double * P, int dimP, double * n, int dimN);
	KVDriver(set_const *, string fname);
	virtual ~KVDriver();

	set_const* getC() const {
		return C;
	}

	void setC(set_const* c) {
		C = c;
	}

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
	double PofE(double);
	double EofN(double);
	double PofN(double);
	double EofP(double);
	double NofP(double);
	double NofE(double);
	int lookFor(double * src, int dim_src, double what);
	void set(double * E, int dimE, double * P, int dimP, double * n, int dimN);
private:
	set_const * C;
	std::string name;
	string fname;
};

#endif /* KVDRIVER_H_ */
