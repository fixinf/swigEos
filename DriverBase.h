/*
 * TOVDriver.h
 *
 *  Created on: 22 июня 2014 г.
 *      Author: const
 */
#include <string>
#include <fstream>
#ifndef TOVDRIVER_H_
#define TOVDRIVER_H_


using namespace std;
class DriverBase {
public:
	DriverBase();
	DriverBase(string fname);
	virtual ~DriverBase();
	string fname;
	virtual void readEos();
	double * E;
	double * P;
	double * n;
	int count;
	virtual double PofE(double E);
	virtual double PofN(double n);
	virtual double EofN(double n);
	virtual double EofP(double P);
};

#endif /* TOVDRIVER_H_ */
