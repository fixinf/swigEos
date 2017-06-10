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
	virtual double NofP(double N);
	virtual double NofE(double N);
	virtual double EofP(double N) = 0;
	virtual double EofN(double N) = 0;
	virtual double PofN(double N) = 0;
	// virtual double PofE(double E) = 0;
	double * lastNstar;
	double * lastRstar;
	double * lastMstar;
	double * lastPstar;
	double * lastEstar;
	int nSize;
	void getLastN(double * N, int dimN);
	void getLastR(double * N, int dimN);
	void getLastM(double * N, int dimN);
	void getLastP(double * N, int dimN);
	void getLastE(double * N, int dimN);
};

#endif /* TOVDRIVER_H_ */
