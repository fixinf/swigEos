/*
 * TOVDriver.cpp
 *
 *  Created on: 22 июня 2014 г.
 *      Author: const
 */

#include "DriverBase.h"

DriverBase::DriverBase() {
	this->fname = "";
	this->lastNstar = 0;
	this->lastRstar = 0;
	this->lastMstar = 0;
	this->lastEstar = 0;
	this->lastPstar = 0;
}

DriverBase::DriverBase(string fname) {
	this->fname = fname;
	readEos();
}

DriverBase::~DriverBase() {

}

void DriverBase::readEos(){

}

double DriverBase::PofE(double E){
	return 0.0;
}

double DriverBase::PofN(double n) {
	return 0;
}

double DriverBase::EofN(double n) {
	return 0.0;
}

double DriverBase::EofP(double P){
	return 0.0;
}

double DriverBase::NofE(double E){
	return 0.;
}

double DriverBase::NofP(double N){
	return 0.;
}

void DriverBase::getLastN(double * N, int dimN){
	if (this->lastNstar){
		for (int i = 0; i < dimN; i++){
			N[i] = this->lastNstar[i];
		}
	}
	else{
		printf("Last N not set! \n");
		return;
	}
}
void DriverBase::getLastR(double * N, int dimN){
	if (this->lastRstar){
		for (int i = 0; i < dimN; i++){
			N[i] = this->lastRstar[i];
		}
	}
	else{
		printf("Last R not set! \n");
		return;
	}
}

void DriverBase::getLastM(double * N, int dimN){
	if (this->lastMstar){
		for (int i = 0; i < dimN; i++){
			N[i] = this->lastMstar[i];
		}
	}
	else{
		printf("Last M not set! \n");
		return;
	}
}

void DriverBase::getLastE(double * N, int dimN){
	if (this->lastEstar){
		for (int i = 0; i < dimN; i++){
			N[i] = this->lastEstar[i];
		}
	}
	else{
		printf("Last E not set! \n");
		return;
	}
}

void DriverBase::getLastP(double * N, int dimN){
	if (this->lastPstar){
		for (int i = 0; i < dimN; i++){
			N[i] = this->lastPstar[i];
		}
	}
	else{
		printf("Last P not set! \n");
		return;
	}
}
