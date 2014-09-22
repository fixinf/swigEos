/*
 * TOVDriver.cpp
 *
 *  Created on: 22 июня 2014 г.
 *      Author: const
 */

#include "DriverBase.h"

DriverBase::DriverBase() {
	this->fname = "";
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
