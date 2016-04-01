//
// Created by const on 8/17/15.
//

#ifndef EOSWRAP_SCDELTA_H
#define EOSWRAP_SCDELTA_H


#include "setconst.h"
#include "KVOR.h"
#include "Walecka.h"
#include "KVORmod2.h"
#include "MKVOR.h"

class SCDelta: virtual public set_const {
public:
    SCDelta(){
        this->setDeltaConstants(2, 0);
        double alpha[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        this->set_hs_alpha(alpha, 9);
        this->set_hs_z(alpha, 9);
        this->offset = 8;
    }

    int setDeltaConstants(int, int);

    void setDeltaRho(double * X, int dimX){
        for (int i = 0; i < dimX; i++){
            this->X_r[i+offset] = X[i];
        }
    }
    void setDeltaOmega(double * X, int dimX){
        for (int i = 0; i < dimX; i++){
            this->X_o[i+offset] = X[i];
        }
    }
    void setDeltaSigma(double * X, int dimX){
        for (int i = 0; i < dimX; i++){
            this->X_s[i+offset] = X[i];
        }
    }
    int offset;
    int setDeltaOnlyConstants();
};

class KVOR_d : public KVOR, public SCDelta{
public:
    KVOR_d() : KVOR(), SCDelta(){
    }
};

class Walecka_d: public Walecka, public SCDelta{
public:
    Walecka_d() : Walecka(), SCDelta(){
        Cs = 246;
        Co = 156.3;
        b = 1.8e-3;
        c = 2.87e-4;
        n0 = 163./160*n0;
        f0 = .22;
    }
};

class MKVOR_d: public KVOR_mod2, public SCDelta{
public:
    MKVOR_d() : KVOR_mod2(), SCDelta(){

    }
};

class KVORcut_d: public KVOR_cut, public SCDelta{
public:
    KVORcut_d() : KVOR_cut(), SCDelta(){

    }
};

class KVORcut_sigma_d: public KVOR_cut_sigma, public SCDelta{
public:
    KVORcut_sigma_d() : KVOR_cut_sigma(), SCDelta(){

    }
};

class MKVOR_delta: public MKVOR, public SCDelta{
public:
    MKVOR_delta(): MKVOR(), SCDelta(){

    }
};

class MKVOR2: public MKVOR, public SCDelta{
public:
    MKVOR2() : MKVOR(), SCDelta(){
        this->fcut_om = 100500.;
        this->fcut_rho = 100500.;
        this->acut_rho = 1;
        this->acut_om = 1;

        this->bcut_rho = 1;
        this->bcut_om = 1;
    }

    double eta_o(double f){
        return MKVOR::eta_o(f) * 0.5 * acut_om * (1
                - pow(tanh(bcut_om * (f - fcut_om)), 1));
    }

    double eta_r(double f){
        return MKVOR::eta_r(f) * 0.5 * acut_rho * (1
                - pow(tanh(bcut_rho * (f - fcut_rho)), 1));
    }

    double fcut_om;
    double bcut_om;
    double acut_om;

    double fcut_rho;
    double acut_rho;
    double bcut_rho;
};
#endif //EOSWRAP_SCDELTA_H
