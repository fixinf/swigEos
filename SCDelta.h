//
// Created by const on 8/17/15.
//

#ifndef EOSWRAP_SCDELTA_H
#define EOSWRAP_SCDELTA_H


#include "setconst.h"
#include "KVOR.h"
#include "Walecka.h"
#include "KVORmod2.h"

class SCDelta: virtual public set_const {
public:
    SCDelta(){
        this->setDeltaConstants(2, 0);
        double alpha[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        this->set_hs_alpha(alpha, 9);
        this->set_hs_z(alpha, 9);
    }

    int setDeltaConstants(int, int);

    void setDeltaRho(double * X, int dimX){
        for (int i = 0; i < dimX; i++){
            this->X_r[i+8] = X[i];
        }
    }
    void setDeltaOmega(double * X, int dimX){
        for (int i = 0; i < dimX; i++){
            this->X_o[i+8] = X[i];
        }
    }
    void setDeltaSigma(double * X, int dimX){
        for (int i = 0; i < dimX; i++){
            this->X_s[i+8] = X[i];
        }
    }
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

#endif //EOSWRAP_SCDELTA_H
