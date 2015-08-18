//
// Created by const on 8/17/15.
//

#ifndef EOSWRAP_SCDELTA_H
#define EOSWRAP_SCDELTA_H


#include "setconst.h"
#include "KVOR.h"

class SCDelta: virtual public set_const {
public:
    SCDelta(){
        this->setDeltaConstants(2, 0);
        double alpha[8] = {0, 0, 0, 0, 0, 0, 0, 0};
        this->set_hs_alpha(alpha, 8);
        this->set_hs_z(alpha, 8);
    }

    int setDeltaConstants(int, int);
};

class KVOR_d : public SCDelta, public KVOR{
public:
    KVOR_d() : KVOR(), SCDelta(){
    }
};


#endif //EOSWRAP_SCDELTA_H
