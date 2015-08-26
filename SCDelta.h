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
};

class KVOR_d : public SCDelta, public KVOR{
public:
    KVOR_d() : SCDelta(), KVOR(){
    }
};

class Walecka_d: public SCDelta, public Walecka{
public:
    Walecka_d() : Walecka(), SCDelta(){
        double mn = 938 / 197.33;
        Cs = 9.927 * mn * mn;
        Co = 4.82 * mn * mn;
        Cr = 4.791 * mn * mn;
        b = 0.008659;
        c = -0.002421;
        f0 = .22;
    }
};

class MKVOR_d: public SCDelta, public KVOR_mod2{
public:
    MKVOR_d() : KVOR_mod2(), SCDelta(){

    }
};

#endif //EOSWRAP_SCDELTA_H
