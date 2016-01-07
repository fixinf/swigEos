//
// Created by const on 8/17/15.
//

#include "SCDelta.h"

int SCDelta::setDeltaConstants(int hyper, int delta) {
    this->SetHyperConstants(hyper);
    this->X_o.push_back(1);
    this->X_r.push_back(1);
    this->X_s.push_back(1);
    this->X_p.push_back(0);
    this->T.push_back(-1.5);
    this->M.push_back(1232./135.);
    this->Q.push_back(-1);
    this->S.push_back(1.5);

    this->X_o.push_back(1);
    this->X_r.push_back(1);
    this->X_s.push_back(1);
    this->X_p.push_back(0);
    this->T.push_back(-0.5);
    this->M.push_back(1232./135.);
    this->Q.push_back(0);
    this->S.push_back(1.5);

    this->X_o.push_back(1);
    this->X_r.push_back(1);
    this->X_s.push_back(1);
    this->X_p.push_back(0);
    this->T.push_back(0.5);
    this->M.push_back(1232./135.);
    this->Q.push_back(1);
    this->S.push_back(1.5);

    this->X_o.push_back(1);
    this->X_r.push_back(1);
    this->X_s.push_back(1);
    this->X_p.push_back(0);
    this->T.push_back(1.5);
    this->M.push_back(1232./135.);
    this->Q.push_back(2);
    this->S.push_back(1.5);
}

int SCDelta::setDeltaOnlyConstants() {

    this->offset = 2;
    this->X_o.clear();
    this->X_r.clear();
    this->X_s.clear();
    this->X_p.clear();
    this->T.clear();
    this->M.clear();
    this->Q.clear();
    this->S.clear();

    this->X_o.push_back(1);
    this->X_r.push_back(1);
    this->X_s.push_back(1);
    this->X_p.push_back(0);
    this->T.push_back(-0.5);
    this->M.push_back(938./135.);
    this->Q.push_back(0);
    this->S.push_back(0.5);

    this->X_o.push_back(1);
    this->X_r.push_back(1);
    this->X_s.push_back(1);
    this->X_p.push_back(0);
    this->T.push_back(0.5);
    this->M.push_back(938./135.);
    this->Q.push_back(1);
    this->S.push_back(0.5);

    this->X_o.push_back(1);
    this->X_r.push_back(1);
    this->X_s.push_back(1);
    this->X_p.push_back(0);
    this->T.push_back(-1.5);
    this->M.push_back(1232./135.);
    this->Q.push_back(-1);
    this->S.push_back(1.5);

    this->X_o.push_back(1);
    this->X_r.push_back(1);
    this->X_s.push_back(1);
    this->X_p.push_back(0);
    this->T.push_back(-0.5);
    this->M.push_back(1232./135.);
    this->Q.push_back(0);
    this->S.push_back(1.5);

    this->X_o.push_back(1);
    this->X_r.push_back(1);
    this->X_s.push_back(1);
    this->X_p.push_back(0);
    this->T.push_back(0.5);
    this->M.push_back(1232./135.);
    this->Q.push_back(1);
    this->S.push_back(1.5);

    this->X_o.push_back(1);
    this->X_r.push_back(1);
    this->X_s.push_back(1);
    this->X_p.push_back(0);
    this->T.push_back(1.5);
    this->M.push_back(1232./135.);
    this->Q.push_back(2);
    this->S.push_back(1.5);
}