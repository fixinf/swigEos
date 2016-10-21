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
#include "InterpolatedScalings.h"

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

class InterScDelta: public InterpolatedScalings, public SCDelta{
  public:
  InterScDelta() : InterpolatedScalings(), SCDelta(){
  
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

class MKVOR_pole: public MKVOR2{
public:
    MKVOR_pole(): MKVOR2(){
    }
        double eta_r(double f){
            return MKVOR::eta_r(f) + acut_rho / pow(fabs(f - fcut_rho), bcut_rho);
        }
    
};

class MKVOR3_tail : public MKVOR2{
  public:
  MKVOR3_tail() : MKVOR2(){
    tail_mult_rho = 1.;
    tail_mult_om = 1.;
  }

  double eta_r(double f){
    return 0.5 * MKVOR::eta_r(f) * (1 + tanh(bcut_rho * (fcut_rho - f))) + 
           0.5 * tail_mult_rho/(pow(f, acut_rho)) * (1 - tanh(bcut_rho * (fcut_rho - f)));
  }


  double eta_o(double f){
    return 0.5 * MKVOR::eta_o(f) * (1 + tanh(bcut_om * (fcut_om - f))) + 
           0.5 * tail_mult_om/(pow(f, acut_om)) * (1 - tanh(bcut_om * (fcut_om - f)));
  }
  
  double tail_mult_rho;
  double tail_mult_om;

};

class MKVOR_tail1: public MKVOR2{
public:
	MKVOR_tail1(): MKVOR2(){
		tail_mult_om = 1;
		tail_mult_rho = 1;
	}
	double eta_r(double f){
		return 0.5 * MKVOR::eta_r(f) * (1 + tanh(bcut_rho * (fcut_rho - f))) +
		           0.5 * tail_mult_rho/(1 + pow(fabs(f)/fcut_rho, acut_rho)) * (1 - tanh(bcut_rho * (fcut_rho - f)));
	}

	double eta_o(double f){
	    return 0.5 * MKVOR::eta_o(f) * (1 + tanh(bcut_om * (fcut_om - f))) +
	           0.5 * tail_mult_om/(1 + pow(fabs(f)/fcut_om, acut_om)) * (1 - tanh(bcut_om * (fcut_om - f)));
	  }

	double tail_mult_rho;
	double tail_mult_om;
};

class MKVOR_tail3: public MKVOR2{
public:
	MKVOR_tail3(): MKVOR2(){
		tail_mult_om = 1;
		tail_mult_rho = 1;
	}
	double eta_r(double f){
		return 0.5 * MKVOR::eta_r(f) * (1 + tanh(bcut_rho * pow(fcut_rho - f, 3))) +
		           0.5 * tail_mult_rho/(1 + pow(fabs(f)/fcut_rho, acut_rho)) * (1 - tanh(bcut_rho * pow(fcut_rho - f, 3)));
	}

	double eta_o(double f){
	    return 0.5 * MKVOR::eta_o(f) * (1 + tanh(bcut_om * pow(fcut_om - f, 3))) +
	           0.5 * tail_mult_om/(1 + pow(fabs(f)/fcut_om, acut_om)) * (1 - tanh(bcut_om * pow(fcut_om - f, 3)));
	  }

	double tail_mult_rho;
	double tail_mult_om;
};

class MKVOR_tail_poly: public MKVOR2{
public:
	MKVOR_tail_poly(): MKVOR2(){
		tail_mult_om = 1;
	}
	double eta_r(double f){
		if (f < fcut_rho){
			return MKVOR::eta_r(f);
		}
		else{
			return eta_r_new(f);
		}
	}

	double eta_r_new(double f){
		return 1./(acut_rho + bcut_rho * (f - fcut_rho) + c_cut_rho * pow(f - fcut_rho, 2) +
							d_cut_rho * pow(f-fcut_rho, 3));
	}

	double eta_o(double f){
	    return 0.5 * MKVOR::eta_o(f) * (1 + tanh(bcut_om * (fcut_om - f))) +
	           0.5 * tail_mult_om/(1 + pow(fabs(f)/fcut_om, acut_om)) * (1 - tanh(bcut_om *(fcut_om - f)));
	  }

	double c_cut_rho;
	double d_cut_rho;
	double tail_mult_om;
};

class MKVOR_tail_poly_exp: public MKVOR2{
public:
	MKVOR_tail_poly_exp(): MKVOR2(){
		tail_mult_om = 1;
	}
	double eta_r(double f){
		if (f < fcut_rho){
			return MKVOR::eta_r(f);
		}
		else{
			return eta_r_new(f);
		}
	}

	double eta_r_new(double f){
		return 1./(acut_rho + bcut_rho * (f - fcut_rho) + c_cut_rho * pow(f - fcut_rho, 2) +
							e_cut_rho * exp(d_cut_rho * pow(f-fcut_rho, 3)));
	}

	double eta_o(double f){
	    return 0.5 * MKVOR::eta_o(f) * (1 + tanh(bcut_om * (fcut_om - f))) +
	           0.5 * tail_mult_om/(1 + pow(fabs(f)/fcut_om, acut_om)) * (1 - tanh(bcut_om *(fcut_om - f)));
	  }

	double c_cut_rho;
	double d_cut_rho;
	double e_cut_rho;
	double tail_mult_om;
};

class MKVOR_tail_poly4: public MKVOR2{
public:
	MKVOR_tail_poly4(): MKVOR2(){
		tail_mult_om = 1;
	}
	double eta_r(double f){
		if (f < fcut_rho){
			return MKVOR::eta_r(f);
		}
		else{
			return eta_r_new(f);
		}
	}

	double eta_r_new(double f){
		return 1./(acut_rho + bcut_rho * (f/fcut_rho - 1) +
							c_cut_rho * pow(f/fcut_rho - 1, 2) +
							d_cut_rho * pow(f/fcut_rho - 1, 3) +
							e_cut_rho * pow(f/fcut_rho - 1, 4));
	}

	double eta_o(double f){
	    return 0.5 * MKVOR::eta_o(f) * (1 + tanh(bcut_om * (fcut_om - f))) +
	           0.5 * tail_mult_om/(1 + pow(fabs(f)/fcut_om, acut_om)) * (1 - tanh(bcut_om *(fcut_om - f)));
	  }

	double c_cut_rho;
	double d_cut_rho;
	double e_cut_rho;
	double tail_mult_om;
};
#endif //EOSWRAP_SCDELTA_H
