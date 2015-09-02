//
// Created by const on 8/26/15.
//

#ifndef EOSWRAP_BARYON_H
#define EOSWRAP_BARYON_H

#include <vector>

typedef std::vector<double> vec;

class Baryon {
public:
    Baryon(double m, double q, double t, double s, double xs, double xo,
           double xr, double xp, double xsp){
        this->m = m;
        this->q = q;
        this->t = t;
        this->s = s;
        this->xs = xs;
        this->xo = xo;
        this->xp = xp;
        this->xsp = xsp;
    }
    double m, q, t, s, xs, xo, xr, xp, xsp;
};


#endif //EOSWRAP_BARYON_H
