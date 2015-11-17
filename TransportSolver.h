//
// Created by Paul VANHAESEBROUCK on 11/11/2015.
//

#ifndef STABLE_FLUIDS_TRANSPORTSOLVER_H
#define STABLE_FLUIDS_TRANSPORTSOLVER_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_interp2d.h>

#include "main.h"

class TransportSolver {

    double h, xa[N0], ya[N1], *u[NDIM];
    gsl_interp2d *interp;

    gsl_interp_accel *xacc, *yacc;

public:
    TransportSolver(const double h);
    void setU(gsl_vector *U[NDIM]);
    ~TransportSolver();

    void evaluateU(const double coord[NDIM], double u_interp[NDIM]);
    void evaluateUDeriv(const double coord[NDIM], double u_deriv[NDIM]);

    void transport(gsl_vector *S1, const gsl_vector *S0, const double dt);
};


#endif //STABLE_FLUIDS_TRANSPORTSOLVER_H
