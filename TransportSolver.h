//
// Created by Paul VANHAESEBROUCK on 11/11/2015.
//

#ifndef STABLE_FLUIDS_TRANSPORTSOLVER_H
#define STABLE_FLUIDS_TRANSPORTSOLVER_H

#include <gsl/gsl_vector.h>

class TransportSolver {

    gsl_vector **U;

    void evaluateU(const double coord[], double u[]);

public:
    TransportSolver(gsl_vector *U[]);
    ~TransportSolver();

    void transport(gsl_vector *S1, const gsl_vector *S0, const double dt);
};


#endif //STABLE_FLUIDS_TRANSPORTSOLVER_H
