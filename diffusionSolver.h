//
// Created by Paul VANHAESEBROUCK on 10/11/2015.
//

#ifndef STABLE_FLUID_DIFFUSIONSOLVER_H
#define STABLE_FLUID_DIFFUSIONSOLVER_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_spmatrix.h>

class DiffusionSolver {

    const double coeff, h, dt;
    gsl_spmatrix *C;

public:

    DiffusionSolver(const double coeff, const double h, const double dt);
    ~DiffusionSolver();

    void diffuse(gsl_vector *U1, const gsl_vector *U0);
};


#endif //STABLE_FLUID_DIFFUSIONSOLVER_H
