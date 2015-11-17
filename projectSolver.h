//
// Created by Paul VANHAESEBROUCK on 17/11/2015.
//

#ifndef STABLE_FLUIDS_PROJECTSOLVER_H
#define STABLE_FLUIDS_PROJECTSOLVER_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_spmatrix.h>

#include "main.h"

class ProjectSolver {

    const double h;
    gsl_spmatrix *UV, *C;

public:

    ProjectSolver(const double h);

    ~ProjectSolver();

    void div(gsl_vector *UV, gsl_vector *U);
    void project(gsl_vector *U1[NDIM], gsl_vector *U0[NDIM]);
};


#endif //STABLE_FLUIDS_PROJECTSOLVER_H
