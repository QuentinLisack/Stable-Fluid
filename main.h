//
// Created by Paul VANHAESEBROUCK on 10/11/2015.
//

#ifndef STABLE_FLUIDS_MAIN_H
#define STABLE_FLUIDS_MAIN_H

#include <stddef.h>
#include <gsl/gsl_vector.h>

#include "DiffusionSolver.h"
#include "TransportSolver.h"

/*
 * Be careful to check everything here if you modify it.
 */

const size_t NDIM = 2, LOG2N = 2;
//8 // POWER OF 2 FOR SPACE DIMENSION
const size_t N = (1 << LOG2N);
const size_t M = N;
#define _at(i, j) (((i)<<LOG2N) + (j)) // FAST CONVERSION TO VECTOR INDEX


void Vstep(gsl_vector *U1[NDIM], gsl_vector *U0[NDIM], const double visc, gsl_vector *F[NDIM], const double dt,
           const DiffusionSolver &diffSolver, const TransportSolver &transpSover);

void Sstep(gsl_vector *S1, gsl_vector *S0, const double kS, const double aS, gsl_vector *U1[NDIM],
           const gsl_vector *source, const double dt, const DiffusionSolver &diffSolver,
           const TransportSolver &transpSover);

#endif //STABLE_FLUIDS_MAIN_H
