//
// Created by Paul VANHAESEBROUCK on 10/11/2015.
//

#ifndef STABLE_FLUIDS_MAIN_H
#define STABLE_FLUIDS_MAIN_H

#include <stddef.h>
#include <gsl/gsl_vector.h>

class TransportSolver;
class DiffusionSolver;

/*
 * Be careful to check everything here if you modify it.
 */

const size_t NDIM = 2, LOG2N = 8;
//8 // POWER OF 2 FOR SPACE DIMENSION
const size_t N0 = (1 << LOG2N);
const size_t N1 = N0;
const size_t N_TOT = N0 * N1;
#define _at(x, y) (((y)<<LOG2N) + (x)) // FAST CONVERSION TO VECTOR INDEX


void Vstep(gsl_vector *U1[NDIM], // Vector to update
           gsl_vector *U0[NDIM],
           gsl_vector *F[NDIM], // Force vector
           const double dt,
           DiffusionSolver &diffSolver,
           TransportSolver &transpSolver);

void Sstep(gsl_vector *S1, // Vector to update
           gsl_vector *S0,
           const double aS,
           const gsl_vector *source, // Source vector
           const double dt,
           DiffusionSolver &diffSolver,
           TransportSolver &transpSolver);
		   
void addForce(gsl_vector *U, const gsl_vector *F);

void addSource(gsl_vector *S, const gsl_vector *source);

void dissipate(gsl_vector *S, const double dt, const double a);

// For debug purposes...

void printArray(double array[NDIM]);
void printVector(gsl_vector *vec);

#endif //STABLE_FLUIDS_MAIN_H
