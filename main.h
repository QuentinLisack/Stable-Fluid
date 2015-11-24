//
// Created by Paul VANHAESEBROUCK on 10/11/2015.
//

#ifndef STABLE_FLUIDS_MAIN_H
#define STABLE_FLUIDS_MAIN_H

#include <gsl/gsl_vector.h>

class TransportSolver;
class DiffusionSolver;
class ProjectSolver;

/*
 * Be careful to check everything here if you modify it.
 */

const size_t NDIM = 2, LOG2N = 8;
const size_t NS = 3;
//8 // POWER OF 2 FOR SPACE DIMENSION
const size_t N0 = (1 << LOG2N);
const size_t N1 = N0;
const size_t N_TOT = N0 * N1;
#define _at(x, y) (((y)<<LOG2N) + (x)) // FAST CONVERSION TO VECTOR INDEX


void Vstep(gsl_vector *U1[],
           gsl_vector *U0[],
           gsl_vector *F[],
           const double dt,
           DiffusionSolver &diffSolver,
           TransportSolver &transpSolver,
           ProjectSolver &projectSolver);

void Sstep(gsl_vector *S1,
           gsl_vector *S0,
           const double aS,
           gsl_vector *source,
           const double dt,
           DiffusionSolver &diffSolver,
           TransportSolver &transpSolver);
		   
inline void addForce(gsl_vector *U, const gsl_vector *F);

inline void addSource(gsl_vector *S, const gsl_vector *source);

inline void dissipate(gsl_vector *S, const double dt, const double a);

// For debug purposes...

void printArray(double array[NDIM]);
void printVector(gsl_vector *vec);

#endif //STABLE_FLUIDS_MAIN_H
