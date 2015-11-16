//
// Created by Paul VANHAESEBROUCK on 10/11/2015.
//

#include "main.h"

#include <iostream>


int main(int argc, char **argv) {

    bool simulating = true;
    const double visc = 1., kS = 1., aS = 1., h = 1., dt = 0.1;

    DiffusionSolver uDiffSolver(visc, dt), sDiffSolver(kS, dt);

    // Vectors allocation
    gsl_vector *U0[NDIM], *U1[NDIM], *S0, *S1,
            *F[NDIM], *Ssource = gsl_vector_calloc(M * N);

    for (size_t i = 0; i < NDIM; i++) {

        // calloc initialises with zeros
        U0[i] = gsl_vector_calloc(M * N);
        U1[i] = gsl_vector_calloc(M * N);

        F[i] = gsl_vector_calloc(M * N);
    }

    S0 = gsl_vector_calloc(M * N);
    S1 = gsl_vector_calloc(M * N);

    // Main loop
//    while (simulating) {

    // Swap vectors
    for (size_t i = 0; i < NDIM; i++) {
        gsl_vector *temp = U0[i];
        U0[i] = U1[i];
        U1[i] = temp;
    }
    gsl_vector *temp = S0;
    S0 = S1;
    S1 = temp;

    TransportSolver transpSolver(U0);

    Vstep(U1, U0, visc, F, dt, uDiffSolver, transpSolver);
    Sstep(S1, S0, kS, aS, U1, Ssource, dt, sDiffSolver, transpSolver);

    // TODO: Read data and display it, update forces and sources
    // ...
//    }

    return 0;
}


/*
 * Step functions
 * Create a separate class (eg: Fluid) if you want to put them outside of this file
 */

void Vstep(gsl_vector *U1[NDIM], // Vector to update
           gsl_vector *U0[NDIM],
           const double visc,
           gsl_vector *F[NDIM], // Force vector
           const double dt,
           const DiffusionSolver &diffSolver,
           const TransportSolver &transpSover) {

    size_t i;
    for (i = 0; i < NDIM; i++)
       addForce(U0[i], F[i]); // works if F is a velocity
    for (i = 0; i < NDIM; i++)
        transpSolver.transport(U1[i], U0[i], dt);
    for (i = 0; i < NDIM; i++)
        uDiffSolver.diffuse(U1[i], U0[i]);
	
	//TODO
		//Project !


}

// function to add forces to velocity
void addForce(gsl_vector *U, gsl_vector *F){
	gsl_vector_add(U, F);
}

void addSources(gsl_vector *S, gsl_vector source, const double dt){
	gsl_vector_scale(source, dt);
	gsl_vector_add(S, source);
	// this function can be used only if we update the source between each iteration. 
	// else we do not need to scale the source (just do it in the initialisation)
}

void dissipate(gsl_vector *S, const double dt, const double a){
	gsl_vector_scale(S, 1 / (1 + aS * dt));
}

void Sstep(gsl_vector *S1, // Vector to update
           gsl_vector *S0,
           const double kS,
           const double aS,
           gsl_vector *U1[NDIM],
           const gsl_vector *source, // Source vector
           const double dt,
           const DiffusionSolver &diffSolver,
           const TransportSolver &transpSover) {

    // TODO
	size_t i;
	addForce(S1, source, dt);
	transpSolver.transport(S1, S0, dt);
	sDiffSolver.diffuse(S1, S0);	
	dissipate(S1, dt, aS);
