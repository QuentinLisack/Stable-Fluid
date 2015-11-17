//
// Created by Paul VANHAESEBROUCK on 10/11/2015.
//

#include "main.h"

#include "TransportSolver.h"
#include "DiffusionSolver.h"
#include "ProjectSolver.h"

#include "image.h"

#include <iostream>
#include <opencv2/highgui.hpp>
#include <time.h>

//vars used for the user interaction
int xf;
int yf;
time_t start;

int main(int argc, char **argv) {

	start = time(0);
	
    bool simulating = true;
    const double visc = 50.0, kS = 5.0, aS = 1e-3, h = 1.0, dt = 1.0;
    size_t x, y, d;

    DiffusionSolver uDiffSolver(visc, h, dt), sDiffSolver(kS, h, dt);
    TransportSolver ts(h);
    ProjectSolver projectSolver(h);


    // Vectors allocation
    gsl_vector *U0[NDIM], *U1[NDIM], *S0, *S1,
            *F[NDIM], *Ssource = gsl_vector_calloc(N_TOT);

    for (d = 0; d < NDIM; d++) {

        // calloc initialises with zeros
        U0[d] = gsl_vector_calloc(N_TOT);
        U1[d] = gsl_vector_calloc(N_TOT);

        F[d] = gsl_vector_calloc(N_TOT);
    }

    S0 = gsl_vector_calloc(N_TOT);
    S1 = gsl_vector_calloc(N_TOT);

	xf = -1;
	yf = -1;

    // Fictive inputs
    double radius = 50, x0 = 128, y0 = 128;
    for (double a = 0; a < 2 * M_PI; a += M_PI_4 / radius) {
        gsl_vector_set(F[0], _at(size_t(round(x0 + radius * cos(a))), size_t(round(y0 + radius * sin(a)))), - radius * sin(a));
        gsl_vector_set(F[1], _at(size_t(round(x0 + radius * cos(a))), size_t(round(y0 + radius * sin(a)))),   radius * cos(a));
    }

    gsl_vector_set(Ssource, _at(128, 165), 5);


    Image<double> result(N0, N1, CV_64F);

	// callbacks for the user interaction
	namedWindow("Result", 1);
	setMouseCallback("Result", CallBackFuncForce, F);
	setMouseCallback("Result", CallBackFuncSource, Ssource);
	
    // Main loop
    while (simulating) {

        // Swap vectors
        for (d = 0; d < NDIM; d++) {
            gsl_vector *temp = U0[d];
            U0[d] = U1[d];
            U1[d] = temp;
        }
//        gsl_vector *temp = S0;
//        S0 = S1;
//        S1 = temp;

        Vstep(U1, U0, F, dt, uDiffSolver, ts, projectSolver);
        Sstep(S1, S0, aS, Ssource, dt, sDiffSolver, ts);

        // TEST : keep it since we update the ource and the force with the user interaction
        gsl_vector_scale(Ssource, 0.8);
        gsl_vector_scale(F[0], 0.8);
        gsl_vector_scale(F[1], 0.8);

        // TODO: Read data and display it
        // ...
        for (y = 0; y < N1; y++) for (x = 0; x < N0; x++)
                result(x, y) = gsl_vector_get(S0, _at(x, y));

        imshow("Result", result.greyImage());
        waitKey(1);
    }

    return 0;
}

void CallBackFuncForce(int event, int x, int y, int flags, void* F)
{
	double forceCoeff = 1;
	double t = difftime(time(0), start);
     if ( event == EVENT_MOUSEMOVE && flags == EVENT_FLAG_CTRLKEY)
     {
		 if(xf >= 0 && yf >= 0){
          gsl_vector_set(F[0], _at(x, y), forceCoeff*(x-xf)/t);
		  gsl_vector_set(F[0], _at(x, y), forceCoeff*(y-yf)/t);
		 }
		 xf = x;
		 yf = y;
	}else if(event == EVENT_MOUSEMOVE){
		xf = -1;
		yf = -1;
		start = time(0);
	}
}

void CallBackFuncSource(int event, int x, int y, int flags, void* Ssource)
{
     if (event == EVENT_LBUTTONDOWN)
     {
		 gsl_vector_set(Ssource, _at(x, y), 5);
	}
}

/*
 * Step functions
 * Create a separate class (eg: Fluid) if you want to put them outside of this file
 */

void Vstep(gsl_vector *U1[], gsl_vector *U0[], gsl_vector *F[], const double dt, DiffusionSolver &diffSolver,
           TransportSolver &transpSolver, ProjectSolver &projectSolver) {

    size_t d;
    for (d = 0; d < NDIM; d++)
        addForce(U0[d], F[d]); // works if F is a velocity

    transpSolver.setU(U0);

    for (d = 0; d < NDIM; d++)
        transpSolver.transport(U1[d], U0[d], dt);
    for (d = 0; d < NDIM; d++)
        diffSolver.diffuse(U0[d], U1[d]);

    projectSolver.project(U1, U0);
}

void Sstep(gsl_vector *S1, // Vector to update
           gsl_vector *S0,
           const double aS,
           const gsl_vector *source, // Source vector
           const double dt,
           DiffusionSolver &diffSolver,
           TransportSolver &transpSolver) {

    addSource(S0, source);
    transpSolver.transport(S1, S0, dt);
    diffSolver.diffuse(S0, S1);
    dissipate(S0, dt, aS);
}


// function to add forces to velocity
void addForce(gsl_vector *U, const gsl_vector *F){
    gsl_vector_add(U, F);
}

void addSource(gsl_vector *S, const gsl_vector *source){
    gsl_vector_add(S, source);
}

void dissipate(gsl_vector *S, const double dt, const double a){
    gsl_vector_scale(S, 1 / (1 + a * dt));
}


// For debug purposes...

void printArray(double array[NDIM]) {
    for (size_t d = 0; d < NDIM; d++)
        std::cout << array[d] << " ";
    std::cout << std::endl;
}

void printVector(gsl_vector *vec) {
    size_t x, y;
    for (y = 0; y < N1; y++) {
        for (x = 0; x < N0; x++)
            std::cout << gsl_vector_get(vec, _at(x, y)) << "\t";
        std::cout << std::endl;
    }

    std::cout << std::endl;
}
