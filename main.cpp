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

#include <boost/thread.hpp>
#include <chrono>

#define WINDOW_NAME "Result"

//vars used for the user interaction


typedef struct {
    gsl_vector **F;
    gsl_vector *S;
    boost::mutex F_mutex, S_mutex;
} Data;

volatile bool simulating = true;

const double forceCoeff = 100.0, sourceCoeff = 5.0;

int xf = -1, yf = -1;
chrono::system_clock::time_point start;

void mouseCallback(int event, int x, int y, int flags, void *ptr) {

    Data *data = (Data*) ptr;
    gsl_vector **F = data->F, *Ssource = data->S;

    if (event == EVENT_MOUSEMOVE && flags & EVENT_FLAG_CTRLKEY) {

        // Measure time elapsed since last event
        chrono::system_clock::time_point now = chrono::system_clock::now();
        double t = (chrono::duration_cast<chrono::duration<double> >(now - start)).count();
        start = now;

        // Add velocity to F
        if (xf >= 0 && yf >= 0 && t > 0.0) {
            data->F_mutex.lock();
            gsl_vector_set(F[0], (size_t) _at(x, y), forceCoeff * (x - xf) / t);
            gsl_vector_set(F[1], (size_t) _at(x, y), forceCoeff * (y - yf) / t);
            data->F_mutex.unlock();
        }

        xf = x;
        yf = y;
    } else {
        xf = -1;
        yf = -1;
    }

    if (event == EVENT_LBUTTONDOWN || (event == EVENT_MOUSEMOVE && flags & EVENT_FLAG_LBUTTON)) {

        // Add source to Ssource
        data->S_mutex.lock();
        gsl_vector_set(Ssource, (size_t) _at(x, y), sourceCoeff);
        data->S_mutex.unlock();
    }
}

void getInput(Data *data, gsl_vector **F, gsl_vector *Ssource) {
    data->F_mutex.lock();
    for (size_t d = 0; d < NDIM; d++) {
            gsl_vector_add(F[d], data->F[d]);
            gsl_vector_set_zero(data->F[d]);
        }
    data->F_mutex.unlock();
    data->S_mutex.lock();
    gsl_vector_add(Ssource, data->S);
    gsl_vector_set_zero(data->S);
    data->S_mutex.unlock();
}

void task(Data *data) {

    gsl_vector *F[NDIM], *Ssource;

    // Space parameters
    const double visc = 25.0, kS = 5.0, aS = 1e-3, h = 1.0, dt = 1.0;
    size_t x, y, d;

    DiffusionSolver uDiffSolver(visc, h, dt), sDiffSolver(kS, h, dt);
    TransportSolver ts(h);
    ProjectSolver projectSolver(h);

    // Vectors allocation
    gsl_vector *U0[NDIM], *U1[NDIM], *S0, *S1;

    for (d = 0; d < NDIM; d++) {

        // calloc initialises with zeros
        U0[d] = gsl_vector_calloc(N_TOT);
        U1[d] = gsl_vector_calloc(N_TOT);
        F[d] = gsl_vector_calloc(N_TOT);
    }

    S0 = gsl_vector_calloc(N_TOT);
    S1 = gsl_vector_calloc(N_TOT);
    Ssource = gsl_vector_calloc(N_TOT);

    // Output image
    Image<double> result(N0, N1, CV_64F);

    // Main loop
    while (simulating) {

        // Swap velocity vectors
        for (d = 0; d < NDIM; d++) {
            gsl_vector *temp = U0[d];
            U0[d] = U1[d];
            U1[d] = temp;
        }

        getInput(data, F, Ssource);

        Vstep(U1, U0, F, dt, uDiffSolver, ts, projectSolver);
        Sstep(S1, S0, aS, Ssource, dt, sDiffSolver, ts);

        // TEST : keep it since we update the source and the force with the user interaction
        gsl_vector_scale(Ssource, 0.8);
        gsl_vector_scale(F[0], 0.8);
        gsl_vector_scale(F[1], 0.8);

        // TODO: Improve result rendering

        for (y = 0; y < N1; y++) for (x = 0; x < N0; x++)
                result(x, y) = gsl_vector_get(S0, _at(x, y));

        imshow(WINDOW_NAME, result.greyImage());
        waitKey(1);
    }
}


int main(int argc, char **argv) {

    // Input vectors
    gsl_vector *F[NDIM], *Ssource;

    Ssource = gsl_vector_calloc(N_TOT);
    for (size_t d = 0; d < NDIM; d++)
        F[d] = gsl_vector_calloc(N_TOT);

    // Data structure shared between threads
    Data data = {F, Ssource};

    // callbacks for the user interaction
    namedWindow(WINDOW_NAME);
    setMouseCallback(WINDOW_NAME, mouseCallback, &data);

    // Start computing thread
    boost::thread t(task, &data);

    // Stop when user presses Escape
    while(waitKey(0) != 27);

    simulating = false;

    // End computing thread
    t.join();

    return 0;
}


/*
 * Step functions
 * Create a separate class (eg: Fluid) if you want to put them outside of this file
 */

void Vstep(gsl_vector *U1[],
           gsl_vector *U0[],
           gsl_vector *F[],
           const double dt,
           DiffusionSolver &diffSolver,
           TransportSolver &transpSolver,
           ProjectSolver &projectSolver) {

    size_t d;
    for (d = 0; d < NDIM; d++)
        addForce(U0[d], F[d]); // works if F is a velocity
    transpSolver.setU(U0); // Update transpSolver with this velocity
    for (d = 0; d < NDIM; d++)
        transpSolver.transport(U1[d], U0[d], dt);
    for (d = 0; d < NDIM; d++)
        diffSolver.diffuse(U0[d], U1[d]);

    projectSolver.project(U1, U0);
}

void Sstep(gsl_vector *S1, // Vector to update
           gsl_vector *S0,
           const double aS,
           gsl_vector *source, // Source vector
           const double dt,
           DiffusionSolver &diffSolver,
           TransportSolver &transpSolver) {

    addSource(S0, source);
    transpSolver.transport(S1, S0, dt);
    diffSolver.diffuse(S0, S1);
    dissipate(S0, dt, aS);
}


// function to add forces to velocity
inline void addForce(gsl_vector *U, const gsl_vector *F) {
    gsl_vector_add(U, F);
}

inline void addSource(gsl_vector *S, const gsl_vector *source) {
    gsl_vector_add(S, source);
}

inline void dissipate(gsl_vector *S, const double dt, const double a) {
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
