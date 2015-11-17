//
// Created by Paul VANHAESEBROUCK on 17/11/2015.
//

#include <gsl/gsl_splinalg.h>
#include "ProjectSolver.h"

ProjectSolver::ProjectSolver(const double h) : h(h) {

    gsl_spmatrix *A = gsl_spmatrix_alloc(N_TOT, N_TOT); /* triplet format */

    /* construct the sparse matrix for the finite difference equation */

    for (size_t y = 0; y < N1; y++) for (size_t x = 0; x < N0; x++) {

            gsl_spmatrix_set(A, _at(x, y), _at(x, y), - 4);
            if (x > 0)     gsl_spmatrix_set(A, _at(x, y), _at(x - 1, y    ), 1);
            if (y > 0)     gsl_spmatrix_set(A, _at(x, y), _at(x    , y - 1), 1);
            if (y < N1-1)  gsl_spmatrix_set(A, _at(x, y), _at(x    , y + 1), 1);
            if (x < N0-1)  gsl_spmatrix_set(A, _at(x, y), _at(x + 1, y    ), 1);
        }

    /* convert to compressed column format */

    C = gsl_spmatrix_compcol(A);

    gsl_spmatrix_free(A);
}


void ProjectSolver::project(gsl_vector **U1, gsl_vector **U0) {

    const double tol = 1.0e-1; /* solution relative tolerance */
    const size_t max_iter = 5; /* maximum iterations */
    const gsl_splinalg_itersolve_type *T = gsl_splinalg_itersolve_gmres;

    gsl_vector *UV = gsl_vector_alloc(N_TOT);
    gsl_vector *u = gsl_vector_alloc(N_TOT);

    for (size_t d = 0; d < NDIM; d++) {

        div(UV, U0[d]);

        gsl_splinalg_itersolve *work = gsl_splinalg_itersolve_alloc(T, N_TOT, 0);
        size_t iter = 0;
        int status;

        /* initial guess u = U1 */
        gsl_vector_memcpy(u, U0[d]);

        /* solve the system A u = f */
        do
            status = gsl_splinalg_itersolve_iterate(C, UV, tol, u, work);
        while (status == GSL_CONTINUE && ++iter < max_iter);

        /* output solution */

        for (size_t y = 0; y < N1; y++) for (size_t x = 0; x < N0; x++) {
                double temp = gsl_vector_get(U0[d], _at(x, y));

                if (d == 0 && x > 0 && x < N0-1) temp -= (gsl_vector_get(u, _at(x + 1, y    )) - gsl_vector_get(u, _at(x - 1, y    ))) * 0.5 / h;
                if (d == 1 && y > 0 && y < N1-1) temp -= (gsl_vector_get(u, _at(x    , y + 1)) - gsl_vector_get(u, _at(x    , y - 1))) * 0.5 / h;

                gsl_vector_set(U1[d], _at(x, y), temp);
            }

        gsl_splinalg_itersolve_free(work);
    }

    gsl_vector_free(u);
    gsl_vector_free(UV);
}


void ProjectSolver::div(gsl_vector *UV, gsl_vector *U) {

    for (size_t y = 0; y < N1; y++) for (size_t x = 0; x < N0; x++) {

            double temp = 0;
            if (x > 0 && x < N0-1) temp += (gsl_vector_get(U, _at(x + 1, y    )) - gsl_vector_get(U, _at(x - 1, y    ))) * h * 0.5;
            if (y > 0 && y < N1-1) temp += (gsl_vector_get(U, _at(x    , y + 1)) - gsl_vector_get(U, _at(x    , y - 1))) * h * 0.5;

            gsl_vector_set(UV, _at(x, y), temp);
        }
}

ProjectSolver::~ProjectSolver() {
    gsl_spmatrix_free(C);
}

