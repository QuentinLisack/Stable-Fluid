//
// Created by Paul VANHAESEBROUCK on 10/11/2015.
//

#include "DiffusionSolver.h"

#include "main.h"

#include <gsl/gsl_splinalg.h>
#include <gsl/gsl_blas.h>


DiffusionSolver::DiffusionSolver(const double coeff, const double h, const double dt): coeff(coeff), h(h), dt(dt) {

    const double K = - coeff * dt / (h * h);
    gsl_spmatrix *A = gsl_spmatrix_alloc(N_TOT, N_TOT); /* triplet format */

    /* construct the sparse matrix for the finite difference equation */

    for (size_t y = 0; y < N1; y++) for (size_t x = 0; x < N0; x++) {

            gsl_spmatrix_set(A, _at(x, y), _at(x, y), 1 - 4 * K);
            if (x > 0)     gsl_spmatrix_set(A, _at(x, y), _at(x - 1, y    ), K);
            if (y > 0)     gsl_spmatrix_set(A, _at(x, y), _at(x    , y - 1), K);
            if (y < N1-1)  gsl_spmatrix_set(A, _at(x, y), _at(x    , y + 1), K);
            if (x < N0-1)  gsl_spmatrix_set(A, _at(x, y), _at(x + 1, y    ), K);
        }

    /* convert to compressed column format */

    C = gsl_spmatrix_compcol(A);

    gsl_spmatrix_free(A);
}

void DiffusionSolver::diffuse(gsl_vector *U1, const gsl_vector *U0) {

    const double tol = 1.0e-1; /* solution relative tolerance */
    const size_t max_iter = 5; /* maximum iterations */
    const gsl_splinalg_itersolve_type *T = gsl_splinalg_itersolve_gmres;
    gsl_splinalg_itersolve *work = gsl_splinalg_itersolve_alloc(T, N_TOT, 0);
    size_t iter = 0;
    int status;

    /* initial guess u = U1 */
    gsl_blas_dcopy(U0, U1);

    /* solve the system A u = f */
    do
        status = gsl_splinalg_itersolve_iterate(C, U0, tol, U1, work);
    while (status == GSL_CONTINUE && ++iter < max_iter);


    gsl_splinalg_itersolve_free(work);
}

DiffusionSolver::~DiffusionSolver() {
    gsl_spmatrix_free(C);
}
