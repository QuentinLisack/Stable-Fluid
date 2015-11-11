//
// Created by Paul VANHAESEBROUCK on 10/11/2015.
//

#include "DiffusionSolver.h"

#include <gsl/gsl_splinalg.h>

#include "main.h"


DiffusionSolver::DiffusionSolver(const double coeff, const double dt): coeff(coeff), dt(dt) {

    const double K = - coeff * dt;
    gsl_spmatrix *A = gsl_spmatrix_alloc(M*N, M*N); /* triplet format */

    /* construct the sparse matrix for the finite difference equation */

    for (size_t i = 0; i < M; i++) for (size_t j = 0; j < N; j++) {

            gsl_spmatrix_set(A, _at(i, j), _at(i, j), 1 - 4*K);
            if (i > 0)   gsl_spmatrix_set(A, _at(i, j), _at(i-1, j  ), K);
            if (j > 0)   gsl_spmatrix_set(A, _at(i, j), _at(i  , j-1), K);
            if (j < N-1) gsl_spmatrix_set(A, _at(i, j), _at(i  , j+1), K);
            if (i < M-1) gsl_spmatrix_set(A, _at(i, j), _at(i+1, j  ), K);
        }

    /* convert to compressed column format */
    C = gsl_spmatrix_compcol(A);

    gsl_spmatrix_free(A);
}

void DiffusionSolver::diffuse(gsl_vector *U1, const gsl_vector *U0) {

    const double tol = 1.0e-6;  /* solution relative tolerance */
    const size_t max_iter = 10; /* maximum iterations */
    const gsl_splinalg_itersolve_type *T = gsl_splinalg_itersolve_gmres;
    gsl_splinalg_itersolve *work = gsl_splinalg_itersolve_alloc(T, N*M, 0);
    size_t iter = 0;
    double residual;
    int status;
    size_t i;

    /* initial guess u = U1 */
    gsl_vector *u = gsl_vector_alloc(M*N);
    gsl_vector_memcpy(u, U0);

//    gsl_vector_set_zero(u);

    /* solve the system A u = f */
    do
    {
        status = gsl_splinalg_itersolve_iterate(C, U0, tol, u, work);

        /* print out residual norm ||A*u - f|| */
        residual = gsl_splinalg_itersolve_normr(work);
        fprintf(stderr, "iter %zu residual = %.12e\n", iter, residual);

        if (status == GSL_SUCCESS)
            fprintf(stderr, "Converged\n");
    }
    while (status == GSL_CONTINUE && ++iter < max_iter);

    /* output solution */
    gsl_vector_memcpy(U1, u);

    gsl_splinalg_itersolve_free(work);
    gsl_vector_free(u);
}

DiffusionSolver::~DiffusionSolver() {
    gsl_spmatrix_free(C);
}
