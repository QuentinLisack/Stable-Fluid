//
// Created by Paul VANHAESEBROUCK on 11/11/2015.
//

#include "TransportSolver.h"

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector_double.h>
#include <string.h>


TransportSolver::TransportSolver(const double h) : h(h) {

    // Data for U bilinear interpolation
    size_t x, y, d;

    const double *za = NULL;
    for (d = 0; d < NDIM; d++)
        u[d] = new double[N_TOT];

    for (x = 0; x < N0; x++)
        xa[x] = h * (0.5 + x);
    for (y = 0; y < N1; y++)
        ya[y] = h * (0.5 + y);

    const gsl_interp2d_type *T = gsl_interp2d_bilinear;
    interp = gsl_interp2d_alloc(T, N0, N1);

    gsl_interp2d_init(interp, xa, ya, za, N0, N1);

    xacc = gsl_interp_accel_alloc();
    yacc = gsl_interp_accel_alloc();
}

void TransportSolver::setU(gsl_vector *U[NDIM]) {

    size_t x, y, d;

    for (y = 0; y < N1; y++) for (x = 0; x < N0; x++) {
            const size_t idx = _at(x, y);
            for (d = 0; d < NDIM; d++)
                u[d][idx] = - gsl_vector_get(U[d], idx);
        }
}

void TransportSolver::evaluateU(const double coord[NDIM], double u_interp[NDIM]) {

    for (size_t d = 0; d < NDIM; d++) // Works only if NDIM == 2 for now
        u_interp[d] = gsl_interp2d_eval_extrap(interp, xa, ya, u[d], coord[0], coord[1], xacc, yacc);
}

void TransportSolver::evaluateUDeriv(const double coord[NDIM], double u_deriv[NDIM]) {

    // Works only if NDIM == 2 for now
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(u_deriv, 2, 2);
    gsl_matrix *m = &dfdy_mat.matrix;
    gsl_matrix_set(m, 0, 0, gsl_interp2d_eval_deriv_x(interp, xa, ya, u[0], coord[0], coord[1], xacc, yacc));
    gsl_matrix_set(m, 0, 1, gsl_interp2d_eval_deriv_y(interp, xa, ya, u[0], coord[0], coord[1], xacc, yacc));
    gsl_matrix_set(m, 1, 0, gsl_interp2d_eval_deriv_x(interp, xa, ya, u[1], coord[0], coord[1], xacc, yacc));
    gsl_matrix_set(m, 1, 1, gsl_interp2d_eval_deriv_y(interp, xa, ya, u[1], coord[0], coord[1], xacc, yacc));
}


// Differential equation system

int func(double t, const double y[], double f[], void *params) {
    TransportSolver *ts = (TransportSolver *) params;
    ts->evaluateU(y, f);
    return GSL_SUCCESS;
}

int jac(double t, const double y[], double *dfdy,
    double dfdt[], void *params) {
    TransportSolver *ts = (TransportSolver *) params;
    ts->evaluateUDeriv(y, dfdy);
    for (size_t d = 0; d < NDIM; dfdt[d++] = 0.0);
    return GSL_SUCCESS;
}


void TransportSolver::transport(gsl_vector *S1, const gsl_vector *S0, const double dt) {

    size_t x, y;

    gsl_odeiv2_system sys = {func, jac, NDIM, this};

    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4,
                                           dt/10, 1e-2, 1e-4);

    for (y = 0; y < N1; y++) for (x = 0; x < N0; x++) {
            double Y[NDIM] = { h * (0.5 + x), h * (0.5 + y) }, t = 0;

            int status = gsl_odeiv2_driver_apply (d, &t, dt, Y);

            if (status != GSL_SUCCESS)
            {
                printf ("ODE error, return value=%d\n", status);
                abort();
            }

            gsl_vector_set(S1, _at(x, y), gsl_interp2d_eval_extrap(interp, xa, ya, S0->data, Y[0], Y[1], xacc, yacc));

            gsl_odeiv2_driver_reset(d);
        }

    gsl_odeiv2_driver_free (d);
}


TransportSolver::~TransportSolver() {
    gsl_interp2d_free(interp);

    gsl_interp_accel_free(xacc);
    gsl_interp_accel_free(yacc);

    for (size_t d = 0; d < NDIM; d++)
        delete[] u[d];
}
