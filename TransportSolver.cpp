//
// Created by Paul VANHAESEBROUCK on 11/11/2015.
//

#include "TransportSolver.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include "main.h"


TransportSolver::TransportSolver(gsl_vector *U[NDIM]) : U(U) {

    // Interpolate U
    // We can use a spline
    // But now we will just use barycentric coordinates
}

void TransportSolver::evaluateU(const double coord[NDIM], double u[NDIM]) {

    // Find in which square we are

    int o_i = (int) (coord[0] - 0.5), o_j = (int) (coord[1] - 0.5);
    double coeff_i = coord[0] - 1. - o_i, coeff_j = coord[1] - 1. - o_j;

    // Linear interpolation

    u[0] = 0;
    u[1] = 0;

    for (size_t i = 0; i != 2; i++) for (size_t j = 0; j != 2; j++) {
            const double c_i = 0.5 + ((i == 0) ? 1. : -1.) * coeff_i, c_j = 0.5 + ((j == 0) ? 1. : -1.) * coeff_j;
            if (c_i >= 0. && c_i <= 1. && c_j >= 0. && c_j <= 1.) {
                u[0] += c_i * c_j * gsl_vector_get(U[0], _at(o_i + i, o_j + j));
                u[1] += c_i * c_j * gsl_vector_get(U[1], _at(o_i + i, o_j + j));
            }
        }
}



// Differential equation system

//int
//func (double t, const double y[], double f[],
//      void *params)
//{
//    double mu = *(double *)params;
//    f[0] = y[1];
//    f[1] = -y[0] - mu*y[1]*(y[0]*y[0] - 1);
//    return GSL_SUCCESS;
//}
//
//int
//jac (double t, const double y[], double *dfdy,
//     double dfdt[], void *params)
//{
//    double mu = *(double *)params;
//    gsl_matrix_view dfdy_mat
//            = gsl_matrix_view_array (dfdy, 2, 2);
//    gsl_matrix * m = &dfdy_mat.matrix;
//    gsl_matrix_set (m, 0, 0, 0.0);
//    gsl_matrix_set (m, 0, 1, 1.0);
//    gsl_matrix_set (m, 1, 0, -2.0*mu*y[0]*y[1] - 1.0);
//    gsl_matrix_set (m, 1, 1, -mu*(y[0]*y[0] - 1.0));
//    dfdt[0] = 0.0;
//    dfdt[1] = 0.0;
//    return GSL_SUCCESS;
//}


void TransportSolver::transport(gsl_vector *S1, const gsl_vector *S0, const double dt) {

}


TransportSolver::~TransportSolver() {

}
