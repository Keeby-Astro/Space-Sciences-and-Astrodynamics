// ALGORITHM 4.3: OBTAIN THE CLASSICAL EULER ANGLE SEQUENCE
// FROM A DIRECTION COSINE MATRIX

#ifndef DCM_TO_EULER_H
#define DCM_TO_EULER_H

#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <tuple>
#include <vector>
#include "atan2d_0_360.h"

/*
    This function finds the angles of the classical Euler sequence
    R3(gamma)*R1(beta)*R3(alpha) from the direction cosine matrix.

    Q     - direction cosine matrix
    alpha - first angle of the sequence (deg)
    beta  - second angle of the sequence (deg)
    gamma - third angle of the sequence (deg)

    User h-function required: atan2d_0_360
*/

inline std::tuple<double, double, double> dcm_to_euler(const std::vector<std::vector<double>> &Q) {
    double alpha = atan2d_0_360(Q[2][0], -Q[2][1]);
    double beta = std::acos(Q[2][2]) * 180.0 / M_PI;
    double gamma = atan2d_0_360(Q[0][2], Q[1][2]);

    return std::make_tuple(alpha, beta, gamma);
}

#endif // DCM_TO_EULER_H
