#ifndef DCM_TO_YPR_H
#define DCM_TO_YPR_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <tuple>
#include "atan2d_0_360.h"

/*
    This function finds the angles of the yaw-pitch-roll sequence
    R1(gamma)*R2(beta)*R3(alpha) from the direction cosine matrix.

    Q     - direction cosine matrix
    yaw   - yaw angle (deg)
    pitch - pitch angle (deg)
    roll  - roll angle (deg)

    User h-function required: atan2d_0_360
*/
inline std::tuple<double, double, double> dcm_to_ypr(const double Q[3][3]) {

    double yaw = atan2d_0_360(Q[0][1], Q[0][0]);
    double pitch = std::asin(-Q[0][2]) * 180.0 / M_PI;
    double roll = atan2d_0_360(Q[1][2], Q[2][2]);

    return std::make_tuple(yaw, pitch, roll);
}

#endif // DCM_TO_YPR_H
