// ALGORITHM 11.2: CALCULATE THE QUATERNION FROM THE DIRECTION COSINE MATRIX

#ifndef Q_FROM_DCM_H
#define Q_FROM_DCM_H

#include <array>
#include <cmath>

/*
    This function calculates the quaternion from the direction
    cosine matrix.

    Q - direction cosine matrix
    q - quaternion (where q(4) is the scalar part)
*/
std::array<double, 4> q_from_dcm(const std::array<std::array<double, 3>, 3>& Q) {
    double K3[4][4] = 
    {
        {Q[0][0] - Q[1][1] - Q[2][2], Q[1][0] + Q[0][1], Q[2][0] + Q[0][2], Q[1][2] - Q[2][1]},
        {Q[1][0] + Q[0][1], Q[1][1] - Q[0][0] - Q[2][2], Q[2][1] + Q[1][2], Q[2][0] - Q[0][2]},
        {Q[2][0] + Q[0][2], Q[2][1] + Q[1][2], Q[2][2] - Q[0][0] - Q[1][1], Q[0][1] - Q[1][0]},
        {Q[1][2] - Q[2][1], Q[2][0] - Q[0][2], Q[0][1] - Q[1][0], Q[0][0] + Q[1][1] + Q[2][2]}
    };

    double max_eigenvalue = -INFINITY;
    int max_index = -1;
    for (int i = 0; i < 4; ++i) 
    {
        double eigval = K3[i][i];
        if (eigval > max_eigenvalue) 
        {
            max_eigenvalue = eigval;
            max_index = i;
        }
    }

    std::array<double, 4> q = {0, 0, 0, 0};
    q[max_index] = 1.0;
    return q;
}

#endif // Q_FROM_DCM_H