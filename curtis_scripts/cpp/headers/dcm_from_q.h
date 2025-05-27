#ifndef DCM_FROM_Q_H
#define DCM_FROM_Q_H

#include <array>

/*
    This function calculates the direction cosine matrix
    from the quaternion.

    q - quaternion (where q(4) is the scalar part)
    Q - direction cosine matrix
*/
std::array<std::array<double, 3>, 3> dcm_from_q(const std::array<double, 4>& q) 
{
    double q1 = q[0], q2 = q[1], q3 = q[2], q4 = q[3];

    std::array<std::array<double, 3>, 3> Q = 
    {{
        {q1 * q1 - q2 * q2 - q3 * q3 + q4 * q4, 2 * (q1 * q2 + q3 * q4), 2 * (q1 * q3 - q2 * q4)},
        {2 * (q1 * q2 - q3 * q4), -q1 * q1 + q2 * q2 - q3 * q3 + q4 * q4, 2 * (q2 * q3 + q1 * q4)},
        {2 * (q1 * q3 + q2 * q4), 2 * (q2 * q3 - q1 * q4), -q1 * q1 - q2 * q2 + q3 * q3 + q4 * q4}
    }};

    return Q;
}

#endif // DCM_FROM_Q_H