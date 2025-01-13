#ifndef QUAT_ROTATE_H
#define QUAT_ROTATE_H

#include <array>
#include "quatmultiply.h"
#include "quatinv.h"

/*
    quat_rotate rotates a vector by a unit quaternion.
    r = quat_rotate(q,v) calculates the rotated vector r for a
        quaternion q and a vector v.
    q   is a 1-by-4 matrix whose norm must be 1. q(1) is the scalar part
        of the quaternion.
    v   is a 1-by-3 matrix.
    r   is a 1-by-3 matrix.

    The 3-vector v is made into a pure quaternion 4-vector V = [0 v]. r is
    produced by the quaternion product R = q*V*qinv. r = [R(2) R(3) R(4)].
    
    User h-functions required: quatmultiply, quatinv.
*/
std::array<double, 3> quat_rotate(const std::array<double, 4>& q, const std::array<double, 3>& v) 
{
    std::array<double, 4> V = {0, v[0], v[1], v[2]};

    std::array<double, 4> qinv_result = quatinv(q);
    std::array<double, 4> intermediate = quatmultiply(q, V);
    std::array<double, 4> R = quatmultiply(intermediate, qinv_result);

    return {R[1], R[2], R[3]};
}

#endif // QUAT_ROTATE_H
