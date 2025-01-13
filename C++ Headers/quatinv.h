#ifndef QUATINV_H
#define QUATINV_H

#include <array>
#include <vector>
#include <stdexcept>

/*
    Compute the inverse of a quaternion or an array of quaternions,
    replicating MATLAB's quatinv behavior.

    q    - m-by-4 array of quaternions, where each quaternion is [w, x, y, z].
    qinv - m-by-4 array of quaternion inverses.
*/
std::array<double, 4> quatinv(const std::array<double, 4>& q) 
{
    double w = q[0], x = q[1], y = q[2], z = q[3];
    double norm_squared = w * w + x * x + y * y + z * z;

    if (norm_squared == 0) 
    {
        throw std::invalid_argument("Quaternion norm is zero, cannot compute inverse.");
    }

    return {w / norm_squared, -x / norm_squared, -y / norm_squared, -z / norm_squared};
}

std::vector<std::array<double, 4>> quatinv(const std::vector<std::array<double, 4>>& quaternions) 
{
    std::vector<std::array<double, 4>> qinv;
    for (const auto& q : quaternions) 
    {
        qinv.push_back(quatinv(q));
    }

    return qinv;
}

#endif // QUATINV_H