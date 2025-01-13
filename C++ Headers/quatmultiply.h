#ifndef QUATMULTIPLY_H
#define QUATMULTIPLY_H

#include <array>
#include <vector>
#include <stdexcept>

/*
    Perform quaternion multiplication, replicating MATLAB's quatmultiply behavior.

    q1     - m-by-4 array of quaternions, where each quaternion is [w, x, y, z].
    q2     - m-by-4 array of quaternions, where each quaternion is [w, x, y, z].
    qmulti - m-by-4 array of quaternion products.
*/
std::array<double, 4> quatmultiply(const std::array<double, 4>& q1, const std::array<double, 4>& q2) 
{
    // Extract components of q1 and q2
    double w1 = q1[0], x1 = q1[1], y1 = q1[2], z1 = q1[3];
    double w2 = q2[0], x2 = q2[1], y2 = q2[2], z2 = q2[3];

    // Compute the product
    double w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2;
    double x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2;
    double y = w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2;
    double z = w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2;

    return {w, x, y, z};
}

std::vector<std::array<double, 4>> quatmultiply(const std::vector<std::array<double, 4>>& q1, const std::vector<std::array<double, 4>>& q2) 
{
    // Handle broadcasting if the number of rows differs
    if (q1.size() != q2.size()) 
    {
        throw std::invalid_argument("Input quaternion arrays must have the same number of elements.");
    }

    std::vector<std::array<double, 4>> qmulti;
    for (size_t i = 0; i < q1.size(); ++i) 
    {
        qmulti.push_back(quatmultiply(q1[i], q2[i]));
    }

    return qmulti;
}

#endif // QUATMULTIPLY_H