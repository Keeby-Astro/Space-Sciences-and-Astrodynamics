#include <iostream>
#include <array>
#include <cmath>
#include "quat_rotate.h"

void quat_rotate_test() 
{
    // Define test quaternions and vectors
    std::array<double, 4> q1 = {sqrt(2) / 2, 0, 0, sqrt(2) / 2}; // 90-degree Z-rotation
    std::array<double, 4> q2 = {0, 0, 1, 0};                     // 180-degree Y-rotation
    std::array<double, 4> q3 = {1, 0, 0, 0};                     // Identity quaternion
    std::array<double, 4> q4 = {sqrt(2) / 2, sqrt(2) / 2, 0, 0}; // 90-degree X-rotation

    // Vector to rotate
    std::array<double, 3> v1 = {1, 0, 0};
    std::array<double, 3> v2 = {1, 0, 0};
    std::array<double, 3> v3 = {1, 2, 3};
    std::array<double, 3> v4 = {0, 1, 0};

    // Rotated vector
    std::array<double, 3> r1 = quat_rotate(q1, v1);
    std::array<double, 3> r2 = quat_rotate(q2, v2);
    std::array<double, 3> r3 = quat_rotate(q3, v3);
    std::array<double, 3> r4 = quat_rotate(q4, v4);

    // Output results
    std::cout << "90-degree rotation around Z-axis\n";
    std::cout << "Input Vector: [" << v1[0] << ", " << v1[1] << ", " << v1[2] << "]\n";
    std::cout << "Rotated Vector: [" << r1[0] << ", " << r1[1] << ", " << r1[2] << "]\n\n";

    std::cout << "180-degree rotation around Y-axis\n";
    std::cout << "Input Vector: [" << v2[0] << ", " << v2[1] << ", " << v2[2] << "]\n";
    std::cout << "Rotated Vector: [" << r2[0] << ", " << r2[1] << ", " << r2[2] << "]\n\n";

    std::cout << "No rotation\n";
    std::cout << "Input Vector: [" << v3[0] << ", " << v3[1] << ", " << v3[2] << "]\n";
    std::cout << "Rotated Vector: [" << r3[0] << ", " << r3[1] << ", " << r3[2] << "]\n\n";

    std::cout << "90-degree rotation around X-axis\n";
    std::cout << "Input Vector: [" << v4[0] << ", " << v4[1] << ", " << v4[2] << "]\n";
    std::cout << "Rotated Vector: [" << r4[0] << ", " << r4[1] << ", " << r4[2] << "]\n\n";
}

int main() {
    quat_rotate_test();
    return 0;
}
