#include <iostream>
#include <iomanip>
#include <array>
#include "dcm_to_ypr.h"

/*
    This script tests the dcm_to_euler function using various input
    direction cosine matrices and verifies the output yaw, pitch, and roll angles.

    Q     - direction cosine matrix

    User h-function required: dcm_to_ypr
*/

void print_result(const double Q[3][3], double yaw, double pitch, double roll) 
{
    // Results
    std::cout << "Input DCM:\n";
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            std::cout << std::setw(10) << Q[i][j] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "Output Yaw, Pitch, and Roll Angles: yaw = " << yaw
              << ", pitch = " << pitch
              << ", roll = " << roll << " (degrees)\n\n";
}

int main() 
{
    // Input Matrices
    std::array<std::array<std::array<double, 3>, 3>, 5> input_matrices = 
    {{
        {{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}},      // Identity matrix (no rotation)
        {{{0.0, -1.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 0.0, 1.0}}},     // 90-degree rotation about Z-axis
        {{{0.0, 0.0, -1.0}, {0.0, 1.0, 0.0}, {1.0, 0.0, 0.0}}},     // 90-degree rotation about X-axis
        {{{0.0, 1.0, 0.0}, {-1.0, 0.0, 0.0}, {0.0, 0.0, 1.0}}},     // -90-degree rotation about Z-axis
        {{{std::cos(M_PI/4), -std::sin(M_PI/4), 0.0},               // 45-degree rotation about Z-axis
          {std::sin(M_PI/4),  std::cos(M_PI/4), 0.0},
          {0.0,              0.0,              1.0}}}
    }};

    double Q_arr[3][3];
    for (const auto& Q : input_matrices) 
    {
        Q_arr[0][0] = Q[0][0]; Q_arr[0][1] = Q[0][1]; Q_arr[0][2] = Q[0][2];
        Q_arr[1][0] = Q[1][0]; Q_arr[1][1] = Q[1][1]; Q_arr[1][2] = Q[1][2];
        Q_arr[2][0] = Q[2][0]; Q_arr[2][1] = Q[2][1]; Q_arr[2][2] = Q[2][2];

        // Compute Yaw, Pitch, and Roll Angles
        auto [yaw, pitch, roll] = dcm_to_ypr(Q_arr);
        print_result(Q_arr, yaw, pitch, roll);
    }

    return 0;
}
