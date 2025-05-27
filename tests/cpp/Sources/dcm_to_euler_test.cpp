#include <iostream>
#include <vector>
#include <iomanip>
#include "dcm_to_euler.h"

/*
    This script tests the dcm_to_euler function using various input
    direction cosine matrices and verifies the output Euler angles.

    Q     - direction cosine matrix

    User h-function required: dcm_to_euler
*/
void print_matrix(const std::vector<std::vector<double>> &matrix) 
{
    for (const auto &row : matrix) 
    {
        for (const auto &value : row) 
        {
            std::cout << std::setw(10) << value << " ";
        }
        std::cout << "\n";
    }
}

int main() {
    // Input Matrices
    std::vector<std::vector<std::vector<double>>> input_matrices = 
    {
        {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}},              // Identity matrix (no rotation)
        {{0, -1, 0}, {1, 0, 0}, {0, 0, 1}},             // 90-degree rotation about Z-axis
        {{0, 0, -1}, {0, 1, 0}, {1, 0, 0}},             // 90-degree rotation about X-axis
        {{0, 1, 0}, {-1, 0, 0}, {0, 0, 1}},             // -90-degree rotation about Z-axis
        {{std::cos(M_PI / 4), -std::sin(M_PI / 4), 0},  // 45-degree rotation about Z-axis
         {std::sin(M_PI / 4), std::cos(M_PI / 4), 0},
         {0, 0, 1}}
    };

    for (const auto &Q : input_matrices) 
    {
        // Compute Euler Angles
        auto [alpha, beta, gamma] = dcm_to_euler(Q);

        // Results
        std::cout << "Input DCM:\n";
        print_matrix(Q);
        std::cout << "Output Euler Angles: alpha = " << std::fixed << std::setprecision(2)
                  << alpha << ", beta = " << beta << ", gamma = " << gamma << " (degrees)\n\n";
    }

    return 0;
}
