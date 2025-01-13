#include <iostream>
#include <iomanip>
#include "dcm_from_q.h"

/*
    Verifies the direction cosine matrix from known quaternions.

    q - quaternion (where q(4) is the scalar part)
    Q - direction cosine matrix

    User h-functions required: dcm_from_q
*/
void dcm_from_q_test() 
{
    std::array<double, 4> q1 = {1, 0, 0, 0};  // Identity quaternion
    std::array<double, 4> q2 = {0, 1, 0, 0};  // 180-degree rotation about x-axis
    std::array<double, 4> q3 = {0, 0, 1, 0};  // 180-degree rotation about y-axis
    std::array<double, 4> q4 = {0, 0, 0, 1};  // 180-degree rotation about z-axis

    // Calculate the direction cosine matrices using the function
    auto Q1 = dcm_from_q(q1);
    auto Q2 = dcm_from_q(q2);
    auto Q3 = dcm_from_q(q3);
    auto Q4 = dcm_from_q(q4);

    auto display_matrix = [](const std::array<std::array<double, 3>, 3>& Q) 
    {
        for (const auto& row : Q) 
        {
            for (double value : row) {
                std::cout << std::fixed << std::setprecision(6) << value << " ";
            }
            std::cout << "\n";
        }
    };

    // Display the results
    std::cout << "Results for quaternion q1: [1, 0, 0, 0]\n";
    std::cout << "Direction Cosine Matrix:\n";
    display_matrix(Q1);

    std::cout << "Results for quaternion q2: [0, 1, 0, 0]\n";
    std::cout << "Direction Cosine Matrix:\n";
    display_matrix(Q2);

    std::cout << "Results for quaternion q3: [0, 0, 1, 0]\n";
    std::cout << "Direction Cosine Matrix:\n";
    display_matrix(Q3);

    std::cout << "Results for quaternion q4: [0, 0, 0, 1]\n";
    std::cout << "Direction Cosine Matrix:\n";
    display_matrix(Q4);
}

int main() 
{
    dcm_from_q_test();
    return 0;
}
