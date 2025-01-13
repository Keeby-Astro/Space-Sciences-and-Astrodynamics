#include <iostream>
#include <array>
#include "q_from_dcm.h"

/*
    Verifies the calculation of the quaternion from the direction
    cosine matrix using known inputs.

    Q - direction cosine matrix
    q - quaternion (where q(4) is the scalar part)

    User h-functions required: q_from_dcm
*/
void q_from_dcm_test() 
{
    // Define test direction cosine matrices
    std::array<std::array<double, 3>, 3> Q1 = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};   // Identity matrix
    std::array<std::array<double, 3>, 3> Q2 = {{{1, 0, 0}, {0, -1, 0}, {0, 0, -1}}}; // 180-degree rotation about x-axis
    std::array<std::array<double, 3>, 3> Q3 = {{{-1, 0, 0}, {0, 1, 0}, {0, 0, -1}}}; // 180-degree rotation about y-axis
    std::array<std::array<double, 3>, 3> Q4 = {{{-1, 0, 0}, {0, -1, 0}, {0, 0, 1}}}; // 180-degree rotation about z-axis

    // Calculate the quaternions using the function
    std::array<double, 4> q1 = q_from_dcm(Q1);
    std::array<double, 4> q2 = q_from_dcm(Q2);
    std::array<double, 4> q3 = q_from_dcm(Q3);
    std::array<double, 4> q4 = q_from_dcm(Q4);

    // Display the results
    std::cout << "Results for direction cosine matrix Q1:" << std::endl;
    std::cout << "Quaternion: [" << q1[0] << ", " << q1[1] << ", " << q1[2] << ", " << q1[3] << "]\n";

    std::cout << "Results for direction cosine matrix Q2:" << std::endl;
    std::cout << "Quaternion: [" << q2[0] << ", " << q2[1] << ", " << q2[2] << ", " << q2[3] << "]\n";

    std::cout << "Results for direction cosine matrix Q3:" << std::endl;
    std::cout << "Quaternion: [" << q3[0] << ", " << q3[1] << ", " << q3[2] << ", " << q3[3] << "]\n";

    std::cout << "Results for direction cosine matrix Q4:" << std::endl;
    std::cout << "Quaternion: [" << q4[0] << ", " << q4[1] << ", " << q4[2] << ", " << q4[3] << "]\n";
}

int main() 
{
    q_from_dcm_test();
    return 0;
}
