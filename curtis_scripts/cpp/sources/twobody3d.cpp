// ALGORITHM 2.1: NUMERICAL SOLUTION OF THE TWO-BODY PROBLEM
// RELATIVE TO AN INERTIAL FRAME

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "rkf45.h"

// Universal gravitational constant (km^3/kg/s^2)
const double G = 6.67259e-20;

// Masses of the two bodies (kg)
const double m1 = 1.0e26;
const double m2 = 1.0e26;

// Function to compute derivatives (rates)
std::vector<double> rates(double t, const std::vector<double> &y) {
    // Extract positions and velocities
    std::vector<double> R1 = {y[0], y[1], y[2]};
    std::vector<double> R2 = {y[3], y[4], y[5]};

    std::vector<double> V1 = {y[6], y[7], y[8]};
    std::vector<double> V2 = {y[9], y[10], y[11]};

    // Compute the relative distance and accelerations
    double r = std::sqrt(std::pow(R2[0] - R1[0], 2) +
                         std::pow(R2[1] - R1[1], 2) +
                         std::pow(R2[2] - R1[2], 2));

    std::vector<double> A1 = {G * m2 * (R2[0] - R1[0]) / std::pow(r, 3),
                              G * m2 * (R2[1] - R1[1]) / std::pow(r, 3),
                              G * m2 * (R2[2] - R1[2]) / std::pow(r, 3)};

    std::vector<double> A2 = {G * m1 * (R1[0] - R2[0]) / std::pow(r, 3),
                              G * m1 * (R1[1] - R2[1]) / std::pow(r, 3),
                              G * m1 * (R1[2] - R2[2]) / std::pow(r, 3)};

    // Return derivatives: [V1; V2; A1; A2]
    return {V1[0], V1[1], V1[2], V2[0], V2[1], V2[2], A1[0], A1[1], A1[2], A2[0], A2[1], A2[2]};
}

// Function to output and plot the results
void output(const std::vector<double> &t, const std::vector<std::vector<double>> &y) {
    std::cout << "\nResults:\n";

    for (size_t i = 0; i < t.size(); ++i) {
        std::cout << "Time: " << t[i] << " s\n";
        std::cout << "Body 1 Position: (" << y[i][0] << ", " << y[i][1] << ", " << y[i][2] << ") km\n";
        std::cout << "Body 2 Position: (" << y[i][3] << ", " << y[i][4] << ", " << y[i][5] << ") km\n";
        std::cout << "\n";
    }

    // Compute the center of mass positions
    std::vector<double> X1, Y1, Z1, X2, Y2, Z2, XG, YG, ZG;
    for (const auto& state : y) {
        X1.push_back(state[0]);
        Y1.push_back(state[1]);
        Z1.push_back(state[2]);
        X2.push_back(state[3]);
        Y2.push_back(state[4]);
        Z2.push_back(state[5]);
        double xg = (m1 * state[0] + m2 * state[3]) / (m1 + m2);
        double yg = (m1 * state[1] + m2 * state[4]) / (m1 + m2);
        double zg = (m1 * state[2] + m2 * state[5]) / (m1 + m2);
        XG.push_back(xg);
        YG.push_back(yg);
        ZG.push_back(zg);
    }

    // Output the trajectories
    std::cout <<"\n";
    std::cout << "Trajectories:\n";
    std::cout << "Body 1 Trajectory: \n";
    for (size_t i = 0; i < X1.size(); ++i) {
        std::cout << "(" << X1[i] << ", " << Y1[i] << ", " << Z1[i] << ") km\n";
    }
    std::cout << "\nBody 2 Trajectory: \n";
    for (size_t i = 0; i < X2.size(); ++i) {
        std::cout << "(" << X2[i] << ", " << Y2[i] << ", " << Z2[i] << ") km\n";
    }
    std::cout << "\nCenter of Mass Trajectory: \n";
    for (size_t i = 0; i < XG.size(); ++i) {
        std::cout << "(" << XG[i] << ", " << YG[i] << ", " << ZG[i] << ") km\n";
    }
}

int main() {
    // Initial time and final time (s)
    double t0 = 0.0;
    double tf = 480.0;

    // Initial positions (km) and velocities (km/s)
    std::vector<double> R1_0 = {0.0, 0.0, 0.0};
    std::vector<double> R2_0 = {3000.0, 0.0, 0.0};

    std::vector<double> V1_0 = {10.0, 20.0, 30.0};
    std::vector<double> V2_0 = {0.0, 40.0, 0.0};

    // Initial state vector: [R1_0; R2_0; V1_0; V2_0]
    std::vector<double> y0 = {R1_0[0], R1_0[1], R1_0[2],
                              R2_0[0], R2_0[1], R2_0[2],
                              V1_0[0], V1_0[1], V1_0[2],
                              V2_0[0], V2_0[1], V2_0[2]};

    // Solve the system using RKF4(5)
    RKF45 result = rkf45(rates, {t0, tf}, y0, 1.e-8);

    // Output the results
    output(result.tout, result.yout);

    return 0;
}
