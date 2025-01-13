#include <iostream>
#include <vector>
#include <cmath>
#include "rv_from_r0v0_ta.h"

/*
    This program computes the state vector [R,V] from the initial
    state vector [R0,V0] and the change in true anomaly, using the data.

    mu - gravitational parameter (km^3/s^2)
    R0 - the initial position vector (km)
    V0 - the initial velocity vector (km/s)
    r0 - magnitude of R0
    v0 - magnitude of V0
    R  - final position vector (km)
    V  - final velocity vector (km/s)
    r  - magnitude of R
    v  - magnitude of V
    dt - change in true anomaly (degrees)

    User h-functions required: rv_from_r0v0_ta
*/

// Utility functions
double vector_norm(const std::vector<double>& vec) {
    double sum = 0.0;
    for (double val : vec) {
        sum += val * val;
    }
    return std::sqrt(sum);
}

void print_vector(const std::vector<double>& vec, const std::string& label) {
    std::cout << label << " = [";
    for (size_t i = 0; i < vec.size(); ++i) {
        std::cout << vec[i];
        if (i < vec.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]\n";
}

int main() {
    const double mu = 398600.4418;

    //...Input data:
    std::vector<double> R0 = {8182.4, -6865.9, 0.0};
    std::vector<double> V0 = {0.47572, 8.8116, 0.0};
    double dt = 120.0;
    //...End input data

    //...Algorithm 2.3:
    auto [R, V] = rv_from_r0v0_ta(R0, V0, dt, mu);

    double r = vector_norm(R);
    double v = vector_norm(V);
    double r0 = vector_norm(R0);
    double v0 = vector_norm(V0);

    std::cout << "---------------------------------------------------------------\n";
    std::cout << "Initial State Vector:\n";
    print_vector(R0, "r0");
    std::cout << "Magnitude = " << r0 << "\n";
    print_vector(V0, "v0");
    std::cout << "Magnitude = " << v0 << "\n\n";

    std::cout << "State vector after " << dt << " degree change in true anomaly:\n";
    print_vector(R, "r");
    std::cout << "Magnitude = " << r << "\n";
    print_vector(V, "v");
    std::cout << "Magnitude = " << v << "\n";
    std::cout << "---------------------------------------------------------------\n";

    return 0;
}
