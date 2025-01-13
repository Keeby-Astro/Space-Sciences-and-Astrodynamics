#include <iostream>
#include <vector>
#include "fDot_and_gDot_ta.h"

/*
    This script tests the function fDot_and_gDot_ta by passing example
    inputs and verifying the outputs. It computes the time derivatives of
    the Lagrange f and g coefficients for a given position and velocity
    vector, and a change in true anomaly.

    r0      - position vector at time t0 (km)
    v0      - velocity vector at time t0 (km/s)
    dt      - change in true anomaly (degrees)
    mu      - gravitational parameter (km^3/s^2)
    fdot    - time derivative of the Lagrange f coefficient (1/s)
    gdot    - time derivative of the Lagrange g coefficient (dimensionless)

    User h-functions required: fDot_and_gDot_ta
*/
int main() {
    // Inputs
    std::vector<double> r0 = {7000.0, 0.0, 0.0}; // Initial position vector (km)
    std::vector<double> v0 = {0.0, 7.546, 0.0};  // Initial velocity vector (km/s)
    double dt = 30.0;                            // Change in true anomaly (degrees)
    double mu = 398600.0;                        // Gravitational parameter (km^3/s^2)

    double fdot, gdot;

    // Call the function
    fDot_and_gDot_ta(r0, v0, dt, mu, fdot, gdot);

    // Display results
    std::cout << "Results:\n";
    std::cout << "Lagrange f_dot: " << fdot << " (1/s)\n";
    std::cout << "Lagrange g_dot: " << gdot << " (dimensionless)\n";

    return 0;
}
