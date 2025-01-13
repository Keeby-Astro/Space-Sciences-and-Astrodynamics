#include <iostream>
#include <vector>
#include "f_and_g_ta.h"

/*
    The script defines input parameters for r0, v0, dt, and mu, and then
    calculates the f and g coefficients using the f_and_g_ta function.
    It also prints the results for inspection.

    Inputs:
    r0  - Initial position vector (km)
    v0  - Initial velocity vector (km/s)
    dt  - Change in true anomaly (degrees)
    mu  - Gravitational parameter (km^3/s^2)

    Outputs:
    f   - Lagrange f coefficient (dimensionless)
    g   - Lagrange g coefficient (s)

    User h-functions required: f_and_g_ta
*/
int main() {

    double mu = 398600.4418;                     // Gravitational parameter for Earth (km^3/s^2)
    std::vector<double> r0 = {7000.0, 0.0, 0.0}; // Initial position vector (km)
    std::vector<double> v0 = {0.0, 7.5, 0.0};    // Initial velocity vector (km/s)
    double dt = 45.0;                            // Change in true anomaly (degrees)

    double f, g;

    // Call the f_and_g_ta function
    try {
        f_and_g_ta(r0, v0, dt, mu, f, g);

        // Display the results
        std::cout << "Results:\n";
        std::cout << "Lagrange f coefficient: " << f << "\n";
        std::cout << "Lagrange g coefficient: " << g << "\n";
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
