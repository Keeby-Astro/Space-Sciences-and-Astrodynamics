#include <iostream>
#include <iomanip>
#include "fDot_and_gDot.h"

/*
    Outputs the time derivatives of the Lagrange f and g coefficients.

    fdot - time derivative of the Lagrange f coefficient (1/s)
    gdot - time derivative of the Lagrange g coefficient (dimensionless)

    User h-functions required: fDot_and_gDot
*/
int main() 
{
    const double mu = 398600.4418;

    // Inputs
    double x = 1.5;          // Universal anomaly (km^0.5)
    double r = 7000.0;       // Radial position after time t (km)
    double ro = 6878.0;      // Radial position at time t0 (km)
    double a = 1.0 / 8000.0; // Reciprocal of semimajor axis (1/km)

    // Calculate fDot and gDot
    auto [fdot, gdot] = fDot_and_gDot(x, r, ro, a, mu);

    // Display results
    std::cout << "Results:" << std::endl;
    std::cout << "Lagrange f_dot coefficient: " << std::scientific << fdot << " 1/s" << std::endl;
    std::cout << "Lagrange g_dot coefficient: " << std::fixed << gdot << " (dimensionless)" << std::endl;

    return 0;
}
