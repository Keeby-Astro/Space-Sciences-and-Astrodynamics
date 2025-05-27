#ifndef FDOT_AND_GDOT_H
#define FDOT_AND_GDOT_H

#include <cmath>
#include "stumpC.h"
#include "stumpS.h"

/*
    This function calculates the time derivatives of the
    Lagrange f and g coefficients.

    mu   - the gravitational parameter (km^3/s^2)
    a    - reciprocal of the semimajor axis (1/km)
    ro   - the radial position at time to (km)
    t    - the time elapsed since initial state vector (s)
    r    - the radial position after time t (km)
    x    - the universal anomaly after time t (km^0.5)
    fdot - time derivative of the Lagrange f coefficient (1/s)
    gdot - time derivative of the Lagrange g coefficient (dimensionless)

    User h-functions required: stumpC, stumpS
*/
inline std::pair<double, double> fDot_and_gDot(double x, double r, double ro, double a, double mu) 
{
    double z = a * x * x;

    //...Equation 3.69c:
    double fdot = std::sqrt(mu) / (r * ro) * (z * stumpS(z) - 1) * x;

    //...Equation 3.69d:
    double gdot = 1 - (x * x / r) * stumpC(z);

    return {fdot, gdot};
}

#endif // FDOT_AND_GDOT_H
