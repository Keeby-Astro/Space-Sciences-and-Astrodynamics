#ifndef F_AND_G_H
#define F_AND_G_H

#include "stumpC.h"
#include "stumpS.h"
#include <cmath>

/*
    This function calculates the Lagrange f and g coefficients.

    mu - the gravitational parameter (km^3/s^2)
    a - reciprocal of the semimajor axis (1/km)
    ro - the radial position at time to (km)
    t - the time elapsed since ro (s)
    x - the universal anomaly after time t (km^0.5)
    f - the Lagrange f coefficient (dimensionless)
    g - the Lagrange g coefficient (s)

    User h-functions required: stumpC, stumpS
*/
void f_and_g(double x, double t, double ro, double a, double& f, double& g, double mu) 
{
    double z = a * x * x;

    //...Equation 3.69a:
    f = 1 - (x * x / ro) * stumpC(z);

    //...Equation 3.69b:
    g = t - (1.0 / std::sqrt(mu)) * x * x * x * stumpS(z);
}

#endif // F_AND_G_H