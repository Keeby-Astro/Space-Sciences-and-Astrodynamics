// ALGORITHM 3.2: SOLUTION OF KEPLER'S EQUATION FOR THE HYPERBOLA USING NEWTON'S METHOD

#ifndef KEPLER_H_H
#define KEPLER_H_H

#include <cmath>
#include <iostream>

/*
    This function uses Newton's method to solve Kepler's equation
    for the hyperbola e*sinh(F) - F = M for the hyperbolic
    eccentric anomaly, given the eccentricity and the hyperbolic
    mean anomaly.

    F - hyperbolic eccentric anomaly (radians)
    e - eccentricity, passed from the calling program
    M - hyperbolic mean anomaly (radians), passed from the
        calling program

    User h-functions required: none
*/
inline double kepler_H(double e, double M)
{
    //...Set an error tolerance:
    const double error = 1.e-8;

    //...Starting value for F:
    double F = M;

    //...Iterate on Equation 3.45 until F is determined to within
    //...the error tolerance:
    double ratio = 1.0;
    while (std::abs(ratio) > error)
    {
        ratio = (e * std::sinh(F) - F - M) / (e * std::cosh(F) - 1);
        F -= ratio;
    }

    return F;
}

#endif
