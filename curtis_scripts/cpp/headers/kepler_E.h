// ALGORITHM 3.1: SOLUTION OF KEPLER'S EQUATION BY NEWTON'S METHOD

#ifndef KEPLER_E_H
#define KEPLER_E_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#define M_PI 3.14159265358979323846

/*
{
    This function uses Newton's method to solve Kepler's
    equation E - e*sin(E) = M for the eccentric anomaly,
    given the eccentricity and the mean anomaly.

    E  - eccentric anomaly (radians)
    e  - eccentricity, passed from the calling program
    M  - mean anomaly (radians), passed from the calling program
    pi - 3.1415926...

    User h-functions required: none
}
*/
inline double kepler_E(double e, double M) 
{
    //...Set an error tolerance:
    const double error = 1.e-8;

    //...Select a starting value for E:
    double E = (M < M_PI) ? (M + e / 2.0) : (M - e / 2.0);

    //...Iterate on Equation 3.17 until E is determined to within
    //...the error tolerance:
    double ratio = 1.0;
    while (std::abs(ratio) > error) {
        ratio = (E - e * std::sin(E) - M) / (1.0 - e * std::cos(E));
        E -= ratio;
    }

    return E;
}

#endif // KEPLER_E_H
