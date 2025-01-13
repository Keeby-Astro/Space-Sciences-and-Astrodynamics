#include <iostream>
#include "kepler_H.h"

/*
    This program uses Algorithm 3.2 and the data to solve 
    Kepler's equation for the hyperbola.

    e - eccentricity
    M - hyperbolic mean anomaly (dimensionless)
    F - hyperbolic eccentric anomaly (dimensionless)

    User h-function required: kepler_H
*/
int main()
{
    //...Data declaration:
    double e = 2.7696;
    double M = 40.69;

    //...Pass the input data to the function kepler_H, which returns F:
    double F = kepler_H(e, M);

    //...Echo the input data and output to the console:
    std::cout << "----------------------------------------------";
    std::cout << "\n Eccentricity                 = " << e;
    std::cout << "\n Hyperbolic mean anomaly      = " << M << "\n";
    std::cout << "\n Hyperbolic eccentric anomaly = " << F;
    std::cout << "\n----------------------------------------------\n";

    return 0;
}
