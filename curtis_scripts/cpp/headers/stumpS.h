#ifndef STUMPS_H
#define STUMPS_H

#include <iostream>
#include <cmath>

/*
    This function evaluates the Stumpff function S(z) according
    to Equation 3.52.

    z - input argument
    s - value of S(z)

    User h-functions required: none
*/
inline double stumpS(double z) 
{
    if (z > 0) {
        return (std::sqrt(z) - std::sin(std::sqrt(z))) / std::pow(std::sqrt(z), 3);
    } else if (z < 0) {
        return (std::sinh(std::sqrt(-z)) - std::sqrt(-z)) / std::pow(std::sqrt(-z), 3);
    } else {
        return 1.0 / 6.0;
    }
}

#endif // STUMPS_H