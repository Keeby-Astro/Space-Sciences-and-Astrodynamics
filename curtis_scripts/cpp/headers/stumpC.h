#ifndef STUMPC_H
#define STUMPC_H

#include <iostream>
#include <cmath>
/*
    This function evaluates the Stumpff function C(z) according
    to Equation 3.53.

    z - input argument
    c - value of C(z)

    User h-functions required: none
*/
inline double stumpC(double z) 
{
    if (z > 0) 
    {
        return (1.0 - std::cos(std::sqrt(z))) / z;
    } else if (z < 0) 
    {
        return (std::cosh(std::sqrt(-z)) - 1.0) / (-z);
    } else 
    {
        return 1.0 / 2.0;
    }
}

#endif // STUMPC_H