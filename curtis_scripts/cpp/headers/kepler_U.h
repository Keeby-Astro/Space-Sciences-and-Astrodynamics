#ifndef KEPLER_U_H
#define KEPLER_U_H

#include <cmath>
#include <iostream>
#include "stumpC.h"
#include "stumpS.h"

inline double kepler_U(double dt, double ro, double vro, double a, double mu) 
{
    //...Set an error tolerance and a limit on the number of iterations:
    const double error_tolerance = 1.e-8;
    const int nMax = 1000;

    //...Starting value for x:
    double x = std::sqrt(mu) * std::abs(a) * dt;

    //...Iterate on Equation 3.65 until until convergence occurs within
    //...the error tolerance:
    double ratio = 1.0;
    int n = 0;
    while (std::abs(ratio) > error_tolerance && n <= nMax) 
    {
        ++n;
        double z = a * x * x;
        double C = stumpC(z);
        double S = stumpS(z);

        double F = ro * vro / std::sqrt(mu) * x * x * C + (1.0 - a * ro) * x * x * x * S + ro * x - std::sqrt(mu) * dt;
        double dFdx = ro * vro / std::sqrt(mu) * x * (1.0 - a * x * x * S) + (1.0 - a * ro) * x * x * C + ro;

        ratio = F / dFdx;
        x -= ratio;
    }

    //...Deliver a value for x, but report that nMax was reached:
    if (n > nMax) 
    {
        std::cerr << "\n **No. iterations of Kepler equation = " << n;
        std::cerr << "\n F/dFdx = " << ratio << "\n";
    }

    return x;
}

#endif // KEPLER_U_H