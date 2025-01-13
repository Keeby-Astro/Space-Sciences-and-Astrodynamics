#include <iostream>
#include <iomanip>
#include <vector>
#include <functional>
#include "bisect.h"

/*
    This program uses the bisection method to find the three roots of
    Equation 2.204 for the earth-moon system.

    m1  - mass of the earth (kg)
    m2  - mass of the moon (kg)
    r12 - distance from the earth to the moon (km)
    p   - ratio of moon mass to total mass
    xl  - vector containing the low-side estimates of the three roots
    xu  - vector containing the high-side estimates of the three roots
    x   - vector containing the three computed roots

    User h-function required: bisect
    User subfunction requred: fun
*/

double fun(double z, double p)
/*
    This subroutine evaluates the function in Equation 2.204

    z - the dimensionless x - coordinate
    p - defined above
    f - the value of the function
*/
{
    return (1 - p) * (z + p) / std::pow(std::abs(z + p), 3) +
           p * (z + p - 1) / std::pow(std::abs(z + p - 1), 3) - z;
}

void output(double m1, double m2, double r12, const std::vector<double> &roots, double p)
/*
    This function prints out the x coordinates of L1, L2 and L3
    relative to the center of mass.
*/
{
    //...Output to the command window:
    std::cout << "\n\n---------------------------------------------\n";
    std::cout << "\n For\n";
    std::cout << "\n m1  = " << m1 << " kg";
    std::cout << "\n m2  = " << m2 << " kg";
    std::cout << "\n r12 = " << r12 << " km\n";
    std::cout << "\n the 3 colinear Lagrange points (the roots of\n";
    std::cout << " Equation 2.204) are:\n";
    
    std::cout << "\n L3: x = " << std::setw(10) << roots[0] * r12 << " km (f(x3) = " << fun(roots[0], p) << ")";
    std::cout << "\n L1: x = " << std::setw(10) << roots[1] * r12 << " km (f(x1) = " << fun(roots[1], p) << ")";
    std::cout << "\n L2: x = " << std::setw(10) << roots[2] * r12 << " km (f(x2) = " << fun(roots[2], p) << ")";
    std::cout << "\n\n---------------------------------------------\n";
}

int main()
{
    //...Input data:
    const double m1 = 5.974e24;
    const double m2 = 7.348e22;
    const double r12 = 3.844e5;

    std::vector<double> xl = {-1.1, 0.5, 1.0};
    std::vector<double> xu = {-0.9, 1.0, 1.5};
    //...End input data

    double p = m2 / (m1 + m2);
    std::vector<double> roots(3);

    for (size_t i = 0; i < xl.size(); ++i)
    {
        roots[i] = bisect([p](double z) { return fun(z, p); }, xl[i], xu[i]);
    }

    //...Output the results
    output(m1, m2, r12, roots, p);

    return 0;
}