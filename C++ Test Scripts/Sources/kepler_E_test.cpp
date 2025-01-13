#include <iostream>
#include <iomanip>
#include "kepler_E.h"

/*
    This program uses Algorithm 3.1 and the data to solve
    Kepler's equation.

    e - eccentricity
    M - mean anomaly (rad)
    E - eccentric anomaly (rad)

    User h-function required: kepler_E
*/

int main() 
{
    //...Data declaration:
    double e = 0.37255;
    double M = 3.6029;
    //...

    //...Pass the input data to the function kepler_E, which returns E:
    double E = kepler_E(e, M);

    //...Echo the input data and output to the console:
    std::cout << "---------------------------------------------------\n";
    std::cout << " Eccentricity                = " << e << "\n";
    std::cout << " Mean anomaly (radians)      = " << M << "\n";
    std::cout << " Eccentric anomaly (radians) = " << E << "\n";
    std::cout << "---------------------------------------------------\n";

    return 0;
}
