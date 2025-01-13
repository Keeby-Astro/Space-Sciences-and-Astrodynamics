#include <iostream>
#include <iomanip>
#include "kepler_U.h"

/*
    This program uses Algorithm 3.3 and the data
    to solve the universal Kepler's equation.

    mu  - gravitational parameter (km^3/s^2)
    x   - the universal anomaly (km^0.5)
    dt  - time since x = 0 (s)
    ro  - radial position when x = 0 (km)
    vro - radial velocity when x = 0 (km/s)
    a   - semimajor axis (km)

    User h-function required: kepler_U
*/
int main() 
{
    const double mu = 398600.4418; // Gravitational parameter (km^3/s^2)

    //...Data declaration:
    double ro = 10000.0;   // Initial radial coordinate (km)
    double vro = 3.0752;   // Initial radial velocity (km/s)
    double dt = 3600.0;    // Elapsed time (s)
    double a = -19655.0;   // Semimajor axis (km)
    //...

    //...Pass the input data to the function keplfr_U, which returns x
    //...(Universal Kepler's requires the reciprocal of semimajor axis):
    double x = kepler_U(dt, ro, vro, 1.0 / a, mu);

    //...Echo the input data and output the results to the console:
    std::cout << "-----------------------------------------------------\n";
    std::cout << "Initial radial coordinate (km) = " << ro << "\n";
    std::cout << "Initial radial velocity (km/s) = " << vro << "\n";
    std::cout << "Elapsed time (seconds)         = " << dt << "\n";
    std::cout << "Semimajor axis (km)            = " << a << "\n";
    std::cout << "\nUniversal anomaly (km^0.5)     = " << x << "\n";
    std::cout << "-----------------------------------------------------\n";

    return 0;
}
