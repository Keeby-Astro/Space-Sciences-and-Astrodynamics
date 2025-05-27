#define _USE_MATH_DEFINES
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "sv_from_coe.h"
#define M_PI 3.14159265358979323846

/*
    This program uses Algorithm 4.5 to obtain the state vector from
    the orbital elements.

    pi  - 3.1415926...
    deg - factor for converting between degrees and radians
    mu  - gravitational parameter (km^3/s^2)
    coe - orbital elements [h e RA incl w TA a]
          where h = angular momentum (km^2/s)
                e = eccentricity
                RA = right ascension of the ascending node (rad)
                incl = orbit inclination (rad)
                w = argument of perigee (rad)
                TA = true anomaly (rad)
                a = semimajor axis (km)
    r   - position vector (km) in geocentric equatorial frame
    v   - velocity vector (km) in geocentric equatorial frame

    User h-function required: sv_from_coe
*/
int main()
{
    const double deg = M_PI / 180.0;
    double mu = 398600.0;

    // ...Data declaration (angles in degrees):
    double h    = 80000.0;
    double e    = 1.4;
    double RA   = 40.0;
    double incl = 30.0;
    double w    = 60.0;
    double TA   = 30.0;
    // ...

    std::vector<double> coe = {h, e, RA * deg, incl * deg, w * deg, TA * deg};

    // ...Algorithm 4.5 (requires angular elements be in radians):
    auto result = sv_from_coe(coe, mu);
    std::vector<double> r = result.first;
    std::vector<double> v = result.second;

    // ...Echo the input data and output the results to the command window:
    std::cout << "-----------------------------------------------------\n";
    std::cout << " Gravitational parameter (km^3/s^2)  = " << mu << "\n\n";
    std::cout << " Angular momentum (km^2/s)           = " << h << "\n";
    std::cout << " Eccentricity                        = " << e << "\n";
    std::cout << " Right ascension (deg)               = " << RA << "\n";
    std::cout << " Argument of perigee (deg)           = " << w << "\n";
    std::cout << " True anomaly (deg)                  = " << TA << "\n\n";
    std::cout << " State vector:\n";
    std::cout << "   r (km)   = [" 
              << r[0] << " " << r[1] << " " << r[2] << "]\n";
    std::cout << "   v (km/s) = ["
              << v[0] << " " << v[1] << " " << v[2] << "]\n";
    std::cout << "-----------------------------------------------------\n";

    return 0;
}