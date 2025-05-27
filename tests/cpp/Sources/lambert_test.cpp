#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <iomanip>
#include "lambert.h"
#include "coe_from_sv.h"

/*
    This program uses Algorithm 5.2 to solve Lambert's problem for the
    data provided.

    deg - factor for converting between degrees and radians
    pi - 3.1415926...
    mu - gravitational parameter (km^3/s^2)
    r1, r2 - initial and final position vectors (km)
    dt - time between r1 and r2 (s)
    string - = 'pro' if the orbit is prograde
    = 'retro if the orbit is retrograde
    v1, v2 - initial and final velocity vectors (km/s)
    coe - orbital elements [h e RA incl w TA a]
    where h = angular momentum (km^2/s)
    e = eccentricity
    RA = right ascension of the ascending node (rad)
    incl = orbit inclination (rad)
    w = argument of perigee (rad)
    TA = true anomaly (rad)
    a = semimajor axis (km)
    TA1 - Initial true anomaly (rad)
    TA2 - Final true anomaly (rad)
    T - period of an elliptic orbit (s)

    User M-functions required: lambert, coe_from_sv
*/
void lambert_test() 
{
    const double mu = 398600.4418;

    //...Data declaration:
    std::vector<double> r1 = {5000.0, 10000.0, 2100.0};
    std::vector<double> r2 = {-14600.0, 2500.0, 7000.0};
    double dt = 3600.0;
    std::string orbit_type = "pro";
    //...

    //...Algorithm 5.2:
    auto [v1, v2] = lambert(r1, r2, dt, orbit_type, mu);

    //...Algorithm 4.1 (using r1 and v1):
    auto coe1 = coe_from_sv(r1, v1, mu);
    //...Save the initial true anomaly:
    double TA1 = coe1[5];

    //...Algorithm 4.1 (using r2 and v2):
    auto coe2 = coe_from_sv(r2, v2, mu);
    //...Save the final true anomaly:
    double TA2 = coe2[5];

    //...Echo the input data and output the results to the console:
    std::cout << "-----------------------------------------------------\n";
    std::cout << "\nInput data:\n\n";
    std::cout << "   Gravitational parameter (km^3/s^2) = " << mu << "\n";
    std::cout << "   r1 (km)                       = [" << r1[0] << " " << r1[1] << " " << r1[2] << "]\n";
    std::cout << "   r2 (km)                       = [" << r2[0] << " " << r2[1] << " " << r2[2] << "]\n";
    std::cout << "   Elapsed time (s)              = " << dt << "\n";

    std::cout << "\nSolution:\n\n";
    std::cout << "   v1 (km/s)                     = [" << v1[0] << " " << v1[1] << " " << v1[2] << "]\n";
    std::cout << "   v2 (km/s)                     = [" << v2[0] << " " << v2[1] << " " << v2[2] << "]\n";

    std::cout << "\nOrbital elements:\n";
    std::cout << "   Angular momentum (km^2/s)     = " << coe1[0] << "\n";
    std::cout << "   Eccentricity                  = " << coe1[1] << "\n";
    std::cout << "   Inclination (deg)             = " << coe1[3] * (180.0 / M_PI) << "\n";
    std::cout << "   RA of ascending node (deg)    = " << coe1[2] * (180.0 / M_PI) << "\n";
    std::cout << "   Argument of perigee (deg)     = " << coe1[4] * (180.0 / M_PI) << "\n";
    std::cout << "   True anomaly initial (deg)    = " << TA1 * (180.0 / M_PI) << "\n";
    std::cout << "   True anomaly final (deg)      = " << TA2 * (180.0 / M_PI) << "\n";
    std::cout << "   Semimajor axis (km)           = " << coe1[6] << "\n";

    double rp = (coe1[0] * coe1[0]) / mu / (1.0 + coe1[1]);
    std::cout << "   Periapse radius (km)          = " << rp << "\n";
    //...If the orbit is an ellipse, output its period:
    if (coe1[1] < 1.0) {
        double T = 2.0 * M_PI / std::sqrt(mu) * std::pow(coe1[6], 1.5);
        std::cout << "\n   Period:\n";
        std::cout << "     Seconds                 = " << T << "\n";
        std::cout << "     Minutes                 = " << T / 60.0 << "\n";
        std::cout << "     Hours                   = " << T / 3600.0 << "\n";
        std::cout << "     Days                    = " << T / 86400.0 << "\n";
    }
    std::cout << "\n-----------------------------------------------------\n\n";
}

int main() 
{
    lambert_test();
    return 0;
}
