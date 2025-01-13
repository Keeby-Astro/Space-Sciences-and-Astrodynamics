#include <iostream>
#include <iomanip>
#include "coe_from_sv.h"

/*
%{
    This program uses Algorithm 4.2 to obtain the orbital
    elements from the state vector.
    pi  - 3.1415926...
    deg - factor for converting between degrees and radians
    mu  - gravitational parameter (km^3/s^2)
    r   - position vector (km) in the geocentric equatorial frame
    v   - velocity vector (km/s) in the geocentric equatorial frame
    coe - orbital elements [h e RA incl w TA a]
          where h    = angular momentum (km^2/s)
                e    = eccentricity
                RA   = right ascension of the ascending node (rad)
                incl = orbit inclination (rad)
                w    = argument of perigee (rad)
                TA   = true anomaly (rad)
                a    = semimajor axis (km)
    T   - Period of an elliptic orbit (s)

    User h-function required: coe_from_sv
%}
*/
int main() {
    constexpr double pi = 3.141592653589793;
    constexpr double deg = pi / 180.0;
    constexpr double mu = 398600.4418;

    //...Data declaration:
    std::vector<double> R = {-6045.0, -3490.0, 2500.0}; // Position vector (km)
    std::vector<double> V = {-3.457, 6.618, 2.533};     // Velocity vector (km/s)
    //...

    try {
        //...Algorithm 4.2:
        std::vector<double> coe = coe_from_sv(R, V, mu);

        //...Echo the input data and output results to the console:
        std::cout << "-----------------------------------------------------\n";
        std::cout << "Gravitational parameter (km^3/s^2) = " << mu << "\n";
        std::cout << "\nState vector:\n";
        std::cout << "r (km): [" << R[0] << ", " << R[1] << ", " << R[2] << "]\n";
        std::cout << "v (km/s): [" << V[0] << ", " << V[1] << ", " << V[2] << "]\n";
        std::cout << "\nResults:\n";
        std::cout << "Angular momentum (km^2/s): " << coe[0] << "\n";
        std::cout << "Eccentricity: " << coe[1] << "\n";
        std::cout << "Right ascension of ascending node (deg): " << coe[2] / deg << "\n";
        std::cout << "Inclination (deg): " << coe[3] / deg << "\n";
        std::cout << "Argument of perigee (deg): " << coe[4] / deg << "\n";
        std::cout << "True anomaly (deg): " << coe[5] / deg << "\n";
        std::cout << "Semimajor axis (km): " << coe[6] << "\n";

        //...if the orbit is an ellipse, output its period (Equation 2.73):
        if (coe[1] < 1.0) {
            double T = 2 * pi / std::sqrt(mu) * std::pow(coe[6], 1.5);
            std::cout << "\nPeriod:\n";
            std::cout << "Seconds: " << T << "\n";
            std::cout << "Minutes: " << T / 60.0 << "\n";
            std::cout << "Hours: " << T / 3600.0 << "\n";
            std::cout << "Days: " << T / 86400.0 << "\n";
        }

        std::cout << "-----------------------------------------------------\n";
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
    }

    return 0;
}
