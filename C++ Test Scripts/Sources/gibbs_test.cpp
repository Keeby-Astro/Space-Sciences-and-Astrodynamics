#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <array>
#include "gibbs.h"
#include "coe_from_sv.h"

/*
    This program uses Algorithm 5.1 (Gibbs method) and Algorithm 4.2
    to obtain the orbital elements from the data.

    deg        - factor for converting between degrees and radians
    pi         - 3.1415926...
    mu         - gravitational parameter (km^3/s^2)
    r1, r2, r3 - three coplanar geocentric position vectors (km)
    ierr       - 0 if r1, r2, r3 are found to be coplanar
                 1 otherwise
    v2         - the velocity corresponding to r2 (km/s)
    coe        - orbital elements [h e RA incl w TA a]
                 where h = angular momentum (km^2/s)
                 e       = eccentricity
                 RA      = right ascension of the ascending node (rad)
                 incl    = orbit inclination (rad)
                 w       = argument of perigee (rad)
                 TA      = true anomaly (rad)
                 a       = semimajor axis (km)
    T          - period of elliptic orbit (s)

    User h-functions required: gibbs, coe_from_sv
*/
int main() 
{
    #define M_PI 3.14159265358979323846
    const double mu = 398600.4418;
    const double deg = 3.141592653589793 / 180.0;

    //...Data declaration:
    std::array<double, 3> r1 = {-294.32, 4265.1, 5986.7};
    std::array<double, 3> r2 = {-1365.5, 3637.6, 6346.8};
    std::array<double, 3> r3 = {-2940.3, 2473.7, 6555.8};
    //...

    //...Echo the input data to the console:
    std::cout << "-----------------------------------------------------\n";
    std::cout << "\nInput data:\n";
    std::cout << "\nGravitational parameter (km^3/s^2) = " << mu << "\n";
    std::cout << "\nr1 (km) = [" << r1[0] << " " << r1[1] << " " << r1[2] << "]\n";
    std::cout << "r2 (km) = [" << r2[0] << " " << r2[1] << " " << r2[2] << "]\n";
    std::cout << "r3 (km) = [" << r3[0] << " " << r3[1] << " " << r3[2] << "]\n\n";

    //...Algorithm 5.1:
    auto [v2, ierr] = gibbs(r1, r2, r3);

    //...If the vectors r1, r2, r3, are not coplanar, abort:
    if (ierr == 1) 
    {
        std::cerr << "\nThese vectors are not coplanar.\n\n";
        return 1;
    }

    //...Algorithm 4.2:
    std::vector<double> r2_vec(r2.begin(), r2.end());
    std::vector<double> v2_vec(v2.begin(), v2.end());
    std::vector<double> coe = coe_from_sv(r2_vec, v2_vec, mu);

    double h = coe[0];
    double e = coe[1];
    double i = coe[3];
    double RAAN = coe[2];
    double omega = coe[4];
    double TA = coe[5];
    double a = coe[6];

    //...Output the results to the console:
    std::cout << "Solution:\n";
    std::cout << "\nv2 (km/s) = [" << v2[0] << " " << v2[1] << " " << v2[2] << "]\n";
    std::cout << "\nOrbital elements:\n";
    std::cout << "Angular momentum (km^2/s) = " << h << "\n";
    std::cout << "Eccentricity = " << e << "\n";
    std::cout << "Inclination (deg) = " << i / deg << "\n";
    std::cout << "RA of ascending node (deg) = " << RAAN / deg << "\n";
    std::cout << "Argument of perigee (deg) = " << omega / deg << "\n";
    std::cout << "True anomaly (deg) = " << TA / deg << "\n";
    std::cout << "Semimajor axis (km) = " << a << "\n";
    //...If the orbit is an ellipse, output the period:
    if (e < 1.0) 
    {
        double T = 2.0 * M_PI / std::sqrt(mu) * std::pow(a, 1.5);
        std::cout << "\nPeriod:\n";
        std::cout << "Seconds = " << T << "\n";
    }
    std::cout << "-----------------------------------------------------\n";

    return 0;
}
