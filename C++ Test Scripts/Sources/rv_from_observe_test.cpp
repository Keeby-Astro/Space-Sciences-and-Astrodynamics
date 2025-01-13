#include <iostream>
#include <iomanip>
#include "rv_from_observe.h"
#include "coe_from_sv.h"

/*
    This program uses Algorithms 5.4 and 4.2 to obtain the orbital
    elements from the observational data.
    
    deg    - conversion factor between degrees and radians
    pi     - 3.1415926...
    mu     - gravitational parameter (km^3/s^2)

    Re     - equatorial radius of the earth (km)
    f      - earth's flattening factor
    wE     - angular velocity of the earth (rad/s)
    omega  - earth's angular velocity vector (rad/s) in the
             geocentric equatorial frame

    rho    - slant range of object (km)
    rhodot - range rate (km/s)
    A      - azimuth (deg) of object relative to observation site
    Adot   - time rate of change of azimuth (deg/s)
    a      - elevation angle (deg) of object relative to observation site
    adot   - time rate of change of elevation angle (degrees/s)

    theta  - local sidereal time (deg) of tracking site
    phi    - geodetic latitude (deg) of site
    H      - elevation of site (km)

    r      - geocentric equatorial position vector of object (km)
    v      - geocentric equatorial velocity vector of object (km)

    coe    - orbital elements [h e RA incl w TA a]
             where
                 h    = angular momentum (km^2/s)
                 e    = eccentricity
                 RA   = right ascension of the ascending node (rad)
                 incl = inclination of the orbit (rad)
                 w    = argument of perigee (rad)
                 TA   = true anomaly (rad)
                 a    = semimajor axis (km)
    rp     - perigee radius (km)
    T      - period of elliptical orbit (s)
    
    User h-functions required: rv_from_observe, coe_from_sv
*/
int main() 
{
    constexpr double pi = 3.141592653589793;
    constexpr double deg = pi / 180.0;
    double f = 1.0 / 298.256421867;
    double Re = 6378.13655; // km
    double wE = 7.292115e-5; // rad/s
    double mu = 398600.4418; // km^3/s^2

    //...Data declaration:
    double rho = 2551.0;
    double rhodot = 0.0;
    double A = 90.0;
    double Adot = 0.1130;
    double a = 30.0;
    double adot = 0.05651;
    double theta = 300.0;
    double phi = 60.0;
    double H = 0.0;
    //...

    //...Algorithm 5.4:
    auto [r, v] = rv_from_observe(rho, rhodot, A, Adot, a, adot, theta, phi, H, Re, f, wE);

    //...Algorithm 4.2:
    auto coe = coe_from_sv(r, v, mu);

    double h = coe[0];
    double e = coe[1];
    double RA = coe[2];
    double incl = coe[3];
    double w = coe[4];
    double TA = coe[5];
    double a_semimajor = coe[6];

    //...Equation 2.40:
    double rp = h * h / mu / (1 + e);

    //...Echo the input data and output the solution to
    // the command window:
    std::cout << "-----------------------------------------------------\n";
    std::cout << "Input data:\n\n";
    std::cout << "Slant range (km)              = " << rho << "\n";
    std::cout << "Slant range rate (km/s)       = " << rhodot << "\n";
    std::cout << "Azimuth (deg)                 = " << A << "\n";
    std::cout << "Azimuth rate (deg/s)          = " << Adot << "\n";
    std::cout << "Elevation (deg)               = " << a << "\n";
    std::cout << "Elevation rate (deg/s)        = " << adot << "\n";
    std::cout << "Local sidereal time (deg)     = " << theta << "\n";
    std::cout << "Latitude (deg)                = " << phi << "\n";
    std::cout << "Altitude above sea level (km) = " << H << "\n\n";

    std::cout << "Solution:\n\n";
    std::cout << "State vector:\n";
    std::cout << "r (km) = [" << r[0] << ", " << r[1] << ", " << r[2] << "]\n";
    std::cout << "v (km/s) = [" << v[0] << ", " << v[1] << ", " << v[2] << "]\n\n";

    std::cout << "Orbital elements:\n\n";
    std::cout << "Angular momentum (km^2/s)     = " << h << "\n";
    std::cout << "Eccentricity                  = " << e << "\n";
    std::cout << "Inclination (deg)             = " << incl / deg << "\n";
    std::cout << "RA of ascending node (deg)    = " << RA / deg << "\n";
    std::cout << "Argument of perigee (deg)     = " << w / deg << "\n";
    std::cout << "True anomaly (deg)            = " << TA / deg << "\n";
    std::cout << "Semimajor axis (km)           = " << a_semimajor << "\n";
    std::cout << "Perigee radius (km)           = " << rp << "\n\n";
    //...If the orbit is an ellipse, output its period:
    if (e < 1.0) 
    {
        double T = 2 * pi / std::sqrt(mu) * std::pow(a_semimajor, 1.5);
        std::cout << "Period:\n";
        std::cout << "Seconds = " << T << "\n";
        std::cout << "Minutes = " << T / 60.0 << "\n";
        std::cout << "Hours   = " << T / 3600.0 << "\n";
        std::cout << "Days    = " << T / 86400.0 << "\n";
    }

    std::cout << "-----------------------------------------------------\n";
    return 0;
}
