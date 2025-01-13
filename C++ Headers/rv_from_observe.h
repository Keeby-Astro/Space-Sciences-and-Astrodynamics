#ifndef RV_FROM_OBSERVE_H
#define RV_FROM_OBSERVE_H

#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>

/*
    This function calculates the geocentric equatorial position and
    velocity vectors of an object from radar observations of range,
    azimuth, elevation angle and their rates.

    deg     - conversion factor between degrees and radians
    pi      - 3.1415926...
    Re      - equatorial radius of the earth (km)
    f       - earth's flattening factor
    wE      - angular velocity of the earth (rad/s)
    omega   - earth's angular velocity vector (rad/s) in the
              geocentric equatorial frame
    theta   - local sidereal time (degrees) of tracking site
    phi     - geodetic latitude (degrees) of site
    H       - elevation of site (km)
    R       - geocentric equatorial position vector (km) of tracking site
    Rdot    - inertial velocity (km/s) of site
    rho     - slant range of object (km)
    rhodot  - range rate (km/s)
    A       - azimuth (degrees) of object relative to observation site
    Adot    - time rate of change of azimuth (degrees/s)
    a       - elevation angle (degrees) of object relative to observation site
    adot    - time rate of change of elevation angle (degrees/s)
    dec     - topocentric equatorial declination of object (rad)
    decdot  - declination rate (rad/s)
    h       - hour angle of object (rad)
    RA      - topocentric equatorial right ascension of object (rad)
    RAdot   - right ascension rate (rad/s)
    Rho     - unit vector from site to object
    Rhodot  - time rate of change of Rho (1/s)
    r       - geocentric equatorial position vector of object (km)
    v       - geocentric equatorial velocity vector of object (km)

    User h-functions required: none
*/
std::pair<std::vector<double>, std::vector<double>> rv_from_observe(
    double rho, double rhodot, double A, double Adot, double a, double adot,
    double theta, double phi, double H, double Re, double f, double wE)
{
    constexpr double pi = 3.141592653589793;
    const double deg = pi / 180.0;
    std::vector<double> omega = {0.0, 0.0, wE};

    //...Convert angular quantities from degrees to radians:
    A *= deg;
    Adot *= deg;
    a *= deg;
    adot *= deg;
    theta *= deg;
    phi *= deg;

    //...Equation 5.56:
    double denom = std::sqrt(1.0 - (2 * f - f * f) * std::pow(std::sin(phi), 2));
    std::vector<double> R = 
    {
        (Re / denom + H) * std::cos(phi) * std::cos(theta),
        (Re / denom + H) * std::cos(phi) * std::sin(theta),
        (Re * (1 - f) * (1 - f) / denom + H) * std::sin(phi)
    };

    //...Equation 5.66:
    std::vector<double> Rdot = 
    {
        -omega[2] * R[1],
        omega[2] * R[0],
        0.0
    };

    //...Equation 5.83a:
    double dec = std::asin(std::cos(phi) * std::cos(A) * std::cos(a) + std::sin(phi) * std::sin(a));

    //...Equation 5.83b:
    double h = std::acos((std::cos(phi) * std::sin(a) - std::sin(phi) * std::cos(A) * std::cos(a)) / std::cos(dec));
    if (A > 0 && A < pi)
        h = 2 * pi - h;

    //...Equation 5.83c:
    double RA = theta - h;

    //...Equation 5.57:
    std::vector<double> Rho = 
    {
        std::cos(RA) * std::cos(dec),
        std::sin(RA) * std::cos(dec),
        std::sin(dec)
    };

    //...Equation 5.63:
    std::vector<double> r = 
    {
        R[0] + rho * Rho[0],
        R[1] + rho * Rho[1],
        R[2] + rho * Rho[2]
    };

    //...Equation 5.84:
    double decdot = (-Adot * std::cos(phi) * std::sin(A) * std::cos(a) +
                     adot * (std::sin(phi) * std::cos(a) - std::cos(phi) * std::cos(A) * std::sin(a))) /
                    std::cos(dec);

    //...Equation 5.85:
    double RAdot = wE + (Adot * std::cos(A) * std::cos(a) -
                         adot * std::sin(A) * std::sin(a) +
                         decdot * std::sin(A) * std::cos(a) * std::tan(dec)) /
                            (std::cos(phi) * std::sin(a) - std::sin(phi) * std::cos(A) * std::cos(a));

    //...Equation 5.69 and 5.72:
    std::vector<double> Rhodot = 
    {
        -RAdot * std::sin(RA) * std::cos(dec) - decdot * std::cos(RA) * std::sin(dec),
        RAdot * std::cos(RA) * std::cos(dec) - decdot * std::sin(RA) * std::sin(dec),
        decdot * std::cos(dec)
    };

    //...Equation 5.64:
    std::vector<double> v = 
    {
        Rdot[0] + rhodot * Rho[0] + rho * Rhodot[0],
        Rdot[1] + rhodot * Rho[1] + rho * Rhodot[1],
        Rdot[2] + rhodot * Rho[2] + rho * Rhodot[2]
    };

    return {r, v};
}

#endif // RV_FROM_OBSERVE_H
