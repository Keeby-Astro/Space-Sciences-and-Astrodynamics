// ALGORITHM 10.2: CALCULATE THE GEOCENTRIC POSITION OF THE SUN AT A GIVEN EPOCH

#ifndef SOLAR_POSITION_H
#define SOLAR_POSITION_H

#include <cmath>
#include <array>
#define _USE_MATH_DEFINES   
#define M_PI 3.14159265358979323846

/*
    This function calculates the geocentric equatorial position vector
    of the sun, given the julian date.
    
    User h-functions required: None
*/
void solar_position(double jd, double &lamda, double &eps, std::array<double, 3> &r_S) 
{
    //...Astronomical unit (km):
    const double AU = 149597870.691;

    //...Julian days since J2000:
    double n = jd - 2451545.0;

    //...Julian centuries since J2000:
    double cy = n / 36525.0;

    //...Mean anomaly (deg):
    double M = 357.528 + 0.9856003 * n;
    M = fmod(M, 360.0);

    //...Mean longitude (deg):
    double L = 280.460 + 0.98564736 * n;
    L = fmod(L, 360.0);

    //...Apparent ecliptic longitude (deg):
    lamda = L + 1.915 * sin(M * M_PI / 180.0) + 0.020 * sin(2 * M * M_PI / 180.0);
    lamda = fmod(lamda, 360.0);

    //...Obliquity of the ecliptic (deg):
    eps = 23.439 - 0.0000004 * n;

    //...Unit vector from Earth to Sun:
    double u_x = cos(lamda * M_PI / 180.0);
    double u_y = sin(lamda * M_PI / 180.0) * cos(eps * M_PI / 180.0);
    double u_z = sin(lamda * M_PI / 180.0) * sin(eps * M_PI / 180.0);

    //...Distance from Earth to Sun (km):
    double rS = (1.00014 - 0.01671 * cos(M * M_PI / 180.0) - 0.000140 * cos(2 * M * M_PI / 180.0)) * AU;

    //...Geocentric position vector (km):
    r_S = {rS * u_x, rS * u_y, rS * u_z};
}

#endif // SOLAR_POSITION_H