// ALGORITHM 10.4: CALCULATE THE GEOCENTRIC POSITION OF THE MOON AT A GIVEN EPOCH

#ifndef LUNAR_POSITION_H
#define LUNAR_POSITION_H

#include <cmath>
#include <array>
#define _USE_MATH_DEFINES
#define M_PI 3.14159265358979323846

/*
    Calculates the geocentric equatorial position vector of the moon
    given the Julian day.
    
    User h-functions required: None
*/
void lunar_position(double jd, std::array<double, 3> &r_moon) {
    // Earth's radius (km)
    const double RE = 6378.14;

    // ------------------------- implementation -----------------
    //...Time in centuries since J2000:
    double T = (jd - 2451545.0) / 36525.0;

    //...Ecliptic longitude (deg):
    double e_long = 218.32 + 481267.881 * T
        + 6.29 * sin((135.0 + 477198.87 * T) * M_PI / 180.0)
        - 1.27 * sin((259.3 - 413335.36 * T) * M_PI / 180.0)
        + 0.66 * sin((235.7 + 890534.22 * T) * M_PI / 180.0)
        + 0.21 * sin((269.9 + 954397.74 * T) * M_PI / 180.0)
        - 0.19 * sin((357.5 + 35999.05 * T) * M_PI / 180.0)
        - 0.11 * sin((186.5 + 966404.03 * T) * M_PI / 180.0);
    e_long = fmod(e_long, 360.0);

    //...Ecliptic latitude (deg):
    double e_lat = 5.13 * sin((93.3 + 483202.02 * T) * M_PI / 180.0)
        + 0.28 * sin((228.2 + 960400.89 * T) * M_PI / 180.0)
        - 0.28 * sin((318.3 + 6003.15 * T) * M_PI / 180.0)
        - 0.17 * sin((217.6 - 407332.21 * T) * M_PI / 180.0);
    e_lat = fmod(e_lat, 360.0);

    //...Horizontal parallax (deg):
    double h_par = 0.9508
        + 0.0518 * cos((135.0 + 477198.87 * T) * M_PI / 180.0)
        + 0.0095 * cos((259.3 - 413335.36 * T) * M_PI / 180.0)
        + 0.0078 * cos((235.7 + 890534.22 * T) * M_PI / 180.0)
        + 0.0028 * cos((269.9 + 954397.74 * T) * M_PI / 180.0);
    h_par = fmod(h_par, 360.0);

    //...Angle between earth's orbit and its equator (deg):
    double obliquity = 23.439291 - 0.0130042 * T;

    //...Direction cosines of the moon's geocentric equatorial position vector:
    double l = cos(e_lat * M_PI / 180.0) * cos(e_long * M_PI / 180.0);
    double m = cos(obliquity * M_PI / 180.0) * cos(e_lat * M_PI / 180.0) * sin(e_long * M_PI / 180.0)
             - sin(obliquity * M_PI / 180.0) * sin(e_lat * M_PI / 180.0);
    double n = sin(obliquity * M_PI / 180.0) * cos(e_lat * M_PI / 180.0) * sin(e_long * M_PI / 180.0)
             + cos(obliquity * M_PI / 180.0) * sin(e_lat * M_PI / 180.0);

    //...Earth-moon distance (km)
    double dist = RE / sin(h_par * M_PI / 180.0);

    //...Moon's geocentric equatorial position vector (km):
    r_moon = {dist * l, dist * m, dist * n};
}

#endif // LUNAR_POSITION_H