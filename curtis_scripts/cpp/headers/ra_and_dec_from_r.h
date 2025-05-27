// ALGORITHM 4.1: OBTAIN THE RIGHT ASCENSION AND DECLINATION
// FROM THE POSITION VECTOR

#ifndef RA_AND_DEC_FROM_R_H
#define RA_AND_DEC_FROM_R_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <utility>
#define M_PI 3.14159265358979323846

/*
    This function calculates the right ascension and the
    declination from the geocentric equatorial position vector.

    r       - position vector
    l, m, n - direction cosines of r
    ra      - right ascension (degrees)
    dec     - declination (degrees)
*/
std::pair<double, double> ra_and_dec_from_r(const std::vector<double>& r)
{
    double norm_r = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);

    double l = r[0] / norm_r;
    double m = r[1] / norm_r;
    double n = r[2] / norm_r;

    double dec = std::asin(n) * 180.0 / M_PI;

    double ra;
    if (m > 0)
    {
        ra = std::acos(l / std::cos(dec * M_PI / 180.0)) * 180.0 / M_PI;
    }
    else
    {
        ra = 360.0 - std::acos(l / std::cos(dec * M_PI / 180.0)) * 180.0 / M_PI;
    }

    return {ra, dec};
}

#endif // RA_AND_DEC_FROM_R_H