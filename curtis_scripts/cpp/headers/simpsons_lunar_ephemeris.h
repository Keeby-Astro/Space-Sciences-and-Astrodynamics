#ifndef SIMPSONS_LUNAR_EPHEMERIS_H
#define SIMPSONS_LUNAR_EPHEMERIS_H

#include <cmath>
#include <array>

/*
    David G. Simpson, "An Alternative Ephemeris Model for
    On-Board Flight Software Use," Proceedings of the 1999 Flight Mechanics
    Symposium, NASA Goddard Space Flight Center, pp. 175 - 184.
    
    This function computes the state vector of the moon at a given time
    relative to the earth's geocentric equatorial frame using a curve fit
    to JPL's DE200 (1982) ephemeris model.
    
    jd      - julian date (days)
    pos     - position vector (km)
    vel     - velocity vector (km/s)
    a       - matrix of amplitudes (km)
    b       - matrix of frequencies (rad/century)
    c       - matrix of phase angles (rad)
    t       - time in centuries since J2000
    tfac    - no. of seconds in a Julian century (36525 days)

    User h-functions required: None
*/
inline void simpsons_lunar_ephemeris(double jd, std::array<double, 3>& pos, std::array<double, 3>& vel) {
    const double tfac = 36525.0 * 3600.0 * 24.0;
    const double t = (jd - 2451545.0) / 36525.0;

    const double a[3][7] = {
        {383.0, 31.5, 10.6, 6.2, 3.2, 2.3, 0.8},
        {351.0, 28.9, 13.7, 9.7, 5.7, 2.9, 2.1},
        {153.2, 31.5, 12.5, 4.2, 2.5, 3.0, 1.8}
    };

    const double b[3][7] = {
        {8399.685, 70.990, 16728.377, 1185.622, 7143.070, 15613.745, 8467.263},
        {8399.687, 70.997, 8433.466, 16728.380, 1185.667, 7143.058, 15613.755},
        {8399.672, 8433.464, 70.996, 16728.364, 1185.645, 104.881, 8399.116}
    };

    const double c[3][7] = {
        {5.381, 6.169, 1.453, 0.481, 5.017, 0.857, 1.010},
        {3.811, 4.596, 4.766, 6.165, 5.164, 0.300, 5.565},
        {3.807, 1.629, 4.595, 6.162, 5.167, 2.555, 6.248}
    };

    for (int i = 0; i < 3; ++i) {
        pos[i] = 0.0;
        vel[i] = 0.0;
        for (int j = 0; j < 7; ++j) {
            pos[i] += a[i][j] * 1e3 * std::sin(b[i][j] * t + c[i][j]);
            vel[i] += a[i][j] * 1e3 * std::cos(b[i][j] * t + c[i][j]) * b[i][j];
        }
        vel[i] /= tfac;
    }
}

#endif // SIMPSONS_LUNAR_EPHEMERIS_H