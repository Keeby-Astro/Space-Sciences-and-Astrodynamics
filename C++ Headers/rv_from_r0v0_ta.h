#ifndef RV_FROM_R0V0_TA_H
#define RV_FROM_R0V0_TA_H

#include <vector>
#include "f_and_g_ta.h"
#include "fDot_and_gDot_ta.h"
/*
    This function computes the state vector (r,v) from the
    initial state vector (r0,v0) and the change in true anomaly.

    mu - gravitational parameter (km^3/s^2)
    r0 - initial position vector (km)
    v0 - initial velocity vector (km/s)
    dt - change in true anomaly (degrees)
    r  - final position vector (km)
    v  - final velocity vector (km/s)

    User h-functions required: f_and_g_ta, fDot_and_gDot_ta
*/
inline std::pair<std::vector<double>, std::vector<double>> rv_from_r0v0_ta(
    const std::vector<double>& r0,
    const std::vector<double>& v0,
    double dt,
    double mu) {
    double f, g, fdot, gdot;

    //...Compute the f and g functions and their derivatives:
    f_and_g_ta(r0, v0, dt, mu, f, g);
    fDot_and_gDot_ta(r0, v0, dt, mu, fdot, gdot);

    //...Compute the final position and velocity vectors:
    std::vector<double> r(3, 0.0);
    std::vector<double> v(3, 0.0);
    for (size_t i = 0; i < 3; ++i) {
        r[i] = f * r0[i] + g * v0[i];
        v[i] = fdot * r0[i] + gdot * v0[i];
    }

    return {r, v};
}

#endif // RV_FROM_R0V0_TA_H
