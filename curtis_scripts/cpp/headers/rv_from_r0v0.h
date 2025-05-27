// ALGORITHM 3.4: CALCULATION OF THE STATE VECTOR GIVEN THE
// INITIAL STATE VECTOR AND THE TIME LAPSE ΔT

#ifndef RV_FROM_R0V0_H
#define RV_FROM_R0V0_H

#include <vector>
#include <cmath>
#include "kepler_U.h"
#include "f_and_g.h"
#include "fDot_and_gDot.h"

/*
    This function computes the state vector (R,V) from the
    initial state vector (R0,V0) and the elapsed time.

    mu - gravitational parameter (km^3/s^2)
    R0 - initial position vector (km)
    V0 - initial velocity vector (km/s)
    t  - elapsed time (s)
    R  - final position vector (km)
    V  - final velocity vector (km/s)

    User h-functions required: kepler_U, f_and_g, fDot_and_gDot
*/
inline void rv_from_r0v0(const std::vector<double>& R0, const std::vector<double>& V0,
                         double t, std::vector<double>& R, std::vector<double>& V, double mu) 
{

    //...Magnitudes of R0 and V0:
    double r0 = std::sqrt(R0[0] * R0[0] + R0[1] * R0[1] + R0[2] * R0[2]);
    double v0 = std::sqrt(V0[0] * V0[0] + V0[1] * V0[1] + V0[2] * V0[2]);

    //...Initial radial velocity:
    double vr0 = (R0[0] * V0[0] + R0[1] * V0[1] + R0[2] * V0[2]) / r0;

    //...Reciprocal of the semimajor axis (from the energy equation):
    double alpha = 2.0 / r0 - v0 * v0 / mu;

    //...Compute the universal anomaly:
    double x = kepler_U(t, r0, vr0, alpha, mu);

    //...Compute the f and g functions:
    double f, g;
    f_and_g(x, t, r0, alpha, f, g, mu);

    //...Compute the final position vector:
    R[0] = f * R0[0] + g * V0[0];
    R[1] = f * R0[1] + g * V0[1];
    R[2] = f * R0[2] + g * V0[2];

    //...Compute the magnitude of R:
    double r = std::sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2]);

    //...Compute the derivatives of f and g:
    auto [fdot, gdot] = fDot_and_gDot(x, r, r0, alpha, mu);

    //...Compute the final velocity:
    V[0] = fdot * R0[0] + gdot * V0[0];
    V[1] = fdot * R0[1] + gdot * V0[1];
    V[2] = fdot * R0[2] + gdot * V0[2];
}

#endif // RV_FROM_R0V0_H
