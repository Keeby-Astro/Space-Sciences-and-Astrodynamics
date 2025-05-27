#ifndef COE_FROM_SV_H
#define COE_FROM_SV_H

#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>

/*
    This function computes the classical orbital elements (coe)
    from the state vector (R,V) using Algorithm 4.1.

    mu      - gravitational parameter (km^3/s^2)
    R       - position vector in the geocentric equatorial frame (km)
    V       - velocity vector in the geocentric equatorial frame (km)
    r, v    - the magnitudes of R and V
    vr      - radial velocity component (km/s)
    H       - the angular momentum vector (km^2/s)
    h       - the magnitude of H (km^2/s)
    incl    - inclination of the orbit (rad)
    N       - the node line vector (km^2/s)
    n       - the magnitude of N
    cp      - cross product of N and R
    RA      - right ascension of the ascending node (rad)
    E       - eccentricity vector
    e       - eccentricity (magnitude of E)
    eps     - a small number below which the eccentricity is considered
              to be zero
    w       - argument of perigee (rad)
    TA      - true anomaly (rad)
    a       - semimajor axis (km)
    pi      - 3.1415926...
    coe     - vector of orbital elements [h e RA incl w TA a]

    User h-functions required: None
*/

std::vector<double> coe_from_sv(const std::vector<double>& R, const std::vector<double>& V, double mu) 
{
    if (R.size() != 3 || V.size() != 3) 
    {
        throw std::invalid_argument("Input vectors R and V must each have exactly three elements.");
    }

    constexpr double eps = 1.e-10;
    constexpr double pi = 3.141592653589793;

    double r = std::sqrt(R[0]*R[0] + R[1]*R[1] + R[2]*R[2]);
    double v = std::sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2]);
    double vr = (R[0]*V[0] + R[1]*V[1] + R[2]*V[2]) / r;

    std::vector<double> H = 
    {
        R[1]*V[2] - R[2]*V[1],
        R[2]*V[0] - R[0]*V[2],
        R[0]*V[1] - R[1]*V[0]
    };

    double h = std::sqrt(H[0]*H[0] + H[1]*H[1] + H[2]*H[2]);

    //...Equation 4.7:
    double incl = std::acos(H[2] / h);

    //...Equation 4.8:
    std::vector<double> N = {-H[1], H[0], 0.0};
    double n = std::sqrt(N[0]*N[0] + N[1]*N[1]);

    //...Equation 4.9:
    double RA = (n != 0) ? std::acos(N[0] / n) : 0.0;
    if (n != 0 && N[1] < 0) 
    {
        RA = 2 * pi - RA;
    }

    //...Equation 4.10:
    std::vector<double> E = 
    {
        (v*v - mu/r)*R[0] - r*vr*V[0],
        (v*v - mu/r)*R[1] - r*vr*V[1],
        (v*v - mu/r)*R[2] - r*vr*V[2]
    };
    for (auto& e : E) e /= mu;
    double e = std::sqrt(E[0]*E[0] + E[1]*E[1] + E[2]*E[2]);

    //...Equation 4.12 (incorporating the case e = 0):
    double w = 0.0;
    if (n != 0 && e > eps) 
    {
        w = std::acos((N[0]*E[0] + N[1]*E[1]) / (n*e));
        if (E[2] < 0) 
        {
            w = 2 * pi - w;
        }
    }

    //...Equation 4.13a (incorporating the case e = 0):
    double TA = 0.0;
    if (e > eps) 
    {
        TA = std::acos((E[0]*R[0] + E[1]*R[1] + E[2]*R[2]) / (e*r));
        if (vr < 0) 
        {
            TA = 2 * pi - TA;
        }
    } 
    else 
    {
        double cp_z = N[0]*R[1] - N[1]*R[0];
        TA = (cp_z >= 0) ? std::acos((N[0]*R[0] + N[1]*R[1]) / (n*r))
                         : 2 * pi - std::acos((N[0]*R[0] + N[1]*R[1]) / (n*r));
    }

    //...Equation 4.62 (a < 0 for a hyperbola):
    double a = (1 - e*e != 0) ? h*h / mu / (1 - e*e) : 0.0;

    return {h, e, RA, incl, w, TA, a};
}

#endif // COE_FROM_SV_H
