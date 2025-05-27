#ifndef FDOT_AND_GDOT_TA_H
#define FDOT_AND_GDOT_TA_H

#define _USE_MATH_DEFINES
#include <vector>
#include <cmath>

/*
    This function calculates the time derivatives of the Lagrange
    f and g coefficients from the change in true anomaly since time t0.

    mu      - gravitational parameter (km^3/s^2)
    dt      - change in true anomaly (degrees)
    r0      - position vector at time t0 (km)
    v0      - velocity vector at time t0 (km/s)
    h       - angular momentum (km^2/s)
    vr0     - radial component of v0 (km/s)
    fdot    - time derivative of the Lagrange f coefficient (1/s)
    gdot    - time derivative of the Lagrange g coefficient (dimensionless)

    User M-functions required: None
*/

std::vector<double> cross(const std::vector<double>& a, const std::vector<double>& b) {
    return {a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]};
}

double dot(const std::vector<double>& a, const std::vector<double>& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double norm(const std::vector<double>& v) {
    return std::sqrt(dot(v, v));
}

// Function to compute the time derivatives of Lagrange coefficients
void fDot_and_gDot_ta(const std::vector<double>& r0, const std::vector<double>& v0,
                      double dt, double mu, double& fdot, double& gdot) 
{
    #define M_PI 3.14159265358979323846

    double h = norm(cross(r0, v0));
    double vr0 = dot(v0, r0) / norm(r0);
    double r0_norm = norm(r0);
    double dt_rad = dt * M_PI / 180.0;
    double c = std::cos(dt_rad);
    double s = std::sin(dt_rad);

    //...Equations 2.158c & d:
    fdot = mu / h * (vr0 / h * (1 - c) - s / r0_norm);
    gdot = 1 - mu * r0_norm / (h * h) * (1 - c);
}

#endif // FDOT_AND_GDOT_TA_H