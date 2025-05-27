#ifndef F_AND_G_TA_H
#define F_AND_G_TA_H

#define _USE_MATH_DEFINES
#include <vector>
#include <cmath>
#include <stdexcept>
#define M_PI 3.14159265358979323846

/*
    This function calculates the Lagrange f and g coefficients from the
    change in true anomaly since time t0.

    mu  - gravitational parameter (km^3/s^2)
    dt  - change in true anomaly (degrees)
    r0  - position vector at time t0 (km)
    v0  - velocity vector at time t0 (km/s)
    h   - angular momentum (km^2/s)
    vr0 - radial component of v0 (km/s)
    r   - radial position after the change in true anomaly
    f   - the Lagrange f coefficient (dimensionless)
    g   - the Lagrange g coefficient (s)

    User h-functions required: None
*/
void f_and_g_ta(const std::vector<double>& r0, const std::vector<double>& v0, double dt, double mu, double& f, double& g) {
    if (r0.size() != 3 || v0.size() != 3) {
        throw std::invalid_argument("Input vectors r0 and v0 must have exactly 3 elements.");
    }

    std::vector<double> h(3);
    h[0] = r0[1] * v0[2] - r0[2] * v0[1];
    h[1] = r0[2] * v0[0] - r0[0] * v0[2];
    h[2] = r0[0] * v0[1] - r0[1] * v0[0];
    double h_norm = std::sqrt(h[0] * h[0] + h[1] * h[1] + h[2] * h[2]);

    double r0_norm = std::sqrt(r0[0] * r0[0] + r0[1] * r0[1] + r0[2] * r0[2]);
    double vr0 = (v0[0] * r0[0] + v0[1] * r0[1] + v0[2] * r0[2]) / r0_norm;

    double s = std::sin(dt * M_PI / 180.0);
    double c = std::cos(dt * M_PI / 180.0);

    //...Equation 2.152:
    double r = h_norm * h_norm / mu / (1.0 + (h_norm * h_norm / mu / r0_norm - 1.0) * c - h_norm * vr0 * s / mu);

    //...Equations 2.158a & b:
    f = 1.0 - mu * r * (1.0 - c) / (h_norm * h_norm);
    g = r * r0_norm * s / h_norm;
}

#endif // F_AND_G_TA_H