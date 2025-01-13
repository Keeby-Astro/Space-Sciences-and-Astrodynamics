#include <iostream>
#include <array>
#include <vector>
#include <cmath>
#include "rkf45.h"
#include "q_from_dcm.h"
#include "dcm_from_q.h"
#include "dcm_to_euler.h"

/*
  This program numerically integrates Euler's equations of motion
  for the spinning top. The  quaternion is used to obtain the time history
  of the top's orientation.
 
  User h-functions required: rkf45, q_from_dcm, dcm_from_q, dcm_to_euler
  User subfunction required: rates
*/
using namespace std;
vector<double> rates(double t, const vector<double>& f);
void plotit(const vector<double>& t, const vector<vector<double>>& f);

int main() {
    // Data declaration
    const double g = 9.807;  // Acceleration of gravity (m/s^2)
    const double m = 0.5;    // Mass in kg
    const double d = 0.05;   // Distance of center of mass from pivot point (m)
    const double A = 12e-4;  // Moment of inertia about body x (kg-m^2)
    const double B = 12e-4;  // Moment of inertia about body y (kg-m^2)
    const double C = 4.5e-4; // Moment of inertia about body z (kg-m^2)
    const double ws0 = 1000 * 2 * M_PI / 60;  // Spin rate (rad/s)

    double wp0 = 51.93 * 2 * M_PI / 60; // Precession rate (rad/s)
    wp0 = 0;  // Use to obtain a specific figure

    const double wn0 = 0;    // Nutation rate (rad/s)
    const double theta = 60; // Initial nutation angle (degrees)

    // Initial orientation (z-axis direction and x-axis direction defining x-z plane)
    array<double, 3> z = {0, -sin(theta * M_PI / 180), cos(theta * M_PI / 180)};
    array<double, 3> p = {1, 0, 0};

    // Compute orthonormal axes
    array<double, 3> y = {
        z[1] * p[2] - z[2] * p[1],
        z[2] * p[0] - z[0] * p[2],
        z[0] * p[1] - z[1] * p[0]
    };
    array<double, 3> x = {
        y[1] * z[2] - y[2] * z[1],
        y[2] * z[0] - y[0] * z[2],
        y[0] * z[1] - y[1] * z[0]
    };

    // Normalize vectors
    auto norm = [](const array<double, 3>& v) {
        return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    };
    for (auto& val : x) val /= norm(x);
    for (auto& val : y) val /= norm(y);
    for (auto& val : z) val /= norm(z);

    // Initial direction cosine matrix
    array<array<double, 3>, 3> QXx = {x, y, z};

    // Convert QXx to std::vector<std::vector<double>> for dcm_to_euler
    vector<vector<double>> QXx_vector = {
        {QXx[0][0], QXx[0][1], QXx[0][2]},
        {QXx[1][0], QXx[1][1], QXx[1][2]},
        {QXx[2][0], QXx[2][1], QXx[2][2]}
    };

    // Initial Euler angles (degrees) and quaternion
    auto [phi0, theta0, psi0] = dcm_to_euler(QXx_vector);
    array<double, 4> q0 = q_from_dcm(QXx);

    // Initial body-frame angular velocity (rad/s)
    array<double, 3> w0 = {
        wp0 * sin(theta0 * M_PI / 180) * sin(psi0 * M_PI / 180) + wn0 * cos(psi0 * M_PI / 180),
        wp0 * sin(theta0 * M_PI / 180) * cos(psi0 * M_PI / 180) - wn0 * sin(psi0 * M_PI / 180),
        ws0 + wp0 * cos(theta0 * M_PI / 180)
    };

    // Initial conditions
    vector<double> f0 = {q0[0], q0[1], q0[2], q0[3], w0[0], w0[1], w0[2]};

    double t0 = 0;       // Initial time (s)
    double tf = 1.153;   // Final time (s)

    // Solve the equations of motion using RKF45
    RKF45 result = rkf45(rates, {t0, tf}, f0);

    // Plot the results
    plotit(result.tout, result.yout);

    return 0;
}

// Function to compute time derivatives
vector<double> rates(double t, const vector<double>& f) {
    const double g = 9.807, m = 0.5, d = 0.05, A = 12e-4, B = 12e-4, C = 4.5e-4;

    array<double, 4> q = {f[0], f[1], f[2], f[3]};
    array<double, 3> w = {f[4], f[5], f[6]};

    // Normalize quaternion
    double norm_q = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
    for (double& val : q) val /= norm_q;

    // Direction cosine matrix from quaternion
    array<array<double, 3>, 3> Q = dcm_from_q(q);

    // Moments of the weight vector about the pivot point
    array<double, 3> M = {
        Q[0][2] * (-m * g * d),
        Q[1][2] * (m * g * d),
        0
    };

    // Quaternion derivative
    array<double, 4> q_dot = {
        0.5 * ( w[0] * q[3] - w[1] * q[2] + w[2] * q[1]),
        0.5 * ( w[0] * q[2] + w[1] * q[3] - w[2] * q[0]),
        0.5 * (-w[0] * q[1] + w[1] * q[0] + w[2] * q[3]),
        0.5 * (-w[0] * q[0] - w[1] * q[1] - w[2] * q[2])
    };

    // Angular velocity derivatives (Euler's equations)
    array<double, 3> w_dot = {
        (M[0] / A) - ((C - B) * w[1] * w[2] / A),
        (M[1] / B) - ((A - C) * w[2] * w[0] / B),
        (M[2] / C) - ((B - A) * w[0] * w[1] / C)
    };

    // Combine results
    vector<double> dfdt = {q_dot[0], q_dot[1], q_dot[2], q_dot[3], w_dot[0], w_dot[1], w_dot[2]};

    return dfdt;
}