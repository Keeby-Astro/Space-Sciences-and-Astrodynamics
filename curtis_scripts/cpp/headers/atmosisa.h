#ifndef ATMOSISA_H
#define ATMOSISA_H

#include <cmath>
#include <stdexcept>
#include <vector>

/*
    Returns standard atmosphere properties using an approximation
    to the International Standard Atmosphere (ISA) model for a
    given geometric altitude h_m (meters).

    T   - Temperature (K)
    a   - Speed of sound (m/s)
    P   - Pressure (Pa)
    rho - Density (g/m^3)

    User h-function required: none
*/
struct Layer {
    double h_base;
    double h_top;
    double T_base;
    double P_base;
    double L;
};

inline void atmosisa(double h_m, double& T, double& a, double& P, double& rho) {
    const double g0 = 9.80665;  // Gravity (m/s^2)
    const double R = 287.053;   // Specific gas constant (J/(kg*K))
    const double gamma = 1.4;   // Adiabatic index

    // Define ISA layers
    std::vector<Layer> layers = {
        {0.0,       11000.0, 288.15, 101325.0, -0.0065},
        {11000.0,   20000.0, 216.65, 22632.06,  0.0},
        {20000.0,   32000.0, 216.65, 5474.89,   0.0010},
        {32000.0,   47000.0, 228.65, 868.02,    0.0028},
        {47000.0,   51000.0, 270.65, 110.91,    0.0},
        {51000.0,   71000.0, 270.65, 66.94,    -0.0028},
        {71000.0,   84852.0, 214.65, 3.96,     -0.0020}
    };

    // Validate altitude
    if (h_m < 0.0) {
        throw std::invalid_argument("Altitude cannot be negative for this standard model.");
    }
    if (h_m > 84852.0) {
        h_m = 84852.0;
    }

    // Find the appropriate layer
    Layer layer;
    for (const auto& l : layers) {
        if (h_m <= l.h_top) {
            layer = l;
            break;
        }
    }

    // Compute properties
    double h_b = layer.h_base;
    double T_b = layer.T_base;
    double P_b = layer.P_base;
    double L = layer.L;

    double delta_h = h_m - h_b;

    if (std::abs(L) > 1e-10) {
        T = T_b + L * delta_h;
        double exponent = -g0 / (R * L);
        P = P_b * std::pow(T / T_b, exponent);
    } else {
        T = T_b;
        P = P_b * std::exp(-g0 * delta_h / (R * T_b));
    }

    rho = P / (R * T);
    a = std::sqrt(gamma * R * T);
}

#endif // ATMOSISA_H
