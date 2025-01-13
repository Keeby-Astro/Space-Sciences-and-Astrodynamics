#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "sv_from_coe.h"
#include "coe_from_sv.h"
#include "rv_from_r0v0.h"
#define _USE_MATH_DEFINES
#define M_PI 3.14159265358979323846

/*
    This function use Encke's method together with MATLAB's ode45
    to integrate Equation 10.2 for a J2 gravitational perturbation 
    given by Equation 10.30.
    
    User M-functions required: sv_from_coe, coe_from_sv, rv_from_r0v0
    User subfunction required: rates 
*/
//...Conversion factors:
constexpr double hours = 3600.0;        // Hours to seconds
constexpr double days = 24.0 * hours;   // Days to seconds
constexpr double deg = M_PI / 180.0;    // Degrees to radians

//...Constants:
constexpr double mu = 398600.4418;      // Graviational parameter (km^3/s^2)
constexpr double RE = 6378.14;          // Earth's radius (km)    
constexpr double J2 = 1082.63e-6;    

std::vector<double> rates(double t, const std::vector<double>& f, const std::vector<double>& R0, const std::vector<double>& V0, double t0) 
{
    /*
        This function calculates the time rates of Encke's deviation in position
        del_r and velocity del_v.
    */
    std::vector<double> del_r = {f[0], f[1], f[2]}; // Position deviation
    std::vector<double> del_v = {f[3], f[4], f[5]}; // Velocity deviation

    //...Compute the state vector on the osculating orbit at time t
    //   (Equation 12.5) using Algorithm 3.4:
    std::vector<double> Rosc(3), Vosc(3);
    rv_from_r0v0(R0, V0, t - t0, Rosc, Vosc, mu);

    //...Calculate the components of the state vector on the perturbed orbit
    //   and their magnitudes:
    std::vector<double> Rpp = {Rosc[0] + del_r[0], Rosc[1] + del_r[1], Rosc[2] + del_r[2]};
    std::vector<double> Vpp = {Vosc[0] + del_v[0], Vosc[1] + del_v[1], Vosc[2] + del_v[2]};
    double rpp = std::sqrt(Rpp[0] * Rpp[0] + Rpp[1] * Rpp[1] + Rpp[2] * Rpp[2]);

    //...Compute the J2 perturbing acceleration from Equation 12.30:
    double fac = 1.5 * J2 * (mu / (rpp * rpp)) * std::pow(RE / rpp, 2);
    std::vector<double> ap = 
    {
        -fac * (1 - 5 * std::pow(Rpp[2] / rpp, 2)) * (Rpp[0] / rpp),
        -fac * (1 - 5 * std::pow(Rpp[2] / rpp, 2)) * (Rpp[1] / rpp),
        -fac * (3 - 5 * std::pow(Rpp[2] / rpp, 2)) * (Rpp[2] / rpp)
    };

    //...Compute the total perturbing ecceleration from Equation 12.7:  
    double rosc = std::sqrt(Rosc[0] * Rosc[0] + Rosc[1] * Rosc[1] + Rosc[2] * Rosc[2]);
    double F = 1.0 - std::pow(rosc / rpp, 3);
    std::vector<double> del_a = 
    {
        -mu / (rosc * rosc * rosc) * (del_r[0] - F * Rpp[0]) + ap[0],
        -mu / (rosc * rosc * rosc) * (del_r[1] - F * Rpp[1]) + ap[1],
        -mu / (rosc * rosc * rosc) * (del_r[2] - F * Rpp[2]) + ap[2]
    };

    return {del_v[0], del_v[1], del_v[2], del_a[0], del_a[1], del_a[2]};
}

int main() 
{
    //...Initial orbital parameters (given):
    double zp0 = 300.0;                             // Perigee altitude (km)
    double za0 = 3062.0;                            // Apogee altitude (km)
    double RA0 = 45.0 * deg;                        // Right ascension of the node (radians)
    double i0 = 28.0 * deg;                         // Inclination (radians)
    double w0 = 30.0 * deg;                         // Argument of perigee (radians)
    double TA0 = 40.0 * deg;                        // True anomaly (radians)

    //...Initial orbital parameters (inferred):
    double rp0 = RE + zp0;                          // Perigee radius (km)
    double ra0 = RE + za0;                          // Apogee radius (km)
    double e0 = (ra0 - rp0) / (ra0 + rp0);          // Eccentricity
    double a0 = (ra0 + rp0) / 2.0;                  // Semi-major axis (km)
    double h0 = std::sqrt(rp0 * mu * (1.0 + e0));   // Angular momentum (km^2/s)
    double T0 = 2.0*M_PI*std::sqrt(a0*a0*a0 / mu);  // Period (s)

    double t0 = 0.0;                                // Initial time (s)
    double tf = 2 * days * hours;                   // Final time (s)
    //...end Input data

    // Store the initial orbital elements in the array coe0:
    std::vector<double> coe0 = {h0, e0, RA0, i0, w0, TA0};
    auto [R0, V0] = sv_from_coe(coe0, mu);

    double del_t = T0 / 100.0; // Time step for Encke procedure

    //...Begin the Encke integration:
    std::vector<double> del_y0(6, 0.0);         // Initialize the state vector perturbation
    double t = t0;                              // Initialize the time scalar
    std::vector<double> tsave = {t0};           // Initialize the vector of solution times
    std::vector<std::vector<double>> y = {R0};  // Initialize the state vector
    y.back().insert(y.back().end(), V0.begin(), V0.end());

    t += del_t; // First time step

    //   Loop over the time interval [t0, tf] with equal increments del_t:
    while (t <= tf + del_t / 2.0) 
    {
        // Integrate Equation 12.7 over the time increment del_t:
        std::vector<double> f = rates(t, del_y0, R0, V0, t0);
        for (size_t i = 0; i < 6; ++i) 
        {
            del_y0[i] += f[i] * del_t;
        }

        // Compute the osculating state vector at time t: 
        std::vector<double> Rosc(3), Vosc(3);
        rv_from_r0v0(R0, V0, t - t0, Rosc, Vosc, mu);
        for (size_t i = 0; i < 3; ++i) 
        {
            // Rectify:
            R0[i] = Rosc[i] + del_y0[i];
            V0[i] = Vosc[i] + del_y0[i + 3];
        }

        // Prepare for next time step:
        tsave.push_back(t);
        std::vector<double> state = R0;
        state.insert(state.end(), V0.begin(), V0.end());
        y.push_back(state);

        t0 = t;
        t += del_t;
    }
    // End the loop

    //...Print selected osculating elements:
    std::cout << "\nTime (s), Semi-major axis (km), Eccentricity, Inclination (deg), RAAN (deg), Argument of Perigee (deg), True Anomaly (deg)\n";
    for (size_t i = 0; i < tsave.size() && i < 500; ++i) 
    {
        auto [R, V] = std::make_pair(std::vector<double>(y[i].begin(), y[i].begin() + 3), std::vector<double>(y[i].begin() + 3, y[i].end()));
        auto coe = coe_from_sv(R, V, mu);
        std::cout << tsave[i];
        std::cout << ", " << coe[0];  // Semi-major axis
        std::cout << ", " << coe[1];  // Eccentricity
        std::cout << ", " << coe[3] / deg;  // Inclination
        std::cout << ", " << coe[2] / deg;  // RAAN
        std::cout << ", " << coe[4] / deg;  // Argument of Perigee
        std::cout << ", " << coe[5] / deg;  // True Anomaly
        std::cout << "\n";
    }

    return 0;
}