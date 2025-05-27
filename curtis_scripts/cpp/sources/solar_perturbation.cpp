#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include "sv_from_coe.h"
#include "solar_position.h"
#include "rsmooth.h"
#define M_PI 3.14159265358979323846

/*
    This function uses ODE to integrate
    Equations 10.84, the Gauss variational equations, for a solar
    gravitational perturbation.
    
    User H-functions required:  sv_from_coe, solar_position
    User subfunctions required: solveit rates
*/
double deg_to_rad(double deg) 
{
    return deg * M_PI / 180.0;
}

double days_to_seconds(double days) 
{
    return days * 24 * 3600;
}

struct solar
{
    double h, e, RA, i, w, TA;
};

std::array<double, 6> rates(double t, const solar&coe, double mu, double mu3, const std::array<double, 3> &R_S) 
{
    //...The orbital elements at time t:
    double h = coe.h;
    double e = coe.e;
    double RA = coe.RA;
    double i = coe.i;
    double w = coe.w;
    double TA = coe.TA;
    double phi = w + TA;

    //...Obtain the state vector at time t from Algorithm 4.5:
    auto [R, V] = sv_from_coe({h, e, RA, i, w, TA}, mu);

    //...Obtain the unit vectors of the rsw system:
    double r = std::sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2]);
    std::array<double, 3> ur = {R[0] / r, R[1] / r, R[2] / r};

    //...Update the Julian Day:
    double JD = 2454283.0 + t / days_to_seconds(1);

    //...Find and normalize the position vector of the sun:
    double r_S = std::sqrt(R_S[0] * R_S[0] + R_S[1] * R_S[1] + R_S[2] * R_S[2]);

    std::array<double, 3> R_rel = {R_S[0] - R[0], R_S[1] - R[1], R_S[2] - R[2]};
    double r_rel = std::sqrt(R_rel[0] * R_rel[0] + R_rel[1] * R_rel[1] + R_rel[2] * R_rel[2]);

    //...See Appendix F:
    double q = (R[0] * (2 * R_S[0] - R[0]) + R[1] * (2 * R_S[1] - R[1]) + R[2] * (2 * R_S[2] - R[2])) / (r_S * r_S);
    double F = (q * q - 3 * q + 3) * q / (1 + std::pow(1 - q, 1.5));

    //...Perturbation components in the rsw system:
    std::array<double, 3> ap = 
    {
        mu3 / std::pow(r_rel, 3) * (F * R_S[0] - R[0]),
        mu3 / std::pow(r_rel, 3) * (F * R_S[1] - R[1]),
        mu3 / std::pow(r_rel, 3) * (F * R_S[2] - R[2])
    };
    double apr = ap[0] * ur[0] + ap[1] * ur[1] + ap[2] * ur[2];

    //...Gauss variational equations (Equations 12.84):
    double hdot = r * (ap[0] * ur[0] + ap[1] * ur[1] + ap[2] * ur[2]);
    double edot = h / mu * std::sin(TA) * apr;
    double RAdot = r / (h * std::sin(i)) * std::sin(phi);
    double idot = r / h * std::cos(phi);
    double wdot = -h * std::cos(TA) / (mu * e) * apr;
    double TAdot = h / (r * r);

    //...Return rates to ode45 in the array dfdt:  
    return {hdot, edot, RAdot, idot, wdot, TAdot};
}

void solve_orbit(const std::string &type, double a0, double e0, double i0_deg, double w0_deg, double RA0_deg, double TA0_deg, double mu, double mu3, double RE) 
{
    /*
        Calculations and plots common to all of the orbits
    */
    //...Initial orbital parameters (calculated from the given data):
    double i0 = deg_to_rad(i0_deg);
    double w0 = deg_to_rad(w0_deg);
    double RA0 = deg_to_rad(RA0_deg);
    double TA0 = deg_to_rad(TA0_deg);

    double h0 = std::sqrt(mu * a0 * (1 - e0 * e0));
    double T0 = 2 * M_PI / std::sqrt(mu) * std::pow(a0, 1.5);

    //...Store initial orbital elements (from above) in the vector coe0:
    solar coe0 = {h0, e0, RA0, i0, w0, TA0};

    //...Use ODE to integrate the Equations 12.84, the Gauss variational
    //   equations with lunar gravity as the perturbation, from t0 to tf:
    double t0 = 0;
    double tf = days_to_seconds(720);
    std::vector<double> tspan(400);
    for (int i = 0; i < 400; ++i) 
    {
        tspan[i] = t0 + i * (tf - t0) / 399;
    }

    std::cout << "type   t   h   e   RA   i   w   TA" << std::endl;
    for (double t : tspan) 
    {
        double lamda, eps;
        std::array<double, 3> R_S;
        solar_position(t / days_to_seconds(1), lamda, eps, R_S);

        auto dfdt = rates(t, coe0, mu, mu3, R_S);

        // Output results
        std::cout << type << " " << t << " " << dfdt[0] << " " << dfdt[1] << " " << dfdt[2] << " " << dfdt[3] << " " << dfdt[4] << " " << dfdt[5] << std::endl;

        // Update the orbital elements
        coe0.h += dfdt[0] * (tf - t0) / 399;
        coe0.e += dfdt[1] * (tf - t0) / 399;
        coe0.RA += dfdt[2] * (tf - t0) / 399;
        coe0.i += dfdt[3] * (tf - t0) / 399;
        coe0.w += dfdt[4] * (tf - t0) / 399;
        coe0.TA += dfdt[5] * (tf - t0) / 399;
        t0 = t;
    }
}

int main() 
{
    double mu = 398600.0;
    double mu3 = 132712e9;
    double RE = 6378.0;

    solve_orbit("GEO", 42164, 0.0001, 1, 0, 0, 0, mu, mu3, RE);
    solve_orbit("HEO", 26553.4, 0.741, 63.4, 270, 0, 0, mu, mu3, RE);
    solve_orbit("LEO", 6678.136, 0.01, 28.5, 0, 0, 0, mu, mu3, RE);

    return 0;
}
