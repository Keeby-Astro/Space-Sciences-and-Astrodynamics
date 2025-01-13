#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include "rkf45.h"
#include "atmosisa.h"

/*
  This program numerically integrates Equations 13.6 through
  13.8 for a gravity turn trajectory.
 
  h-functions required:      atmosisa
  User h-functions required: rkf45
  User subfunction requred:  rates
*/
const double PI = 3.141592653589793;
const double DEG_TO_RAD = PI / 180.0;
const double G0 = 9.81;                     // Sea-level acceleration of gravity (m/s^2)
const double RE = 6378e3;                   // Radius of the Earth (m)
const double H_SCALE = 7.5e3;               // Density scale height (m)
const double RHO0 = 1.225;                  // Sea-level density of atmosphere (kg/m^3)

const double DIAM = 196.85 / 12.0 * 0.3048; // Vehicle diameter (m)
const double A = PI / 4.0 * DIAM * DIAM;    // Frontal area (m^2)
const double CD = 0.5;                      // Drag coefficient
const double M0 = 149912 * 0.4536;          // Lift-off mass (kg)
const double MASS_RATIO = 7.0;              // Mass ratio
const double T2W = 1.4;                     // Thrust-to-weight ratio
const double ISP = 390.0;                   // Specific impulse (s)

const double M_FINAL = M0 / MASS_RATIO;     // Burnout mass (kg)
const double THRUST = T2W * M0 * G0;        // Rocket thrust (N)
const double M_DOT = THRUST / (ISP * G0);   // Propellant mass flow rate (kg/s)
const double M_PROP = M0 - M_FINAL;         // Propellant mass (kg)
const double T_BURN = M_PROP / M_DOT;       // Burn time (s)
const double H_TURN = 130.0;                // Height at which pitchover begins (m)

std::vector<double> rates(double t, const std::vector<double>& y);
void output_results(const std::vector<double>& t, const std::vector<std::vector<double>>& f);

int main() 
{
    //...Initial conditions
    double v0 = 0.0;                          // Initial velocity (m/s)
    double gamma0 = 89.85 * DEG_TO_RAD;       // Initial flight path angle (rad)
    double x0 = 0.0;                          // Initial downrange distance (m)
    double h0 = 0.0;                          // Initial altitude (m)
    double vD0 = 0.0;                         // Initial velocity loss due to drag (m/s)
    double vG0 = 0.0;                         // Initial velocity loss due to gravity (m/s)

    //...Initial conditions vector:
    std::vector<double> f0 = {v0, gamma0, x0, h0, vD0, vG0};
    std::pair<double, double> tspan = {0.0, T_BURN};

    //...Call to Runge-Kutta numerical integrator 'rkf45'
    //   rkf45 solves the system of equations df/dt = f(t):
    RKF45 result = rkf45(rates, tspan, f0);

    //...Output results
    output_results(result.tout, result.yout);

    return 0;
}

std::vector<double> rates(double t, const std::vector<double>& y) 
{
    /*
        Calculates the time rates dy/dt of the variables y(t) 
        in the equations of motion of a gravity turn trajectory.
    */
    //...Initialize dydt as a column vector:
    std::vector<double> dydt(6, 0.0);

    double v = y[0];         // Velocity
    double gamma = y[1];     // Flight path angle
    double x = y[2];         // Downrange distance
    double h = y[3];         // Altitude
    double vD = y[4];        // Velocity loss due to drag
    double vG = y[5];        // Velocity loss due to gravity

    // Determine thrust and mass
    double T = (t < T_BURN) ? THRUST : 0.0;
    double m = (t < T_BURN) ? (M0 - M_DOT * t) : (M0 - M_DOT * T_BURN);

    // Gravitational acceleration and air density
    double g = G0 / std::pow(1.0 + h / RE, 2);
    double rho = RHO0 * std::exp(-h / H_SCALE);
    double D = 0.5 * rho * v * v * A * CD; // Drag

    //...Define the first derivatives of v, gamma, x, h, vD and vG
    //   ("dot" means time derivative):
    //v_dot = T/m - D/m - g*sin(gamma); % Equation 13.6
    
    //...Start the gravity turn when h = hturn:
    if (h <= H_TURN) 
    {
        dydt[0] = T / m - D / m - g;       // v_dot
        dydt[1] = 0.0;                     // gamma_dot
        dydt[2] = 0.0;                     // x_dot
        dydt[3] = v;                       // h_dot
        dydt[4] = -D / m;                  // vD_dot
        dydt[5] = -g;                      // vG_dot
    } 
    else 
    {
        dydt[0] = T / m - D / m - g * std::sin(gamma);       
        dydt[1] = -(1.0 / v) * (g - v * v / (RE + h)) * std::cos(gamma); // Equation 13.7
        dydt[2] = (RE / (RE + h)) * v * std::cos(gamma);                 // Equation 13.8(1)
        dydt[3] = v * std::sin(gamma);                                   // Equation 13.8(2)
        dydt[4] = -D / m;                                                // Gravity loss rate
        dydt[5] = -g * std::sin(gamma);                                  //  Equation 13.27(1)
    }

    return dydt;
}

void output_results(const std::vector<double>& t, const std::vector<std::vector<double>>& f) {
    std::cout << "\n\n -----------------------------------\n";
    std::cout << "\n Initial flight path angle = " << 89.85 << " deg ";
    std::cout << "\n Pitchover altitude        = " << H_TURN << " m ";
    std::cout << "\n Burn time                 = " << T_BURN << " s ";

    //...Maximum dynamic pressure and corresponding time, speed, altitude and
    //   Mach number:
    double maxQ = 0.0;
    double tQ = 0.0, vQ = 0.0, hQ = 0.0, MQ = 0.0;

    //...Dynamic pressure vs time:
    for (size_t i = 0; i < t.size(); ++i) 
    {
        double h = f[i][3]; // Altitude
        double v = f[i][0]; // Velocity

        double rho = RHO0 * std::exp(-h / H_SCALE);
        double q = 0.5 * rho * v * v;
        if (q > maxQ) 
        {
            maxQ = q;
            tQ = t[i];
            vQ = v;
            hQ = h;

            double T, a, P, rho_temp;
            atmosisa(h, T, a, P, rho_temp);
            MQ = v / a;
        }
    }

    std::cout << "\n Maximum dynamic pressure  = " << maxQ * 9.869e-6 << " atm ";
    std::cout << "\n    Time                   = " << tQ / 60.0 << " min ";
    std::cout << "\n    Speed                  = " << vQ * 1e-3 << " km/s ";
    std::cout << "\n    Altitude               = " << hQ * 1e-3 << " km ";
    std::cout << "\n    Mach Number            = " << MQ << " ";

    std::cout << "\n At burnout:";
    std::cout << "\n    Speed                  = " << f.back()[0] * 1e-3 << " km/s ";
    std::cout << "\n    Flight path angle      = " << f.back()[1] / DEG_TO_RAD << " deg ";
    std::cout << "\n    Altitude               = " << f.back()[3] * 1e-3 << " km ";
    std::cout << "\n    Downrange distance     = " << f.back()[2] * 1e-3 << " km ";
    std::cout << "\n    Drag loss              = " << -f.back()[4] * 1e-3 << " km/s ";
    std::cout << "\n    Gravity loss           = " << -f.back()[5] * 1e-3 << " km/s ";
    std::cout << "\n\n -----------------------------------\n";
}
