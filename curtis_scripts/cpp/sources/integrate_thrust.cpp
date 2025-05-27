#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "rkf45.h"
#include "coe_from_sv.h"
#include "rv_from_r0v0_ta.h"

/*
    This function uses rkf45 to numerically integrate Equation 6.26 during
    the delta-v burn and then find the apogee of the post-burn orbit.

    mu      - gravitational parameter (km^3/s^2)
    RE      - earth radius (km)
    g0      - sea level acceleration of gravity (m/s^2)
    T       - rated thrust of rocket engine (kN)
    Isp     - specific impulse of rocket engine (s)
    m0      - initial spacecraft mass (kg)
    r0      - initial position vector (km)
    v0      - initial velocity vector (km/s)
    t0      - initial time (s)
    t_burn  - rocket motor burn time (s)
    y0      - column vector containing r0, v0 and m0
    t       - column vector of the times at which the solution is found (s)
    y       - a matrix whose elements are:
                columns 1, 2 and 3:
                    The solution for the x, y and z components of the
                    position vector r at the times t
                columns 4, 5 and 6:
                    The solution for the x, y and z components of the
                    velocity vector v at the times t
                column 7:
                    The spacecraft mass m at the times t
    r1      - position vector after the burn (km)
    v1      - velocity vector after the burn (km/s)
    m1      - mass after the burn (kg)
    coe     - orbital elements of the post-burn trajectory
            (h e RA incl w TA a)
    ra      - position vector vector at apogee (km)
    va      - velocity vector at apogee (km)
    rmax    - apogee radius (km)

    User h-functions required: rkf45, coe_from_sv, rv_from_r0v0_ta
    User subfunctions required: rates, output
*/

//...Preliminaries:
const double deg = M_PI / 180.0;
const double mu = 398600.4418;
const double RE = 6378.14; 
const double g0 = 9.807; 
const double T = 10.0; 
const double Isp = 300.0; 

// Function prototypes
std::vector<double> rates(double t, const std::vector<double>& f);
void output(const std::vector<double>& r0, const std::vector<double>& v0, double m0,
            const std::vector<double>& r1, const std::vector<double>& v1, double m1,
            const std::vector<double>& ra, const std::vector<double>& va, double rmax,
            const std::vector<double>& coe, double t_burn);

int main() 
{
    //...Input data:
    std::vector<double> r0 = {RE + 480.0, 0.0, 0.0};
    std::vector<double> v0 = {0.0, 7.7102, 0.0};
    double m0 = 2000.0;
    double t0 = 0.0;
    double t_burn = 261.1127;
    //...end Input data

    //...Integrate the equations of motion over the burn time:
    std::vector<double> y0 = {r0[0], r0[1], r0[2], v0[0], v0[1], v0[2], m0};
    std::pair<double, double> tspan = {t0, t_burn};
    RKF45 result = rkf45(rates, tspan, y0, 1.e-16);
    const auto& y_end = result.yout.back();

    //...Compute the state vector and mass after the burn:
    std::vector<double> r1 = {y_end[0], y_end[1], y_end[2]};
    std::vector<double> v1 = {y_end[3], y_end[4], y_end[5]};
    double m1 = y_end[6];
    std::vector<double> coe = coe_from_sv(r1, v1, mu);
    double e = coe[1];
    double TA = coe[5];
    double a = coe[6];

    //...Find the state vector at apogee of the post-burn trajectory:
    double dtheta = (TA <= M_PI) ? M_PI - TA : 3 * M_PI - TA;
    auto [ra, va] = rv_from_r0v0_ta(r1, v1, dtheta / deg, mu);
    double rmax = sqrt(ra[0] * ra[0] + ra[1] * ra[1] + ra[2] * ra[2]);

    output(r0, v0, m0, r1, v1, m1, ra, va, rmax, coe, t_burn);

    return 0;
}

// Function to compute derivatives (rates of change)
std::vector<double> rates(double t, const std::vector<double>& f) 
{
    double x = f[0], y = f[1], z = f[2];
    double vx = f[3], vy = f[4], vz = f[5];
    double m = f[6];

    double r = sqrt(x * x + y * y + z * z);
    double v = sqrt(vx * vx + vy * vy + vz * vz);

    double ax = -mu * x / (r * r * r) + T / m * vx / v;
    double ay = -mu * y / (r * r * r) + T / m * vy / v;
    double az = -mu * z / (r * r * r) + T / m * vz / v;
    double mdot = -T * 1000.0 / g0 / Isp;

    return {vx, vy, vz, ax, ay, az, mdot};
}

void output(const std::vector<double>& r0, const std::vector<double>& v0, double m0,
            const std::vector<double>& r1, const std::vector<double>& v1, double m1,
            const std::vector<double>& ra, const std::vector<double>& va, double rmax,
            const std::vector<double>& coe, double t_burn) 
{
    std::cout << "\n\n---------------------------------------------------------\n";
    std::cout << "\nBefore ignition:";
    std::cout << "\n  Mass = " << m0 << " kg";
    std::cout << "\n  State vector:";
    std::cout << "\n    r = [" << r0[0] << ", " << r0[1] << ", " << r0[2] << "] (km)";
    std::cout << "\n      Radius = " << sqrt(r0[0] * r0[0] + r0[1] * r0[1] + r0[2] * r0[2]);
    std::cout << "\n    v = [" << v0[0] << ", " << v0[1] << ", " << v0[2] << "] (km/s)";
    std::cout << "\n      Speed = " << sqrt(v0[0] * v0[0] + v0[1] * v0[1] + v0[2] * v0[2]) << "\n";
    std::cout << "\nThrust          = " << T << " kN";
    std::cout << "\nBurn time       = " << t_burn << " s";
    std::cout << "\nMass after burn = " << m1 << " kg\n";
    std::cout << "\nEnd-of-burn-state vector:";
    std::cout << "\n    r = [" << r1[0] << ", " << r1[1] << ", " << r1[2] << "] (km)";
    std::cout << "\n      Radius = " << sqrt(r1[0] * r1[0] + r1[1] * r1[1] + r1[2] * r1[2]);
    std::cout << "\n    v = [" << v1[0] << ", " << v1[1] << ", " << v1[2] << "] (km/s)";
    std::cout << "\n      Speed = " << sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]) << "\n";
    std::cout << "\nPost-burn trajectory:";
    std::cout << "\n    Eccentricity = " << coe[1];
    std::cout << "\n  Semimajor axis = " << coe[6] << " km";
    std::cout << "\n  Apogee state vector:";
    std::cout << "\n    r = [" << ra[0] << ", " << ra[1] << ", " << ra[2] << "] (km)";
    std::cout << "\n      Radius = " << rmax;
    std::cout << "\n    v = [" << va[0] << ", " << va[1] << ", " << va[2] << "] (km/s)";
    std::cout << "\n      Speed = " << sqrt(va[0] * va[0] + va[1] * va[1] + va[2] * va[2]);
    std::cout << "\n\n---------------------------------------------------------\n\n";
}
