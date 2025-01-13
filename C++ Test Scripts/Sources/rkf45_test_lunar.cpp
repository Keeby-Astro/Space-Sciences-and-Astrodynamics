#include <iostream>
#include <vector>

#define _USE_MATH_DEFINES
#include <cmath>
#include <iomanip>
#include "rkf45.h"

/*
    This program uses the Runge-Kutta-Fehlberg 4(5) method to solve the
    earth-moon restricted three-body problem (Equations 2.192a and 2.192b)
    for the trajectory of a spacecraft.

    The numerical integration is done in the external function 'rkf45',
    which uses the subfunction 'rates' herein to compute the derivatives.

    days      - converts days to seconds
    G         - universal graviational constant (km^3/kg/s^2)
    rmoon     - radius of the moon (km)
    rearth    - radius of the earth (km)
    r12       - distance from center of earth to center of moon (km)
    m1,m2     - masses of the earth and of the moon, respectively (kg)
    M         - total mass of the restricted 3-body system (kg)
    mu        - gravitational parameter of earth-moon system (km^3/s^2)
    mu1,mu2   - gravitational parameters of the earth and of the moon,
                respectively (km^3/s^2)
    pi_1,pi_2 - ratios of the earth mass and the moon mass, respectively,
                to the total earth-moon mass
    W         - angular velocity of moon around the earth (rad/s)
    x1,x2     - x-coordinates of the earth and of the moon, respectively,
                relative to the earth-moon barycenter (km)
    d0        - initial altitude of spacecraft (km)
    phi       - polar azimuth coordinate (degrees) of the spacecraft
                measured positive counterclockwise from the earth-moon line
    v0        - initial speed of spacecraft relative to rotating earth-moon
                system (km/s)
    gamma     - initial flight path angle (degrees)
    r0        - intial radial distance of spacecraft from the earth (km)
    x,y       - x and y coordinates of spacecraft in rotating earth-moon
                system (km)
    vx,vy     - x and y components of spacecraft velocity relative to
                rotating earth-moon system (km/s)
    f0        - column vector containing the initial valus of x, y, vx and vy
    t0,tf     - initial time and final times (s)
    t         - column vector of times at which the solution was computed
    f         - a matrix whose columns are:
                column 1: solution for x at the times in t
                column 2: solution for y at the times in t
                column 3: solution for vx at the times in t
                column 4: solution for vy at the times in t
    xf,yf     - x and y coordinates of spacecraft in rotating earth-moon
                system at tf
    vxf, vyf  - x and y components of spacecraft velocity relative to
                rotating earth-moon system at tf
    df        - distance from surface of the moon at tf
    vf        - relative speed at tf

    User h-functions required: rkf45
    User subfunctions required: rates, circle
*/

// Constants
const double days = 24 * 3600;            // Convert days to seconds
const double G = 6.6742e-20;              // Universal gravitational constant (km^3/kg/s^2)
const double rmoon = 1737;                // Radius of the Moon (km)
const double rearth = 6378;               // Radius of the Earth (km)
const double r12 = 384400;                // Earth-Moon distance (km)
const double m1 = 5.974e24;               // Mass of Earth (kg)
const double m2 = 7.348e22;               // Mass of Moon (kg)

const double M = m1 + m2;                 // Total mass of the Earth-Moon system (kg)
const double pi_1 = m1 / M;               // Mass ratio of Earth
const double pi_2 = m2 / M;               // Mass ratio of Moon

const double mu1 = 398600;                // Gravitational parameter of Earth (km^3/s^2)
const double mu2 = 4903.02;               // Gravitational parameter of Moon (km^3/s^2)
const double mu = mu1 + mu2;              // Combined gravitational parameter (km^3/s^2)

const double W = std::sqrt(mu / std::pow(r12, 3)); // Angular velocity (rad/s)
const double x1 = -pi_2 * r12;                     // X-coordinate of Earth
const double x2 = pi_1 * r12;                      // X-coordinate of Moon
#define M_PI 3.14159265358979323846

std::vector<double> rates(double t, const std::vector<double>& f) 
/*
    This subfunction calculates the components of the relative acceleration
    for the restricted 3-body problem, using Equations 2.192a and 2.192b

    ax,ay - x and y components of relative acceleration (km/s^2)
    r1    - spacecraft distance from the earth (km)
    r2    - spacecraft distance from the moon (km)
    f     - column vector containing x, y, vx and vy at time t
    fdt   - column vector containing vx, vy, ax and ay at time t

    All other variables are defined above.

    User h-functions required: none
*/
{
    double x = f[0];
    double y = f[1];
    double vx = f[2];
    double vy = f[3];

    double r1 = std::sqrt(std::pow(x + pi_2 * r12, 2) + std::pow(y, 2));
    double r2 = std::sqrt(std::pow(x - pi_1 * r12, 2) + std::pow(y, 2));

    double ax = 2 * W * vy + std::pow(W, 2) * x - mu1 * (x - x1) / std::pow(r1, 3) - mu2 * (x - x2) / std::pow(r2, 3);
    double ay = -2 * W * vx + std::pow(W, 2) * y - (mu1 / std::pow(r1, 3) + mu2 / std::pow(r2, 3)) * y;

    return {vx, vy, ax, ay};
}

std::vector<std::pair<double, double>> circle(double xc, double yc, double radius)
/*
    This subfunction calculates the coordinates of points spaced
    0.1 degree apart around the circumference of a circle

    x,y    - x and y coordinates of a point on the circumference
    xc,yc  - x and y coordinates of the center of the circle
    radius - radius of the circle
    xy     - an array containing the x coordinates in column 1 and the
             y coordinates in column 2

    User h-functions required: none
*/ 
{
    std::vector<std::pair<double, double>> points;
    for (double angle = 0; angle <= 360; angle += 0.1) {
        double x = xc + radius * std::cos(angle * M_PI / 180.0);
        double y = yc + radius * std::sin(angle * M_PI / 180.0);
        points.emplace_back(x, y);
    }
    return points;
}

// Main function
int main() 
{
    // Input data
    double d0 = 200;                       // Initial altitude of spacecraft (km)
    double phi = -90;                      // Initial polar azimuth angle (degrees)
    double v0 = 10.9148;                   // Initial speed (km/s)
    double gamma = 20;                     // Initial flight path angle (degrees)
    double t0 = 0;                         // Initial time (s)
    double tf = 3.16689 * days;            // Final time (s)

    double r0 = rearth + d0;               // Initial radial distance from Earth (km)
    double x = r0 * std::cos(phi * M_PI / 180.0) + x1;
    double y = r0 * std::sin(phi * M_PI / 180.0);

    double vx = v0 * (std::sin(gamma * M_PI / 180.0) * std::cos(phi * M_PI / 180.0) - 
                      std::cos(gamma * M_PI / 180.0) * std::sin(phi * M_PI / 180.0));
    double vy = v0 * (std::sin(gamma * M_PI / 180.0) * std::sin(phi * M_PI / 180.0) + 
                      std::cos(gamma * M_PI / 180.0) * std::cos(phi * M_PI / 180.0));

    std::vector<double> f0 = {x, y, vx, vy};

    // Compute trajectory
    RKF45 result = rkf45(rates, {t0, tf}, f0, 1e-8);

    // Output results
    double xf = result.yout.back()[0];
    double yf = result.yout.back()[1];
    double vxf = result.yout.back()[2];
    double vyf = result.yout.back()[3];

    double df = std::sqrt(std::pow(xf - x2, 2) + std::pow(yf, 2)) - rmoon;
    double vf = std::sqrt(std::pow(vxf, 2) + std::pow(vyf, 2));

    std::cout << "-----------------------------------------------------------\n";
    std::cout << " Lunar trajectory using the restricted\n";
    std::cout << " three-body equations.\n";
    std::cout << "\n Initial Earth altitude (km)         = " << d0;
    std::cout << "\n Initial angle between radial\n and earth-moon line (degrees)       = " << phi;
    std::cout << "\n Initial flight path angle (degrees) = " << gamma;
    std::cout << "\n Flight time (days)                  = " << tf / days;
    std::cout << "\n Final distance from the moon (km)   = " << df;
    std::cout << "\n Final relative speed (km/s)         = " << vf;
    std::cout << "\n-----------------------------------------------------------\n";

    return 0;
}
