#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "rkf45.h"

// Universal constants
const double HOURS_TO_SECONDS = 3600.0;
const double G = 6.6742e-20;

// Input data
const double EARTH_MASS = 5.974e24;
const double EARTH_RADIUS = 6378.0;
const double SPACECRAFT_MASS = 1000.0;
const std::vector<double> INITIAL_POSITION = {8000.0, 0.0, 6000.0};
const std::vector<double> INITIAL_VELOCITY = {0.0, 7.0, 0.0};
const double START_TIME = 0.0;
const double END_TIME = 4.0 * HOURS_TO_SECONDS;

std::vector<double> rates(double t, const std::vector<double> &f)
{
    // Extract components from the state vector
    double x = f[0];
    double y = f[1];
    double z = f[2];
    double vx = f[3];
    double vy = f[4];
    double vz = f[5];

    // Compute the magnitude of the position vector
    double r = std::sqrt(x * x + y * y + z * z);

    // Compute the gravitational parameter
    static double mu = G * (EARTH_MASS + SPACECRAFT_MASS);

    // Compute accelerations
    double ax = -mu * x / (r * r * r);
    double ay = -mu * y / (r * r * r);
    double az = -mu * z / (r * r * r);

    // Return derivatives
    return {vx, vy, vz, ax, ay, az};
}

void output(const std::vector<double> &t, const std::vector<std::vector<double>> &y)
{
    std::vector<double> r(t.size());
    for (size_t i = 0; i < t.size(); ++i)
    {
        r[i] = std::sqrt(y[i][0] * y[i][0] + y[i][1] * y[i][1] + y[i][2] * y[i][2]);
    }

    auto max_it = std::max_element(r.begin(), r.end());
    auto min_it = std::min_element(r.begin(), r.end());

    double rmax = *max_it;
    double rmin = *min_it;
    size_t imax = std::distance(r.begin(), max_it);
    size_t imin = std::distance(r.begin(), min_it);

    double v_at_rmax = std::sqrt(y[imax][3] * y[imax][3] + y[imax][4] * y[imax][4] + y[imax][5] * y[imax][5]);
    double v_at_rmin = std::sqrt(y[imin][3] * y[imin][3] + y[imin][4] * y[imin][4] + y[imin][5] * y[imin][5]);

    std::cout << "\n--------------------------------------------------------------\n";
    std::cout << "\n Earth Orbit\n";
    std::cout << "\n The initial position is [" << INITIAL_POSITION[0] << ", " << INITIAL_POSITION[1] << ", " << INITIAL_POSITION[2] << "] (km).";
    std::cout << "\n Magnitude = " << std::sqrt(INITIAL_POSITION[0] * INITIAL_POSITION[0] + INITIAL_POSITION[1] * INITIAL_POSITION[1] + INITIAL_POSITION[2] * INITIAL_POSITION[2]) << " km\n";
    std::cout << "\n The initial velocity is [" << INITIAL_VELOCITY[0] << ", " << INITIAL_VELOCITY[1] << ", " << INITIAL_VELOCITY[2] << "] (km/s).";
    std::cout << "\n Magnitude = " << std::sqrt(INITIAL_VELOCITY[0] * INITIAL_VELOCITY[0] + INITIAL_VELOCITY[1] * INITIAL_VELOCITY[1] + INITIAL_VELOCITY[2] * INITIAL_VELOCITY[2]) << " km/s\n";
    std::cout << "\n Initial time = " << START_TIME / HOURS_TO_SECONDS << " h.\n Final time = " << END_TIME / HOURS_TO_SECONDS << " h.\n";
    std::cout << "\n The minimum altitude is " << rmin - EARTH_RADIUS << " km at time = " << t[imin] / HOURS_TO_SECONDS << " h.";
    std::cout << "\n The speed at that point is " << v_at_rmin << " km/s.\n";
    std::cout << "\n The maximum altitude is " << rmax - EARTH_RADIUS << " km at time = " << t[imax] / HOURS_TO_SECONDS << " h.";
    std::cout << "\n The speed at that point is " << v_at_rmax << " km/s\n";
    std::cout << "\n--------------------------------------------------------------\n\n";

    std::cout << "Orbit coordinates:\n";
    for (size_t i = 0; i < t.size(); ++i)
    {
        std::cout << "t = " << t[i] / HOURS_TO_SECONDS << " h: ("
                << y[i][0] << ", " << y[i][1] << ", " << y[i][2] << ")\n";
    }

    std::cout << "\nStarting point (o): (" << y[0][0] << ", " << y[0][1] << ", " << y[0][2] << ")\n";
    std::cout << "Final point (f): (" << y.back()[0] << ", " << y.back()[1] << ", " << y.back()[2] << ")\n";
}

int main()
{
    // Initial state vector
    std::vector<double> y0 = {INITIAL_POSITION[0], INITIAL_POSITION[1], INITIAL_POSITION[2], INITIAL_VELOCITY[0], INITIAL_VELOCITY[1], INITIAL_VELOCITY[2]};

    // Time span
    std::pair<double, double> tspan(START_TIME, END_TIME);

    // Solve using RKF45
    RKF45 result = rkf45(rates, tspan, y0, 1.e-8);

    // Output results
    output(result.tout, result.yout);

    return 0;
}
