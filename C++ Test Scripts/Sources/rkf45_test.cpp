/*
{
rkf45_test.cpp

    This program uses RKF4(5) with adaptive step size control
    to solve the differential equation

        x'' + mu/x^2 = 0

    The numerical integration is done by the function 'rkf45' (included via
    "rkf45.cpp"), which calls the subfunction 'rates' herein to compute the
    derivatives.

    x     - displacement (km)
    '     - shorthand for d/dt
    t     - time (s)
    mu    - = go*RE^2 (km^3/s^2), where go is the sea level gravitational
              acceleration and RE is the radius of the earth
    x0    - initial value of x
    v0    - initial value of the velocity (x')
    y0    - column vector containing x0 and v0
    t0    - initial time
    tf    - final time
    tspan - a row vector with components t0 and tf
    t     - column vector of the times at which the solution is found
    f     - a matrix whose columns are:
            column 1: solution for x at the times in t
            column 2: solution for x' at the times in t

    User M-function required: rkf45 (in "rkf45.cpp")
    User subfunction required: rates
}
*/

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "rkf45.h"

/*
{
    This function calculates first and second time derivatives of x
    governed by the equation of two-body rectilinear motion.

    x'' + mu/x^2 = 0

    Dx   - velocity x'
    D2x  - acceleration x''
    f    - column vector containing x and Dx at time t
    dfdt - column vector containing Dx and D2x at time t

    User M-functions required: none
}
*/
std::vector<double> rates(double t, const std::vector<double> &f)
{
    double x = f[0];
    double Dx = f[1];

    static double mu = 398600.0;
    double D2x = -mu / (x * x);

    // Return dfdt
    std::vector<double> dfdt(2);
    dfdt[0] = Dx;
    dfdt[1] = D2x;
    return dfdt;
}

void plotit(const std::vector<double> &t,
            const std::vector<std::vector<double>> &f,
            double minutes)
{
    // Print the results
    std::cout << "\n------------------------------------------\n";
    std::cout << "Position vs time (km vs. minutes):\n";
    std::cout << "   time(min)     position(km)\n";
    std::cout << "------------------------------------------\n";
    for (size_t i = 0; i < t.size(); ++i)
    {
        double timeInMin = t[i] / minutes;
        double xVal = f[i][0];
        std::cout << std::setw(12) << timeInMin << "    "
                  << std::setw(12) << xVal << "\n";
    }

    std::cout << "\n------------------------------------------\n";
    std::cout << "Velocity vs time (km/s vs. minutes):\n";
    std::cout << "   time(min)     velocity(km/s)\n";
    std::cout << "------------------------------------------\n";
    for (size_t i = 0; i < t.size(); ++i)
    {
        double timeInMin = t[i] / minutes;
        double vVal = f[i][1];
        std::cout << std::setw(12) << timeInMin << "    "
                  << std::setw(12) << vVal << "\n";
    }
}

int main()
{
    // mu      = 398600;
    double minutes = 60.0; // Conversion from minutes to seconds

    std::vector<double> y0{6500.0, 7.8};

    double t0 = 0.0;
    double tf = 70.0 * minutes;

    std::pair<double, double> tspan(t0, tf);

    RKF45 result = rkf45(rates, tspan, y0, 1.e-8);

    plotit(result.tout, result.yout, minutes);

    return 0;
}