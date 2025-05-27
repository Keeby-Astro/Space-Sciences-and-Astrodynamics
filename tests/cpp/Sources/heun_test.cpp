#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "heun.h"

/*
    This program uses Heun's method with two different time steps to solve
    for and plot the response of a damped single degree of freedom
    spring-mass system to a sinusoidal forcing function, represented by

    x'' + 2*z*wn*x' + wn^2*x = (Fo/m)*sin(w*t)

    The numerical integration is done in the external function 'heun',
    which uses the subfunction 'rates' herein to compute the derivatives.

    x     - displacement (m)
    '     - shorthand for d/dt
    t     - time (s)
    wn    - natural circular frequency (radians/s)
    z     - damping factor
    Fo    - amplitude of the sinusoidal forcing function (N)
    m     - mass (kg)
    w     - forcing frequency (radians/s)
    t0    - initial time (s)
    tf    - final time (s)
    h     - uniform time step (s)
    tspan - row vector containing t0 and tf
    x0    - value of x at t0 (m)
    Dx0   - value of dx/dt at t0 (m/s)
    f0    - column vector containing x0 and Dx0
    t     - column vector of times at which the solution was computed
    f     - a matrix whose columns are:
            column 1: solution for x at the times in t
            column 2: solution for x' at the times in t

    User h-functions required: heun
    User subfunctions required: rates
*/

std::vector<double> rates(double t, const std::vector<double>& f) \
{
    //...System properties:
    static const double m = 1.0;
    static const double z = 0.03;
    static const double wn = 1.0;
    static const double Fo = 1.0;
    static const double w = 0.4 * wn;

    //...Initial conditions:
    double x = f[0];
    double Dx = f[1];
    double D2x = (Fo / m) * std::sin(w * t) - 2 * z * wn * Dx - wn * wn * x;

    return {Dx, D2x};
}

void plot(const std::vector<double>& t1, const std::vector<std::vector<double>>& f1,
          const std::vector<double>& t2, const std::vector<std::vector<double>>& f2) {
    std::cout << "time(s)\t x(m) (h = 1.0)\t x(m) (h = 0.1)\n";
    for (size_t i = 0; i < t1.size() && i < t2.size(); ++i) {
        std::cout << t1[i] << "\t" << f1[i][0] << "\t\t" << f2[i][0] << "\n";
    }
}

int main() 
{
    //...Time range:
    double t0 = 0.0;
    double tf = 110.0;
    std::pair<double, double> tspan(t0, tf);

    std::vector<double> f0 = {0.0, 0.0};

    //...Calculate and plot the solution for h = 1.0:
    double h1 = 1.0;
    HeunResult result1 = heun(rates, tspan, f0, h1);

    //...Calculate and plot the solution for h = 0.1:
    double h2 = 0.1;
    HeunResult result2 = heun(rates, tspan, f0, h2);

    // Output the results
    plot(result1.tout, result1.yout, result2.tout, result2.yout);

    return 0;
}
