#include <iostream>
#include <vector>
#include <cmath>
#include "rk1_4.h"

/*
    This function uses the RK1 through RK4 methods with two
    different time steps each to solve for and outputs the response
    of a damped single degree of freedom spring-mass system to
    a sinusoidal forcing function, represented by

    x'' 2*z*wn*x' + wn^2*x = (Fo/m)*sin(w*t)
    
    The numerical integration is done by the external
    function 'rk1_4', which uses the subfunction 'rates'
    herein to compute the derivatives.

    This function also prints the exact solution for comparison.

    x           - displacement (m)
    '           - shorthand for d/dt
    t           - time (s)
    wn          - natural circular frequency (radians/s)
    z           - damping factor
    wd          - damped natural frequency
    Fo          - amplitude of the sinusoidal forcing function (N)
    m           - mass (kg)
    w           - forcing frequency (radians/s)
    t0          - initial time (s)
    tf          - final time (s)
    h           - uniform time step (s)
    tspan       - a row vector containing t0 and tf
    x0          - value of x at t0 (m)
    x_dot0      - value of dx/dt at t0 (m/s)
    f0          - column vector containing x0 and x_dot0
    rk          - = 1 for RK1; = 2 for RK2; = 3 for RK3; = 4 for RK4
    t           - solution times for the exact solution
    t1, ...,t4  - solution times for RK1,...,RK4 for smaller
    t11,...,t41 - solution times for RK1,...,RK4 for larger h
    f1, ...,f4  - solution vectors for RK1,...,RK4 for smaller h
    f11,...,f41 - solution vectors for RK1,...,RK4 for larger h

    User h-functions required: rk1_4
    User subfunctions required: rates
*/

std::vector<double> rates(double t, const std::vector<double> &f)
{
    //...Input data:
    double m = 1.0;
    double z = 0.03;
    double wn = 1.0;
    double Fo = 1.0;
    double w = 0.4 * wn;

    double x = f[0];
    double Dx = f[1];
    double D2x = (Fo / m) * std::sin(w * t) - 2 * z * wn * Dx - wn * wn * x;
    //...End input data

    return {Dx, D2x};
}

void outputResults(const RKResult &result) {
    for (size_t i = 0; i < result.tout.size(); ++i) {
        std::cout << "Time: " << result.tout[i] << ", X: " << result.yout[i][0]
                  << ", DX: " << result.yout[i][1] << std::endl;
    }
}

int main() {
    double t0 = 0.0;
    double tf = 110.0;
    std::vector<double> y0 = {0.0, 0.0};
    std::pair<double, double> tspan = {t0, tf};

    //...Solve using RK1 through RK4, using the same and a larger
    // time step for each method:
    std::vector<int> methods = {1, 2, 3, 4};
    std::vector<double> stepSizes = {0.01, 0.1, 0.5, 1.0, 2.0};

    for (int rk : methods) {
        std::cout << "\nRK" << rk << " Method:\n";
        for (double h : stepSizes) {
            std::cout << "\nStep size: " << h << std::endl;
            RKResult result = RK1_4::integrate(rates, tspan, y0, h, rk);
            outputResults(result);
        }
    }

    return 0;
}