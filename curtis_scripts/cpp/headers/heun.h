// ALGORITHM 1.2: NUMERICAL INTEGRATION BY HEUN'S PREDICTOR-CORRECTOR METHOD

#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <limits>

struct HeunResult 
{
    std::vector<double> tout;
    std::vector<std::vector<double>> yout;
};

/*
    This function uses the predictor-corrector method to integrate a system
    of first-order differential equations dy/dt = f(t,y).

    y               - column vector of solutions
    f               - column vector of the derivatives dy/dt
    ode_function    - handle for the user M-function in which the derivatives
                      f are computed
    t               - time
    t0              - initial time
    tf              - final time
    tspan           - the vector [t0 tf] giving the time interval for the
                      solution
    h               - time step
    y0              - column vector of initial values of the vector y
    tout            - column vector of the times at which y was evaluated
    yout            - a matrix, each row of which contains the components of y
                      evaluated at the correponding time in tout
    feval           - a built-in MATLAB function which executes 'ode_function'
                      at the arguments t and y
    tol             - Maximum allowable relative error for determining
                      convergence of the corrector
    itermax         - maximum allowable number of iterations for corrector
                      convergence
    iter            - iteration number in the corrector convergence loop
    t1              - time at the beginning of a time step
    y1              - value of y at the beginning of a time step
    f1              - derivative of y at the beginning of a time step
    f2              - derivative of y at the end of a time step
    favg            - average of f1 and f2
    y2p             - predicted value of y at the end of a time step
    y2              - corrected value of y at the end of a time step
    err             - maximum relative error (for all components) between y2p
                      and y2 for given iteration
    eps             - unit roundoff error (the smallest number for which
                      1 + eps > 1). Used to avoid a zero denominator.
*/

HeunResult heun(
    std::function<std::vector<double>(double, const std::vector<double>&)> ode_function,
    const std::pair<double, double>& tspan,
    const std::vector<double>& y0,
    double h
) 
{
    const double tol = 1.e-6;
    const int itermax = 100;

    double t0 = tspan.first;
    double tf = tspan.second;

    double t = t0;
    std::vector<double> y = y0;

    HeunResult result;
    result.tout.push_back(t);
    result.yout.push_back(y);

    while (t < tf) 
    {
        h = std::min(h, tf - t);
        double t1 = t;
        std::vector<double> y1 = y;
        std::vector<double> f1 = ode_function(t1, y1);

        std::vector<double> y2 = y1;
        for (size_t i = 0; i < y1.size(); ++i) 
        {
            y2[i] += f1[i] * h;
        }

        double t2 = t1 + h;
        double err = tol + 1;
        int iter = 0;

        while (err > tol && iter <= itermax) 
        {
            std::vector<double> y2p = y2;
            std::vector<double> f2 = ode_function(t2, y2p);

            std::vector<double> favg(y1.size());
            for (size_t i = 0; i < y1.size(); ++i) 
            {
                favg[i] = (f1[i] + f2[i]) / 2;
                y2[i] = y1[i] + favg[i] * h;
            }

            err = 0.0;
            for (size_t i = 0; i < y1.size(); ++i) 
            {
                err = std::max(err, std::fabs((y2[i] - y2p[i]) / (y2[i] + std::numeric_limits<double>::epsilon())));
            }
            ++iter;
        }

        if (iter > itermax) 
        {
            std::cerr << "\n Maximum number of iterations exceeded at time = " << t << " in function 'heun'.\n";
            break;
        }

        t += h;
        y = y2;

        result.tout.push_back(t); // adds t to the bottom of the column vector tout
        result.yout.push_back(y); // adds y to the bottom of the matrix yout
    }

    return result;
}
